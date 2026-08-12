      SUBROUTINE read_input_mg
! ---
         USE MD_parameter, ONLY: ndim, maxit, mdf_matrix, ip_inter, ndom, ipar, &
                                 isemi, crit, teta, teta_p, alpha
         USE MD_MG_index, ONLY: nlevel, ncycle, mxnbne, maxit_1, ip_nmax, &
                                iter_mg, n_GC, isth, AR_hi, crit_1, relax, &
                                iter_max, nlevel_N, n1_min, n2_min, ioplv, ip_lev, &
                                isol_mg, id_GS_sym, itergs, icheb, icase_MG, &
                                crit_bcg_mg, ihybrid, isol_start
         USE MD_MPI, ONLY: myrank, myrankt
         USE MD_MG_Global_C, ONLY: i_dir, nlv_glo, nlv_glomax, igather
         USE MD_MPI_MG, ONLY: icommu, iGS, nGS, iallocate_c
         USE MD_OpenMP
         USE Zparam, ONLY:  ndim_cupid => ndim
         USE Zbicg, ONLY: eps_mg => eps_bicg
         USE Zcore, ONLY: myrank_mg => myrank, np_mg => np
!
         IMPLICIT NONE
! = = = = = = = = = = = = = = = = = = = = = = = !
!       control value for PMG solver            !
!                                               !
!       crit: criteria convergence              !
!       crit = res/res0                         !
!       mdf_matrix =1 => divided by diagonal    !
!     n_GC=1 =>global domain on coarsest level  !
!     ipar=1 => pardiso for coarsest level      !
!     isth=0:GS smoothing, isth=1:CG smoothing  !
!      iVcy = 0: new V-cycle:                   !
!       only 1 iter. on second level            !
!     iter_mg: No. iteration for smoothing      !
!     nlevel: no. level for PMG                 !
! crit_1:criteria convergence for the coarsest  !
!   levle in case of ipar = 0 (using CG)        !
!    teta: using for semi-coarsing and          !
!      for interpolation
! = = = = = = = = = = = = = = = = = = = = = = = !
! FOR EACH CASE, USING SUITABLE nlevle and teta !
! - - - - - - - - - - - - - - - - - - - - - - - !

!
         INTEGER i, ncycle_pre, iu, ios
         CHARACTER :: smothing*3
         REAL(8) ::  crit_pre

!
         NAMELIST /MG_method/ isol_mg, icase_MG, mdf_matrix, isol_start
         NAMELIST /MG_level/ nlevel, nlv_glomax, ioplv
         NAMELIST /MG_coarsening/ teta, isemi, AR_hi
         NAMELIST /MG_smoothing/ smothing, itergs, id_GS_sym, relax, icheb
         NAMELIST /MG_interpolation/ teta_p, ip_nmax, ip_inter, ip_lev, alpha
         NAMELIST /MG_coarsest/ n_GC, i_dir, ipar, iGS, nGS, iallocate_c, crit_1
         NAMELIST /MG_MPI/ icommu, igather
         NAMELIST /MG_OpenMP/ ihybrid, nthre
         NAMELIST /MG_Precond/ ncycle_pre, crit_pre
         NAMELIST /MG_more_option/ mxnbne

!
         myrank = myrank_mg
         myrankt = myrank_mg
         ndim = ndim_cupid
         crit = eps_mg
         ndom = np_mg
!
         OPEN (newunit=iu, file='mg.in', status='old', action='read', iostat=ios)
         IF (ios /= 0) THEN
            WRITE (*, *) 'read_input_mg: cannot open mg.in, rank', myrank
            STOP
         ENDIF
!
         READ (iu, nml=MG_method)
         READ (iu, nml=MG_level)
         READ (iu, nml=MG_coarsening)
         READ (iu, nml=MG_smoothing)
         READ (iu, nml=MG_interpolation)
         READ (iu, nml=MG_coarsest)
         READ (iu, nml=MG_MPI)
         READ (iu, nml=MG_OpenMP)
         READ (iu, nml=MG_Precond)
         READ (iu, nml=MG_more_option)

         CLOSE (iu)

!

! set default parameters - - - - - - - -

         maxit = 1000
         maxit_1 = 1000
         ncycle = 500
!      mxnbne = 300
         iter_mg = itergs(1)          ! not use
         iter_max = iter_mg            ! not use

         IF (ioplv .EQ. 1) THEN
!          IF(ndim.EQ.2) THEN
!              n1_min = 100
!          ELSE
            n1_min = 100
!          ENDIF

            n2_min = 5

            nlevel = 20

! set nvl_glo:
            IF (ndom .LE. 10) THEN
               nlv_glomax = 0
            ELSEIF (ndom .LE. 50) THEN
               nlv_glomax = 1
            ELSEIF (ndom .LE. 1000) THEN
               nlv_glomax = 2
            ELSE
               nlv_glomax = 3
            END IF
!
         END IF

! set nvl_glo:
!      IF(ndom.LE.10) THEN
!          nlv_glomax = 0
!      ELSEIF(ndom.LE.30) THEN
!          nlv_glomax = 1
!      ELSEIF(ndom.LE.60) THEN
!          nlv_glomax = 2
!      ELSE
!          nlv_glomax = 3
!      ENDIF

!
         nlevel = nlevel - nlv_glomax
         nlevel_N = nlevel
         nlv_glo = nlv_glomax
!
         IF (n_GC .EQ. 0) THEN
            i_dir = 0
            nlv_glo = 0
         END IF
!
! for itergs:
         DO i = 2, 20
            IF (itergs(i) == 0) THEN
               itergs(i) = itergs(i - 1)
            END IF
         END DO

! for smothing scheme:

         IF (smothing == 'GAS') THEN
            isth = 0
         ELSEIF (smothing == 'JAC') THEN
            isth = 1
         ELSEIF (smothing == 'ILU') THEN
            isth = 2
         ELSEIF (smothing == 'POL') THEN
            isth = 3
         ELSEIF (smothing == 'COG') THEN
            isth = 4
         ELSE
            WRITE (*, *) 'need to provide smothing approach'
            WRITE (999, *) 'need to provide smothing approach'
            PAUSE
            STOP
         END IF

! MG as preconditioner

         IF (isol_mg .LE. 0) THEN
            ncycle = ncycle_pre

            crit_bcg_mg = crit
            crit = crit_pre

            IF ((ncycle .GT. 2) .OR. (crit_1 .LT. 1.d-3)) THEN
               WRITE (*, *) 'no good for MG pre, check ncycle_pre, crit_pre'
            END IF
         END IF

!
         RETURN
      END
