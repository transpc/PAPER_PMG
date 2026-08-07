!
      SUBROUTINE cupid_solvers_vector(diag,src,sol,iswitch_factor)
!  
!     CG solver based on CSR format for vectors, such as vl_n, vg_n
!     Used only for SMAC=3 (or implicit calculation)
!
      USE Zinterface
      USE Zmpi          , ONLY: ncell_fp,maxmt,maxmt_pad,maxmt_lu0,maxmt_lu1, &
                                au,ia_a,ju_a,ja_a,iend,                       &
                                diag_lu,alu0,alu1,ja0,ja1,ia0,ia1,            &
                                ap,iap,jap,jaa,jaar,iaa,ngroup,nbgroup,       &
                                maxmt2,maxmt_r,ia_r,ja_r,ju_r,                &
                                perm_r,permi_r,index_r,                       &
                                levt,lev_typet
      USE Zzone         , ONLY: ncell_fluid,ncell_fluid_pad,ncell_fluid_padv
      USE Zcore         , ONLY: np
      USE Zimplicit     , ONLY: eps_imp_mom,max_iter_mom
      USE Zparam        , ONLY: ndim
!
      IMPLICIT NONE
!
!.....Input
      REAL(8),DIMENSION(ncell_fluid) :: diag
      REAL(8),DIMENSION(ncell_fluid_pad,ndim) :: src
!.....Output
      REAL(8),DIMENSION(ncell_fp,ndim) :: sol
!.....Local variables
      INTEGER :: i,ix
      INTEGER :: izone=0
!.....Local arrrays
      INTEGER,DIMENSION(ndim) :: iswitch_factor
      REAL(8),DIMENSION(ncell_fp) :: sol0
!
      IF(levt.eq.0) THEN
         IF(lev_typet.eq.0) THEN
!
!...........CG solver ilu0 with METIS,CUTHILL reordering
!
            CALL ilupc11(ncell_fluid,maxmt,maxmt_lu0,maxmt_lu1, &
                         diag,au,ia_a,ja_a,ju_a,iend,           &
                         diag_lu,alu0,alu1,ia0,ia1,ja0,ja1)
         ELSE
!
!...........CG solver ilu0 with METIS,CUTHILL reordering
!
            CALL ilupc11p(ncell_fluid,maxmt,maxmt_r,maxmt_lu0,maxmt_lu1, &
                          diag,au,ia_a,ja_a,iend,                        &
                          ia_r,ja_r,ju_r,                                &
                          diag_lu,alu0,alu1,ia0,ia1,ja0,ja1,             &
                          perm_r,permi_r)
         ENDIF
      ELSEIF(levt.gt.0) THEN
         IF(np.eq.1) THEN
            IF(lev_typet.eq.0) THEN
!
!..............CG solver ilup
!
               CALL factor(ncell_fluid,maxmt,maxmt_r,maxmt_lu0,maxmt_lu1, &
                           diag,au,ia_a,ja_a,iend,                        &
                           ia_r,ja_r,ju_r,                                &
                           diag_lu,alu0,alu1,ia0,ia1,ja0,ja1)
            ELSE
!
!..............CG solver ilup with METIS,CUTHILL reordering
!
               CALL factorp(ncell_fluid,maxmt,maxmt2,maxmt_r,maxmt_lu0,maxmt_lu1, &
                            diag,au,ia_a,ja_a,iend,                               &
                            ia_r,ja_r,ju_r,                                       &
                            diag_lu,alu0,alu1,ia0,ia1,ja0,ja1,                    &
                            perm_r,permi_r,index_r)
            ENDIF
         ELSE
            WRITE(*,*) 'cupid_solvers_vector: ilup solver not yet implemented for MPI'
            CALL finalize_mpi
            STOP
         ENDIF
      ELSEIF(levt.eq.-1) THEN
         IF(np.eq.1) THEN
            IF(lev_typet.eq.0) THEN
!
!...........Direct solver no reordering
!
               IF(ndim.eq.2) THEN
                  CALL factor_solvev2(ncell_fluid,ncell_fluid_pad,maxmt,maxmt_r,maxmt_lu1, &
                                      diag,au,ia_a,ja_a,iend,                              &
                                      ia_r,ja_r,ju_r,                                      &
                                      ncell_fp,src,sol)
               ELSE
                  CALL factor_solvev3(ncell_fluid,ncell_fluid_pad,maxmt,maxmt_r,maxmt_lu1, &
                                      diag,au,ia_a,ja_a,iend,                              &
                                      ia_r,ja_r,ju_r,                                      &
                                      ncell_fp,src,sol)
               ENDIF
            ELSE
!
!..............Direct solver with METIS,CUTHILL ordering save LU
!
               IF(ndim.eq.2) THEN
                  CALL factor_solve0pv2(ncell_fluid,ncell_fluid_pad,maxmt,maxmt2,maxmt_r,maxmt_lu1, &
                                        diag,au,ia_a,ja_a,iend,                                     &
                                        ia_r,ja_r,ju_r,                                             &
                                        ncell_fp,src,sol,perm_r,permi_r,index_r)
               ELSE
                  CALL factor_solve0pv3(ncell_fluid,ncell_fluid_pad,maxmt,maxmt2,maxmt_r,maxmt_lu1, &
                                        diag,au,ia_a,ja_a,iend,                                     &
                                        ia_r,ja_r,ju_r,                                             &
                                        ncell_fp,src,sol,perm_r,permi_r,index_r)
               ENDIF
            ENDIF
         ELSE
            WRITE(*,*) 'cupid_solvers_vector: Direct solver not yet implemented for MPI'
            CALL finalize_mpi
            STOP
         ENDIF
      ELSEIF(levt.eq.-2) THEN
         IF(np.eq.1) THEN
            IF(ndim.eq.2) THEN
               CALL factor_solve_tri02(ncell_fluid,ncell_fluid_pad,maxmt, &
                                       ia_a,au,                           &
                                       ncell_fp,src,sol)
            ELSE
               CALL factor_solve_tri03(ncell_fluid,ncell_fluid_pad,maxmt, &
                                       ia_a,au,                           &
                                       ncell_fp,src,sol)
            ENDIF
         ELSE
            WRITE(*,*) 'cupid_solvers_vector: Tridiagonal solvet not yet implemented for MPI'
            CALL finalize_mpi
            STOP
         ENDIF
      ENDIF
!
      IF(levt.ge.0) THEN
!
!........Copy a array in block format for vector processing
!
         CALL copy_a_vector(ncell_fluid,maxmt,maxmt_pad,   &
                            ia_a,au,                       &
                            iap,ap,jaa,iaa,ngroup,nbgroup)
!
!........Apply ILU pre-conditioned bi-conjugate gradient solver
!
         DO ix=1,ndim
            IF(np.eq.1 .and. iswitch_factor(ix).eq.1) CYCLE
            CALL pbcg_ilu(eps_imp_mom,max_iter_mom,ncell_fluid,ncell_fluid_padv,maxmt_pad,maxmt_lu0,maxmt_lu1, &
                          diag,                                                                                &
                          diag_lu,alu0,alu1,ia0,ia1,ja0,ja1,                                                   &
                          iap,jap,ap,jaar,iaa,ngroup,nbgroup,                                                  &
                          lev_typet,perm_r,                                                                    &
                          ncell_fp,src(1,ix),sol0,izone)
            DO i=1,ncell_fp
               sol(i,ix)=sol0(i)
            ENDDO
         ENDDO 
      ENDIF
!
      END SUBROUTINE cupid_solvers_vector
