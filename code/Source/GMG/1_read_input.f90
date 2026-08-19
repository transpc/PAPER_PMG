      SUBROUTINE read_input_mg
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
!  PMG 입력 (클리닝 P5 재작성)
!
!  고정 구성(현 조합 전용 — 코드 상수, mg.in 으로 변경 불가):
!    isol_mg=-2(BiCGSTAB+MG 예조건자), POL(Chebyshev)+Gershgorin λ상계,
!    mdf_matrix=1, n_GC=1, i_dir=1(직접해), icommu=2, igather=1,
!    ioplv=1(레벨 자동), 예조건자 = V-cycle 1회
!
!  mg.in (선택 — 없으면 전부 기본값):
!    &MG_tuning  teta, teta_p, alpha, itergs, icheb, ip_nmax, ip_inter, ip_lev
!    &MG_options il1_gs, isetup_comm, nthre
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
         USE MD_parameter, ONLY: ndim, maxit, mdf_matrix, ip_inter, ndom, &
                                 crit, teta, teta_p, alpha
         USE MD_MG_index, ONLY: nlevel, ncycle, mxnbne, maxit_1, ip_nmax, &
                                iter_mg, n_GC, crit_1, &
                                iter_max, nlevel_N, n1_min, n2_min, ioplv, ip_lev, &
                                itergs, icheb, icase_MG, &
                                crit_bcg_mg, il1_gs, isetup_comm
         USE MD_MPI, ONLY: myrank, myrankt
         USE MD_OpenMP, ONLY: nthre
         USE MD_MG_Global_C, ONLY: i_dir, nlv_glo, nlv_glomax, igather
         USE MD_MPI_MG, ONLY: icommu
         USE Zparam, ONLY: ndim_cupid => ndim
         USE Zbicg, ONLY: eps_mg => eps_bicg
         USE Zcore, ONLY: myrank_mg => myrank, np_mg => np
!
         IMPLICIT NONE
!
         INTEGER i, iu, ios
!
         NAMELIST /MG_tuning/ teta, teta_p, alpha, itergs, icheb, &
                              ip_nmax, ip_inter, ip_lev
         NAMELIST /MG_options/ il1_gs, isetup_comm, nthre
!
         myrank = myrank_mg
         myrankt = myrank_mg
         ndim = ndim_cupid
         ndom = np_mg
!
! ---- 고정 구성 ----
         mdf_matrix = 1
         icase_MG = 2
         n_GC = 1
         i_dir = 1
         icommu = 2
         igather = 1
         ioplv = 1
         n1_min = 100
         n2_min = 5
         mxnbne = 100
         maxit = 1000
         maxit_1 = 1000
         crit_1 = 1.d-1
! ---- 예조건자 모드 (구 isol_mg=-2 경로의 net 효과) ----
         crit_bcg_mg = eps_mg       ! BiCGSTAB 수렴 판정
         crit = 1.d-1               ! V-cycle 내부 판정 (구 crit_pre)
         ncycle = 1                 ! 예조건자 적용당 V-cycle 수 (구 ncycle_pre)
! ---- 튜닝 기본값 ----
         teta = 0.6d0
         teta_p = 0.65d0
         alpha = 0.005d0
         itergs = 0
         itergs(1) = 1
         itergs(2) = 1
         itergs(3) = 2
         icheb = 0
         icheb(1) = 2               ! Chebyshev 반복수 (2~4)
         ip_nmax = 4
         ip_inter = 1
         ip_lev = 1
! ---- 옵션 기본값 ----
         il1_gs = 0                 ! (1) l1-보정 코어스 GS (Baker et al. 2011)
         isetup_comm = 0            ! (1) 셋업 분배 MPI 통신 모드
         nthre = 1                  ! OpenMP 스레드 수
!
! ---- mg.in (선택) : 그룹 순서 무관, 부재 허용 ----
         OPEN (newunit=iu, file='mg.in', status='old', action='read', iostat=ios)
         IF (ios == 0) THEN
            READ (iu, nml=MG_tuning, iostat=ios)
            REWIND (iu)
            READ (iu, nml=MG_options, iostat=ios)
            CLOSE (iu)
         ELSE
            IF (myrank == 0) WRITE (*, *) 'read_input_mg: mg.in not found - using defaults'
         END IF
!
! ---- 레벨 자동 선택 (ioplv=1 로직) ----
         nlevel = 20
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
         nlevel = nlevel - nlv_glomax
         nlevel_N = nlevel
         nlv_glo = nlv_glomax
!
         iter_mg = itergs(1)
         iter_max = iter_mg
! itergs: 미지정 레벨은 직전 값 상속
         DO i = 2, 20
            IF (itergs(i) == 0) THEN
               itergs(i) = itergs(i - 1)
            END IF
         END DO
!
         RETURN
      END
