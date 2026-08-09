!=======================================================================
!     driver_pmg — 합성 3D Poisson 으로 PMG(BiCGSTAB+MG 예조건) 단독 검증
!     (LOOP C007)
!
!     프로덕션(read_grid.f90/pressure_solve.f90)과 동일한 호출 체인:
!       read_input_mg → subdomain_infor_MG → read_mesh_MPI
!       → Prep_fine_P → Prep_MG_GarL → assemble_FVM → SOLVE_GMG
!
!     사용법:  driver_pmg [nx ny nz aspect]   (기본 24 24 24 1.0)
!       - aspect > 1 : z 방향 격자 늘림 → 이방성 문제 (semi-coarsening 검증)
!       - 실행 디렉토리에 mg.in 필요. MG_tmp/ 는 체인이 생성.
!
!     판정 (LOOP.md §2 L2):
!       - 독립 residual  ‖b−A·u‖/‖b‖ ≤ 10*crit  (솔버 자기보고에 의존 안 함)
!       - 제작해(manufactured solution) 오차  ‖u−u_ex‖/‖u_ex‖  보고
!       - PASS → exit 0, FAIL → exit 1.  its 는 fort.501 에서 러너가 파싱.
!=======================================================================
      PROGRAM driver_pmg
!
      USE Zcore,        ONLY: np_z => np, myrank_z => myrank
      USE Zparam,       ONLY: ndim_z => ndim, ns
      USE Zbicg,        ONLY: eps_bicg
      USE Zcoord1,      ONLY: xloc_tmp
      USE Zmpi,         ONLY: celem
      USE md_geometry,  ONLY: nelem_mg, num_neigh_mg, neigh_mg, coord
      USE MD_parameter, ONLY: nf_max
      USE MD_MPI,       ONLY: nintf
      USE MD_matrix,    ONLY: nnz, ia, ja, au, u, b
!
      IMPLICIT NONE
      INCLUDE 'mpif.h'
!
      INTEGER :: ierr, nx, ny, nz, ncell, i, j, k, c, kk, irow, jcol
      INTEGER :: nargs
      REAL(8) :: aspect, dist2, coef, vol, diag_const, bn, rn, en, uexn
      REAL(8) :: rel_res, rel_err, d(3)
      CHARACTER(32) :: arg
      REAL(8), ALLOCATABLE :: diag_my(:), au_my(:), src(:), uex(:), r(:)
      REAL(8), ALLOCATABLE :: mycoord(:,:)   ! coord 스냅샷 (Prep_fine_P 가 해제하므로)
!
!.....MPI + CUPID 측 최소 상태 (프로덕션에서 본체가 하는 일)
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank_z, ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, np_z, ierr)
      IF (np_z .NE. 1) STOP 'driver_pmg: serial only (np=1) at C007'
!
      ndim_z   = 3
      ns       = 1
      eps_bicg = 1.d-8          ! read_input_mg: crit = eps_bicg → crit_bcg_mg
!
!.....합성 격자 파라미터
      nx = 24; ny = 24; nz = 24; aspect = 1.d0
      nargs = COMMAND_ARGUMENT_COUNT()
      IF (nargs .GE. 3) THEN
         CALL GET_COMMAND_ARGUMENT(1, arg); READ (arg, *) nx
         CALL GET_COMMAND_ARGUMENT(2, arg); READ (arg, *) ny
         CALL GET_COMMAND_ARGUMENT(3, arg); READ (arg, *) nz
      END IF
      IF (nargs .GE. 4) THEN
         CALL GET_COMMAND_ARGUMENT(4, arg); READ (arg, *) aspect
      END IF
      ncell = nx*ny*nz
      WRITE (*, '(A,3I5,A,F8.2,A,I9)') ' C007 grid:', nx, ny, nz,          &
             '  aspect=', aspect, '  ncell=', ncell
!
!.....GMG 셋업 입력 채우기 (프로덕션에선 read_grid/subdomain_info_ser 가 채움)
!     - 7점 스텐실 인접성, 좌표(z 는 aspect 로 늘림), 단일 도메인
      nelem_mg = ncell
      nf_max   = 6
      ALLOCATE (celem(ncell));  celem = 1
      ALLOCATE (xloc_tmp(ncell, 3))
      ALLOCATE (num_neigh_mg(ncell), neigh_mg(nf_max, ncell))
      num_neigh_mg = 0
      neigh_mg     = 0
      DO k = 1, nz
      DO j = 1, ny
      DO i = 1, nx
         c = i + (j-1)*nx + (k-1)*nx*ny
         xloc_tmp(c, 1) = (DBLE(i) - 0.5d0)
         xloc_tmp(c, 2) = (DBLE(j) - 0.5d0)
         xloc_tmp(c, 3) = (DBLE(k) - 0.5d0)*aspect
         kk = 0
         IF (i .GT. 1)  THEN; kk = kk+1; neigh_mg(kk, c) = c - 1;     END IF
         IF (i .LT. nx) THEN; kk = kk+1; neigh_mg(kk, c) = c + 1;     END IF
         IF (j .GT. 1)  THEN; kk = kk+1; neigh_mg(kk, c) = c - nx;    END IF
         IF (j .LT. ny) THEN; kk = kk+1; neigh_mg(kk, c) = c + nx;    END IF
         IF (k .GT. 1)  THEN; kk = kk+1; neigh_mg(kk, c) = c - nx*ny; END IF
         IF (k .LT. nz) THEN; kk = kk+1; neigh_mg(kk, c) = c + nx*ny; END IF
         num_neigh_mg(c) = kk
      END DO
      END DO
      END DO
!
!.....프로덕션 초기화 체인 (read_grid.f90 순서 그대로)
      CALL read_input_mg
      CALL subdomain_infor_MG
      CALL read_mesh_MPI
!.....GMG 로컬 순서의 좌표 스냅샷 — Prep_fine_P 가 coord 를 해제함 (3_Prep_fine_P.f90:130)
      ALLOCATE (mycoord(3, ncell))
      mycoord(:, 1:ncell) = coord(:, 1:ncell)
      CALL Prep_fine_P
      CALL Prep_MG_GarL
!
      WRITE (*, '(A,I9,A,I10)') ' C007 setup done: nintf=', nintf,         &
             '  nnz(off-diag)=', nnz
      IF (nintf .NE. ncell) STOP 'driver_pmg: nintf /= ncell (mapping?)'
!
!.....행렬 값 생성 — 셋업 후 GMG 자신의 ia/ja/coord 기준 (재배열에 무관)
!     GMG CSR 은 대각을 au 안에 포함 (ju 가 대각 위치; mt_amux/resi_normP 는
!     별도 diag 항 없이 full-row 로 A·x 계산). assemble_FVM 의 diag 인자는
!     스무더 스케일링용 (diagt = 1/diag).
!     off-diag 계수: -vol/dist^2 (축정렬 균일 늘림 격자에서 -area/dist 와 동치)
!     diag: Dirichlet ghost(=0) 포함 전면 합 = 2cx+2cy+2cz (상수) → 정칙
      vol = aspect                              ! dx=dy=1, dz=aspect
      diag_const = 2.d0*aspect + 2.d0*aspect + 2.d0/aspect
      ALLOCATE (diag_my(ncell), au_my(nnz), src(ncell), uex(ncell), r(ncell))
      diag_my = diag_const
      DO irow = 1, ncell
         DO kk = ia(irow), ia(irow+1) - 1
            jcol = ja(kk)
            IF (jcol .EQ. irow) THEN
               au_my(kk) = diag_const
            ELSE
               d(:) = mycoord(:, irow) - mycoord(:, jcol)
               dist2 = d(1)*d(1) + d(2)*d(2) + d(3)*d(3)
               au_my(kk) = -vol/dist2
            END IF
         END DO
      END DO
!
!.....제작해:  u_ex 결정적 지정 → b = A·u_ex  (해가 알려진 문제, full-row)
      DO irow = 1, ncell
         uex(irow) = 1.d0 + SIN(0.5d0*DBLE(irow))
      END DO
      DO irow = 1, ncell
         src(irow) = 0.d0
         DO kk = ia(irow), ia(irow+1) - 1
            src(irow) = src(irow) + au_my(kk)*uex(ja(kk))
         END DO
      END DO
!
!.....솔브 (pressure_solve.f90 의 MG 경로 그대로)
      u = 0.d0
      CALL assemble_FVM(1, nnz, src, au_my, diag_my)
      CALL SOLVE_GMG(1)
!
!.....독립 검증 — 드라이버 보유 행렬로 residual 직접 계산
      bn = 0.d0; rn = 0.d0; en = 0.d0; uexn = 0.d0
      DO irow = 1, ncell
         r(irow) = src(irow)
         DO kk = ia(irow), ia(irow+1) - 1
            r(irow) = r(irow) - au_my(kk)*u(ja(kk))
         END DO
         bn   = bn   + src(irow)*src(irow)
         rn   = rn   + r(irow)*r(irow)
         en   = en   + (u(irow) - uex(irow))**2
         uexn = uexn + uex(irow)*uex(irow)
      END DO
      rel_res = SQRT(rn/bn)
      rel_err = SQRT(en/uexn)
!
      WRITE (*, '(A,ES12.4,A,ES12.4)') ' C007-RESULT  rel_res=', rel_res,  &
             '  rel_err=', rel_err
      IF (rel_res .LE. 10.d0*eps_bicg) THEN
         WRITE (*, '(A)') ' C007-VERDICT PASS'
         CALL MPI_FINALIZE(ierr)
      ELSE
         WRITE (*, '(A)') ' C007-VERDICT FAIL'
         CALL MPI_FINALIZE(ierr)
         CALL EXIT(1)
      END IF
!
      END PROGRAM driver_pmg
