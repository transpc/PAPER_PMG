!=======================================================================
!     driver_hypre — PMG vs hypre(BiCGSTAB+BoomerAMG) A/B 비교 드라이버
!     (2026-08-20, hypre v2.33.0 / IJ 인터페이스, Fortran 바인딩)
!
!     driver_pmg 와 동일한 셋업 체인으로 동일한 fine CSR 행렬을 만든 뒤,
!     한 프로세스 안에서 두 솔버를 순차 실행한다:
!       (A) PMG:   assemble_FVM → SOLVE_GMG  (프로덕션 경로 그대로)
!       (B) hypre: IJMatrix/IJVector 조립 → BiCGSTAB + BoomerAMG 예조건
!
!     사용법 (driver_pmg 와 동일):
!       (1) 합성:  driver_hypre [nx ny nz aspect [ipart]]   (기본 24 24 24 1.0)
!       (2) 재생:  driver_hypre replay <golden_dir> <step>  (np=1 전용)
!       공통: 실행 디렉토리에 mg.in 필요.
!
!     공정 비교 조건:
!       - 두 솔버 모두 zero 초기추정의 잔차방정식 형태로 통일:
!         PMG 는 u0 상태에서 그대로 풀고(프로덕션 동일), hypre 는
!         A·δ = r0 (= b − A·u0), δ0 = 0 을 풀어 u = u0 + δ.
!         → hypre BiCGSTAB 의 수렴 판정 ‖r‖/‖b‖ ≤ tol 이 PMG 의
!           ‖r‖/‖r0‖ ≤ eps (6_solver_pbcg_mg.f90:263) 와 정확히 일치.
!       - tol = eps_bicg = 1e-8 (somaFlow.in 의 eps_bicg 와 동일)
!       - 시간 비교의 비대칭 주의: PMG 의 per-solve 셋업(RAP 등)은
!         SOLVE_GMG(icase=1) 내부에서 수행 → pmg_t 에 포함.
!         계층 위상 셋업(Prep_*)은 1회성 → prep_t 로 별도 보고.
!         hypre 는 hyp_tset(IJ 조립+AMG setup) / hyp_tsol(반복) 분리 보고.
!
!     출력 (러너 파싱용, rank 0):
!       RESULT <tag> pmg_its= .. pmg_t= .. pmg_res= .. hyp_its= .. \
!              hyp_tset= .. hyp_tsol= .. hyp_res= .. reldiff= ..
!       pmg_its 는 fort.501 재판독(rank0, 동일 유닛 REWIND)으로 회수 — 실패 시 -1.
!       *_res 는 독립 잔차 게이트 ‖b−A·u‖/(10·eps·‖r0‖) — 자기보고에 의존하지 않음.
!       reldiff = ‖u_hypre − u_pmg‖_2/‖u_pmg‖_2 (전역) — 두 해의 일치도.
!
!     판정: 합성 = 두 솔버 res 게이트 ≤ 1 / 재생 = hypre 자기보고 relres ≤ tol
!       (재생의 PMG res 게이트는 재귀 잔차 드리프트로 보고 전용 — LOG C009-r2)
!       PASS → exit 0, FAIL → exit 1
!=======================================================================
      PROGRAM driver_hypre
!
      USE Zcore,        ONLY: np_z => np, myrank_z => myrank
      USE Zparam,       ONLY: ndim_z => ndim, ns
      USE Zbicg,        ONLY: eps_bicg
      USE Zcoord1,      ONLY: xloc_tmp
      USE Zmpi,         ONLY: celem
      USE md_geometry,  ONLY: nelem_mg, num_neigh_mg, neigh_mg, coord,   &
                              nnode, nnodegl
      USE MD_parameter, ONLY: nf_max
      USE MD_MPI,       ONLY: nintf
      USE MD_matrix,    ONLY: nnz, ia, ja, u
!
      IMPLICIT NONE
      INCLUDE 'mpif.h'
!.....clock_gettime(CLOCK_MONOTONIC) 바인딩 — wtime() 참조 (계측 전용)
      INTERFACE
         FUNCTION clock_gettime_c(clkid, tp) BIND(C, NAME='clock_gettime')
         USE ISO_C_BINDING, ONLY: C_INT, C_LONG
         INTEGER(C_INT), VALUE :: clkid
         INTEGER(C_LONG) :: tp(2)
         INTEGER(C_INT) :: clock_gettime_c
         END FUNCTION clock_gettime_c
      END INTERFACE
!
      INTEGER :: ierr, nx, ny, nz, ncell, i, j, k, c, kk, irow, jcol
      INTEGER :: nargs, istep, iu, iver, idum(5), nfail, nbad
      INTEGER :: ipart, px, py, pz, ba, bb, bc
      INTEGER :: ioffs, its_p, its_h
      REAL(8) :: aspect, dist2, vol, diag_const, d(3)
      REAL(8) :: r0n, t0, prep_t, pmg_t, hset_t, hsol_t, hasm_t, uscale
      REAL(8) :: ser_t, fan_t, rap_t
      REAL(8) :: res_p, res_h, relres_h, reldiff, upn, dn, gtmp(2)
      LOGICAL :: replay
      CHARACTER(32)  :: arg
      CHARACTER(256) :: gdir, fn
      REAL(8), ALLOCATABLE :: diag_my(:), au_my(:), src(:), uex(:)
      REAL(8), ALLOCATABLE :: mycoord(:,:)
      REAL(8), ALLOCATABLE :: src2(:), diag2(:), au2(:)
      REAL(8), ALLOCATABLE :: gidr(:), upmg(:), uh1(:), uh2(:), rh(:)
      INTEGER, ALLOCATABLE :: gcol(:), grows(:)
!.....hypre 핸들 (셋업 재사용을 위해 프로그램 스코프에 보존)
      INTEGER(8) :: hA, hb, hx, hparA, hparb, hparx, hsol, hpre
      LOGICAL    :: hyp_alive = .FALSE.
      INTEGER, PARAMETER :: HYPRE_PARCSR = 5555
!.....hypre AMG 설정 (환경변수 주입 — hypre_config)
      INTEGER :: cf_relax, cf_sweeps, cf_agg, cf_print
      REAL(8) :: cf_strong
      CHARACTER(64) :: cfg_name
!
!.....MPI + CUPID 측 최소 상태 (driver_pmg 와 동일)
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank_z, ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, np_z, ierr)
!
      ndim_z   = 3
      ns       = 1
      eps_bicg = 1.d-8
      CALL hypre_config()
!
!.....모드 판정 + 입력 구성 (driver_pmg 의 셋업부 그대로)
      nargs = COMMAND_ARGUMENT_COUNT()
      replay = .FALSE.
      IF (nargs .GE. 1) THEN
         CALL GET_COMMAND_ARGUMENT(1, arg)
         replay = (TRIM(arg) .EQ. 'replay')
      END IF
!
      IF (replay) THEN
         IF (nargs .LT. 3) STOP 'usage: driver_hypre replay <golden_dir> <step>'
         IF (np_z .NE. 1) STOP 'replay: np=1 전용 (골든이 np1 채취)'
         CALL GET_COMMAND_ARGUMENT(2, gdir)
         CALL GET_COMMAND_ARGUMENT(3, arg); READ (arg, *) istep
!
         WRITE (fn, '(A,A)') TRIM(gdir), '/setup_r0.bin'
         OPEN (NEWUNIT=iu, FILE=fn, FORM='unformatted', ACCESS='stream',  &
               STATUS='old', ACTION='read')
         READ (iu) iver, idum(1), nelem_mg, nf_max
         IF (idum(1) .NE. 3) STOP 'replay: ndim/=3 in setup dump'
         ncell = nelem_mg
         ALLOCATE (num_neigh_mg(ncell), neigh_mg(nf_max, ncell))
         ALLOCATE (xloc_tmp(ncell, 3), celem(ncell))
         READ (iu) num_neigh_mg(1:ncell)
         READ (iu) neigh_mg(1:nf_max, 1:ncell)
         READ (iu) xloc_tmp(1:ncell, 1:3)
         READ (iu) celem(1:ncell)
         CLOSE (iu)
         WRITE (*, '(A,I9,A,I3,A,I6)') ' hypre-AB replay: ncell=', ncell, &
                '  nf_max=', nf_max, '  step=', istep
      ELSE
!.....합성: 7점 스텐실 Poisson (driver_pmg C007 과 동일 구성)
         nx = 24; ny = 24; nz = 24; aspect = 1.d0
         IF (nargs .GE. 3) THEN
            CALL GET_COMMAND_ARGUMENT(1, arg); READ (arg, *) nx
            CALL GET_COMMAND_ARGUMENT(2, arg); READ (arg, *) ny
            CALL GET_COMMAND_ARGUMENT(3, arg); READ (arg, *) nz
         END IF
         IF (nargs .GE. 4) THEN
            CALL GET_COMMAND_ARGUMENT(4, arg); READ (arg, *) aspect
         END IF
         ipart = 0
         IF (nargs .GE. 5) THEN
            CALL GET_COMMAND_ARGUMENT(5, arg); READ (arg, *) ipart
         END IF
         px = 1; py = 1; pz = np_z
         IF (ipart .EQ. 1) THEN
            kk = HUGE(kk)
            DO ba = 1, np_z
               IF (MOD(np_z, ba) .NE. 0) CYCLE
               DO bb = ba, np_z/ba
                  IF (MOD(np_z/ba, bb) .NE. 0) CYCLE
                  bc = np_z/(ba*bb)
                  IF (bc .LT. bb) CYCLE
                  IF (bc - ba .LT. kk) THEN
                     kk = bc - ba; px = ba; py = bb; pz = bc
                  END IF
               END DO
            END DO
         END IF
         ncell = nx*ny*nz
         IF (myrank_z .EQ. 0)                                            &
            WRITE (*, '(A,3I5,A,F8.2,A,I9,A,I4)') ' hypre-AB grid:', nx, &
                ny, nz, '  aspect=', aspect, '  ncell=', ncell, '  np=', np_z
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
            IF (ipart .EQ. 1) THEN
               ba = ((i-1)*px)/nx; bb = ((j-1)*py)/ny; bc = ((k-1)*pz)/nz
               celem(c) = 1 + ba + bb*px + bc*px*py
            ELSE
               celem(c) = 1 + ((k-1)*np_z)/nz
            END IF
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
      END IF
!
!.....프로덕션 초기화 체인 (driver_pmg 와 동일) — Prep 구간은 1회성 셋업으로 계시
      CALL read_input_mg
!.....세 구간을 나눠 계측 (np 확장성 분석용, 2026-08-21):
!       ser_t  = subdomain_infor_MG — 전역 격자 조대화·P/R 구축. rank 0 단독
!                실행이라 np 가 늘어도 줄지 않는다 (read_grid.f90:249 동일 구조).
!       fan_t  = MG_tmp 파일 fan-out 판독 (rank 별 자기 몫)
!       prep_t = Prep_fine_P + Prep_MG_GarL (분산)
      t0 = wtime()
      IF (myrank_z .EQ. 0) CALL subdomain_infor_MG
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      ser_t = wtime() - t0
      t0 = wtime()
      CALL read_mesh_MPI
      fan_t = wtime() - t0
      IF (.NOT. replay) THEN
         ALLOCATE (mycoord(3, nnodegl))
         mycoord(:, 1:nnodegl) = coord(:, 1:nnodegl)
      END IF
      t0 = wtime()
      CALL Prep_fine_P
      CALL Prep_MG_GarL
      prep_t = wtime() - t0
!
      CALL MPI_ALLREDUCE(nintf, i, 1, MPI_INTEGER, MPI_SUM,              &
                         MPI_COMM_WORLD, ierr)
      IF (myrank_z .EQ. 0)                                               &
         WRITE (*, '(A,I9,A,I10,A,F9.3)') ' setup done: ncell=', i,      &
                '  nnz_loc=', nnz, '  prep_t=', prep_t
      IF (i .NE. ncell) STOP 'driver_hypre: sum(nintf) /= ncell (mapping?)'
!
!.....전역 행 번호: rank 별 내부행 연속 소유 [ioffs+1, ioffs+nintf]
!     고스트 열의 전역 번호는 gid 실수 벡터를 communicate_s 로 교환해 획득
!     (분할 형상 무관 — hypre IJ 의 연속 소유 요구를 구성적으로 충족)
      ioffs = 0
      CALL MPI_EXSCAN(nintf, ioffs, 1, MPI_INTEGER, MPI_SUM,             &
                      MPI_COMM_WORLD, ierr)
      IF (myrank_z .EQ. 0) ioffs = 0
      ALLOCATE (gidr(nnode))
      gidr = 0.d0
      DO i = 1, nintf
         gidr(i) = DBLE(ioffs + i)
      END DO
      IF (np_z .GT. 1) CALL communicate_s(gidr)
      ALLOCATE (gcol(nnz), grows(nintf))
      nbad = 0
      DO irow = 1, nintf
         grows(irow) = ioffs + irow
         DO kk = ia(irow), ia(irow+1) - 1
            IF (gidr(ja(kk)) .LT. 0.5d0) nbad = nbad + 1
            gcol(kk) = NINT(gidr(ja(kk)))
         END DO
      END DO
      IF (nbad .GT. 0) STOP 'driver_hypre: ghost 전역 id 미해결 (communicate?)'
!
      nfail = 0
      rap_t = 0.d0
      ALLOCATE (upmg(nintf), uh1(nintf), uh2(nintf), rh(nintf))
!
      IF (replay) THEN
!=====================================================================
!     재생: k1 (u0=0) → k2 (u0 = 각자의 k1 해; 행렬 k1 것 재사용, 프로덕션 동일)
!=====================================================================
         ALLOCATE (src(ncell), diag_my(ncell), au_my(nnz))
         CALL load_pre(gdir, istep, 1, ncell, nnz, src, au_my, diag_my)
!
!........(A) PMG k1
         u = 0.d0
         CALL assemble_FVM(1, nnz, src, au_my, diag_my)
         r0n = norm_res(nintf, src, au_my, uh1, .TRUE.)   ! x=0 → ‖b‖
         t0 = wtime()
         CALL SOLVE_GMG(1)
         pmg_t = wtime() - t0
         CALL negchk(pmg_t, t0)
         its_p = its_from_501()
         upmg(1:nintf) = u(1:nintf)
         res_p = norm_res(nintf, src, au_my, upmg, .FALSE.)
         uscale = set_uscale()          ! fid 정규화 = k1 장 스케일 (driver_pmg 동일)
!........RAP 비용 분리: 같은 계를 icase=2 (RAP 생략)로 다시 풀어 차이를 취한다.
!        u 를 0 으로 되돌려 반복수가 1회차와 같아지게 한 뒤, k2 단계가
!        프로덕션과 같은 초기추정(=k1 해)을 쓰도록 복원한다 (재생은 np=1 전용).
         u = 0.d0
         t0 = wtime()
         CALL SOLVE_GMG(2)
         rap_t = pmg_t - (wtime() - t0)
         u = 0.d0
         u(1:nintf) = upmg(1:nintf)
!........(B) hypre k1 — 셋업(IJ 조립 + AMG) + 솔브
         CALL hypre_setup(au_my, hset_t)
         CALL hypre_apply(src, uh1, its_h, relres_h, hsol_t)
         res_h = norm_res(nintf, src, au_my, uh1, .FALSE.)
         CALL report('k1', r0n)
         IF (relres_h .GT. eps_bicg) nfail = nfail + 1
!
!........(A) PMG k2 — u 는 PMG k1 해 유지 (icase=2: b 만 교체)
         ALLOCATE (src2(ncell), diag2(ncell), au2(nnz))
         CALL load_pre(gdir, istep, 2, ncell, nnz, src2, au2, diag2)
         CALL assemble_FVM(2, nnz, src2, au2, diag2)
         r0n = norm_res(nintf, src2, au_my, upmg, .FALSE.)  ! r0 = b2 − A·u_k1
         t0 = wtime()
         CALL SOLVE_GMG(2)
         pmg_t = wtime() - t0
         CALL negchk(pmg_t, t0)
         its_p = its_from_501()
         upmg(1:nintf) = u(1:nintf)
         res_p = norm_res(nintf, src2, au_my, upmg, .FALSE.)
!........(B) hypre k2 — 잔차방정식: A·δ = b2 − A·u_h,k1, δ0=0, u = u_h,k1 + δ
!        행렬 불변이므로 AMG 계층 재사용 (PMG 의 icase=2 가 RAP 을 건너뛰는 것과 대응)
         CALL resid_vec(nintf, src2, au_my, uh1, rh)
         hset_t = 0.d0
         hasm_t = 0.d0
         CALL hypre_apply(rh, uh2, its_h, relres_h, hsol_t)
         uh1(1:nintf) = uh1(1:nintf) + uh2(1:nintf)
         res_h = norm_res(nintf, src2, au_my, uh1, .FALSE.)
         CALL report('k2', r0n)
         IF (relres_h .GT. eps_bicg) nfail = nfail + 1
!
      ELSE
!=====================================================================
!     합성: 제작해 기반 — 두 솔버 모두 u0=0, 독립 잔차 게이트로 판정
!=====================================================================
         vol = aspect
         diag_const = 2.d0*aspect + 2.d0*aspect + 2.d0/aspect
         ALLOCATE (diag_my(nintf), au_my(nnz), src(nintf), uex(nnodegl))
         diag_my = diag_const
         DO irow = 1, nnode
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
         DO c = 1, nnodegl
            i = NINT(mycoord(1, c) + 0.5d0)
            j = NINT(mycoord(2, c) + 0.5d0)
            k = NINT(mycoord(3, c)/aspect + 0.5d0)
            irow = i + (j-1)*nx + (k-1)*nx*ny
            uex(c) = 1.d0 + SIN(0.5d0*DBLE(irow))
         END DO
         DO irow = 1, nintf
            src(irow) = 0.d0
            DO kk = ia(irow), ia(irow+1) - 1
               src(irow) = src(irow) + au_my(kk)*uex(ja(kk))
            END DO
         END DO
!
!........(A) PMG
         u = 0.d0
         CALL assemble_FVM(1, nnz, src, au_my, diag_my)
         r0n = norm_res(nintf, src, au_my, uh1, .TRUE.)
         t0 = wtime()
         CALL SOLVE_GMG(1)
         pmg_t = wtime() - t0
         CALL negchk(pmg_t, t0)
         its_p = its_from_501()
         upmg(1:nintf) = u(1:nintf)
         res_p = norm_res(nintf, src, au_my, upmg, .FALSE.)
         uscale = set_uscale()
!........동일 계를 icase=2 (RAP 생략)로 다시 풀어 RAP 비용을 분리.
!        u 를 0 으로 되돌려 반복수가 1회차와 같아지게 한다 → 차이 = RAP.
         u = 0.d0
         t0 = wtime()
         CALL SOLVE_GMG(2)
         rap_t = pmg_t - (wtime() - t0)
!........(B) hypre
         CALL hypre_setup(au_my, hset_t)
         CALL hypre_apply(src, uh1, its_h, relres_h, hsol_t)
         res_h = norm_res(nintf, src, au_my, uh1, .FALSE.)
         CALL report('syn', r0n)
         IF (res_p/(10.d0*eps_bicg*r0n + 1.d-300) .GT. 1.d0) nfail = nfail + 1
         IF (res_h/(10.d0*eps_bicg*r0n + 1.d-300) .GT. 1.d0) nfail = nfail + 1
!
!........제작해 대비 오차 (참고 보고)
         upn = 0.d0; dn = 0.d0
         DO i = 1, nintf
            upn = upn + (upmg(i) - uex(i))**2
            dn  = dn  + (uh1(i)  - uex(i))**2
         END DO
         gtmp(1) = upn; gtmp(2) = dn
         CALL MPI_ALLREDUCE(MPI_IN_PLACE, gtmp, 2, MPI_DOUBLE_PRECISION, &
                            MPI_SUM, MPI_COMM_WORLD, ierr)
         IF (myrank_z .EQ. 0)                                            &
            WRITE (*, '(A,ES11.3,A,ES11.3)') ' err-vs-exact  pmg=',      &
                SQRT(gtmp(1)), '  hypre=', SQRT(gtmp(2))
      END IF
!
      CALL hypre_free()
      IF (nfail .EQ. 0) THEN
         IF (myrank_z .EQ. 0) WRITE (*, '(A)') ' VERDICT PASS'
         CALL MPI_FINALIZE(ierr)
      ELSE
         IF (myrank_z .EQ. 0) WRITE (*, '(A)') ' VERDICT FAIL'
         CALL MPI_FINALIZE(ierr)
         CALL EXIT(1)
      END IF
!
      CONTAINS
!
!-----------------------------------------------------------------------
      SUBROUTINE load_pre(gdir, istep, isite, n, nnz_expect, src, auv, dg)
      CHARACTER(*), INTENT(IN) :: gdir
      INTEGER, INTENT(IN)  :: istep, isite, n, nnz_expect
      REAL(8), INTENT(OUT) :: src(n), auv(nnz_expect), dg(n)
      INTEGER :: iu, hdr(6)
      CHARACTER(256) :: fn
      WRITE (fn, '(A,A,I0,A,I0,A)') TRIM(gdir), '/s', istep, '_k',       &
             isite, '_r0_c1.pre'
      OPEN (NEWUNIT=iu, FILE=fn, FORM='unformatted', ACCESS='stream',    &
            STATUS='old', ACTION='read')
      READ (iu) hdr
      IF (hdr(5) .NE. n)          STOP 'load_pre: n mismatch'
      IF (hdr(6) .NE. nnz_expect) STOP 'load_pre: nnz mismatch'
      READ (iu) src
      READ (iu) auv
      READ (iu) dg
      CLOSE (iu)
      END SUBROUTINE load_pre
!
!-----------------------------------------------------------------------
      SUBROUTINE resid_vec(n, rhs, auv, x, r)
!     r = rhs − A·x (full-row; x 는 내부행만 — 고스트는 교환해 사용)
      INTEGER, INTENT(IN)  :: n
      REAL(8), INTENT(IN)  :: rhs(n), auv(nnz), x(n)
      REAL(8), INTENT(OUT) :: r(n)
      INTEGER :: irow, kk
      REAL(8), ALLOCATABLE :: xt(:)
      ALLOCATE (xt(nnode))
      xt = 0.d0
      xt(1:n) = x(1:n)
      IF (np_z .GT. 1) CALL communicate_s(xt)
      DO irow = 1, n
         r(irow) = rhs(irow)
         DO kk = ia(irow), ia(irow+1) - 1
            r(irow) = r(irow) - auv(kk)*xt(ja(kk))
         END DO
      END DO
      DEALLOCATE (xt)
      END SUBROUTINE resid_vec
!
!-----------------------------------------------------------------------
      REAL(8) FUNCTION norm_res(n, rhs, auv, x, xzero)
!     ‖rhs − A·x‖_2 (전역). xzero=.TRUE. 면 x 를 0 으로 간주 (‖rhs‖).
      INTEGER, INTENT(IN) :: n
      REAL(8), INTENT(IN) :: rhs(n), auv(nnz), x(n)
      LOGICAL, INTENT(IN) :: xzero
      INTEGER :: i, ierr
      REAL(8) :: rn, rng
      REAL(8), ALLOCATABLE :: r(:)
      ALLOCATE (r(n))
      IF (xzero) THEN
         r(1:n) = rhs(1:n)
      ELSE
         CALL resid_vec(n, rhs, auv, x, r)
      END IF
      rn = 0.d0
      DO i = 1, n
         rn = rn + r(i)*r(i)
      END DO
      CALL MPI_ALLREDUCE(rn, rng, 1, MPI_DOUBLE_PRECISION, MPI_SUM,      &
                         MPI_COMM_WORLD, ierr)
      norm_res = SQRT(rng)
      DEALLOCATE (r)
      END FUNCTION norm_res
!
!-----------------------------------------------------------------------
      INTEGER FUNCTION its_from_501()
!     solve_pbcg_mg 가 rank0 의 unit 501 에 누적한 'its ro/ro0' 의 마지막
!     줄을 같은 유닛 REWIND 재판독으로 회수. EOF 상태를 남기지 않도록
!     2-패스(줄 수 세기 → 정확히 그 수만큼 재판독)로 파일 위치를 마지막
!     레코드 뒤에 두어 이후 솔브의 append 를 보존한다. 실패 시 -1.
      INTEGER :: itv, ios, nrec, i
      REAL(8) :: dum
      its_from_501 = -1
      IF (myrank_z .NE. 0) RETURN
      REWIND (501, IOSTAT=ios)
      IF (ios .NE. 0) RETURN
      nrec = 0
      DO
         READ (501, *, IOSTAT=ios) itv, dum
         IF (ios .NE. 0) EXIT
         nrec = nrec + 1
      END DO
      REWIND (501, IOSTAT=ios)
      DO i = 1, nrec
         READ (501, *, IOSTAT=ios) itv, dum
         IF (ios .NE. 0) RETURN
      END DO
      IF (nrec .GT. 0) its_from_501 = itv
      END FUNCTION its_from_501
!
!-----------------------------------------------------------------------
      REAL(8) FUNCTION wtime()
!     진짜 단조 시계 — clock_gettime(CLOCK_MONOTONIC) 직접 바인딩.
!     주의 (2026-08-21 확정): MPI_WTIME 도, ifort 의 SYSTEM_CLOCK 도
!     이 환경에선 벽시계(CLOCK_REALTIME)다 — SYSTEM_CLOCK 이 돌려준
!     count 가 epoch 기준 1.787e15 µs 로 확인됐다. WSL2 는 호스트 시각
!     동기화 때 REALTIME 이 역행하며, 실제로 pmg_t=-14.4 s 같은 음수
!     경과시간이 관측됐다. 계측은 반드시 MONOTONIC 이어야 한다.
      INTEGER(8) :: tv(2)
      INTEGER :: rc
      rc = clock_gettime_c(1, tv)          ! 1 = CLOCK_MONOTONIC (Linux)
      wtime = DBLE(tv(1)) + 1.d-9*DBLE(tv(2))
      END FUNCTION wtime
!
!-----------------------------------------------------------------------
      SUBROUTINE negchk(dt, tstart)
!     음수 경과시간 진단 (2026-08-21): 시계는 단조임을 확인했으므로
!     (clock_gettime 20초 조밀 샘플링 역행 0회, SYSTEM_CLOCK µs 해상도)
!     음수가 나오면 SOLVE_GMG 실행 중 호출자 프레임의 tstart 가
!     덮였다는 뜻 — 메모리 손상 신호다. 원시값을 남겨 판별한다.
      REAL(8), INTENT(IN) :: dt, tstart
      IF (dt .GE. 0.d0) RETURN
      WRITE (*, '(A,ES24.16,A,ES24.16,A,ES14.6)')                        &
         ' !! NEGATIVE elapsed: t0=', tstart, ' t1=', tstart+dt,         &
         ' dt=', dt
      END SUBROUTINE negchk
!
!-----------------------------------------------------------------------
      SUBROUTINE hypre_config()
!     AMG 설정을 환경변수로 주입 (재빌드 없이 스윕하기 위함).
!     기본값은 "hypre 기본 설정" — SetRelaxType 을 호출하지 않으면
!     hypre 2.33 기본 스무더 ℓ1-GS forward/backward(13/14, 최조대는 가우스
!     소거 9) 가 유지된다. 이게 분할-강건(np-robust) 설정이라 np 스케일
!     비교의 정당한 상대다. HYPRE_RELAX=6 을 주면 hybrid sym-GS/SSOR
!     (비-ℓ1, 랭크 수에 따라 품질 저하) 로 바꿔 그 영향을 측정할 수 있다.
!     조대화 10(HMIS)·보간 6(ext+i)·P_max 4·강연결 0.25 는 hypre 기본값과
!     동일하므로 명시하지 않는다 (기본값 추종).
      CHARACTER(64) :: ev
      INTEGER :: ln, ios
      cf_relax  = -1        ! <0 = hypre 기본(13/14) 유지
      cf_sweeps =  1
      cf_agg    =  0
      cf_strong = -1.d0     ! <0 = hypre 기본(0.25) 유지
      cfg_name  = 'default'
      CALL GET_ENVIRONMENT_VARIABLE('HYPRE_RELAX', ev, ln)
      IF (ln .GT. 0) READ (ev, *, IOSTAT=ios) cf_relax
      CALL GET_ENVIRONMENT_VARIABLE('HYPRE_SWEEPS', ev, ln)
      IF (ln .GT. 0) READ (ev, *, IOSTAT=ios) cf_sweeps
      CALL GET_ENVIRONMENT_VARIABLE('HYPRE_AGG', ev, ln)
      IF (ln .GT. 0) READ (ev, *, IOSTAT=ios) cf_agg
      CALL GET_ENVIRONMENT_VARIABLE('HYPRE_STRONG', ev, ln)
      IF (ln .GT. 0) READ (ev, *, IOSTAT=ios) cf_strong
      CALL GET_ENVIRONMENT_VARIABLE('HYPRE_CFG', ev, ln)
      IF (ln .GT. 0) cfg_name = ev
!.....HYPRE_PRINT=3 이면 BoomerAMG 셋업이 계층 표(레벨·격자/연산자 복잡도)를
!     stdout 에 찍는다 — 셋업 비용 차이의 원인 분석용
      cf_print = 0
      CALL GET_ENVIRONMENT_VARIABLE('HYPRE_PRINT', ev, ln)
      IF (ln .GT. 0) READ (ev, *, IOSTAT=ios) cf_print
      END SUBROUTINE hypre_config
!
!-----------------------------------------------------------------------
      SUBROUTINE hypre_setup(auv, t_setup)
!     전역계 (행 grows, 열 gcol, 값 auv 의 내부행 구간) 를 hypre IJ 로
!     조립하고 BiCGSTAB + BoomerAMG(V(1,1)) 의 셋업까지 수행.
!     여기까지가 PMG 의 SOLVE_GMG(icase=1) 내부 RAP 셋업에 대응한다 —
!     행렬이 안 바뀌는 후속 솔브(k2)는 hypre_apply 만 다시 부르면 된다.
      REAL(8), INTENT(IN)  :: auv(nnz)
      REAL(8), INTENT(OUT) :: t_setup
      INTEGER :: irow, n1, one, ierr, rows1(1)
      REAL(8) :: tw
      REAL(8), ALLOCATABLE :: z0(:)
!
      tw = wtime()
      one = 1
!.....행렬 조립 (행 단위)
      CALL HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ioffs+1, ioffs+nintf,    &
                                ioffs+1, ioffs+nintf, hA, ierr)
      CALL HYPRE_IJMatrixSetObjectType(hA, HYPRE_PARCSR, ierr)
      CALL HYPRE_IJMatrixInitialize(hA, ierr)
      DO irow = 1, nintf
         n1 = ia(irow+1) - ia(irow)
         rows1(1) = ioffs + irow
         CALL HYPRE_IJMatrixSetValues(hA, one, n1, rows1,                &
                 gcol(ia(irow):ia(irow+1)-1), auv(ia(irow):ia(irow+1)-1), ierr)
      END DO
      CALL HYPRE_IJMatrixAssemble(hA, ierr)
      CALL HYPRE_IJMatrixGetObject(hA, hparA, ierr)
!.....하네스 고유 비용: CUPID CSR → hypre IJ 변환·조립. 네이티브 hypre 앱이면
!     처음부터 IJ 로 조립하므로 이 몫은 비교에서 빼고 봐야 공정하다.
      hasm_t = wtime() - tw
!.....벡터 b, x (값은 hypre_apply 에서 매 솔브 갱신)
      ALLOCATE (z0(nintf)); z0 = 0.d0
      CALL HYPRE_IJVectorCreate(MPI_COMM_WORLD, ioffs+1, ioffs+nintf, hb, ierr)
      CALL HYPRE_IJVectorSetObjectType(hb, HYPRE_PARCSR, ierr)
      CALL HYPRE_IJVectorInitialize(hb, ierr)
      CALL HYPRE_IJVectorSetValues(hb, nintf, grows, z0, ierr)
      CALL HYPRE_IJVectorAssemble(hb, ierr)
      CALL HYPRE_IJVectorGetObject(hb, hparb, ierr)
!
      CALL HYPRE_IJVectorCreate(MPI_COMM_WORLD, ioffs+1, ioffs+nintf, hx, ierr)
      CALL HYPRE_IJVectorSetObjectType(hx, HYPRE_PARCSR, ierr)
      CALL HYPRE_IJVectorInitialize(hx, ierr)
      CALL HYPRE_IJVectorSetValues(hx, nintf, grows, z0, ierr)
      CALL HYPRE_IJVectorAssemble(hx, ierr)
      CALL HYPRE_IJVectorGetObject(hx, hparx, ierr)
      DEALLOCATE (z0)
!
!.....BoomerAMG 예조건자 — Tol(0)+MaxIter(1) 로 "적용당 V-cycle 1회" 고정
!     (기본값 max_iter=20/tol=1e-6 이면 예조건자가 아니라 완전 솔브가 됨)
      CALL HYPRE_BoomerAMGCreate(hpre, ierr)
      CALL HYPRE_BoomerAMGSetTol(hpre, 0.d0, ierr)
      CALL HYPRE_BoomerAMGSetMaxIter(hpre, 1, ierr)
      CALL HYPRE_BoomerAMGSetNumSweeps(hpre, cf_sweeps, ierr)
      CALL HYPRE_BoomerAMGSetPrintLevel(hpre, cf_print, ierr)
      IF (cf_relax  .GE. 0)    CALL HYPRE_BoomerAMGSetRelaxType(hpre, cf_relax, ierr)
      IF (cf_strong .GT. 0.d0) CALL HYPRE_BoomerAMGSetStrongThrshld(hpre, cf_strong, ierr)
      IF (cf_agg    .GT. 0)    CALL HYPRE_BoomerAMGSetAggNumLevels(hpre, cf_agg, ierr)
!
!.....BiCGSTAB (수렴 판정: ‖r‖/‖b‖ ≤ tol — x0=0 이므로 PMG 기준과 동일)
      CALL HYPRE_ParCSRBiCGSTABCreate(MPI_COMM_WORLD, hsol, ierr)
      CALL HYPRE_ParCSRBiCGSTABSetTol(hsol, eps_bicg, ierr)
      CALL HYPRE_ParCSRBiCGSTABSetMaxIter(hsol, 1000, ierr)
      CALL HYPRE_ParCSRBiCGSTABSetPrintLev(hsol, 0, ierr)
      CALL HYPRE_ParCSRBiCGSTABSetPrecond(hsol, 2, hpre, ierr)
      CALL HYPRE_ParCSRBiCGSTABSetup(hsol, hparA, hparb, hparx, ierr)
      hyp_alive = .TRUE.
      t_setup = wtime() - tw
      END SUBROUTINE hypre_setup
!
!-----------------------------------------------------------------------
      SUBROUTINE hypre_apply(rhs, x, its, relres, t_solve)
!     셋업된 계층을 그대로 쓰고 b 만 교체해 x0=0 에서 푼다.
      REAL(8), INTENT(IN)  :: rhs(nintf)
      REAL(8), INTENT(OUT) :: x(nintf), relres, t_solve
      INTEGER, INTENT(OUT) :: its
      INTEGER :: ierr
      REAL(8) :: tw
      x(1:nintf) = 0.d0
      CALL HYPRE_IJVectorSetValues(hb, nintf, grows, rhs, ierr)
      CALL HYPRE_IJVectorSetValues(hx, nintf, grows, x, ierr)
      tw = wtime()
      CALL HYPRE_ParCSRBiCGSTABSolve(hsol, hparA, hparb, hparx, ierr)
      t_solve = wtime() - tw
      CALL HYPRE_ParCSRBiCGSTABGetNumIter(hsol, its, ierr)
      CALL HYPRE_ParCSRBiCGSTABGetFinalRel(hsol, relres, ierr)
      CALL HYPRE_IJVectorGetValues(hx, nintf, grows, x, ierr)
      END SUBROUTINE hypre_apply
!
!-----------------------------------------------------------------------
      SUBROUTINE hypre_free()
      INTEGER :: ierr
      IF (.NOT. hyp_alive) RETURN
      CALL HYPRE_ParCSRBiCGSTABDestroy(hsol, ierr)
      CALL HYPRE_BoomerAMGDestroy(hpre, ierr)
      CALL HYPRE_IJMatrixDestroy(hA, ierr)
      CALL HYPRE_IJVectorDestroy(hb, ierr)
      CALL HYPRE_IJVectorDestroy(hx, ierr)
      hyp_alive = .FALSE.
      END SUBROUTINE hypre_free
!
!-----------------------------------------------------------------------
      SUBROUTINE report(tag, r0n)
!     RESULT 라인 (러너 파싱 대상) — uh1 vs upmg 일치도 포함
!     reldiff 는 상수 영공간 사영 후 값: 재생(압력) 행렬은 전-Neumann 계라
!     해가 상수만큼 자유 — 두 솔버가 주입하는 상수 성분이 달라 원시 차는
!     O(1) 이 될 수 있다 (iSMR s10 k2 에서 관측). 물리적으로 유의미한
!     비교는 mean(Δ) 제거 후. 원시 차는 reldiff_raw 로 병기.
      CHARACTER(*), INTENT(IN) :: tag
      REAL(8), INTENT(IN) :: r0n
      INTEGER :: i, ierr, ng
      REAL(8) :: en, un, dsum, dbar, g(3), en0, gate_p, gate_h, rdraw, fidm
      en = 0.d0; un = 0.d0; dsum = 0.d0
      DO i = 1, nintf
         en   = en   + (uh1(i) - upmg(i))**2
         un   = un   + upmg(i)*upmg(i)
         dsum = dsum + (uh1(i) - upmg(i))
      END DO
      g(1) = en; g(2) = un; g(3) = dsum
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, g, 3, MPI_DOUBLE_PRECISION,       &
                         MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE(nintf, ng, 1, MPI_INTEGER, MPI_SUM,             &
                         MPI_COMM_WORLD, ierr)
      dbar  = g(3)/DBLE(ng)
      en0   = MAX(g(1) - DBLE(ng)*dbar*dbar, 0.d0)   ! Σ(Δ−mean)² = ΣΔ²−N·mean²
      rdraw   = SQRT(g(1))/(SQRT(g(2)) + 1.d-300)
      reldiff = SQRT(en0)/(SQRT(g(2)) + 1.d-300)
!.....fid: mean(Δ) 제거 후 max|Δ|/uscale — k2 처럼 자기 노름이 ~0 인 보정
!     솔브에서 reldiff 가 부풀려지는 문제의 보정 (driver_pmg 의 fidelity 와
!     동일 철학: 유의미 스케일은 k1 압력장, LOG C009-r2)
      fidm = 0.d0
      DO i = 1, nintf
         fidm = MAX(fidm, ABS(uh1(i) - upmg(i) - dbar))
      END DO
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, fidm, 1, MPI_DOUBLE_PRECISION,    &
                         MPI_MAX, MPI_COMM_WORLD, ierr)
      fidm = fidm/uscale
      gate_p = res_p/(10.d0*eps_bicg*r0n + 1.d-300)
      gate_h = res_h/(10.d0*eps_bicg*r0n + 1.d-300)
      IF (myrank_z .EQ. 0) THEN
         WRITE (*, '(3A,I5,A,F9.3,A,ES10.2,A,I5,A,F9.3,A,F9.3,A,ES10.2,A,ES10.2,A,ES10.2,A,F9.3,A,F9.3,A,F9.3,A,F9.3,A,F9.3,A)') &
             ' RESULT ', tag, ' pmg_its=', its_p, ' pmg_t=', pmg_t,      &
             ' pmg_res=', gate_p, ' hyp_its=', its_h, ' hyp_tset=',      &
             hset_t, ' hyp_tsol=', hsol_t, ' hyp_res=', gate_h,          &
             ' fid=', fidm, ' reldiff=', reldiff, ' hyp_tasm=', hasm_t,  &
             ' pmg_tprep=', prep_t, ' pmg_tser=', ser_t,          &
             ' pmg_tfan=', fan_t, ' pmg_trap=', rap_t,               &
             ' cfg='//TRIM(cfg_name)
         IF (fidm .GT. 1.d-4)                                            &
            WRITE (*, '(3A)') ' WARN[', tag, '] fid > 1e-4 (두 해 괴리 큼)'
      END IF
      END SUBROUTINE report
!
!-----------------------------------------------------------------------
      REAL(8) FUNCTION set_uscale()
!     fid 정규화 스케일 = 1 + 전역 max|u_pmg,k1| (driver_pmg 의 uscale 동일)
      INTEGER :: i, ierr
      REAL(8) :: um
      um = 0.d0
      DO i = 1, nintf
         um = MAX(um, ABS(upmg(i)))
      END DO
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, um, 1, MPI_DOUBLE_PRECISION,      &
                         MPI_MAX, MPI_COMM_WORLD, ierr)
      set_uscale = 1.d0 + um
      END FUNCTION set_uscale
!
      END PROGRAM driver_hypre
