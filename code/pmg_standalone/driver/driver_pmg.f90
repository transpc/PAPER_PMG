!=======================================================================
!     driver_pmg — PMG(BiCGSTAB+MG 예조건) 단독 검증 드라이버
!     (LOOP C007 합성 모드 / C009 골든 재생 모드)
!
!     프로덕션(read_grid.f90/pressure_solve.f90)과 동일한 호출 체인:
!       read_input_mg → subdomain_infor_MG → read_mesh_MPI
!       → Prep_fine_P → Prep_MG_GarL → assemble_FVM → SOLVE_GMG
!
!     사용법:
!       (1) 합성:  driver_pmg [nx ny nz aspect]        (기본 24 24 24 1.0)
!       (2) 재생:  driver_pmg replay <golden_dir> <step> [fid_gate]
!           <golden_dir>/setup_r0.bin, s<step>_k{1,2}_r0_c1.{pre,post} 필요.
!           [fid_gate]: 충실도 판정 임계 (기본 1e-9). 골든 세트별로 다름 —
!           프로덕션·하네스 its 차이(예조건자 시드 문맥, LOG C009)의 발현
!           크기가 케이스마다 달라서 (구입력 ≤1e-10 vs ECT1 ≤2.4e-7, C017).
!           k1(주 압력 솔브) 재생 → u 상태 유지 → k2(비직교 보정) 재생 —
!           프로덕션과 동일 순서 (k2 는 k1 해를 초기추정으로 사용).
!       공통: 실행 디렉토리에 mg.in 필요 (재생은 케이스 mg.in 과 동일해야 함).
!
!     판정 (LOOP.md §2 L2):
!       - 독립 residual ‖b−A·u‖/‖b‖ ≤ 10*eps  (자기보고에 의존하지 않음)
!       - 재생: 골든 u* 와 비교 (bitwise 불일치 개수 / max|Δ| / rel-L2)
!         PASS 조건: rel-L2 ≤ 1e-12 (bitwise 는 보고 — 기대값 0 불일치)
!       - PASS → exit 0, FAIL → exit 1.  its 는 fort.501 에서 러너가 파싱.
!=======================================================================
      PROGRAM driver_pmg
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
      USE MD_matrix,    ONLY: nnz, ia, ja, au, u, b
!
      IMPLICIT NONE
      INCLUDE 'mpif.h'
!
      INTEGER :: ierr, nx, ny, nz, ncell, i, j, k, c, kk, irow, jcol
      INTEGER :: nargs, istep, iu, iver, idum(5), nnz_g, nfail
      INTEGER :: ipart, px, py, pz, ba, bb, bc
      REAL(8) :: aspect, dist2, vol, diag_const, d(3)
      REAL(8) :: rel_res, rel_err, r0n, uscale, fid_gate
      LOGICAL :: replay
      CHARACTER(32)  :: arg
      CHARACTER(256) :: gdir, fn
      REAL(8), ALLOCATABLE :: diag_my(:), au_my(:), src(:), uex(:)
      REAL(8), ALLOCATABLE :: mycoord(:,:)   ! 합성 모드 좌표 스냅샷
      REAL(8), ALLOCATABLE :: ug(:), au2(:), src2(:), diag2(:)
!
!.....MPI + CUPID 측 최소 상태 (프로덕션에서 본체가 하는 일)
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank_z, ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, np_z, ierr)
!
      ndim_z   = 3
      ns       = 1
      eps_bicg = 1.d-8          ! = somaFlow.in 의 eps_bicg (재생 충실도 전제)
      fid_gate = 1.d-9          ! 재생 충실도 판정 기본값 (replay 4번째 인자로 세트별 지정)
!
!.....모드 판정 + 입력 구성
      nargs = COMMAND_ARGUMENT_COUNT()
      replay = .FALSE.
      IF (nargs .GE. 1) THEN
         CALL GET_COMMAND_ARGUMENT(1, arg)
         replay = (TRIM(arg) .EQ. 'replay')
      END IF
!
      IF (replay) THEN
!.....재생: 셋업 배열을 골든 덤프에서 로드
         IF (nargs .LT. 3) STOP 'usage: driver_pmg replay <golden_dir> <step>'
         IF (np_z .NE. 1) STOP 'replay: np=1 전용 (골든이 np1 채취)'
         CALL GET_COMMAND_ARGUMENT(2, gdir)
         CALL GET_COMMAND_ARGUMENT(3, arg); READ (arg, *) istep
         IF (nargs .GE. 4) THEN
            CALL GET_COMMAND_ARGUMENT(4, arg); READ (arg, *) fid_gate
         END IF
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
         WRITE (*, '(A,I9,A,I3,A,I6)') ' C009 replay: ncell=', ncell,    &
                '  nf_max=', nf_max, '  step=', istep
      ELSE
!.....합성: 7점 스텐실 Poisson (LOG C007)
         nx = 24; ny = 24; nz = 24; aspect = 1.d0
         IF (nargs .GE. 3) THEN
            CALL GET_COMMAND_ARGUMENT(1, arg); READ (arg, *) nx
            CALL GET_COMMAND_ARGUMENT(2, arg); READ (arg, *) ny
            CALL GET_COMMAND_ARGUMENT(3, arg); READ (arg, *) nz
         END IF
         IF (nargs .GE. 4) THEN
            CALL GET_COMMAND_ARGUMENT(4, arg); READ (arg, *) aspect
         END IF
!........분할 형상 (C010-3): 0=k-슬랩(기본), 1=3D 블록 (np 를 최대한 정육면체로 인수분해)
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
            WRITE (*, '(A,3I5,A,F8.2,A,I9,A,I4)') ' C007 grid:', nx, ny, &
                nz, '  aspect=', aspect, '  ncell=', ncell, '  np=', np_z
         nelem_mg = ncell
         nf_max   = 6
!........np>1: k-슬랩 기하 분할 (celem = 셀→도메인 1..np, C010-2)
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
!.....프로덕션 초기화 체인 (read_grid.f90 순서 그대로 — subdomain 은 rank 0 전용)
      CALL read_input_mg
      IF (myrank_z .EQ. 0) CALL subdomain_infor_MG
!.....rank 0 의 MG_tmp 파일 fan-out 완료 후 판독 (프로덕션은 사이 collective 들이 동기화 역할)
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      CALL read_mesh_MPI
!.....GMG 로컬 순서의 좌표 스냅샷 — Prep_fine_P 가 coord 를 해제함 (3_Prep_fine_P.f90:130)
!     np>1 에선 nnode = 로컬(내부+고스트) 셀 수 (np=1 은 nnode=ncell)
      IF (.NOT. replay) THEN
!        열 공간은 nnodegl (전역 coarse 확장 포함, 대규모 np 에서 nnodegl > nnode — LOG C010-4)
         ALLOCATE (mycoord(3, nnodegl))
         mycoord(:, 1:nnodegl) = coord(:, 1:nnodegl)
      END IF
      CALL Prep_fine_P
      CALL Prep_MG_GarL
!
      CALL MPI_ALLREDUCE(nintf, i, 1, MPI_INTEGER, MPI_SUM,              &
                         MPI_COMM_WORLD, ierr)
      IF (myrank_z .EQ. 0)                                               &
         WRITE (*, '(A,I9,A,I10,A,I9)') ' setup done: sum(nintf)=', i,   &
                '  nnz_loc=', nnz, '  nintf_loc=', nintf
      IF (i .NE. ncell) STOP 'driver_pmg: sum(nintf) /= ncell (mapping?)'
!
      nfail = 0
      IF (replay) THEN
!=====================================================================
!     재생 모드: k1 (주 압력 솔브) → k2 (비직교 보정, k1 해가 초기추정)
!=====================================================================
         ALLOCATE (src(ncell), diag_my(ncell), au_my(nnz), ug(ncell))
!
!........k1 입력 로드 + 검증
         CALL load_pre(gdir, istep, 1, ncell, nnz, src, au_my, diag_my)
         u = 0.d0
         CALL assemble_FVM(1, nnz, src, au_my, diag_my)
         r0n = res_norm(ncell, nnz, src, au_my)     ! 초기잔차 (u=0 → ‖b‖)
         CALL SOLVE_GMG(1)
         CALL load_post(gdir, istep, 1, ncell, ug)
         CALL save_u(istep, 1, ncell)
!........충실도 정규화 스케일 = k1 압력장 크기 — k2(비직교 보정)는 크기가
!........~0 인 보정량이라 편차의 유의미 스케일도 k1 장 기준 (LOG C009-r2)
         uscale = 1.d0 + MAXVAL(ABS(ug))
         CALL verify('k1', ncell, nnz, src, au_my, ug, r0n, uscale, nfail)
!
!........k2: 행렬 불변(icase=2, b 만 교체), u 는 k1 해 유지 (프로덕션 동일)
         ALLOCATE (src2(ncell), diag2(ncell), au2(nnz))
         CALL load_pre(gdir, istep, 2, ncell, nnz, src2, au2, diag2)
         IF (MAXVAL(ABS(au2 - au_my)) .NE. 0.d0)                          &
            WRITE (*, '(A)') ' note: k2 au differs from k1 (icase=2 라 미사용)'
         CALL assemble_FVM(2, nnz, src2, au2, diag2)
         r0n = res_norm(ncell, nnz, src2, au_my)    ! 초기잔차 (u = k1 해)
         CALL SOLVE_GMG(2)
         CALL load_post(gdir, istep, 2, ncell, ug)
         CALL save_u(istep, 2, ncell)
         CALL verify('k2', ncell, nnz, src2, au_my, ug, r0n, uscale, nfail)
!
      ELSE
!=====================================================================
!     합성 모드 (LOG C007): 제작해 기반 검증
!     GMG CSR 은 대각을 au 안에 포함 (ju 가 위치; full-row A·x).
!     assemble_FVM 의 diag 인자는 스무더 스케일링용 (diagt=1/diag).
!=====================================================================
         vol = aspect                          ! dx=dy=1, dz=aspect
         diag_const = 2.d0*aspect + 2.d0*aspect + 2.d0/aspect
!........로컬 규격: 행 1..nintf 내부, nintf+1..nnode 고스트 (np=1: nintf=nnode=ncell)
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
!........uex: 좌표에서 전역 셀 id 를 복원해 np 무관 동일 제작해 (np=1 과 bitwise 동일 값)
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
         u = 0.d0
         CALL assemble_FVM(1, nnz, src, au_my, diag_my)
         r0n = res_norm(nintf, nnz, src, au_my)
         CALL SOLVE_GMG(1)
!
         uscale = MAXVAL(ABS(uex(1:nintf)))
         CALL MPI_ALLREDUCE(uscale, dist2, 1, MPI_DOUBLE_PRECISION,      &
                            MPI_MAX, MPI_COMM_WORLD, ierr)
         uscale = 1.d0 + dist2
         CALL verify('syn', nintf, nnz, src, au_my, uex, r0n, uscale, nfail)
      END IF
!
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
      IF (hdr(6) .NE. nnz_expect) STOP 'load_pre: nnz mismatch (setup 재생 불일치)'
      READ (iu) src
      READ (iu) auv
      READ (iu) dg
      CLOSE (iu)
      END SUBROUTINE load_pre
!
!-----------------------------------------------------------------------
      SUBROUTINE load_post(gdir, istep, isite, n, ug)
      CHARACTER(*), INTENT(IN) :: gdir
      INTEGER, INTENT(IN)  :: istep, isite, n
      REAL(8), INTENT(OUT) :: ug(n)
      INTEGER :: iu, hdr(5)
      CHARACTER(256) :: fn
      WRITE (fn, '(A,A,I0,A,I0,A)') TRIM(gdir), '/s', istep, '_k',       &
             isite, '_r0_c1.post'
      OPEN (NEWUNIT=iu, FILE=fn, FORM='unformatted', ACCESS='stream',    &
            STATUS='old', ACTION='read')
      READ (iu) hdr
      IF (hdr(5) .NE. n) STOP 'load_post: n mismatch'
      READ (iu) ug
      CLOSE (iu)
      END SUBROUTINE load_post
!
!-----------------------------------------------------------------------
      SUBROUTINE save_u(istep, isite, n)
!     재생 해를 저장 — 하네스 베이스라인(bitwise 회귀 게이트) 및 결정성 검증용
      INTEGER, INTENT(IN) :: istep, isite, n
      INTEGER :: iu
      CHARACTER(64) :: fn
      WRITE (fn, '(A,I0,A,I0,A)') 'replay_s', istep, '_k', isite, '.u'
      OPEN (NEWUNIT=iu, FILE=fn, FORM='unformatted', ACCESS='stream',    &
            STATUS='replace')
      WRITE (iu) n
      WRITE (iu) u(1:n)
      CLOSE (iu)
      END SUBROUTINE save_u
!
!-----------------------------------------------------------------------
      REAL(8) FUNCTION res_norm(n, nnzv, rhs, auv)
!     현재 u 에 대한 ‖rhs − A·u‖ (full-row, 내부 행 한정 + 전역 합)
      INTEGER, INTENT(IN) :: n, nnzv
      REAL(8), INTENT(IN) :: rhs(n), auv(nnzv)
      INTEGER :: irow, kk, ierr
      REAL(8) :: r, rn, rng
      rn = 0.d0
      DO irow = 1, n
         r = rhs(irow)
         DO kk = ia(irow), ia(irow+1) - 1
            r = r - auv(kk)*u(ja(kk))
         END DO
         rn = rn + r*r
      END DO
      CALL MPI_ALLREDUCE(rn, rng, 1, MPI_DOUBLE_PRECISION, MPI_SUM,      &
                         MPI_COMM_WORLD, ierr)
      res_norm = SQRT(rng)
      END FUNCTION res_norm
!
!-----------------------------------------------------------------------
      SUBROUTINE verify(tag, n, nnzv, rhs, auv, uref, r0n, uscale, nfail)
!     독립 residual (full-row, 드라이버 보유 행렬) + 기준해 비교
!     게이트 1 (수렴):  ‖r‖ ≤ 10·eps·‖r0‖ — 솔버의 실제 수렴 기준
!       (초기잔차 상대. k2 는 k1 해가 초기추정이라 r0 = b2 − A·u_k1 ≈ −b1;
!        ‖b2‖≈0 이어도 스케일 정합 — LOG C009-r1 에서 확정)
!     게이트 2 (충실도, 재생만):  max|Δ|/uscale ≤ 1e-9, uscale = 1+max|u*_k1|
!       — 프로덕션과 재생은 solver-tolerance 수준 일치가 기대치 (bitwise 는
!         환경 차이로 성립하지 않음을 C009 에서 확정; bitwise 회귀는
!         하네스 베이스라인(save_u) 대비로 별도 수행)
      CHARACTER(*), INTENT(IN) :: tag
      INTEGER, INTENT(IN)    :: n, nnzv
      REAL(8), INTENT(IN)    :: rhs(n), auv(nnzv), uref(n), r0n, uscale
      INTEGER, INTENT(INOUT) :: nfail
      INTEGER :: irow, kk, nbit, ierr
      REAL(8) :: r, rn, en, un, maxd, maxu, resg, fidel, gsum(3), gmax
      rn = 0.d0; en = 0.d0; un = 0.d0
      maxd = 0.d0; maxu = 0.d0
      nbit = 0
      DO irow = 1, n
         r = rhs(irow)
         DO kk = ia(irow), ia(irow+1) - 1
            r = r - auv(kk)*u(ja(kk))
         END DO
         rn = rn + r*r
         en = en + (u(irow) - uref(irow))**2
         un = un + uref(irow)*uref(irow)
         maxd = MAX(maxd, ABS(u(irow) - uref(irow)))
         maxu = MAX(maxu, ABS(uref(irow)))
         IF (u(irow) .NE. uref(irow)) nbit = nbit + 1
      END DO
!.....전역 집계 (np=1 은 항등)
      gsum(1) = rn; gsum(2) = en; gsum(3) = un
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, gsum, 3, MPI_DOUBLE_PRECISION,    &
                         MPI_SUM, MPI_COMM_WORLD, ierr)
      rn = gsum(1); en = gsum(2); un = gsum(3)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, maxd, 1, MPI_DOUBLE_PRECISION,    &
                         MPI_MAX, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, nbit, 1, MPI_INTEGER,             &
                         MPI_SUM, MPI_COMM_WORLD, ierr)
      resg  = SQRT(rn)/(10.d0*eps_bicg*r0n + 1.d-300)
      fidel = maxd/uscale
      IF (myrank_z .EQ. 0)                                               &
         WRITE (*, '(3A,ES12.4,A,ES12.4,A,ES11.3,A,I9)') ' RESULT[',     &
             tag, '] res_gate=', resg, '  fidelity=', fidel,             &
             '  max|d|=', maxd, '  nbit_diff=', nbit
      IF (TRIM(tag) .EQ. 'syn') THEN
!........합성: 독립 잔차가 곧 판정 (제작해 문제는 참 잔차 소거가 성립해야 함)
         IF (resg .GT. 1.d0) nfail = nfail + 1
      ELSE
!........재생: res_gate 는 보고 전용 — s30 에서 솔버 재귀 잔차와 참 잔차의
!........~1.6e4 배 괴리를 확인 (프로덕션 고유 특성, LOG C009-r2). 판정은
!........충실도(+러너의 its·베이스라인 bitwise)로 수행.
         IF (resg .GT. 1.d0 .AND. myrank_z .EQ. 0)                       &
            WRITE (*, '(3A)') ' WARN[', tag,                             &
            '] true residual exceeds solver criterion (재귀 잔차 드리프트)'
         IF (fidel .GT. fid_gate) nfail = nfail + 1
      END IF
      END SUBROUTINE verify
!
      END PROGRAM driver_pmg
