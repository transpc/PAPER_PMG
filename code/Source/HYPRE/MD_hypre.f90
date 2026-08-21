! ====================================================================!
!  MD_hypre — CUPID 압력계를 hypre(BiCGSTAB + BoomerAMG)로 푸는 경로의
!             상태·파라미터·내부 루틴 일체
!
!  위치: PMG(SOLVE_GMG)와 같은 층. assemble_FVM 이 채운 fine CSR
!        (MD_matrix 의 ia/ja/au/b/u, MD_MPI 의 nintf)을 그대로 받는다.
!        진입점은 SOLVE_HYPRE(icase) — SOLVE_GMG(icase)와 같은 시그니처.
!
!  구조 (PMG 와 대칭):
!    외부 반복 : BiCGSTAB           <-> solve_pbcg_mg
!    예조건자  : BoomerAMG V(1,1)   <-> SOLVER_NEW (ncycle=1)
!                Tol(0)+MaxIter(1) 로 "적용당 V-cycle 1회" 고정.
!                (hypre 기본 max_iter=20/tol=1e-6 이면 예조건자가 아니라
!                 완전 솔브가 되어 비교가 성립하지 않는다)
!
!  파라미터는 현재 전부 하드와이어(아래 PARAMETER 블록). 추후 mg.in 계열
!  파일 입력으로 이관 예정 — 그때 이 블록만 읽기 대상으로 바꾸면 된다.
!  단 수렴 tol 은 예외로 Zbicg 의 eps_bicg 를 그대로 쓴다. PMG 의
!  crit_bcg_mg 도 같은 변수에서 오므로(1_read_input.f90) 두 솔버가 한 개의
!  tol 을 공유하는 것이 코드로 보장된다 — 비교 공정성의 전제.
!---------------------------------------------------------------------+
      MODULE MD_hypre
!
      IMPLICIT NONE
      SAVE
!
!.....BoomerAMG 설정 (하드와이어). 음수는 "hypre 기본값 유지" 를 뜻하며
!     해당 Set 루틴을 아예 호출하지 않는다.
!       relax  : <0 = hypre 2.33 기본 l1-GS forward/backward(13/14),
!                최조대는 가우스소거(9). 분할-강건이라 np 스케일 비교의
!                정당한 상대다. 6 이면 hybrid sym-GS/SSOR(비-l1).
!       strong : <0 = hypre 기본 0.25
!       조대화 10(HMIS)·보간 6(ext+i)·P_max 4 는 hypre 기본값과 같으므로
!       명시하지 않는다 (기본값 추종).
      INTEGER, PARAMETER :: hyp_sweeps = 1        ! 예조건자 레벨당 스무딩 횟수
      INTEGER, PARAMETER :: hyp_relax  = -1       ! <0 = 기본(13/14)
      REAL(8), PARAMETER :: hyp_strong = -1.d0    ! <0 = 기본(0.25)
      INTEGER, PARAMETER :: hyp_agg    = 0        ! aggressive coarsening 레벨수
      INTEGER, PARAMETER :: hyp_print  = 0        ! 3 이면 AMG 계층표 출력
      INTEGER, PARAMETER :: hyp_maxit  = 1000     ! BiCGSTAB 최대 반복
!
!.....hypre 상수 (hypre_IJ 오브젝트 타입)
      INTEGER, PARAMETER :: HYPRE_PARCSR = 5555
!
!.....hypre 핸들 (셋업 재사용을 위해 모듈 스코프에 보존)
      INTEGER(8) :: hA, hb, hx, hparA, hparb, hparx, hsol, hpre
      LOGICAL    :: hyp_alive  = .FALSE.   ! 위 핸들이 살아 있는가
      LOGICAL    :: gnum_ready = .FALSE.   ! 전역 행번호가 준비되었는가
!
!.....전역 행 번호 (rank 별 내부행 연속 소유 [ioffs+1, ioffs+nintf])
      INTEGER :: ioffs
      INTEGER, ALLOCATABLE :: gcol(:), grows(:)
!
!.....잔차방정식 작업 배열 (hypre_gnum 에서 1회 할당)
      REAL(8), ALLOCATABLE :: hyp_wrk(:)   ! (nnode) r0 = b - A*u0
      REAL(8), ALLOCATABLE :: hyp_dlt(:)   ! (nintf) 보정해 d
!
!.....직전 솔브 계측 (러너/논문 계측용 — 소비처는 SOLVE_HYPRE)
      REAL(8) :: hyp_t_asm = 0.d0   ! IJ 조립
      REAL(8) :: hyp_t_amg = 0.d0   ! BoomerAMG + BiCGSTAB 셋업
      REAL(8) :: hyp_t_sol = 0.d0   ! BiCGSTAB 반복
      REAL(8) :: hyp_relres = 0.d0  ! hypre 자기보고 상대잔차
      INTEGER :: hyp_its = 0        ! BiCGSTAB 반복수
!
      CONTAINS
!
!-----------------------------------------------------------------------
      SUBROUTINE hypre_gnum()
!     전역 행 번호 구성 (최초 1회).
!     rank 는 내부행 nintf 개를 [ioffs+1, ioffs+nintf] 로 연속 소유한다
!     (hypre IJ 의 연속 소유 요구). 고스트 열의 전역 번호는 gid 를 실수
!     벡터에 실어 PMG 헤일로 교환(communicate_s)으로 획득하므로 분할
!     형상에 의존하지 않는다.
!     전제: Prep_fine_P 까지 끝나 ia/ja 가 확정되어 있을 것.
!
      USE MD_matrix,    ONLY: nnz, ia, ja
      USE MD_MPI,       ONLY: nintf, myrank
      USE MD_geometry,  ONLY: nnode
      USE MD_parameter, ONLY: ndom
!
      IMPLICIT NONE
!DEC$IF defined (mpi_flag)
      INCLUDE 'mpif.h'
      INTEGER :: ierr
!DEC$ENDIF
      INTEGER :: i, irow, kk, nbad
      CHARACTER :: report_text*100
      REAL(8), ALLOCATABLE :: gidr(:)
!
      ioffs = 0
!DEC$IF defined (mpi_flag)
      CALL MPI_EXSCAN(nintf, ioffs, 1, MPI_INTEGER, MPI_SUM,             &
                      MPI_COMM_WORLD, ierr)
      IF (myrank .EQ. 0) ioffs = 0
!DEC$ENDIF
!
      ALLOCATE (gidr(nnode))
      gidr = 0.d0
      DO i = 1, nintf
         gidr(i) = DBLE(ioffs + i)
      END DO
      IF (ndom .GT. 1) CALL communicate_s(gidr)
!
      ALLOCATE (gcol(nnz), grows(nintf))
      nbad = 0
      DO irow = 1, nintf
         grows(irow) = ioffs + irow
         DO kk = ia(irow), ia(irow+1) - 1
            IF (gidr(ja(kk)) .LT. 0.5d0) nbad = nbad + 1
            gcol(kk) = NINT(gidr(ja(kk)))
         END DO
      END DO
      DEALLOCATE (gidr)
!
      IF (nbad .GT. 0) THEN
         WRITE (*, *) 'HYPRE error: ghost global id unresolved, rank', myrank, ' n=', nbad
         report_text = 'HYPRE: ghost global id unresolved (communicate_s?)'
         CALL STOP_MPI(report_text)
      END IF
!
!.....잔차방정식용 작업 배열 (매 솔브 할당을 피해 상주)
      ALLOCATE (hyp_wrk(nnode), hyp_dlt(MAX(nintf,1)))
      hyp_wrk = 0.d0
      hyp_dlt = 0.d0
!
      gnum_ready = .TRUE.
      END SUBROUTINE hypre_gnum
!
!-----------------------------------------------------------------------
      SUBROUTINE hypre_build()
!     현재 au 로 IJ 행렬을 조립하고 BiCGSTAB + BoomerAMG 셋업까지 수행.
!     PMG 의 SOLVE_GMG(icase=1) 내부 RAP 셋업에 대응한다 — 행렬이 안 바뀌는
!     후속 솔브는 hypre_apply 만 다시 부르면 된다.
!
!     재진입 안전: 기존 핸들이 살아 있으면 먼저 해제한다. 프로덕션은 매
!     타임스텝 호출되므로 이 가드가 없으면 IJ 오브젝트와 AMG 계층이 계속
!     쌓인다.
!
      USE MD_matrix,   ONLY: nnz, ia, au
      USE MD_MPI,      ONLY: nintf
      USE Zbicg,       ONLY: eps_bicg
!
      IMPLICIT NONE
!DEC$IF defined (mpi_flag)
      INCLUDE 'mpif.h'
!DEC$ENDIF
      INTEGER :: irow, n1, one, ierr, rows1(1)
      REAL(8) :: tw
      REAL(8), ALLOCATABLE :: z0(:)
!
      CALL hypre_free()          ! 재진입 가드 (살아 있지 않으면 즉시 반환)
!
      one = 1
      tw  = hyp_wtime()
!
!.....행렬 조립 (행 단위)
      CALL HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ioffs+1, ioffs+nintf,    &
                                ioffs+1, ioffs+nintf, hA, ierr)
      CALL HYPRE_IJMatrixSetObjectType(hA, HYPRE_PARCSR, ierr)
      CALL HYPRE_IJMatrixInitialize(hA, ierr)
      DO irow = 1, nintf
         n1 = ia(irow+1) - ia(irow)
         rows1(1) = ioffs + irow
         CALL HYPRE_IJMatrixSetValues(hA, one, n1, rows1,                &
                 gcol(ia(irow):ia(irow+1)-1), au(ia(irow):ia(irow+1)-1), ierr)
      END DO
      CALL HYPRE_IJMatrixAssemble(hA, ierr)
      CALL HYPRE_IJMatrixGetObject(hA, hparA, ierr)
!
!.....벡터 b, x (값은 hypre_apply 에서 매 솔브 갱신)
      ALLOCATE (z0(MAX(nintf,1)))
      z0 = 0.d0
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
      hyp_t_asm = hyp_wtime() - tw
      tw = hyp_wtime()
!
!.....BoomerAMG 예조건자 — Tol(0)+MaxIter(1) 로 "적용당 V-cycle 1회" 고정
      CALL HYPRE_BoomerAMGCreate(hpre, ierr)
      CALL HYPRE_BoomerAMGSetTol(hpre, 0.d0, ierr)
      CALL HYPRE_BoomerAMGSetMaxIter(hpre, 1, ierr)
      CALL HYPRE_BoomerAMGSetNumSweeps(hpre, hyp_sweeps, ierr)
      CALL HYPRE_BoomerAMGSetPrintLevel(hpre, hyp_print, ierr)
      IF (hyp_relax  .GE. 0)    CALL HYPRE_BoomerAMGSetRelaxType(hpre, hyp_relax, ierr)
      IF (hyp_strong .GT. 0.d0) CALL HYPRE_BoomerAMGSetStrongThrshld(hpre, hyp_strong, ierr)
      IF (hyp_agg    .GT. 0)    CALL HYPRE_BoomerAMGSetAggNumLevels(hpre, hyp_agg, ierr)
!
!.....BiCGSTAB (수렴 판정: ||r||/||b|| <= tol. 잔차방정식 형태로 호출하므로
!     ||b|| = ||r0|| 이 되어 PMG 의 ro/ro0 <= eps 와 분모가 일치한다)
      CALL HYPRE_ParCSRBiCGSTABCreate(MPI_COMM_WORLD, hsol, ierr)
      CALL HYPRE_ParCSRBiCGSTABSetTol(hsol, eps_bicg, ierr)
      CALL HYPRE_ParCSRBiCGSTABSetMaxIter(hsol, hyp_maxit, ierr)
      CALL HYPRE_ParCSRBiCGSTABSetPrintLev(hsol, 0, ierr)
      CALL HYPRE_ParCSRBiCGSTABSetPrecond(hsol, 2, hpre, ierr)
      CALL HYPRE_ParCSRBiCGSTABSetup(hsol, hparA, hparb, hparx, ierr)
!
      hyp_alive = .TRUE.
      hyp_t_amg = hyp_wtime() - tw
      END SUBROUTINE hypre_build
!
!-----------------------------------------------------------------------
      SUBROUTINE hypre_apply(rhs, x)
!     A*x = rhs 를 x0=0 에서 푼다 (rhs, x 는 내부행 nintf 개).
!     반복수·상대잔차는 모듈 변수(hyp_its, hyp_relres)로 회수.
!
      USE MD_MPI, ONLY: nintf
!
      IMPLICIT NONE
      REAL(8), INTENT(IN)  :: rhs(*)
      REAL(8), INTENT(OUT) :: x(*)
      INTEGER :: ierr
      REAL(8) :: tw
!
      x(1:nintf) = 0.d0
      CALL HYPRE_IJVectorSetValues(hb, nintf, grows, rhs, ierr)
      CALL HYPRE_IJVectorSetValues(hx, nintf, grows, x,   ierr)
!
      tw = hyp_wtime()
      CALL HYPRE_ParCSRBiCGSTABSolve(hsol, hparA, hparb, hparx, ierr)
      hyp_t_sol = hyp_wtime() - tw
!
      CALL HYPRE_ParCSRBiCGSTABGetNumIter(hsol, hyp_its, ierr)
      CALL HYPRE_ParCSRBiCGSTABGetFinalRel(hsol, hyp_relres, ierr)
      CALL HYPRE_IJVectorGetValues(hx, nintf, grows, x, ierr)
      END SUBROUTINE hypre_apply
!
!-----------------------------------------------------------------------
      SUBROUTINE hypre_free()
!     핸들 해제. 살아 있지 않으면 아무것도 하지 않는다 (재진입/종료 공용).
      IMPLICIT NONE
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
      REAL(8) FUNCTION hyp_wtime()
!     단조 벽시계. MPI 빌드에서는 MPI_Wtime, 아니면 SYSTEM_CLOCK.
      IMPLICIT NONE
!DEC$IF defined (mpi_flag)
      INCLUDE 'mpif.h'
      hyp_wtime = MPI_Wtime()
!DEC$ELSE
      INTEGER(8) :: cnt, rate
      CALL SYSTEM_CLOCK(COUNT=cnt, COUNT_RATE=rate)
      hyp_wtime = DBLE(cnt)/DBLE(rate)
!DEC$ENDIF
      END FUNCTION hyp_wtime
!
      END MODULE MD_hypre
