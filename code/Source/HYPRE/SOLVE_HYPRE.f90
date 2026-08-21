! ====================================================================!
!  SOLVE_HYPRE(icase) — 압력계를 hypre(BiCGSTAB + BoomerAMG)로 푼다.
!
!  SOLVE_GMG(icase) 와 같은 시그니처·같은 계약:
!    입력  : assemble_FVM 이 채운 MD_matrix 의 ia/ja/au/b, 초기추정 u
!    출력  : u (내부행 + 고스트 갱신까지)
!    icase = 1 : 행렬이 바뀌었다 → IJ 재조립 + AMG 셋업 후 솔브
!            2 : 행렬 불변, b 만 교체 → 셋업 재사용하고 솔브만
!            (PMG 의 icase 의미와 동일. stiffness_MG/RAP 재계산 여부에 대응)
!
!  잔차방정식 형태로 통일한다:
!      A*d = b - A*u0,  d0 = 0,  u = u0 + d
!  이유 세 가지.
!   (1) hypre BiCGSTAB 의 판정은 ||r||/||b|| <= tol 인데, 이 형태에서는
!       ||b|| = ||r0|| 이 되어 PMG 의 ro/ro0 <= eps 와 분모가 정확히 일치.
!   (2) 프로덕션 비직교 보정 루프(icase=2)는 u 를 0 으로 리셋하지 않고
!       직전 해에서 출발한다 — 그 warm start 를 그대로 흡수한다.
!   (3) icase=1 은 호출부가 u=0 을 넣으므로 r0 = b 로 자동 축약된다.
! ---------------------------------------------------------------------+
      SUBROUTINE SOLVE_HYPRE(icase)
!
      USE MD_hypre
      USE MD_matrix,    ONLY: nnz, ia, ja, au, u, b
      USE MD_MPI,       ONLY: nintf, myrank
      USE MD_geometry,  ONLY: nnode
      USE MD_parameter, ONLY: ndom
      USE Zbicg,        ONLY: eps_bicg
!
      IMPLICIT NONE
!
      INTEGER(4), INTENT(IN) :: icase
      INTEGER :: i
      CHARACTER :: report_text*100
!
!.....전역 행번호는 최초 1회만 (Prep_fine_P 이후면 언제든 유효)
      IF (.NOT. gnum_ready) CALL hypre_gnum()
!
      IF (icase .EQ. 1) THEN
         CALL hypre_build()
      ELSE
!........셋업 재사용 — 계측상 이번 스텝의 셋업 비용은 0
         hyp_t_asm = 0.d0
         hyp_t_amg = 0.d0
      END IF
!
      IF (.NOT. hyp_alive) THEN
         report_text = 'HYPRE: icase=2 called before any icase=1 setup'
         CALL STOP_MPI(report_text)
      END IF
!
!.....r0 = b - A*u0   (u0 의 고스트가 유효해야 amux0P 가 성립)
      IF (ndom .GT. 1) CALL communicate_s(u)
      CALL amux0P(nintf, nnode, nnz, u, hyp_wrk, au, ja, ia)
      DO i = 1, nintf
         hyp_wrk(i) = b(i) - hyp_wrk(i)
      END DO
!
!.....A*d = r0,  d0 = 0
      CALL hypre_apply(hyp_wrk, hyp_dlt)
!
!.....u = u0 + d,  고스트 갱신 (호출부의 grad_press 가 고스트를 읽는다)
      DO i = 1, nintf
         u(i) = u(i) + hyp_dlt(i)
      END DO
      IF (ndom .GT. 1) CALL communicate_s(u)
!
!.....솔브당 1행 기록 (rank0, fort.502) — PMG 가 fort.501 에 남기는
!     'its ro/ro0' 과 대칭. 계측 소비처: 논문용 타임스텝별 비용 곡선.
!       icase  its  relres  t_asm  t_amg  t_sol
      IF (myrank .EQ. 0)                                                 &
         WRITE (502, *) icase, hyp_its, hyp_relres,                      &
                        hyp_t_asm, hyp_t_amg, hyp_t_sol
!
!.....미수렴은 에러 채널(콘솔, rank0)로 — PMG 의 PBCG_MG 초과 보고와 동격
      IF (hyp_relres .GT. eps_bicg .AND. myrank .EQ. 0) THEN
         WRITE (*, '(A,I6,A,ES12.4,A,ES12.4)')                           &
             ' HYPRE error: BiCGSTAB not converged  its=', hyp_its,      &
             '  relres=', hyp_relres, '  tol=', eps_bicg
      END IF
!
      RETURN
      END SUBROUTINE SOLVE_HYPRE
