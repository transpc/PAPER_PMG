! ====================================================================!
!  HYPRE_FINALIZE — hypre 핸들 해제 (계산 종료 시 1회)
!
!  MPI_FINALIZE 보다 반드시 먼저 불려야 한다 (hypre 오브젝트 파괴가
!  통신자를 참조). cupid_main 의 반환 직전에서 호출한다.
!  셋업이 없으면(= hypre 경로를 쓰지 않은 실행) 아무 일도 하지 않는다.
! ---------------------------------------------------------------------+
      SUBROUTINE HYPRE_FINALIZE()
!
      USE MD_hypre, ONLY: hypre_free
!
      IMPLICIT NONE
!
      CALL hypre_free()
!
      RETURN
      END SUBROUTINE HYPRE_FINALIZE
