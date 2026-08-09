!     직렬(np=1) 전용 communicate 스텁 (LOOP C006)
!     프로덕션 semantics: 서브도메인 경계(고스트) 교환. np=1 에선 이웃이 없어 항등 연산.
!     MPI 하네스 확장 시(C009 이후) 06_MPI/communicate.f90 실물 + 모듈 클로저로 교체할 것.
      SUBROUTINE communicate(a,izone)
      IMPLICIT NONE
      INTEGER :: izone
      REAL(8) :: a(*)
      END SUBROUTINE communicate
