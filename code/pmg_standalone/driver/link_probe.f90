!     링크 클로저 프로브 (LOOP C006) — 전체 오브젝트를 링크해 미정의 외부 심볼을 노출시키는 용도
      PROGRAM link_probe
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER ierr
      CALL MPI_INIT(ierr)
      WRITE(*,*) 'link_probe: OK'
      CALL MPI_FINALIZE(ierr)
      END PROGRAM link_probe
