      MODULE Zturbzeq
! 
      IMPLICIT NONE
      SAVE
!
      INTEGER,Allocatable::cell_Hindex(:)
      REAL(8) tlengs
      REAL(8),Allocatable::tleng(:),vorticity(:),ChHeight(:,:)
!.....input
      CHARACTER(30) s_turb_zero
!
      END MODULE Zturbzeq