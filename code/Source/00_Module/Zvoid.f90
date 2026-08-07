      MODULE Zvoid
! 
      IMPLICIT NONE
      SAVE
!
      REAL(8) gradient_void
!     daldx not needed anymore was used previously in pressure_solve
!     REAL(8),Allocatable::dagdx(:,:),daldx(:,:),gamma_void(:)
      REAL(8),Allocatable::dagdx(:,:),gamma_void(:)
!
      END MODULE Zvoid
