!
      SUBROUTINE udf_free_surface_tracking
!
!     Define free surface cell criteria with gamma
!
      USE Zvoid   , ONLY: gradient_void
!
      IMPLICIT NONE
!
      gradient_void = 0.4d0
!
      RETURN
      END SUBROUTINE udf_free_surface_tracking
!
