!
      SUBROUTINE udfn_weight(weight,x,x2,x1)
!
!     User-defined weight for heat transfer coefficients
!     (only when "i_weight" is used)
!
      IMPLICIT NONE
!
      REAL(8) weight,x,x2,x1
!
      weight=DMIN1(1.0d0,DMAX1(0.0d0,DEXP(-8.0d0*(x-x1)/(x2-x1))))
!      
      RETURN
      END SUBROUTINE udfn_weight