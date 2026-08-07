!
      SUBROUTINE fanning_factor(rey, rough, diameter, f)
!
!     This routine calculates fanning factor used for wall drag coefficient.
!
      IMPLICIT NONE
!
!      INTEGER i
!      
      REAL(8) a, b, c, f
      REAL(8) rey,rough,diameter !pik 
!
!     calculate fanning factor
!
      c=(7.0/rey)**0.9+0.27*DMIN1(DMAX1(rough/diameter,1.0e-9),0.02)
      a=(2.475*DLOG(c))**16
      b=(37530.0/rey)**16
      f=2.0*((8.0/rey)**12+1.0/(a+b)**(3.0/2.0))**(1.0/12.0)
!
      RETURN
      END SUBROUTINE fanning_factor
