!
      SUBROUTINE v_tangent(velt,xn,v)
!
!     This routine calculates tangential velocities.
!       
      USE Zparam, ONLY: ndim
!
      IMPLICIT NONE
!
!     input
      REAL(8) xn(ndim),v(ndim)
!     output
      REAL(8) velt
!     local variables
      REAL(8) vnorm
! 
!.....velt=DMAX1(velt,1.d-20) was considered.
!      
!
      IF(ndim.eq.2) THEN
         vnorm=xn(1)*v(1)+xn(2)*v(2)
         velt=(v(1)-xn(1)*vnorm)**2+(v(2)-xn(2)*vnorm)**2
         velt=SQRT(velt)
      ELSE
         vnorm=xn(1)*v(1)+xn(2)*v(2)+xn(3)*v(3)
         velt=(v(1)-xn(1)*vnorm)**2+(v(2)-xn(2)*vnorm)**2+(v(3)-xn(3)*vnorm)**2
         velt=SQRT(velt)
      ENDIF
!
      RETURN
      END  SUBROUTINE v_tangent
