!
      SUBROUTINE cub_interp(x,x1,x2,y1,y2,dydx1,dydx2,c,dcdx)
!
!     Interploation using cubic polynomial
!
      IMPLICIT NONE
!      
      REAL(8) c,dcdx,x,x1,x2,y1,y2,dydx1,dydx2
      REAL(8) a1,a2,a3,a4,dx,dy,ddy,dxx,dxxx,xx1,xx,xxx
!
      dx=x2-x1
      dy=y2-y1
      ddy=dydx2-dydx1
      dxx=dx*dx
      dxxx=dxx*dx
      xx1=x1*x1
      xx=x*x
      xxx=xx*x
!
      a1=ddy/dxx+2.0d0*(dydx1/dxx-dy/dxxx)
      a2=(ddy/dx-3.0d0*(x2+x1)*a1)/2.0d0
      a3=dydx1-3.0d0*a1*xx1-2.0d0*a2*x1
      a4=y1-a1*xx1*x1-a2*xx1-a3*x1
!
      c=a1*xxx+a2*xx+a3*x+a4
      dcdx=3.0d0*xx*a1+2.0d0*x*a2+a3
!
      RETURN
      END SUBROUTINE cub_interp
!