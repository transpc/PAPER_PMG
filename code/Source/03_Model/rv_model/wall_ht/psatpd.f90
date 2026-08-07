!
      SUBROUTINE psatpd_cupid(t,press,presdt,itype,err) 
!
!  Calculate saturation pressure (press) and dpdt (presdt) for a given
!  temperature (t) if itype=1.
!  Calculate saturation temperature (t) and dpdt (presdt) for a given
!  pressure (press) if itype=2.
!
      IMPLICIT none 
!
      REAL(8) t,press,presdt 
      INTEGER itype 
      LOGICAL err 
!
!  Local variables.
      INTEGER i 
      LOGICAL istop 
      REAL(8) a1,a2,a3,crp,dt,k(9),p,pxxy(5),tcinv,theta,theta1,txxy(5) 
!
      DATA k/-7.691234564d0,-26.08023696d0,-168.1706546d0,6.423285504d1,&
      -1.189646225d2,4.167117320d0,2.097506760d1,1.d9,6.d0/             
      DATA pxxy/1.3d5,9.3d5,2.6d6,5.2d6,9.4d6/ 
      DATA txxy/326.0d0,415.0d0,475.0d0,520.0d0,560.0d0/ 
      DATA crp/22120000.0d0/,tcinv/1.544878727d-3/ 
!
!
!     Calculate saturation pressure and dpdt for the given temperature.
      IF(itype.ne.2)then 
         IF(t.lt.273.16d0.or.t.gt.647.3d0) GOTO 201 
         theta=t*tcinv 
         theta1=1.0d0-theta 
         a1=theta1*(k(1)+theta1*(k(2)+theta1*(k(3)+theta1*(k(4)+k(5)*   &
            theta1))))
         a2=1.0d0/(theta*(1.0d0+theta1*(k(6)+k(7)*theta1))) 
         a3=1.0d0/(k(8)*theta1*theta1+k(9)) 
         press=crp*DEXP(a1*a2-theta1*a3) 
         presdt=-(a2*(k(1)+theta1*(2.0d0*k(2)+theta1*(3.0d0*k(3)+theta1*&
                (4.0d0*k(4)+5.0d0*k(5)*theta1)))+a1*(1.0d0+k(6)*(theta1-&
                theta)+k(7)*theta1*(1.0d0-theta * 3.0d0))*a2)+(k(8)*    &
                theta1*theta1-k(9))*a3**2)*press*tcinv
         err=.false. 
         RETURN 
      ELSE 
!
!        Find saturation temperature and dpdt for a given pressure.
!
         IF(press.gt.crp) GOTO 201 
         DO i=1,5 
            IF(press.lt.pxxy(i))then 
               t=txxy(i) 
               GOTO 106 
            ENDIF 
         END DO 
         t=614.0d0 
  106    istop=.false. 
         DO i=1,15 
            theta=t*tcinv 
            theta1=1.0d0-theta 
            a1=theta1*(k(1)+theta1*(k(2)+theta1*(k(3)+theta1*(k(4)+k(5)*&
               theta1))))
            a2=1.0d0/(theta*(1.0d0+theta1*(k(6)+k(7)*theta1))) 
            a3=1.0d0/(k(8)*theta1*theta1+k(9)) 
            p=crp*DEXP(a1*a2-theta1*a3) 
            presdt=-(a2*(k(1)+theta1*(2.0d0*k(2)+theta1*(3.0d0*k(3)+    &
                   theta1*(4.0d0*k(4)+5.0d0*k(5)*theta1)))+a1*(1.0d0+   &
                   k(6)*(theta1-theta)+k(7)*theta1*(1.0d0-theta*3.0d0))*&
                   a2)+(k(8)*theta1*theta1-k(9))*a3**2)*p*tcinv
            IF(istop) GOTO 134 
            dt=(press-p)/presdt 
            t=DMAX1(DMIN1(t+dt,647.3d0),273.16d0) 
            IF(DABS(dt).lt.t * 1.0d-5)istop=.true. 
         END DO 
  134    err=.false. 
         RETURN 
      ENDIF 
  201 err=.true. 
      RETURN 
      END SUBROUTINE psatpd_cupid