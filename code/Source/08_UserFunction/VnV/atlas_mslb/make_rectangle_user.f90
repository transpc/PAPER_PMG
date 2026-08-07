      SUBROUTINE make_rectangle_user(x,y,node)
!
      IMPLICIT NONE
!
      INTEGER i,j,node(4)
      INTEGER sortopt
      REAL (8) x(4),y(4),xmid,ymid,theta(4),storet,storen,storex,storey
      xmid=(x(1)+x(2)+x(3)+x(4))/4.0d0
      ymid=(y(1)+y(2)+y(3)+y(4))/4.0d0
      x(:)=x(:)-xmid
      y(:)=y(:)-ymid
!.....make 4 angles
      DO i=1,4
         CALL calc_angle_radian(x(i),y(i),theta(i))  
      ENDDO      
!.....sort 4 angles
       DO j = 1, 10000
        sortopt = 0
        DO i = 1, 3
         IF(theta(i) .gt. theta(i+1) )then
           storet = theta(i)
           storen = node(i)
           storex = x(i)
           storey = y(i)
           
           theta(i) = theta(i+1)
           node(i) = node(i+1)
           x(i) = x(i+1)
           y(i) = y(i+1)
           
           theta(i+1) = storet
           node(i+1) = storen
           x(i+1)=storex
           y(i+1)=storey

           sortopt = sortopt + 1
         ENDIF
        ENDDO !i

        IF(sortopt .eq. 0) exit
      ENDDO !j
      x(:)=x(:)+xmid
      y(:)=y(:)+ymid

       RETURN
       ENDSUBROUTINE make_rectangle_user
!---------------------------------------------------------------------           
      SUBROUTINE calc_angle_degree(x,y,theta)
!
      USE Zparam         , ONLY: pi
!
      IMPLICIT NONE
!
      REAL (8) x,y,r,theta
!      
!.....calculate angle from positive x-axis
!
       r=(x**2.d0+y**2.d0)**0.5d0
       IF(r.ne.0.d0)then
          theta=x/r
          theta=dacos(theta) !theta=0~180
       ELSE
          PRINT *,'r is 0 in make_rectangle_point!'
          PAUSE
          STOP
       ENDIF               
       IF(y.gt.0.d0)then
          theta=theta !0~180 ==> 0-180      
       ELSEIF(y.le.0.d0)then
          theta=-theta !-0 ~ -180
          theta=theta+2.0d0*pi ! 360~180 
       ELSE
          PRINT *,'theta cannot be defined in make_rectangle_point!'
          STOP
       ENDIF
       theta=theta*180/pi
!
       RETURN
       ENDSUBROUTINE calc_angle_degree   
!---------------------------------------------------------------------           
      SUBROUTINE calc_angle_radian(x,y,theta)
!
      USE Zparam         , ONLY: pi
!
      IMPLICIT NONE
!
      REAL (8) x,y,r,theta
!      
!.....calculate angle from positive x-axis
!
       r=(x**2.d0+y**2.d0)**0.5d0
       IF(r.ne.0.d0)then
          theta=x/r
          theta=dacos(theta) !theta=0~180
       ELSE
          PRINT *,'r is 0 in calc_angle_radian!'
          PAUSE
          STOP
       ENDIF               
       IF(y.gt.0.d0)then
          theta=theta !0~180 ==> 0-180      
       ELSEIF(y.le.0.d0)then
          theta=-theta !-0 ~ -180
          theta=theta+2.0d0*pi ! 360~180 
       ELSE
          PRINT *,'theta cannot be defined in calc_angle_radian!'
          STOP
       ENDIF
!      
!       IF(x.ne.0.0d0)write(*,"(5(1pe12.3))")x/r,DCOS(theta),x,y,theta
!       IF(x.ne.0.0d0)write(97,"(5(1pe12.3))")x/r,DCOS(theta),x,y,theta       
! 
       RETURN
       ENDSUBROUTINE calc_angle_radian          
