!
      SUBROUTINE dambreaking_front_user(time)
!
!     Dam breaking problem only
!   
      USE VOL_DATA              
      USE Zcore    , ONLY: np,myrank
      USE Zcoord1  , ONLY: xloc
      USE Zzone    , ONLY: ncell_fluid
!
      IMPLICIT NONE
!
      INTEGER i,myix 
!
      LOGICAL, SAVE::INITIAL
!      
      REAL(8) time,time_non 
      REAL(8) alphag_crit,front
      REAL(8),SAVE::front_old      
!      
      PARAMETER(myix=1,alphag_crit=0.99d0)
!      
      DATA INITIAL/.TRUE./
!
      front=-1.d-10
      IF(myrank.eq.0.and.INITIAL)THEN
         front_old=front
         OPEN(40,file='vft8_vv.dat')
!         WRITE(40,5)
         INITIAL=.FALSE.
      ENDIF
!
      DO i=1,ncell_fluid
         IF(cell%alphag(i).lt.alphag_crit)THEN
            IF(front.lt.xloc(i,myix)) front = xloc(i,myix)
         ENDIF
      ENDDO
!
      IF(np.gt.1) CALL allreducei_max_r1(front)
!
      IF(myrank.eq.0)THEN
         time_non=time*19.8090882d0
         front=front/0.05d0
         IF(front_old.eq.front)RETURN         
         WRITE(40,10) time_non,front
         front_old=front
      ENDIF
!
!5     FORMAT('variables = non-dimensional time,non-dimensional x-front')
10    FORMAT(2f15.10)
!    
      RETURN
      END SUBROUTINE dambreaking_front_user
