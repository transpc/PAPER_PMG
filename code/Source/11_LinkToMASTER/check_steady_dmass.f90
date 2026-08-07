!
      SUBROUTINE check_steady_dmass(dmass_max)
!
      USE Zcore        , ONLY: np,myrank
      USE MASTER4      , ONLY: dmass_pass,mas_dmass,count_dmass,dmass     
!      
      IMPLICIT NONE
      INCLUDE '../10_LinkToMARS/c3com.h' !n_junleg,mflow_junle
!      
      INTEGER i
      INTEGER,SAVE :: count=0
      LOGICAL,SAVE :: initial=.TRUE.
!      
      REAL(8) dmass_max
      REAL(8),ALLOCATABLE,SAVE::mflow_junleg_old(:)
!      
      IF(initial)THEN
         initial=.FALSE.
         ALLOCATE(dmass(i3n_junleg))
         IF(i3n_junleg.gt.20)THEN
           IF(myrank.eq.0)WRITE(*,"(11x,a)")'i3n_junleg should be less than 20!!!'
           PAUSE
           STOP
         ENDIF
         ALLOCATE(mflow_junleg_old(i3n_junleg))
         mflow_junleg_old(:)=c3mflow_junleg(1:i3n_junleg)   
      ENDIF
      dmass_max=0.0d0
      DO i=1,i3n_junleg
          dmass(i)=ABS(c3mflow_junleg(i)-mflow_junleg_old(i))
          dmass_max=MAX(dmass_max,dmass(i))
      ENDDO            
      IF(np.gt.1) CALL allreducei_max_r1(dmass_max)
!
      IF(dmass_max.lt.mas_dmass)THEN
         count=count+1
      ELSE
         count=0
      ENDIF   
      IF(count.ge.count_dmass)dmass_pass=1
!
      mflow_junleg_old(:)=c3mflow_junleg(1:i3n_junleg)       
!
      CALL check_steady_dp
!
      END SUBROUTINE check_steady_dmass
