     SUBROUTINE print_mars 
!
      USE Zconst1     , ONLY: cplmaster,cplmars
      USE Zcore       , only: myrank
      USE c3com_cupid , ONLY: nvols_mars
      USE Zmars       , only: time_mars
      USE MASTER4     , ONLY: dpres       
!     
      IMPLICIT NONE
!      
!DEC$IF defined (MCC)
      INCLUDE 'c3com.h'  
!DEC$ELSEIF defined (MCC_DLL)      
      INCLUDE 'c3com.h'  
      !dec$ attributes dllexport :: print_mars      
!DEC$ENDIF      
!      
      INTEGER i  
!
      LOGICAL,SAVE:: initial
!     
      REAL(8),SAVE:: time_interval,print_time
      INTEGER,SAVE:: nprn
!
      DATA initial/.true./
!
      time_mars=c3time_sys
!     
      IF (initial) then
         time_interval=1.0d0
         print_time=0.0d0
         nprn=0
      ENDIF
!      
      IF(i3marsin.ne.0)RETURN  
!
!.....Print mars power ratio
!
      IF(time_mars.gt.print_time.and.cplmars.gt.0)then
         print_time=time_mars+time_interval
         !nprn=nprn+1
         !IF(nprn.eq.10)THEN
         !   nprn=0
         IF(myrank.eq.0)WRITE(*,"(11x,a,3(1pe11.3))")'time,power ratio,kfactor=',time_mars,c3rktpow_ctl_val,c3kfactor_ctl_val !92-171          
         ! ENDIF  
!
!.....Print dp      
         CALL check_steady_dp
         CALL print_ss_status_dp         
      ENDIF 
!
      initial=.false.      
!
      RETURN
      END SUBROUTINE print_mars 
            
      
          
