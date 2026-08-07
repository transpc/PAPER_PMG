!
      SUBROUTINE change_tend_mars_cupid
!
      USE Ztimecon        , ONLY: ctrl_opt,ictrl,nctrl,time,t_end,t_end_ctrl
      USE Zconst1         , ONLY: cplmars
      USE Zcore           , ONLY: myrank
      USE MASTER4         , ONLY: dmaster_pass,mas_delay
      USE Zmars           , ONLY: time_mars
!     
      IMPLICIT NONE
!      
      INCLUDE 'c3com.h'
!      
      LOGICAL,SAVE::initial
!      
      INTEGER master_pass,i
!      
      REAL(8) time_end     
!      
      DATA initial/.TRUE./
!      
      IF(dmaster_pass.eq.0)RETURN
!
      IF(cplmars.ge.2.and.cplmars.le.3)THEN      
!      
         IF(initial)THEN
            initial=.FALSE.
!      
!...........set the end time      
            time_end=time_mars+mas_delay
!
!...........set the end of time of MARS
            c3tend=time_end
!         
!...........set the end time of CUPID
            t_end=time_end
!            
            IF(myrank.eq.0)THEN
               WRITE(*,"(a,11f12.6)")'t_end=',t_end
               write(662,"(a,11f12.6)")'t_end=',t_end
               WRITE(97,"(a,11f12.6)")'t_end=',t_end
            ENDIF   
            DO i=ictrl,nctrl
               t_end_ctrl(i)=time_end
               IF(myrank.eq.0)THEN
                  WRITE(*,"(a,1i5,1f12.6)")'i,t_end_ctrl(i)=',i,t_end_ctrl(i)
                  write(662,"(a,1i5,1f12.6)")'i,t_end_ctrl(i)=',i,t_end_ctrl(i)
                  WRITE(97,"(a,1i5,1f12.6)")'i,t_end_ctrl(i)=',i,t_end_ctrl(i)
               ENDIF   
            ENDDO 
         ENDIF
!      
      ENDIF   
!
      RETURN    
      ENDSUBROUTINE change_tend_mars_cupid
                  