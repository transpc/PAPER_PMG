!
      SUBROUTINE tend_mars   
!
      USE Zcore ,ONLY:myrank,nrank
!
      USE TSTP_CT
!
      IMPLICIT NONE
!      
      INCLUDE 'contrl.h' !timehy
!DEC$IF defined (MCC_DLL)
      !DEC$ ATTRIBUTES DLLIMPORT:: /c3com/,/c3com_dll/
!DEC$ENDIF
      INCLUDE 'c3com.h'
!      
      INTEGER ncard,icard,i
      LOGICAL,SAVE :: initial
      DATA initial/.true./
!
      c3time_sys=timehy
!
      IF(initial)THEN 
         IF(c3tend.gt.0.0d0)THEN
            initial=.FALSE.
            IF(myrank.eq.0)THEN            
               ncard=tstp_hd%curclm(2)
               icard=1+tstp_hd%curctl(2)    
               DO i=icard,ncard
                     tstp_da(i)%tspend=c3tend
                     WRITE(*,"(a,1i5,1f12.6)")'i,tstp_da(i)%tspend=',i,tstp_da(i)%tspend  
                     write(662,"(a,1i5,1f12.6)")'i,tstp_da(i)%tspend=',i,tstp_da(i)%tspend  
                     !WRITE(97,"(a,1i5,1f12.6)")'i,tstp_da(i)%tspend=',i,tstp_da(i)%tspend  
               ENDDO   
            ENDIF
            !IF(nrank.gt.1)CALL broadcast_i(ncard,1) 
            !IF(nrank.gt.1)CALL broadcast_r(tstp_da(:)%tspend,ncard) 
         ENDIF 
      ENDIF
!
      CALL check_trip_cntrl 
!     
      RETURN
      ENDSUBROUTINE tend_mars   
 !------------------------------------------------------------------------------------------ 
      SUBROUTINE print_indta_1st_correction(iflag)  
!
      USE Zmars_index      
!      
      IMPLICIT NONE
!            
 !DEC$IF defined (MCC_DLL)
      !DEC$ ATTRIBUTES DLLIMPORT:: /c3com/,/c3com_dll/
!DEC$ENDIF
      INCLUDE 'c3com.h'
!
      INTEGER :: iflag
      REAL(8),SAVE::c3RPV3d_ctl_val_init,c3p_tmdpvol2nd_init(2),c3t_tmdpvol2nd_init(2)
!
      IF(iflag.eq.0)THEN
         c3RPV3d_ctl_val_init=c3RPV3d_ctl_val
         c3p_tmdpvol2nd_init(:)=c3p_tmdpvol2nd(:)
         c3t_tmdpvol2nd_init(:)=c3t_tmdpvol2nd(:)
         RETURN
      ENDIF   
!
      RPV3Ddp_ctl_num=RPV3Ddp_ctl_num*100+20500000
      kfactor_ctl_num=kfactor_ctl_num*100+20500000
!      
      IF(i3run_mode.eq.2.or.i3run_mode.eq.4)THEN
         OPEN(97,file='indta_1st_correction.i',status='replace')
         WRITE(97,"(a)")'*'
         WRITE(97,"(a)")'*before correction of K factor'
         WRITE(97,"(a)")'*'
         WRITE(97,"(a,1i8,a,1f12.1)")'*',RPV3Ddp_ctl_num,' dpcon constant ',c3RPV3d_ctl_val_init
         WRITE(97,"(a,1i8,a)")'*',kfactor_ctl_num,  ' kfactor function 1.0 0.0 1 0 '
         WRITE(97,"(a,1i8,a)")'*',kfactor_ctl_num+1,' cntrlvar 846 900 '
         WRITE(97,"(a)")'*'
         WRITE(97,"(a)")'*after correction of K factor'
         WRITE(97,"(a)")'*'
         WRITE(97,"(1i8,a,1f12.1)")RPV3Ddp_ctl_num,' dpcon constant ',c3RPV3d_ctl_val   
         WRITE(97,"(1i8,a,1f12.6)")kfactor_ctl_num,' kfactor constant ',c3kfactor_ctl_val         
         WRITE(97,"(a)")'*'
         CLOSE(97)  
      ENDIF
!
      number_tmdpvol2nd(:)=number_tmdpvol2nd(:)/1000000
      number_tmdpvol2nd(:)=number_tmdpvol2nd(:)*10000
      number_tmdpvol2nd(:)=number_tmdpvol2nd(:)+201
!      
      IF(i3run_mode.eq.3.or.i3run_mode.eq.4)THEN      
         OPEN(97,file='indta_2nd_correction.i',status='replace')
         WRITE(97,"(a)")'*'
         WRITE(97,"(a)")'*before correction of pressure at tmdpvol'
         WRITE(97,"(a)")'*'
         WRITE(97,"(a,1i7,a,1e15.8,1f12.2)")'*',number_tmdpvol2nd(1),'  0.0 ',c3p_tmdpvol2nd_init(2),c3t_tmdpvol2nd_init(2)         
         WRITE(97,"(a,1i7,a,1e15.8,1f12.2)")'*',number_tmdpvol2nd(2),'  0.0 ',c3p_tmdpvol2nd_init(2),c3t_tmdpvol2nd_init(2)         
         WRITE(97,"(a)")'*'
         WRITE(97,"(a)")'*after correction of pressure at tmdpvol'
         WRITE(97,"(a)")'*'
         WRITE(97,"(1i7,a,1e15.8,1f12.2)")number_tmdpvol2nd(1),'  0.0 ',c3p_tmdpvol2nd(1),c3t_tmdpvol2nd(1)
         WRITE(97,"(1i7,a,1e15.8,1f12.2)")number_tmdpvol2nd(2),'  0.0 ',c3p_tmdpvol2nd(2),c3t_tmdpvol2nd(2)         
         WRITE(97,"(a)")'*'
         CLOSE(97)     
      ENDIF
!
      RETURN
      ENDSUBROUTINE print_indta_1st_correction 
   
    
