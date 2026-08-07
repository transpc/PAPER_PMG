      SUBROUTINE check_trip_cntrl
!
      USE TRP_BLK
      USE VOL_DAT 
      USE JUN_DAT 
      USE VOL_DAT     
      USE CON_VARC    !ctl_da(:)%var_n        
!
      USE Zcore       ,ONLY: myrank,nrank
      USE Zmars_index      
!      
      IMPLICIT NONE  
!
!DEC$IF defined (MCC_DLL)
      !DEC$ ATTRIBUTES DLLIMPORT:: /c3com/,/c3com_dll/
!DEC$ENDIF      

      INCLUDE 'c3com.h' !c3marsin,c3cplmaster
!
      INTEGER i,j,k
      LOGICAL,SAVE:: initial,initial_vv,initial_tmdp
      DATA initial,initial_vv,initial_tmdp /.true.,.true.,.true./
!
      IF(initial)THEN
         initial=.false.
!!!         IF(myrank.eq.0)CALL print_logo
         IF(myrank.eq.0)CALL get_mars_index    
      ENDIF !initial  
!
!.....Rx trip
!
      IF(i3marsin.ne.0)RETURN
!         
      IF(initial_vv)THEN
         initial_vv=.false.
                  
         rx_trip_mars=0
         rcp_trip_mars=0
         mslb_trip_mars=0
         degba_trip_mars=0
         rod_trip_mars=0
         
         idx_rx_trip_mars=0
         idx_rcp_trip_mars=0
         idx_degba_trip_mars=0
         idx_mslb_trip_mars=0
         idx_rktpow_ctl_num=0
         idx_kfactor_ctl_num=0
         idx_RPV3Ddp_ctl_num=0
         
         c3RPV3d_ctl_val=-1.0d0
!            
!........check control             
         IF(myrank.eq.0)THEN
            DO i=1,ctl_hd%nCompnt(2) 
               IF(ctl_da(i)%CompntNum(2).eq.rktpow_ctl_num)THEN
                  idx_rktpow_ctl_num=i
                  IF(myrank.eq.0)WRITE(*,"(11x,a,2i5)")'power ratio control number, index=',ctl_da(i)%CompntNum(2),idx_rktpow_ctl_num
                  IF(myrank.eq.0)write(662,"(11x,a,2i5)")'power ratio control number, index=',ctl_da(i)%CompntNum(2),idx_rktpow_ctl_num
               ENDIF
               IF(ctl_da(i)%CompntNum(2).eq.kfactor_ctl_num)THEN
                  idx_kfactor_ctl_num=i
                  IF(myrank.eq.0)WRITE(*,"(11x,a,2i5)")'K facotr control number, index=',ctl_da(i)%CompntNum(2),idx_kfactor_ctl_num
                  IF(myrank.eq.0)write(662,"(11x,a,2i5)")'K factor control number, index=',ctl_da(i)%CompntNum(2),idx_kfactor_ctl_num
               ENDIF
               IF(ctl_da(i)%CompntNum(2).eq.RPV3Ddp_ctl_num)THEN
                  idx_RPV3Ddp_ctl_num=i
                  IF(myrank.eq.0)WRITE(*,"(11x,a,2i5)")'RPV3D dpcon control number, index=',ctl_da(i)%CompntNum(2),idx_RPV3Ddp_ctl_num
                  IF(myrank.eq.0)write(662,"(11x,a,2i5)")'RPV3D dpcon control number, index=',ctl_da(i)%CompntNum(2),idx_RPV3Ddp_ctl_num
               ENDIF
            ENDDO   
         ENDIF
         IF(nrank.gt.1)CALL broadcast_i(idx_rktpow_ctl_num,1)              
         IF(nrank.gt.1)CALL broadcast_i(idx_kfactor_ctl_num,1)   
         IF(nrank.gt.1)CALL broadcast_i(idx_RPV3Ddp_ctl_num,1)            
!           
!........check trip    
         IF(myrank.eq.0)THEN                 
            DO i=1,trp_hd%nVarTrip(2)+trp_hd%nLogicTrip(2) 
               IF(trp_da(i)%number(2).eq.rx_trip_num)THEN !reactor trip
                  idx_rx_trip_mars=i
                  IF(myrank.eq.0)WRITE(*,"(11x,a,2i5)")'reactor trip number, index=',trp_da(i)%number(2),idx_rx_trip_mars
                  IF(myrank.eq.0)write(662,"(11x,a,2i5)")'reactor trip number, index=',trp_da(i)%number(2),idx_rx_trip_mars
               ENDIF  
               IF(trp_da(i)%number(2).eq.rod_trip_num)THEN !rod trip
                  idx_rod_trip_mars=i
                  IF(myrank.eq.0)WRITE(*,"(11x,a,2i5)")'rod trip number, index=',trp_da(i)%number(2),idx_rod_trip_mars
                  IF(myrank.eq.0)write(662,"(11x,a,2i5)")'rod trip number, index=',trp_da(i)%number(2),idx_rod_trip_mars
               ENDIF  
               IF(trp_da(i)%number(2).eq.rcp_trip_num)THEN !rcp trip
                  idx_rcp_trip_mars=i
                  IF(myrank.eq.0)WRITE(*,"(11x,a,2i5)")'rcp trip number, index=',trp_da(i)%number(2),idx_rcp_trip_mars
                  IF(myrank.eq.0)write(662,"(11x,a,2i5)")'rcp trip number, index=',trp_da(i)%number(2),idx_rcp_trip_mars
               ENDIF  
               IF(trp_da(i)%number(2).eq.degba_trip_num)THEN !degba trip
                  idx_degba_trip_mars=i
                  IF(myrank.eq.0)WRITE(*,"(11x,a,2i5)")'degba trip number, index=',trp_da(i)%number(2),idx_degba_trip_mars
                  IF(myrank.eq.0)write(662,"(11x,a,2i5)")'degba trip number, index=',trp_da(i)%number(2),idx_degba_trip_mars
               ENDIF  
               IF(trp_da(i)%number(2).eq.mslb_trip_num)THEN !mslb trip
                  idx_mslb_trip_mars=i
                  IF(myrank.eq.0)WRITE(*,"(11x,a,2i5)")'mslb trip number, index=',trp_da(i)%number(2),idx_mslb_trip_mars
                  IF(myrank.eq.0)write(662,"(11x,a,2i5)")'mslb trip number, index=',trp_da(i)%number(2),idx_mslb_trip_mars
               ENDIF  
!!!               IF(trp_da(i)%number(2).eq.event_trip_num)THEN !event trip
!!!                  idx_event_trip_mars=i
!!!                  IF(myrank.eq.0)WRITE(*,"(11x,a,2i5)")'mslb trip number, index=',trp_da(i)%number(2),idx_event_trip_mars
!!!                  IF(myrank.eq.0)write(662,"(11x,a,2i5)")'mslb trip number, index=',trp_da(i)%number(2),idx_event_trip_mars
!!!               ENDIF 
!!!               i=idx_event_trip_mars 
!!!               trp_da(i)%tripcon=INT(c3time_sys+event_delay_time+1.0d0)
            ENDDO 
         ENDIF
         IF(nrank.gt.1)CALL broadcast_i(idx_rx_trip_mars,1)
         IF(nrank.gt.1)CALL broadcast_i(idx_rod_trip_mars,1)
         IF(nrank.gt.1)CALL broadcast_i(idx_rcp_trip_mars,1)
         IF(nrank.gt.1)CALL broadcast_i(idx_degba_trip_mars,1)
         IF(nrank.gt.1)CALL broadcast_i(idx_mslb_trip_mars,1)
!            
         IF(idx_rx_trip_mars.eq.0)THEN
             IF(myrank.eq.0)WRITE(*,"(11x,a)")'## lack of reactor trip number in print_mars!'            
             IF(myrank.eq.0)write(662,"(11x,a)")'## lack of reactor trip number in print_mars!'            
!!!                PAUSE
         ENDIF
         IF(idx_rod_trip_mars.eq.0)THEN
             IF(myrank.eq.0)WRITE(*,"(11x,a)")'## lack of rod trip number in print_mars!'            
             IF(myrank.eq.0)write(662,"(11x,a)")'## lack of rod trip number in print_mars!'            
!!!                PAUSE
         ENDIF
         IF(idx_rktpow_ctl_num.eq.0)THEN
             IF(myrank.eq.0)WRITE(*,"(11x,a)")'## lack of reactor power control number in print_mars!'            
             IF(myrank.eq.0)write(662,"(11x,a)")'## lack of reactor power control number in print_mars!'            
!!!                PAUSE
         ENDIF
         IF(idx_kfactor_ctl_num.eq.0)THEN
             IF(myrank.eq.0)WRITE(*,"(11x,a)")'## lack of K factor control number in print_mars!'            
             IF(myrank.eq.0)write(662,"(11x,a)")'## lack of K factor control number in print_mars!'            
!!!                PAUSE
         ENDIF
         IF(idx_RPV3Ddp_ctl_num.eq.0)THEN
             IF(myrank.eq.0)WRITE(*,"(11x,a)")'## lack of K factor control number in print_mars!'            
             IF(myrank.eq.0)write(662,"(11x,a)")'## lack of K factor control number in print_mars!'            
             c3RPV3d_ctl_val=0.0d0
!!!                PAUSE
         ELSE    
            IF(myrank.eq.0)c3RPV3d_ctl_val=ctl_da(idx_RPV3Ddp_ctl_num)%Var_n
         ENDIF         
!          
      ENDIF !initial_vv   
!
!.....assign mass flow rates and pressures at legs
!
      IF(myrank.eq.0)THEN
         DO k=1,i3n_junleg
            c3mflow_junleg(k)=j_da(index_junleg(k))%mflowj
         ENDDO  
         DO k=1,i3n_volleg
             c3p_volleg(k)=v_da(index_volleg(k))%p
             c3tl_volleg(k)=v_da(index_volleg(k))%tf
         ENDDO 
      ENDIF         
!!!IF(nrank.gt.1)CALL broadcast_r(mflow_junleg,n_junleg)
!!!IF(nrank.gt.1)CALL broadcast_r(p_volleg,n_volleg)
!!!IF(nrank.gt.1)CALL broadcast_r(tl_volleg,n_volleg)
!
!........put the value of a control variable
!
      IF(myrank.eq.0)THEN
         IF(idx_rktpow_ctl_num.gt.0)THEN 
            IF(i3cplmaster.gt.0)ctl_da(idx_rktpow_ctl_num)%Var_n=c3rktpow_ctl_val
         ENDIF
         IF(idx_kfactor_ctl_num.gt.0)THEN   
            c3kfactor_ctl_val=ctl_da(idx_kfactor_ctl_num)%Var_n
         ELSE
            c3kfactor_ctl_val=-1.0d0
         ENDIF
         IF(idx_RPV3Ddp_ctl_num.gt.0 .and. c3RPV3d_ctl_val.gt.0.0d0)THEN   
            ctl_da(idx_RPV3Ddp_ctl_num)%Var_n=c3RPV3d_ctl_val    
         ENDIF
      ENDIF   
      IF(myrank.eq.0)THEN
         IF(initial_tmdp)THEN
            initial_tmdp=.FALSE.
            DO k=1,i3n_tmdpvol2nd
                c3p_tmdpvol2nd(k)=v_da(index_tmdpvol2nd(k))%p
                c3t_tmdpvol2nd(k)=v_da(index_tmdpvol2nd(k))%tf
            ENDDO
            CALL print_indta_1st_correction(0)
         ELSE   
            DO k=1,i3n_tmdpvol2nd
                v_da(index_tmdpvol2nd(k))%p=c3p_tmdpvol2nd(k)
            ENDDO 
         ENDIF
      ENDIF   
      IF(nrank.gt.1)CALL broadcast_r(c3kfactor_ctl_val,1) !mcc-mpi         
!
!.....get trip         
!
      IF(myrank.eq.0)THEN
         IF(rx_trip_mars.lt.1.and.idx_rx_trip_mars.gt.0)THEN
            IF(myrank.eq.0.and.trp_da(idx_rx_trip_mars)%TripTime.ne.-1.0d0)then
               rx_trip_mars=1
               WRITE(*,"(11x,a,2e12.5,1i5)")'## time,Rx trip=',c3time_sys,trp_da(idx_rx_trip_mars)%TripTime,trp_da(idx_rx_trip_mars)%number(2)
               write(662,"(11x,a,2e12.5,1i5)")'## time,Rx trip=',c3time_sys,trp_da(idx_rx_trip_mars)%TripTime,trp_da(idx_rx_trip_mars)%number(2)  
            ENDIF
         ENDIF   
         IF(rod_trip_mars.lt.1.and.idx_rod_trip_mars.gt.0)THEN
            IF(myrank.eq.0.and.trp_da(idx_rod_trip_mars)%TripTime.ne.-1.0d0)then
               rod_trip_mars=1
               WRITE(*,"(11x,a,2e12.5,1i5)")'## time,Rod trip=',c3time_sys,trp_da(idx_rod_trip_mars)%TripTime,trp_da(idx_rod_trip_mars)%number(2)
               write(662,"(11x,a,2e12.5,1i5)")'## time,Rod trip=',c3time_sys,trp_da(idx_rod_trip_mars)%TripTime,trp_da(idx_rod_trip_mars)%number(2)  
            ENDIF
         ENDIF
         IF(rcp_trip_mars.lt.1.and.idx_rcp_trip_mars.gt.0)THEN
            IF(myrank.eq.0.and.trp_da(idx_rcp_trip_mars)%TripTime.ne.-1.0d0)then
               rcp_trip_mars=1
               WRITE(*,"(11x,a,2e12.5,1i5)")'## time,RCP trip=',c3time_sys,trp_da(idx_rcp_trip_mars)%TripTime,trp_da(idx_rcp_trip_mars)%number(2)
               write(662,"(11x,a,2e12.5,1i5)")'## time,RCP trip=',c3time_sys,trp_da(idx_rcp_trip_mars)%TripTime,trp_da(idx_rcp_trip_mars)%number(2)  
             ENDIF
         ENDIF 
         IF(degba_trip_mars.lt.1.and.idx_degba_trip_mars.gt.0)THEN
            IF(myrank.eq.0.and.trp_da(idx_degba_trip_mars)%TripTime.ne.-1.0d0)then
               degba_trip_mars=1
               WRITE(*,"(11x,a,2e12.5,1i5)")'## time, degba trip=',c3time_sys,trp_da(idx_degba_trip_mars)%TripTime,trp_da(idx_degba_trip_mars)%number(2)
               write(662,"(11x,a,2e12.5,1i5)")'## time, degba trip=',c3time_sys,trp_da(idx_degba_trip_mars)%TripTime,trp_da(idx_degba_trip_mars)%number(2)  
            ENDIF
         ENDIF   
         IF(mslb_trip_mars.lt.1.and.idx_mslb_trip_mars.gt.0)THEN
            IF(myrank.eq.0.and.trp_da(idx_mslb_trip_mars)%TripTime.ne.-1.0d0)then
               mslb_trip_mars=1
               WRITE(*,"(11x,a,2e12.5,1i5)")'## time, mslb trip=',c3time_sys,trp_da(idx_mslb_trip_mars)%TripTime,trp_da(idx_mslb_trip_mars)%number(2)
               write(662,"(11x,a,2e12.5,1i5)")'## time, mslb trip=',c3time_sys,trp_da(idx_mslb_trip_mars)%TripTime,trp_da(idx_mslb_trip_mars)%number(2)  
            ENDIF
         ENDIF  
      ENDIF
      IF(nrank.gt.1)CALL broadcast_i(rx_trip_mars,1)
      IF(nrank.gt.1)CALL broadcast_i(rod_trip_mars,1)
      IF(nrank.gt.1)CALL broadcast_i(rcp_trip_mars,1)
      IF(nrank.gt.1)CALL broadcast_i(degba_trip_mars,1)
      IF(nrank.gt.1)CALL broadcast_i(mslb_trip_mars,1)
!
         i3rx_trip=rx_trip_mars
         i3rcp_trip=rcp_trip_mars
         i3mslb_trip=mslb_trip_mars
         i3degba_trip=degba_trip_mars
         i3rod_trip=rod_trip_mars
!           
      RETURN
      ENDSUBROUTINE check_trip_cntrl
      
           