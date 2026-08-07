!
      MODULE Zmars_index
      IMPLICIT NONE
      SAVE
      INTEGER,ALLOCATABLE::mars_vol_index(:),mars_jun_index(:),  &
                             mars_vol_jun_from(:),mars_vol_jun_to(:) 
	   INTEGER ntrip_mars,nctl_mars
	   INTEGER rx_trip_mars,rcp_trip_mars,degba_trip_mars,mslb_trip_mars,rod_trip_mars
	   INTEGER idx_rx_trip_mars,idx_rcp_trip_mars,idx_degba_trip_mars,idx_mslb_trip_mars,idx_rod_trip_mars
	   INTEGER rx_trip_num,rcp_trip_num,mslb_trip_num,degba_trip_num,rod_trip_num   
	   INTEGER rktpow_ctl_num,idx_rktpow_ctl_num,kfactor_ctl_num,idx_kfactor_ctl_num
      INTEGER RPV3Ddp_ctl_num,idx_RPV3Ddp_ctl_num
	   INTEGER,ALLOCATABLE:: number_trip(:),number_ctl(:)                             
    	INTEGER,ALLOCATABLE:: number_volleg(:),index_volleg(:)
   	INTEGER,ALLOCATABLE:: number_junleg(:),index_junleg(:)
      INTEGER,ALLOCATABLE:: number_tmdpvol2nd(:),index_tmdpvol2nd(:)    
!
      ENDMODULE Zmars_index
!------------------------------------------------------------------------------------  
!
      SUBROUTINE read_cupid_mars
!
      USE Zcore        ,ONLY: myrank
      USE Zmars_index      
!      
      IMPLICIT NONE
!     
!DEC$IF defined (MCC_DLL)
      !DEC$ ATTRIBUTES DLLIMPORT:: /c3com/,/c3com_dll/
!DEC$ENDIF      
      
      INCLUDE 'c3com.h' !c3marin
!
      INTEGER i,nout    
!
      rx_trip_mars=0
      rcp_trip_mars=0
      degba_trip_mars=0  
!                  
      OPEN(812,file='cupid_mars.in',status='old',iostat=i3marsin) 
      IF(i3marsin.eq.0)THEN
          READ(812,*)ntrip_mars
          ALLOCATE(number_trip(ntrip_mars))
          IF(ntrip_mars.gt.10)THEN
             IF(myrank.eq.0)WRITE(*,"(a)")'ntrip_mars should be less 10 !'
             PAUSE
             STOP
          ENDIF
          DO i=1,ntrip_mars
             READ(812,*)number_trip(i)
          ENDDO         
          rod_trip_num=number_trip(1)
          rx_trip_num=number_trip(2)
          rcp_trip_num=number_trip(3)
          mslb_trip_num=number_trip(4)
          degba_trip_num=number_trip(5)
          READ(812,*)nctl_mars
          ALLOCATE(number_ctl(nctl_mars))
          IF(nctl_mars.gt.10)THEN
             IF(myrank.eq.0)WRITE(*,"(a)")'number_ctl should be less 10 !'
             PAUSE
             STOP
          ENDIF
          DO i=1,nctl_mars
             READ(812,*)number_ctl(i)
          ENDDO            
          rktpow_ctl_num=number_ctl(1)
          kfactor_ctl_num=number_ctl(2)
          RPV3Ddp_ctl_num=number_ctl(3)
          c3rktpow_ctl_val=1.0d0
          c3kfactor_ctl_val=0.0d0
          READ(812,*)i3n_junleg
          IF(i3n_junleg.gt.10)THEN
             IF(myrank.eq.0)WRITE(*,"(a)")'n_junleg should be less than 10 !'
             PAUSE
             STOP
          ENDIF
          ALLOCATE(number_junleg(i3n_junleg*2))
          ALLOCATE(index_junleg(i3n_junleg*2))
!          ALLOCATE(c3mflow_junleg(i3n_junleg*2))
          DO i=1,i3n_junleg
             READ(812,*)number_junleg(i),number_junleg(i+i3n_junleg)
          ENDDO  
          IF(i3cplmars.ge.2)THEN
             DO i=1,i3n_junleg
                number_junleg(i)=number_junleg(i+i3n_junleg)
             ENDDO   
          ENDIF 
!          
          READ(812,*)i3n_volleg    
          IF(i3n_volleg.gt.10)THEN
             IF(myrank.eq.0)WRITE(*,"(a)")'n_volleg should be less than 10 !'
             PAUSE
             STOP
          ENDIF                
          ALLOCATE(number_volleg(i3n_volleg*2))
          ALLOCATE(index_volleg(i3n_volleg*2))
!          ALLOCATE(c3p_volleg(i3n_volleg*2),c3tl_volleg(i3n_volleg*2))
          DO i=1,i3n_volleg
             READ(812,*)number_volleg(i),number_volleg(i+i3n_volleg)
          ENDDO            
          IF(i3cplmars.le.1)THEN
             DO i=1,i3n_volleg
                number_volleg(i+i3n_volleg)=number_volleg(i)
             ENDDO   
          ENDIF           
          i3n_volleg=i3n_volleg*2
!
          IF(i3cplmars.eq.2)THEN
             READ(812,*)i3n_tmdpvol2nd    
             IF(i3n_tmdpvol2nd.gt.2)THEN
                IF(myrank.eq.0)WRITE(*,"(a)")'n_tmdpvol2nd should be less than 2 !'
                PAUSE
                STOP
             ENDIF                
             ALLOCATE(number_tmdpvol2nd(i3n_tmdpvol2nd))
             ALLOCATE(index_tmdpvol2nd(i3n_tmdpvol2nd))
             DO i=1,i3n_tmdpvol2nd
                READ(812,*)number_tmdpvol2nd(i)
             ENDDO  
             READ(812,*)c3run_mode3_dur,c3run_mode3_dpcri
             READ(812,*)c3run_mode2_dur,c3run_mode2_dpcri
          ELSE
             c3run_mode3_dur=0.0d0
             i3n_tmdpvol2nd=0
          ENDIF
!          
          CLOSE(812)
      ELSE
         IF(myrank.eq.0)WRITE(*,"(11x,a)")'No trip information due to lack of <cupid_mars.in>.'
         IF(myrank.eq.0)WRITE(97,"(11x,a)")'No trip information due to lack of <cupid_mars.in>.' 
      ENDIF
!      
      IF(i3cplmars.gt.0)THEN
         IF(myrank.eq.0)WRITE(*,"(11x,a,  a,a)")'****************','***',  '*********' 
         IF(myrank.eq.0)WRITE(*,"(11x,a,  a,a)")'***             ','   ',  '      ***' 
         IF(myrank.eq.0)WRITE(*,"(11x,a,1i3,a)")'*** CUPID_MARS= ',i3cplmars,' STEP ***' 
         IF(myrank.eq.0)WRITE(*,"(11x,a,  a,a)")'***             ','   ',  '      ***' 
         IF(myrank.eq.0)WRITE(*,"(11x,a,  a,a)")'****************','***',  '*********'
       
         IF(myrank.eq.0)WRITE(97,"(11x,a,  a,a)")'****************','***',  '*********' 
         IF(myrank.eq.0)WRITE(97,"(11x,a,  a,a)")'***             ','   ',  '      ***' 
         IF(myrank.eq.0)WRITE(97,"(11x,a,1i3,a)")'*** CUPID_MARS= ',i3cplmars,' STEP ***' 
         IF(myrank.eq.0)WRITE(97,"(11x,a,  a,a)")'***             ','   ',  '      ***' 
         IF(myrank.eq.0)WRITE(97,"(11x,a,  a,a)")'****************','***',  '*********'      
      ELSEIF(i3cplmars.eq.0)THEN
         IF(myrank.eq.0)WRITE(*,"(11x,a,  a,a)")'****************','***',  '*********' 
         IF(myrank.eq.0)WRITE(*,"(11x,a,  a,a)")'***             ','   ',  '      ***' 
         IF(myrank.eq.0)WRITE(*,"(11x,a,1i3,a)")'***       MARS= ',i3cplmars+1,' STEP ***' 
         IF(myrank.eq.0)WRITE(*,"(11x,a,  a,a)")'***             ','   ',  '      ***' 
         IF(myrank.eq.0)WRITE(*,"(11x,a,  a,a)")'****************','***',  '*********'
       
         IF(myrank.eq.0)WRITE(97,"(11x,a,  a,a)")'****************','***',  '*********' 
         IF(myrank.eq.0)WRITE(97,"(11x,a,  a,a)")'***             ','   ',  '      ***' 
         IF(myrank.eq.0)WRITE(97,"(11x,a,1i3,a)")'***       MARS= ',i3cplmars+1,' STEP ***' 
         IF(myrank.eq.0)WRITE(97,"(11x,a,  a,a)")'***             ','   ',  '      ***' 
         IF(myrank.eq.0)WRITE(97,"(11x,a,  a,a)")'****************','***',  '*********'           
      ENDIF   
!      
      RETURN
      ENDSUBROUTINE read_cupid_mars
     
