!
      SUBROUTINE CM_Control_Reactivity
!    
!.....Control reactivity 
!      
!PPM : 노심 전체의 평균 보론 농도 (PPM)
!NCB : 제어봉의 종류 개수 
!CBNAM : 제어봉의 이름 (MASTER4 입력에서 확인할 수 있다)
!ZCB : 각 제어봉의 위치 (노심 하부 반사체로부터의 높이, cm)
!
! v(t)=206.62t-150.33t2+61.661t3-11.55t4
! V(t)=381.0cm -[ 103.31t2-50.11t3+15.415t4-2.31t5 ] >30.0cm
!
      USE MASTER4      ,ONLY: ncb_mas,cbnam_mas,zcb_mas
      USE MASTER4      ,ONLY: ncb,cbnam,zcb
      USE Zmars        ,ONLY: time_scram 
      USE Ztimecon     ,ONLY: time
      USE Zconst1      ,ONLY: vv_prob
      USE Zcore        ,ONLY: myrank
      USE MASTER4      ,ONLY: ZCB,NCB,CBNAM,mas_rx_trip
      USE Zconst1      ,ONLY: restart
      USE Zuserdefined ,ONLY: user_iary
      USE Zio_unit     , ONLY: unit_log
!     
      IMPLICIT NONE
!
      INCLUDE '../10_LinkToMARS/c3com.h' 
!     
      INTEGER i
!
      LOGICAL,SAVE :: initial=.TRUE.
      LOGICAL,SAVE :: initial_scram=.TRUE.
!
!.....Save reactivity at first
!     
      IF(initial)THEN
         initial=.FALSE.
!
!........CEA information (should be EXACTLY matched with MASTER input, MAS_INP)
!
         NCB=11
         ALLOCATE(ZCB(NCB),CBNAM(NCB))
         CBNAM=(/"R1  ", &
                 "S   ", &
                 "R2  ", &
                 "R3  ", &
                 "B   ", &
                 "B12 ", &
                 "R4  ", &
                 "P   ", &
                 "R5  ", &
                 "A   ", &
                 "B11 "   /)
         ZCB=381.0d0
!
!........Assign into local array (ncb_mas,cbnam_mas,zcb_mas)
!
         ncb_mas=NCB
         ALLOCATE(cbnam_mas(ncb_mas))
         ALLOCATE(zcb_mas(ncb_mas))
         DO i=1,ncb_mas
            cbnam_mas(i)=TRIM(CBNAM(I))
            zcb_mas(i)=ZCB(I)
         ENDDO
         IF(restart.ne.0)mas_rx_trip=user_iary(31)         
      ENDIF
!
!.....Setpoint criteria (121%)
!      
     IF(mas_rx_trip.eq.1)i3rod_trip=1
     user_iary(31)=mas_rx_trip
!     
     IF(i3rod_trip.eq.0)THEN
        IF(vv_prob.eq.'opr1000_rv')THEN
           IF(c3rktpow_ctl_val.gt.121.0d0)THEN
              i3rod_trip=1
              IF(myrank.eq.0)WRITE(*,"(11x,a,1i3,1f12.6)")'========> i3rod_trip at time=',i3rod_trip,time
              IF(myrank.eq.0)WRITE(unit_log,"(11x,a,1i3,1f12.6)")'========> i3rod_trip at time=',i3rod_trip,time
           ENDIF   
        ENDIF   
        IF(vv_prob.eq.'apr1400_rv')THEN
           IF(c3rktpow_ctl_val.gt.121.0d0)THEN
              i3rod_trip=1
              IF(myrank.eq.0)WRITE(*,"(11x,a,1i3,1f12.6)")'========> i3rod_trip at time=',i3rod_trip,time
              IF(myrank.eq.0)WRITE(unit_log,"(11x,a,1i3,1f12.6)")'========> i3rod_trip at time=',i3rod_trip,time
           ENDIF   
        ENDIF   
     ENDIF   
!
!....Get reactor trip signal
!
!     rx_trip_signal=reactor_trip_function(rx_trip_signal)
!     IF(rx_trip_sginal.eq.1)rod_trip_mars=1
!
!....Define scram
!
      IF(i3rod_trip.eq.0)THEN
!     
         RETURN
!         
      ELSE 
!             
!.......set reactor trip time of master
         IF(initial_scram.eq..true.)THEN
            initial_scram=.FALSE.
            time_scram=time
         ENDIF
!
!........set control rods as scram
         CALL udfn_scram_actuation
!        
      ENDIF    
!         
      END SUBROUTINE CM_Control_Reactivity
!      
!------------------------------------------------------------------------------------
!
      SUBROUTINE udfn_scram_actuation !should be ch
!    
!.....scramming for a time, all rods to control or to shut down reactor are inserted
!      
      USE Zcore        , ONLY: myrank
      USE Ztimecon     , ONLY: time
      USE master4      , ONLY: ncb_mas,zcb_mas,time_cri
      USE Zconst1      , ONLY: restart
      USE Zuserdefined , ONLY: user_rary     
      USE Zio_unit     , ONLY: unit_log
!     
      IMPLICIT NONE
!     
      INTEGER i  
!
      REAL(8) cm_cri_distance
!
      LOGICAL,SAVE :: initial   =.TRUE.
      REAL(8),SAVE :: cri_loc=381.d0
      REAL(8),SAVE :: cri_loc_init=381.d0
      REAL(8),SAVE :: cri_loc_min=30.d0
!
      IF(initial)THEN
!     
         initial=.FALSE.
         time_cri=-1.0d0
         IF(restart.ne.0)time_cri=user_rary(31)
         IF(time_cri.lt.0.0d0)THEN
            IF(myrank.eq.0)WRITE(*,*)'1time_cri,user_rary(31)=',time_cri,user_rary(31)
            time_cri=time
            user_rary(31)=time
            IF(myrank.eq.0)WRITE(*,*)'2time_cri,user_rary(31)=',time_cri,user_rary(31)
         ENDIF   
         cri_loc_init=MAXVAL(zcb_mas(1:ncb_mas))
         cri_loc_min=0.0d0
         IF(myrank.eq.0)THEN
            WRITE(*,"(11x,a)")'==============================================='
            WRITE(*,"(11x,a,1f8.1,a)")'==reactor scram of MASTER4 at ',time_cri,'s'
            WRITE(*,"(11x,a)")'==============================================='
            WRITE(unit_log,"(11x,a)")'==============================================='
            WRITE(unit_log,"(11x,a,1f8.1,a)")'==reactor scram of MASTER4 at ',time_cri,'s'
            WRITE(unit_log,"(11x,a)")'==============================================='
         ENDIF
!
     ENDIF    
!
         IF(cri_loc.gt.cri_loc_min)THEN
            cri_loc=cri_loc_init-cm_cri_distance(time-time_cri)         
            cri_loc=DMAX1(cri_loc_min,cri_loc)
            DO i=1,ncb_mas
               IF(zcb_mas(i).gt.cri_loc)THEN
                  zcb_mas(i)=cri_loc
               ENDIF
            ENDDO
            !IF(myrank.eq.0.and.MOD(itim,10).eq.0)THEN
            IF(myrank.eq.0)THEN
               WRITE(*,"(11x,a,1f9.3,2f6.1)")'==time(s),CR Loc.(cm),CR D.(cm)=',time-time_cri,cri_loc,cm_cri_distance(time-time_cri)  
               WRITE(unit_log,"(11x,a,1f9.3,2f6.1)")'==time(s),CR Loc.(cm),CR D.(cm)=',time-time_cri,cri_loc,cm_cri_distance(time-time_cri)  
            ENDIF
         ENDIF  
!     
      RETURN   
      END SUBROUTINE udfn_scram_actuation
!      
!------------------------------------------------------------------------------------
!     
      FUNCTION cm_cri_distance(time_cri)     
!     
      IMPLICIT NONE
      INTEGER i, i1
      REAL(8) cm_cri_distance
      REAL(8) time_cri
      REAL(8) cm_cri_time(15),cm_cri_position(15)
      DATA cm_cri_time/0.0,0.2,0.4,0.6,0.8,&
                       1.0,1.2,1.4,1.6,1.8,&
                       2.0,2.2,2.4,2.6,1.d10/
      DATA cm_cri_position/0.00,3.81,14.44,50.52,236.68,&
                           304.50,338.29,355.01,364.27,370.34,&
                           375.21,377.19,379.10,381.00,381.00/ 
!     DATA cm_cri_position/0.0   ,  3.5 , 13.3 ,46.5 ,218.0,& 
!                           280.5 ,311.7 ,327.1 ,335.6 ,341.2,& 
!                           345.7 ,347.5 ,349.2 ,351.0 ,351.0/ 
      DO i=1,14
         IF(time_cri.ge.cm_cri_time(i).and.time_cri.lt.cm_cri_time(i+1))THEN
            i1=i
         ENDIF
      ENDDO  
      i=i1
      cm_cri_distance=(cm_cri_position(i+1)-cm_cri_position(i))/(cm_cri_time(i+1)-cm_cri_time(i))*(time_cri-cm_cri_time(i))+cm_cri_position(i)
!     
      END FUNCTION cm_cri_distance
