      SUBROUTINE udfn_set_qvol_liq_atlas
!
!     volumetric heat source in porous cells
!
      USE VOL_DATA     
      USE Zmpi      , ONLY: ncell_fp
      USE Zzone     , ONLY: ncell_fluid
      USE Zcore     , ONLY: np,myrank
      USE Ztimecon  , ONLY: time,itim
      USE Zparam    , ONLY: pi         
      USE Zconst1   , ONLY: vv_prob,cplmars
      USE Zcoord1   , ONLY: xloc
      USE Zcoord3   , ONLY: volp    
      USE Zqvol     , ONLY: qvol_liq
      USE Zmars     , ONLY: rx_trip_signal, time_scram
      USE Zio_unit  , ONLY: unit_log
      
!
      IMPLICIT NONE
!      
      INTEGER(4) i     
      INTEGER ::reactor_trip_function
!
      LOGICAL, SAVE :: initial_atlas,initial
      DATA initial_atlas,initial/.true.,.true./           
!      
      REAL(8) :: x,y,z,r,hfactor,heat_factor,time_after_scram
      REAL(8) :: qvolsum,ql
      REAL(8), SAVE :: corevol
      REAL(8) :: zh,zl,zc,zi,zr,rdc,tpower            
      REAL(8),ALLOCATABLE,SAVE :: zfactor(:)      


     
      REAL(8) :: tscram(15),tfactor
      DATA tscram/0,5,10,20,40,50,100, 200,300,400, &
                  500,1000,2500,5005,10000/
      REAL(8), SAVE :: time_ramp

      REAL(8),ALLOCATABLE,SAVE::qvol_liq_init(:)
!      
      IF(initial) THEN
         ALLOCATE(qvol_liq_init(ncell_fp))
         initial=.false.
      ENDIF
      IF(vv_prob.eq.'atlas_mc_porous')THEN
          zh=-0.771d0
          zl=-2.676d0
          zc=(zh+zl)/2.0d0
          zi=(zh-zl)/11         
          rdc=0.168d0
          !tpower=1.566e6
          tpower=1.633777e6 !kdh-consider heat loss 
      ELSEIF(vv_prob.eq.'pwr_mc_poro')THEN
          zh=-1.1d0    !OLD
          zl=-4.9d0
          zc=(zh+zl)/2.0d0
          zi=(zh-zl)/11         
          rdc=1.994d0
          tpower=3.983e9 !rktpow in mars input !apr1400=4000MW, opr1400=3000MW      
      ELSEIF(vv_prob.eq.'apr1400_mc_poro')THEN
          zh=-5.4d0     !NEW 
          zl=-9.22d0
          zc=(zh+zl)/2.0d0
          zi=(zh-zl)/11         
          rdc=1.77d0
          tpower=3.983e9          
!      ELSEIF(vv_prob.eq.'opr1000_mc')THEN
!          zh=4.53d0 !3.81
!          zl=0.72d0
!          rdc=1.60d0
!          zc=(zh+zl)/2.0d0
!          zi=(zh-zl)/11             
!          tpower=2806.03e6 !rktpow in mars input !apr1400=4000MW, opr1400=3000MW          
      ENDIF
!                             
     IF(initial_atlas.eq..true.)then
         time_ramp=time           
         rx_trip_signal=0
         !####
         IF(time.gt.5110.0d0)then
           IF(myrank.eq.0)THEN
              WRITE(unit_log,*)'########## rx_trip_signal and time_scram are set manually!'
              WRITE(*,*)'########## rx_trip_signal and time_scram are set manually!'
           ENDIF   
           rx_trip_signal=1
           time_scram=5110.0d0
         ENDIF
         !####  
         corevol=0.0d0
         DO i=1,ncell_fluid
            x=xloc(i,1)
            y=xloc(i,2)
            z=xloc(i,3)
            r=DSQRT(x*x+y*y)
            IF(r.le.rdc.and.z.le.zh.and.z.ge.zl)then
               corevol=corevol+volp(i)
            ENDIF
         ENDDO
         IF(np.gt.1)then
            CALL allreducei_r1(corevol)
         ENDIF   
         qvolsum=0.0d0
         ALLOCATE(zfactor(ncell_fluid))
         zfactor(:)=0.0d0  
         DO i=1,ncell_fluid
            x=xloc(i,1)
            y=xloc(i,2)
            z=xloc(i,3)
            r=DSQRT(x*x+y*y)
            zr=DABS(z-zc)
            IF(zr.le.zi/2.d0)THEN
               zfactor(i)=1.466d0
            ELSEIF(zr.le.zi/2.d0+zi)THEN
               zfactor(i)=1.414d0   
            ELSEIF(zr.le.zi/2.d0+2*zi)THEN
               zfactor(i)=1.264d0
            ELSEIF(zr.le.zi/2.d0+3*zi)THEN
               zfactor(i)=1.024d0
            ELSEIF(zr.le.zi/2.d0+4*zi)THEN
               zfactor(i)=0.713d0
            ELSEIF(zr.le.zi/2.d0+5*zi)THEN
               zfactor(i)=0.352d0
            ELSE
               zfactor(i)=0.0d0
            ENDIF             
            qvol_liq(i)=0.0d0
            qvol_liq_init(i)=0.0d0
            IF(r.le.rdc.and.z.le.zh.and.z.ge.zl)then
               qvol_liq(i)=tpower/corevol            
               qvolsum=qvolsum+qvol_liq(i)*volp(i)*zfactor(i)
               qvol_liq_init(i)=qvol_liq(i)
            ENDIF
         ENDDO   
         IF(np.gt.1)then
            CALL allreducei_r1(qvolsum)
         ENDIF                   
         DO i=1,ncell_fluid
             x=xloc(i,1)
             y=xloc(i,2)
             z=xloc(i,3)
             r=DSQRT(x*x+y*y)
             qvol_liq(i)=0.0d0
             qvol_liq_init(i)=0.0d0
             IF(r.le.rdc.and.z.le.zh.and.z.ge.zl)then
               qvol_liq(i)=tpower/corevol*(tpower/qvolsum)       
               qvolsum=qvolsum+qvol_liq(i)*volp(i)*zfactor(i)
               qvol_liq_init(i)=qvol_liq(i)
             ENDIF
          ENDDO           
         IF(myrank.eq.0)then
            WRITE(*,"(11x,a)")'##Confirm heater region in udfn_set_qvol_liquid_atlas!'
            WRITE(*,"(11x,a,1pe17.5,a)")'##qvol_liq_sum=',qvolsum,'J'         
            WRITE(*,"(11x,a,1pe17.5,a)")'##core_liquid_volume=',corevol,'m3'        
               DO i=1,15
               time_scram=-5.0d0+tscram(i)
               hfactor=heat_factor(time_scram)  
            ENDDO
         ENDIF  
     ENDIF
!
!.....reactor trip and decay power
!
     
      IF(rx_trip_signal.eq.0)then
         hfactor=1.0d0  
         rx_trip_signal=reactor_trip_function(rx_trip_signal)
         IF(rx_trip_signal.eq.1)then
            CALL broadcast_i1(rx_trip_signal)
            IF(vv_prob.eq.'atlas_mc_porous')THEN
                time_scram=time+12.07d0
            ELSEIF(vv_prob.eq.'pwr_mc_poro'.or.vv_prob.eq.'apr1400_mc_poro'.or.vv_prob.eq.'opr1000_mc_poro')THEN
               time_scram=time 
            ENDIF         
            IF(myrank.eq.0)then
                WRITE(unit_log,"(11x,a)")'##==============================='
                WRITE(unit_log,"(11x,a)")'##==Reactor is tripped in CUPID=='
                WRITE(unit_log,"(11x,a,e17.5)")'##==time_scram=',time_scram
                WRITE(unit_log,"(11x,a)")'##==============================='
            ENDIF
         ENDIF
      ELSEIF(rx_trip_signal.eq.1)then
         hfactor=1.0d0
         time_after_scram=time-time_scram
         hfactor=heat_factor(time_after_scram)
      ENDIF
!
!.....ramping power
!
      tfactor=1.0d0
      IF(cplmars.eq.2)THEN
         IF(initial_atlas)THEN
            time_ramp=time
            IF(myrank.eq.0)WRITE(unit_log,"(11x,a,1i3,2(1pe14.4))")'cplmars,time,time_ramp=',cplmars,time,time_ramp
            IF(myrank.eq.0)WRITE(*,"(11x,a,1i3,2(1pe14.4))")'cplmars,time,time_ramp=',cplmars,time,time_ramp
         ENDIF
         IF(time.lt.time_ramp+1.0d0)THEN
            tfactor=(time-time_ramp)
         ELSE
            tfactor=1.0d0
         ENDIF
      ENDIF
!      
      DO i=1,ncell_fluid
         qvol_liq(i)=qvol_liq_init(i)*hfactor*zfactor(i)*tfactor
      ENDDO
!
!.....sum qporous_liq,qporous_gas
!     
      ql=0.0d0
      IF(MOD(itim,100).eq.0)THEN
         DO i=1,ncell_fluid
            ql=ql+qvol_liq(i)*volp(i)
         ENDDO
         CALL allreducei_r1(ql)
         IF(myrank.eq.0)WRITE(*,1009)rx_trip_signal,ql,hfactor
         IF(myrank.eq.0)WRITE(unit_log,1009)rx_trip_signal,ql,hfactor
    1009 FORMAT(11x,'--Rxtrip,qliq,hf=',1i3,2(1pe14.4))        
      ENDIF
!      
      initial_atlas=.false.
! 
      RETURN
      ENDSUBROUTINE udfn_set_qvol_liq_atlas     
