      SUBROUTINE udfn_set_qvol_solid_atlas
!
!     volumetric heat source in porous cells
!
      USE VOL_DATA  
      USE SOLID_DATA, ONLY: solid
      USE Zparam    , ONLY: pi            
      USE Zconst1   , ONLY: vv_prob,cplmars
      USE Zcoord1   , ONLY: xloc
      USE Zcoord3   , ONLY: vol,porosity  
      USE Zcore     , ONLY: np,myrank
      USE Zmpi      , ONLY: ncell_ps
      USE Znum_cell , ONLY: n_fluid
      USE Zqvol     , ONLY: qvol_ice_solid,qporous_gas,qporous_liq
      USE Ztimecon  , ONLY: time,itim
      USE Zzone     , ONLY: ncell_cond
      USE Zmars     , ONLY: rx_trip_signal, time_scram
      USE Zio_unit  , ONLY: unit_log
!
      IMPLICIT NONE
!      
      INTEGER :: i,ii

      INTEGER ::reactor_trip_function
!
      LOGICAL, SAVE :: initial_atlas,initial
      DATA initial_atlas,initial/.true.,.true./           
!      
      REAL(8) :: x,y,z,r,hfactor,heat_factor,time_after_scram,tfactor
      REAL(8) :: qvolsum
      REAL(8),SAVE :: corevol,coremas
      REAL(8) :: zh,zl,zc,zi,zr,rdc,tpower           
      REAL(8),ALLOCATABLE:: zfactor(:)
      SAVE zfactor
      REAL(8) :: qporl,qporg,qpor
      REAL(8) :: porous_s !mcc-pik
      REAL(8) :: tscram(15)
      DATA tscram/0,5,10,20,40,50,100, 200,300,400, &
                  500,1000,2500,5005,10000/
      REAL(8), SAVE :: time_ramp
!
      REAL(8),ALLOCATABLE :: qvol_ice_solid_init(:),volps(:)
      SAVE qvol_ice_solid_init
      SAVE volps
      IF(initial) THEN
         ALLOCATE(qvol_ice_solid_init(ncell_ps),volps(ncell_ps))
         initial=.false.
      ENDIF
!
!.....core bottom,top,center and radius, total power
!
      IF(vv_prob.eq.'atlas_mc_porous')THEN
          zh=-0.771d0
          zl=-2.676d0
          zc=(zh+zl)/2.0d0
          zi=(zh-zl)/11         
          rdc=0.168d0
          !tpower=1.566e6
          tpower=1.633777e6 !kdh-consider heat loss  
      ELSEIF(vv_prob.eq.'pwr_mc_poro')THEN
          zh=-1.1d0      !OLD
          zl=-4.9d0
          zc=(zh+zl)/2.0d0
          zi=(zh-zl)/11    
          rdc=1.994d0
          tpower=3.983e9    
      ELSEIF(vv_prob.eq.'apr1400_mc_poro')THEN
          zh=-5.4d0     !NEW 
          zl=-9.22d0
          zc=(zh+zl)/2.0d0
          zi=(zh-zl)/11         
          rdc=1.77d0
          tpower=3.983e9          
      ELSEIF(vv_prob.eq.'opr1000_mc_poro')THEN
          zh=4.0d0
          zl=0.2d0
          rdc=1.7526d0
          zc=(zh+zl)/2.0d0
          zi=(zh-zl)/11             
          tpower=2806.03d6 !rktpow in mars input !apr1400=4000MW, opr1400=3000MW
      ENDIF
!                                
!.....initialize power and axial power shape
!                                
      IF(initial_atlas.eq..true.)then
         initial_atlas=.false.
         rx_trip_signal=0
         IF(time.gt.5110.0d0)then
           IF(myrank.eq.0)THEN
              WRITE(unit_log,*)'########## rx_trip_signal and time_scram are set manually!'
              WRITE(*,*)'########## rx_trip_signal and time_scram are set manually!'
           ENDIF   
           rx_trip_signal=1
           time_scram=5110.0d0
         ENDIF
         corevol=0.0d0
         coremas=0.0d0
         time_ramp=time
!
!........core volume and mass
!               
         DO i=1,ncell_cond
            ii=n_fluid(i)
            porous_s=1.0d0-porosity(ii)
            volps(i)=vol(ii)*porous_s
            x=xloc(ii,1)
            y=xloc(ii,2)
            z=xloc(ii,3)
            r=DSQRT(x*x+y*y)
            IF(r.le.rdc.and.z.le.zh.and.z.ge.zl)then
               corevol=corevol+volps(i)
               coremas=coremas+volps(i)*solid%rhocps(i)/500.d0               
            ENDIF
         ENDDO
         IF(np.gt.1)then
             CALL allreducei_r1(corevol)
             CALL allreducei_r1(coremas)             
         ENDIF   
!
!........assign initial power and axial power shape(zfactor)
!         
         qvolsum=0.0d0
         qvol_ice_solid(:)=0.0d0
         qvol_ice_solid_init(:)=0.0d0
         ALLOCATE(zfactor(ncell_cond))
         zfactor(:)=1.0d0  
         DO i=1,ncell_cond
            ii=n_fluid(i)
            x=xloc(ii,1)
            y=xloc(ii,2)
            z=xloc(ii,3)
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
            IF(r.le.rdc.and.z.le.zh.and.z.ge.zl)then
               qvol_ice_solid(i)=tpower/corevol            
               qvolsum=qvolsum+qvol_ice_solid(i)*volps(i)*zfactor(i)
               qvol_ice_solid_init(i)=qvol_ice_solid(i) !must be saved
            ENDIF
         ENDDO   
         IF(np.gt.1) CALL allreducei_r1(qvolsum)
         IF(myrank.eq.0)then
            WRITE(*,"(11x,a)")'Confirm heater region in udfn_set_qvol_solid_atlas!'
            WRITE(*,"(11x,a,1pe17.5,a)")'qvolsolid_sum=',qvolsum/1.0d6,'MJ'         
            WRITE(*,"(11x,a,1pe17.5,a)")'core solid volume=',corevol,'m3'   
            WRITE(*,"(11x,a,1pe17.5,a)")'core solid mass=',coremas/1000.d0,'ton' 
            WRITE(unit_log,"(11x,a)")'Confirm heater region in udfn_set_qvol_solid_atlas!'
            WRITE(unit_log,"(11x,a,1pe17.5,a)")'core power=',qvolsum/1.0d6,'MJ'         
            WRITE(unit_log,"(11x,a,1pe17.5,a)")'core solid volume=',corevol,'m3'   
            WRITE(unit_log,"(11x,a,1pe17.5,a)")'core solid mass=',coremas/1000.d0,'ton'             
            WRITE(unit_log,*)'----------print hfactor----------'
!
!...........print time-dependent power shape
!            
            DO i=1,15
               time_scram=-5.0d0+tscram(i)
               hfactor=heat_factor(time_scram)  
               IF(time_scram.ge.5000.0d0)EXIT           
               WRITE(unit_log,1000)time_scram,hfactor
          1000 FORMAT(10x,'time_cram,hfactor=',2(1pe17.5))
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
            IF(vv_prob.eq.'atlas_mc_porous')THEN
               time_scram=time+12.07d0
            ELSEIF(vv_prob.eq.'pwr_mc_poro'.or.vv_prob.eq.'apr1400_mc_poro'.or.vv_prob.eq.'opr1000_mc_poro')THEN
               time_scram=time !apr1400-rx-trip
            ENDIF         
            IF(myrank.eq.0)then
                WRITE(*,"(11x,a)")         '======================================='
                WRITE(*,"(11x,a)")         '========Reactor is tripped in CUPID===='
                WRITE(*,"(11x,a,1e17.5,a)")'=====Time_scram:',time_scram,"======="
                WRITE(*,"(11x,a)")         '======================================='
                WRITE(unit_log,"(11x,a)")         '======================================='
                WRITE(unit_log,"(11x,a)")         '========Reactor is tripped in CUPID===='
                WRITE(unit_log,"(11x,a,1e17.5,a)")'=====Time_scram:',time_scram,"======="
                WRITE(unit_log,"(11x,a)")         '======================================='
             ENDIF
         ENDIF
      ELSEIF(rx_trip_signal.eq.1)then
         hfactor=1.0d0
         time_after_scram=time-time_scram
         hfactor=heat_factor(time_after_scram)
      ENDIF
!
!.....time-dependent power
!
      tfactor=1.0d0
      IF(cplmars.eq.2)THEN
         IF(time.lt.time_ramp+1.0d0)THEN
            tfactor=(time-time_ramp)
         ELSE
      tfactor=1.0d0
         ENDIF
      ENDIF  
      DO i=1,ncell_cond
         ii=n_fluid(i)
         qvol_ice_solid(i)=qvol_ice_solid_init(i)*hfactor*zfactor(i)*tfactor
      ENDDO
!
!.....sum qporous_liq,qporous_gas
!     
      qporl=0.0d0
      qporg=0.0d0
      qpor=0.0d0
      IF(MOD(itim,100).eq.0)THEN
         DO i=1,ncell_cond
            ii=n_fluid(i)
            qporl=qporl+qporous_liq(ii)
            qporg=qporg+qporous_gas(ii)     
            qpor=qpor+qvol_ice_solid(i)*volps(i)
         ENDDO
         IF(np.gt.1)THEN
            CALL allreducei_r1(qporl)
            CALL allreducei_r1(qporg)
            CALL allreducei_r1(qpor)
         ENDIF
         IF(myrank.eq.0)WRITE(*,1009)rx_trip_signal,qpor,qporl,qporg,hfactor
         IF(myrank.eq.0)WRITE(unit_log,1009)rx_trip_signal,qpor,qporl,qporg,hfactor
    1009 FORMAT(11x,'--Rxtrip,qsol,qliq,qgas,hf=',1i3,4(1pe10.3))        
      ENDIF 
!
      RETURN
!     
      ENDSUBROUTINE udfn_set_qvol_solid_atlas           
                
