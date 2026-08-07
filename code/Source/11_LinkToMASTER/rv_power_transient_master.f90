!
      SUBROUTINE rv_power_transient_master                  
!      
      USE Zqvol        , ONLY: qporous_gas,qporous_liq,qporous_gamma,qrv_gas,qrv_liq,qrv_gamma       
      USE Ztimecon     , ONLY: time
      USE Zcore        , ONLY: np,myrank 
      USE zconst1      , ONLY: cplmars,restart
      USE Zrv_model    , ONLY: rv_ht_w,rv_ht_str      
!      USE Zuserdefined , ONLY: udfl_rv_wall_control
      USE Zzone        , ONLY: ncell_fluid         
      USE Zrv_hts_2d   , ONLY: power_2d 
!      USE Zmars        , ONLY: c3rktpow_ctl_val 
      USE MASTER4      , ONLY: power_master,ppm_mas,PPCT_MASTER
!     USE MASTER4      , ONLY: mas_wait,mas_interval,mas_delay,mas_interval_delay,mas_init_duration
      USE MASTER4      , ONLY: mas_wait,mas_interval,mas_interval_delay,mas_init_duration
      USE MASTER4      , ONLY: mas_dtemp_opt,mas_dpower_opt,mas_dmass_opt
      USE MASTER4      , ONLY: dtemp_pass,dpower_pass,dmass_pass,dmaster_pass
      USE Zio_unit     , ONLY: unit_log
!
      IMPLICIT NONE
!      
      INCLUDE '../10_LinkToMARS/c3com.h' !c3rktpow_ctl_val
!
!.....Local variables
      INTEGER i 
      INTEGER,SAVE:: i_opt,rv_ht_w_old,rv_ht_str_old
      INTEGER,SAVE:: unit_power=87
      LOGICAL,SAVE:: initial=.true.
      LOGICAL,SAVE:: initial_cplmars0=.true.
!     LOGICAL,SAVE:: initial_cplmars3=.true.
      LOGICAL,SAVE:: initial_cplmars4=.true.
      REAL(8),SAVE:: dpower_max=0.0d0
      REAL(8),SAVE:: dmpower_max=0.0d0
      REAL(8),SAVE:: dmass_max=0.0d0
      REAL(8),SAVE:: dtemp_max=0.0d0
      REAL(8),SAVE:: power_max,mpower_max      
      REAL(8),SAVE:: power_master_init
      REAL(8),SAVE:: call_time,call_time_interval,call_time_initial
      REAL(8) qporous_sum,qrv_sum
!.....Local arrays
      REAL(8) :: tmp(2)
!      
!.....initialize power and master control option
!
      IF(initial)THEN  
         initial=.false.
!                  
!........transfer power from master to cupid                  
         CALL rv_power_master_to_cupid 
!         
!........set master control parameter
         power_master_init=power_master      
         IF(cplmars.eq.2)power_master=power_master_init*0.1d0 !see below
         call_time=time+mas_wait
         call_time_interval=mas_interval
         call_time_initial=call_time 
         IF(myrank.eq.0)WRITE(* ,"(11x,a,3e11.4,a)")'--time,master_steady,interval=',time,call_time,call_time_interval, 's'  
         IF(myrank.eq.0)WRITE(unit_log,"(11x,a,3e11.4,a)")'--time,master_steady,interval=',time,call_time,call_time_interval, 's'                     
         dmaster_pass=0
!
!........check qporous and qrv
         qporous_sum=0.0
         qrv_sum=0.0
         DO i=1,ncell_fluid
            qporous_sum=qporous_sum+(qporous_liq(i)+qporous_gas(i)+qporous_gamma(i)) !*aporous(i*volp(i)
            qrv_sum    =qrv_sum    +(qrv_liq(i)+qrv_gas(i)+qrv_gamma(i)) !*aporous(i)*volp(i)
         ENDDO
         IF(np.gt.1) THEN
            tmp(1)=qporous_sum
            tmp(2)= qrv_sum
            CALL allreducei_r(tmp,2)
            qporous_sum=tmp(1)
            qrv_sum    =tmp(2)
         ENDIF
         IF(myrank.eq.0)WRITE(*,"(11x,a,e15.5,3f8.2,a,1i2)")'--1time,power,qporous,qrv,rv_ht_w =',time,power_2d/1.0d6,qporous_sum/1e6,qrv_sum/1e6,' MW',rv_ht_w         
         IF(myrank.eq.0)WRITE(unit_log,"(11x,a,e15.5,3f8.2,a,1i2)")'--1time,power,qporous,qrv,rv_ht_w =',time,power_2d/1.0d6,qporous_sum/1e6,qrv_sum/1e6,' MW',rv_ht_w         
!               
!........skip rv heat structure and wall heat transfer during first 2 seconds.          
!         IF(udfl_rv_wall_control)THEN 
            rv_ht_w_old=rv_ht_w
            rv_ht_str_old=rv_ht_str
            rv_ht_w=0 
            rv_ht_str=0              
            IF(myrank.eq.0)WRITE(*,"(11x,a,2i3)")'--reserve rv wall heat',rv_ht_w,rv_ht_str
            IF(myrank.eq.0)WRITE(unit_log,"(11x,a,2i3)")'--reserve rv wall heat',rv_ht_w,rv_ht_str
!         ENDIF
!            
!........open a power log file
!            
         IF(myrank.eq.0)OPEN(unit_power,file='power_master_reactor.dat',status='replace')  
         IF(myrank.eq.0)WRITE(unit_power,"(a)")'Time(s)    MASTER(MW) CUPID(MW)  Qrod(MW)   Qporous(MW) Mode      Boron(PPM)'          
      ENDIF
!
!.....bypass calling MASTER
!
      IF(time.lt.call_time)THEN
         IF(time.ge.call_time_initial)THEN
!
!........Check the steady state of the core                               
            IF(dmaster_pass.eq.0)THEN          
                IF(mas_dtemp_opt)CALL check_steady_dtl(dtemp_max)
                IF(mas_dmass_opt)CALL check_steady_dmass(dmass_max)
                dmaster_pass=dpower_pass*dtemp_pass*dmass_pass
            ENDIF
!         
!...........t_end=time(dmaster_pass=1)+mas_delay in subroutine change_tend_mars_cupid or change_tend_cupid        
            IF(cplmars.eq.0)THEN 
               IF(restart.eq.0)THEN
                  IF((time-call_time_initial).lt.mas_init_duration)dmaster_pass=0         
               ENDIF   
            ELSEIF(cplmars.ge.2.and.cplmars.le.3)THEN
               IF((time-call_time_initial).lt.mas_init_duration)dmaster_pass=0         
            ENDIF 
!
!...........print variations of temperaure, power, mass flow rate           
            qrv_sum=0.0
            DO i=1,ncell_fluid
               qrv_sum=qrv_sum+(qrv_liq(i)+qrv_gas(i)+qrv_gamma(i)) !*aporous(i)*volp(i)
            ENDDO
            IF(np.gt.1) CALL allreducei_r1(qrv_sum)
            CALL print_ss_status(power_2d,qrv_sum,i_opt,dtemp_max,dpower_max,dmpower_max,dmass_max,dmaster_pass,power_max,mpower_max)
!
!...........Control call_time_interval
            IF(cplmars.eq.0)THEN
               IF(restart.eq.0)THEN
                  IF(dmaster_pass.eq.1)call_time_interval=mas_interval_delay
               ENDIF   
            ELSEIF(cplmars.eq.2)THEN
               IF(dmaster_pass.eq.1)call_time_interval=mas_interval_delay
            ELSEIF(cplmars.eq.3)THEN
               IF(dmaster_pass.eq.1)call_time_interval=mas_interval_delay
            ENDIF 
!                                     
         ENDIF
!         
         RETURN
      ENDIF   
!      
!.....set MASTER mode 
!
!
!.....back to original rv_ht_w & rv_ht_str. 
!      IF(udfl_rv_wall_control)THEN 
         rv_ht_w=rv_ht_w_old 
         rv_ht_str=rv_ht_str_old 
!      ENDIF   
! 
!.....default MASTER mode==1 (initialization)    
      i_opt=1
!         
!.....MASTER/CUPID without MASRS        
      IF(cplmars.eq.0)THEN
!
!........master/cupid 1st step, change MASTER MODE to operation mode              
         IF(restart.eq.0)THEN
            IF(dmaster_pass.eq.1)then
               i_opt=2
            ENDIF
!
!........master/cupid 2nd step, restart MASTER in operation mode              
         ELSE
            IF(initial_cplmars0)THEN
               initial_cplmars0=.FALSE.
               i_opt=3
            ELSE   
               i_opt=2   
            ENDIF 
         ENDIF
!         
!.....MASTER/CUPID/MARS 2nd step, ramp power       
      ELSEIF(cplmars.eq.2)THEN
         IF(power_master.lt.power_master_init)THEN
            power_master=power_master+power_master_init*0.2
            IF(power_master.gt.power_master_init)power_master=power_master_init         
         ENDIF    
!         
!.....MASTER/CUPID/MARS 3rd step, change MASTER MODE to operation mode         
      ELSEIF(cplmars.eq.3)THEN
         IF(dmaster_pass.eq.1)then
            i_opt=2
         ENDIF   
!            IF(initial_cplmars3)THEN !it doesnot work properly.
!               initial_cplmars3=.FALSE.
!               i_opt=4
!            ENDIF            
!         
!.....MASTER/CUPID/MARS 4th step, restart MASTER in operation mode         
      ELSEIF(cplmars.eq.4)THEN
         IF(initial_cplmars4)THEN
            initial_cplmars4=.FALSE.
            i_opt=3
         ELSE   
            i_opt=2
         ENDIF
!
      ENDIF 
!         
!.....call MASTER          
!
      CALL t_masterC(i_opt)
      call_time=call_time+call_time_interval
!                  
!.....transfer power from master to cupid                  
      CALL rv_power_master_to_cupid   
! 
!.....check qporous and qrv       
      qporous_sum=0.0
      qrv_sum=0.0
      DO i=1,ncell_fluid
         qporous_sum=qporous_sum+(qporous_liq(i)+qporous_gas(i)+qporous_gamma(i)) 
         qrv_sum    =qrv_sum    +(qrv_liq(i)+qrv_gas(i)+qrv_gamma(i)) 
      ENDDO
      IF(np.gt.1) THEN
         tmp(1)=qporous_sum
         tmp(2)= qrv_sum
         CALL allreducei_r(tmp,2)
         qporous_sum=tmp(1)
         qrv_sum    =tmp(2)
      ENDIF
      IF(myrank.eq.0)THEN
         WRITE(* ,"(11x,a,2f8.2,a,1i2,1f8.2,a,2f8.2,a)")'4Fuel&MAS Power,Mode,Boron,qporous,qrv=',power_2d/1.0d6,PPCT_MASTER/1.0d6,' MW',i_opt,ppm_mas,' PPM',qporous_sum/1e6,qrv_sum/1e6,' MW'      
         WRITE(unit_log,"(11x,a,2f8.2,a,1i2,1f8.2,a,2f8.2,a)")'4Fuel&MAS Power,Mode,Boron,qporous,qrv=',power_2d/1.0d6,PPCT_MASTER/1.0d6,' MW',i_opt,ppm_mas,' PPM',qporous_sum/1e6,qrv_sum/1e6,' MW'      
         WRITE(unit_power,"(1e12.5,4f11.3,1i3,1f11.3)")Time,PPCT_MASTER/1.0d6,power_2d/1.0d6,qrv_sum/1e6,qporous_sum/1e6,i_opt,ppm_mas
      ENDIF   
!        
!.....mars control value for overpower trip 
      IF(np.gt.1) CALL broadcast_r1(PPCT_MASTER)
      c3rktpow_ctl_val=PPCT_MASTER/power_master*100.0d0                      
!
!.....Check the steady state of the core       
      IF(dmaster_pass.eq.0)THEN          
          IF(mas_dpower_opt) THEN
             CALL check_steady_dpower(dpower_max,power_max)
             IF(DABS((power_2d-qrv_sum)*1.d-6).gt.1.0d-1)dpower_pass=0
          ENDIF   
          !IF(mas_dpower_opt) CALL check_steady_dpower_master(dmpower_max,mpower_max)
          IF(mas_dtemp_opt)  CALL check_steady_dtl(dtemp_max)
          IF(mas_dmass_opt)  CALL check_steady_dmass(dmass_max)
          dmaster_pass=dpower_pass*dtemp_pass*dmass_pass
      ENDIF
!         
!.....t_end=time(dmaster_pass=1)+mas_delay in subroutine change_tend_mars_cupid or change_tend_cupid        
      IF(cplmars.eq.0)THEN 
         IF(restart.eq.0)THEN
            IF((time-call_time_initial).lt.mas_init_duration)dmaster_pass=0         
         ENDIF   
      ELSEIF(cplmars.ge.2.and.cplmars.le.3)THEN
         IF((time-call_time_initial).lt.mas_init_duration)dmaster_pass=0         
      ENDIF 
!
!.....print variations of temperaure, power, mass flow rate           
      CALL print_ss_status(power_2d,qrv_sum,i_opt,dtemp_max,dpower_max,dmpower_max,dmass_max,dmaster_pass,power_max,mpower_max)
!
!.....Control call_time_interval
      IF(cplmars.eq.0)THEN
         IF(restart.eq.0)THEN
            IF(dmaster_pass.eq.1)call_time_interval=mas_interval_delay
         ENDIF   
      ELSEIF(cplmars.eq.2)THEN
         IF(dmaster_pass.eq.1)call_time_interval=mas_interval_delay
      ELSEIF(cplmars.eq.3)THEN
         IF(dmaster_pass.eq.1)call_time_interval=mas_interval_delay
      ENDIF    
!      
      END SUBROUTINE rv_power_transient_master
                 
