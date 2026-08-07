!
      SUBROUTINE read_cupid_master
! 
      USE Zconst1      , ONLY: cplmars
      USE Zcore        , ONLY: myrank,cupid_mars
      USE Ztimecon     , ONLY: time
      USE MASTER4      , ONLY: power_master !,cpc_vopt
      USE MASTER4      , ONLY: mas_wait,mas_interval,mas_delay,mas_interval_delay,mas_init_duration
      USE MASTER4      , ONLY: mas_dtemp,mas_dpower,mas_dmass,mas_dtemp_opt,mas_dpower_opt,mas_dmass_opt
      USE MASTER4      , ONLY: dtemp_pass,dpower_pass,dmass_pass,dmaster_pass,count_dtemp,count_dmass,count_dpower 
      USE MASTER4      , ONLY: time_cri,mas_rx_trip      
      USE Zio_unit     , ONLY: unit_log
!
      IMPLICIT NONE
!
      INCLUDE '../10_LinkToMARS/c3com.h' !i3rod_trip      
!      
      INTEGER tripin
!
!DEC$IF defined (master_flag) 
!DEC$ELSE
      IF(myrank.eq.0)WRITE(*,"(11x,a)")'master_flag should be defined as true. Check compiler option!'
      PAUSE
      STOP 
!DEC$ENDIF  
!
      time_cri=-1.0d0
      mas_rx_trip=0
      i3rod_trip=0
!
      IF(cupid_mars.eq.0)THEN
         cplmars=0 
         time=0.0d0
      ENDIF   
!
      CLOSE(812)
!
      OPEN(812,file='cupid_master.in',status='old',iostat=tripin)
!
      IF(tripin.eq.0)THEN
          IF(myrank.eq.0)WRITE(*,"(11x,a)")'Reading cupid_master.in...'      
          IF(cplmars.eq.4)THEN
             READ(812,*)mas_wait,mas_interval
             mas_delay=0.0d0        
             mas_init_duration=0.0d0
             READ(812,*)power_master !,cpc_vopt
             mas_dtemp_opt=0
             mas_dtemp=0.0d0
             count_dtemp=0
             mas_dmass_opt=0
             mas_dmass=0.0d0 
             count_dmass=0   
             mas_dpower_opt=0
             mas_dpower=0.0
             count_dpower=0
          ELSE
             READ(812,*)mas_wait,mas_interval,mas_delay,mas_interval_delay
             READ(812,*)mas_init_duration
             READ(812,*)power_master !,cpc_vopt
             READ(812,*)mas_dtemp_opt, mas_dtemp, count_dtemp
             READ(812,*)mas_dmass_opt, mas_dmass, count_dmass   
             READ(812,*)mas_dpower_opt,mas_dpower, count_dpower
          ENDIF 
          power_master=power_master*1.d6
          dtemp_pass=1
          dpower_pass=1
          dmass_pass=1      
          IF(mas_dtemp_opt.eq.1)dtemp_pass=0
          IF(mas_dpower_opt.eq.1)dpower_pass=0
          IF(mas_dmass_opt.eq.1)dmass_pass=0
          dmaster_pass=dtemp_pass*dpower_pass*dmass_pass
          CLOSE(812)
      ELSE
         IF(myrank.eq.0)WRITE(*,"(11x,a)")'Input file is missing...: <cupid_master.in>.'
         IF(myrank.eq.0)WRITE(unit_log,"(11x,a)")'Input file is missing...: <cupid_master.in>.' 
         PAUSE
         STOP
      ENDIF
!
      RETURN 
      END SUBROUTINE read_cupid_master