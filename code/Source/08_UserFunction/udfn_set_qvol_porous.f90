!
      SUBROUTINE udfn_set_qvol_porous
!
!     Define volumetric heat source in porous cells
!
      USE VOL_DATA              
      USE Zconst1   , ONLY: vv_prob,cplmaster,cplmars,restart
      USE Zcoord1   , ONLY: xloc,xloc_c
      USE Zcoord3   , ONLY: volp       
      USE Zcore     , ONLY: np,cupid_mars,myrank
      USE Znum_cell , ONLY: n_fluid
      USE Zqvol     , ONLY: qvol_ice_solid,qvol_liq,qvol_gas,q0_liq,q0_ice_solid
      USE Ztimecon  , ONLY: time
      USE Zzone     , ONLY: ncell_fluid,ncell_cond,nzone,nzone_c
      
      !OPR1000 MASTER
      USE MASTER4   , ONLY:power_master,ppm_mas
      USE MASTER4   , ONLY: mas_wait,mas_interval,mas_init_duration
      USE Zio_unit  , ONLY: unit_log
!
      IMPLICIT NONE
!
      INCLUDE '../10_LinkToMARS/c3com.h' 
!      
!.....Local variables
      INTEGER(4) i,k         
      INTEGER ii,izone               
      INTEGER, SAVE:: dtemp_pass
      REAL(8), SAVE :: qvol_sg,q_average,q_tot
      REAL(8), SAVE :: coeff_PASCAL(40)
      ! OPR1000-MASTER
      integer,SAVE:: i_opt
      REAL(8),SAVE:: call_time,call_time_interval,set_time_trn,call_time_initial
!     LOGICAL,SAVE:: initial,initial_cplmars3,initial_rod,final_rod
!     DATA initial,initial_cplmars3,initial_rod,final_rod/.true.,.true.,.true.,.true./
      LOGICAL,SAVE :: initial=.true.
      LOGICAL,SAVE :: initial_cplmars3=.true.
!      
      DATA coeff_PASCAL/1.3177014597787, 1.2561009375486, 1.1945004153184, 1.1410647917041, 1.1012912261697, &
                        1.0615176606352, 1.0381547329128, 1.0485831677785, 1.0590116026443, 1.0694400375100, &
                        1.0789792259918, 1.0873058057683, 1.0956323855449, 1.1039589653214, 1.1111537769729, &
                        1.1173785016602, 1.1236032263475, 1.1299087916152, 1.1272410524635, 1.1206929654548, &
                        1.1141448784460, 1.1076776320177, 1.0821320086256, 1.0562630229122, 1.0303131966184, &
                        1.0044442109049, 0.9831831382717, 0.9623262685402, 0.9414693988088, 0.9190765580506, &
                        0.8658834561773, 0.8126903543040, 0.7594972524307, 0.7137414839501, 0.7039597737272, &
                        0.6941780635043, 0.6843155127010, 0.6887617446205, 0.7271610202889, 0.7655602959573/
!
!.....Halden650_5
!
      IF(vv_prob.eq.'halden650_5')THEN   
         DO i=1,ncell_cond
            qvol_ice_solid(i)=0.0d0            
            IF(time.le.44.6604d0)THEN
               IF(q0_ice_solid(nzone_c(i)).gt.0.0D0)qvol_ice_solid(i)=16.31182d0/17.0d0*q0_ice_solid(nzone_c(i))
            ELSEIF(time.le.151.42521)THEN
               IF(q0_ice_solid(nzone_c(i)).gt.0.0D0)qvol_ice_solid(i)=16.39859/17.0d0*q0_ice_solid(nzone_c(i))  
            ELSEIF(time.le.157.68504)THEN
               IF(q0_ice_solid(nzone_c(i)).gt.0.0D0)qvol_ice_solid(i)=16.81802/17.0d0*q0_ice_solid(nzone_c(i))  
            ELSEIF(time.le.262.7110)THEN
               IF(q0_ice_solid(nzone_c(i)).gt.0.0D0)qvol_ice_solid(i)=16.73124/17.0d0*q0_ice_solid(nzone_c(i))  
            ELSEIF(time.le.334.69901)THEN
               IF(q0_ice_solid(nzone_c(i)).gt.0.0D0)qvol_ice_solid(i)=16.9048/17.0d0*q0_ice_solid(nzone_c(i))  
            ELSEIF(time.le.451.89687)THEN
               IF(q0_ice_solid(nzone_c(i)).gt.0.0D0)qvol_ice_solid(i)=17.0/17.0d0*q0_ice_solid(nzone_c(i))                
            ELSEIF(time.le.456.07008)THEN
               IF(q0_ice_solid(nzone_c(i)).gt.0.0D0)qvol_ice_solid(i)=12.45022/17.0d0*q0_ice_solid(nzone_c(i)) 
            ELSEIF(time.le.456.17008)THEN
               IF(q0_ice_solid(nzone_c(i)).gt.0.0D0)qvol_ice_solid(i)=0.0d0                
            ELSE
               qvol_ice_solid(i)=0.0d0
            ENDIF
         ENDDO
      ENDIF      
!
!.....PAFS-POOL  
!
      IF(vv_prob.eq.'PAFS-POOL')THEN   
         IF(time.le.10000)THEN
            qvol_sg=3.16776e-17*time**5-7.68614e-13*time**4+5.99256e-09*time**3-1.28272e-05*time**2+1.74947e-2*time+1.64534e+02
         ELSE
            qvol_sg=524.8415195
         ENDIF
         IF(np.le.1) THEN 
            q_average = qvol_sg*1000.d0/(40.d0*(1.d0-0.90952d0)*0.112d0*0.2d0*0.20303d0)
            DO i=1,ncell_cond
               ii=n_fluid(i)
               IF(ii.eq.1720) THEN
                  qvol_ice_solid(i)=coeff_PASCAL(1)*q_average*cell%alphal(ii)
               ELSEIF(ii.eq.1719) THEN
                  qvol_ice_solid(i)=coeff_PASCAL(2)*q_average*cell%alphal(ii)                  
               ELSEIF(ii.eq.1698) THEN
                  qvol_ice_solid(i)=coeff_PASCAL(3)*q_average*cell%alphal(ii)                                    
              ELSEIF(ii.eq.1699) THEN
                  qvol_ice_solid(i)=coeff_PASCAL(4)*q_average*cell%alphal(ii)
              ELSEIF(ii.ge.1700.and.ii.le.1714) THEN
                  k=5+(ii-1700)
                  qvol_ice_solid(i)=coeff_PASCAL(k)*q_average*cell%alphal(ii)
              ELSEIF(ii.eq.1718) THEN
                  qvol_ice_solid(i)=coeff_PASCAL(20)*q_average*cell%alphal(ii)                   
              ELSEIF(ii.eq.1717) THEN
                  qvol_ice_solid(i)=coeff_PASCAL(21)*q_average*cell%alphal(ii)                                     
              ELSEIF(ii.eq.1718) THEN
                  qvol_ice_solid(i)=coeff_PASCAL(22)*q_average*cell%alphal(ii)
              ELSEIF(ii.le.1697.and.ii.ge.1681) THEN
                  k=24+(1697-ii)
                  qvol_ice_solid(i)=coeff_PASCAL(k)*q_average*cell%alphal(ii)
              ELSEIF(ii.eq.1716) THEN
                  qvol_ice_solid(i)=coeff_PASCAL(39)*q_average*cell%alphal(ii)
              ELSEIF(ii.eq.1715) THEN
                  qvol_ice_solid(i)=coeff_PASCAL(40)*q_average*cell%alphal(ii)
              ENDIF    
            ENDDO   
         ELSE
            q_average = qvol_sg*1000.d0/(40.d0*(1.d0-0.90952d0)*0.112d0*0.2d0*0.20303d0)
            DO i=1,ncell_cond
               ii=n_fluid(i)
               qvol_ice_solid(i)=0.d0
               IF(nzone_c(i).eq.3) THEN
                  IF    (xloc_c(i,2).gt.2.8d0) THEN
                     qvol_ice_solid(i)=coeff_PASCAL(1)*q_average*cell%alphal(ii)
                  ELSEIF(xloc_c(i,2).gt.2.6d0) THEN
                     qvol_ice_solid(i)=coeff_PASCAL(2)*q_average*cell%alphal(ii)
                  ELSEIF(xloc_c(i,2).gt.2.2d0) THEN
                     qvol_ice_solid(i)=coeff_PASCAL(20)*q_average*cell%alphal(ii)
                  ELSEIF(xloc_c(i,2).gt.2.d0) THEN
                     qvol_ice_solid(i)=coeff_PASCAL(21)*q_average*cell%alphal(ii)
                  ELSEIF(xloc_c(i,2).gt.1.6d0) THEN
                     qvol_ice_solid(i)=coeff_PASCAL(39)*q_average*cell%alphal(ii)
                  ELSE
                     qvol_ice_solid(i)=coeff_PASCAL(40)*q_average*cell%alphal(ii)
                  ENDIF   
               ENDIF   
               IF(nzone_c(i).eq.2.and.xloc_c(i,2).ge.2.2d0) THEN
                  DO izone=1,17
                     IF(xloc_c(i,1).gt.0.203d0*(DFLOAT(izone)).and.xloc_c(i,1).lt.0.203d0*(DFLOAT(izone)+1.d0)) THEN
                        qvol_ice_solid(i)=coeff_PASCAL(3+izone-1)*q_average*cell%alphal(ii)
                        EXIT
                     ENDIF
                  ENDDO   
               ELSEIF(nzone_c(i).eq.2.and.xloc_c(i,2).lt.2.2d0) THEN
                  DO izone=1,17
                     IF(xloc_c(i,1).gt.0.203d0*(DFLOAT(izone)).and.xloc_c(i,1).lt.0.203d0*(DFLOAT(izone)+1.d0)) THEN
                        qvol_ice_solid(i)=coeff_PASCAL(38-izone+1)*q_average*cell%alphal(ii)
                        EXIT
                     ENDIF
                  ENDDO   
               ENDIF                  
            ENDDO
         ENDIF
      ENDIF
!      
!.....MP1/MP2 (CUPID/MASTER-LJR)
!.....VFS9/VFS10 (Porous-LJR)
!
      IF(vv_prob.eq.'CEA_drop'.or.vv_prob.eq.'CEA_ejection'.or.vv_prob.eq.'annul_porous') THEN 
         DO i=1,ncell_cond
            qvol_ice_solid(i)=0.d0
            qvol_ice_solid(i)=q0_ice_solid(nzone_c(i))
         ENDDO
      ENDIF
!
!.....STERN moderator
!      
      IF (vv_prob.eq.'stern')THEN
         DO i=1,ncell_cond
            qvol_ice_solid(i)=q0_ice_solid(nzone_c(i))
         ENDDO
      ENDIF
!
!.....2D_boiling, 3D_boiling
! 
      IF (vv_prob.eq.'2D_boiling'.or.vv_prob.eq.'3D_boiling') THEN
         DO i=1,ncell_fluid
            qvol_gas(i)=0.d0
            IF(time.lt.10)THEN
               qvol_liq(i)=q0_liq(nzone(i))*time/10.d0
            ELSE
               qvol_liq(i)=q0_liq(nzone(i))
            ENDIF
         ENDDO
      ENDIF
!
!.....2D_loca
! 
      IF (vv_prob.eq.'2D_loca') THEN
         DO i=1,ncell_fluid
            IF(xloc(i,1).le.0.3d0 .and. xloc(i,2).le.1.d0) qvol_liq(i)=cell%alphal(i)*1.d6
            IF(xloc(i,1).ge.4.7d0 .and. xloc(i,2).gt.8.d0) qvol_liq(i)=-cell%alphal(i)*1.d6
         ENDDO
      ENDIF     
!
!.....atlas_mslb
! 
      IF (vv_prob.eq.'atlas_mc_porous'.or.vv_prob.eq.'pwr_mc_poro'.or.vv_prob.eq.'apr1400_mc_poro'.or.vv_prob.eq.'opr1000_mc_poro')THEN          
          CALL udfn_set_qvol_solid_atlas
      ENDIF     
!
      IF(vv_prob.eq.'sgp_separator') CALL udfn_sg_heat_source
!
!.....OPR1000 MASTER power
!
      IF(vv_prob.eq.'OPR1000_fullvessel_1x1')THEN
         IF(initial)then
            power_master=2815.d6 
            call_time=time+0.d0
            set_time_trn=250.d0
            initial=.false.
!      
!...........MASTER call
!           i_flag=0 : ALLOCATION MASTER Variables
!           i_flag=1 : Steady Calculation
!           i_flag=2 : Transient Calculation
            IF(cplmaster.gt.0)THEN   
               i_opt=0 
               CALL t_masterC(i_opt)

               if(myrank.eq.0)write(*,*)'i_flag3 in udfn_set_qvol_porous'

               IF(restart.eq.1) CALL t_masterC(3)
            ENDIF
         ENDIF    
         IF(time.lt.call_time)RETURN
      
         IF(cplmaster.gt.0)THEN    
             IF(time.le.set_time_trn)then
               call_time=call_time+.5d0
               i_opt=1
               CALL t_masterC(i_opt)
            ELSEIF(time.gt.set_time_trn)then
               call_time=call_time+0.5d0
               i_opt=2
               CALL t_masterC(i_opt)
            ENDIF
         ENDIF
      ENDIF
!   
      IF(vv_prob.eq.'opr1000_mc_rv'.and.cplmaster.gt.0)THEN
!
!........Initialize parameters to control calling master
         IF(initial)then
            initial=.false.
            power_master=2806.03d6
            set_time_trn=50.d0
            call_time=time+1.d0
            call_time_interval=1.d0
            IF(cupid_mars.eq.1)THEN 
               set_time_trn=mas_init_duration
               call_time=time+mas_wait
               call_time_interval=mas_interval
            ENDIF
            call_time_initial=call_time   
            IF(myrank.eq.0)WRITE(* ,"(11x,a,3e11.4,a)")'--time,master_steady,interval=',time,call_time,call_time_interval, 's'  
            IF(myrank.eq.0)WRITE(unit_log,"(11x,a,3e11.4,a)")'--time,master_steady,interval=',time,call_time,call_time_interval, 's'                     
            dtemp_pass=0
         ENDIF    
!  
!.......call master       
         IF(time.ge.call_time)THEN
            i_opt=1
            IF(cplmars.eq.3.and.(time-call_time_initial).gt.set_time_trn.and.dtemp_pass.eq.1)then !MASTER-OPERATION
               IF(initial_cplmars3)THEN
                  initial_cplmars3=.FALSE.
                  call_time_interval=0.1d0
                  IF(myrank.eq.0)WRITE(* ,"(11x,a,3e12.5,a)")'--time, master_trans,interval=',time,call_time,call_time_interval, 's'  
                  IF(myrank.eq.0)WRITE(unit_log,"(11x,a,3e12.5,a)")'--time, master_trans,interval=',time,call_time,call_time_interval, 's'                     
               ENDIF
               i_opt=2
            ENDIF   
            CALL t_masterC(i_opt)
            call_time=call_time+call_time_interval
            q_tot=0.d0
            DO i=1,ncell_fluid
               q_tot=q_tot+qvol_liq(i)*volp(i)
            ENDDO
            IF(np.gt.1) CALL allreducei_r1(q_tot)
            IF(myrank.eq.0) WRITE(* ,"(11x,a,1f10.2,a,1i3,1f10.3,a)")'--MASTER Power, Mode, Boron=',q_tot/1.d6,' MW',i_opt,ppm_mas,' PPM'  
            IF(myrank.eq.0) WRITE(unit_log,"(11x,a,1f10.2,a,1i3,1f10.3,a)")'--MASTER Power, Mode, Boron=',q_tot/1.d6,' MW',i_opt,ppm_mas,' PPM'  
            c3rktpow_ctl_val=q_tot/power_master*100.d0                   
         ENDIF
!
!........print master power and temperature variance
         CALL judge_steady_th(q_tot,i_opt)

!         
      ENDIF
!   
      END SUBROUTINE udfn_set_qvol_porous
!------------------------------------------------------------------------------------------------------
      SUBROUTINE judge_steady_th(q_tot,i_opt)
      USE VOL_DATA    
      USE Zzone      ,ONLY:ncell_fluid,nzone
      USE Zcore      ,ONLY:np,myrank
      USE Ztimecon   ,ONLY:time
      USE Zconst2    ,ONLY:dt  
      USE MASTER4    ,ONLY:mas_dtemp,dtemp_pass
!      
      IMPLICIT NONE
!      
      INTEGER i,i_opt
      REAL(8) q_tot
      REAL(8) dtemp_max,dtemp
      REAL(8) dtemp_max_core
      REAL(8),ALLOCATABLE,SAVE:: tlo(:)
      LOGICAL,SAVE::initial
      DATA dtemp_max,initial/0.d0,.true./ 
!      
      IF(initial)THEN
         IF(myrank.eq.0)THEN
            OPEN(812,file='mas_power.dat')
            WRITE(812,"(a)")'time   dtemp_max   dtemp_max_core   power   mode   dtemp_criterion'
         ENDIF   
         ALLOCATE(tlo(ncell_fluid))
         tlo(1:ncell_fluid)=cell%tl(1:ncell_fluid)
      ENDIF
!      
      dtemp_max=0.d0
      dtemp_max_core=0.d0
      DO i=1,ncell_fluid
         dtemp=DABS(cell%tl(i)-tlo(i))
         dtemp_max=MAX(dtemp_max,dtemp)
         IF(nzone(i).eq.6)dtemp_max_core=MAX(dtemp_max_core,dtemp)
      ENDDO
      IF(np.gt.1)THEN
         CALL allreducei_max_r1(dtemp_max)
         CALL allreducei_max_r1(dtemp_max_core)
      ENDIF      
!
      IF(dtemp_pass.eq.0)THEN
         IF(.not.initial.and.dtemp_max.lt.mas_dtemp)dtemp_pass=1      
      ENDIF   
!
      IF(myrank.eq.0)WRITE(812,"(5e14.6,2i3)")time,dtemp_max,dtemp_max_core,q_tot/1.d6,dt,i_opt,dtemp_pass      
!      
      tlo(1:ncell_fluid)=cell%tl(1:ncell_fluid)
!
      initial=.FALSE.
!
      END SUBROUTINE judge_steady_th
