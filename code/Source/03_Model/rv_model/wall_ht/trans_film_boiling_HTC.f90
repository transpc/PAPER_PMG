!
      SUBROUTINE trans_film_boiling_HTC(i,mode)
!
!     This routine calculates wall heat transfer coefficient and heat flux for transition boiling and film boiling
!     when reflood option is OFF
!    
      USE VOL_DATA      , ONLY: cell
      USE STM_TBL_cupid , ONLY: st_tbl,             &
                                nt,np,ns,ns2,ndxstd
      USE Zmpi         , ONLY: jperm
      USE Zparam       , ONLY: pi
      USE Ztimecon     , ONLY: time
      USE Zwall_HTC    , ONLY: qflux_t,qflux_l,qflux_g,HTC_t,HTC_tl,HTC_tst,HTC_tg,HTC_tgp, &
                               f_direc,c_direc1,c_direc2,dia_rod,                           &
                               mflux_tota,vfg,                                              &
                               hfg,chf,rho,cond,viscos,cp,beta,sigma,                       &
                               tw,tsat_t,dt_sat
!
      IMPLICIT NONE
!      
      INTEGER i,it,mode
!
      REAL(8) gravity
      REAL(8) prop(36),drho
      REAL(8) ag,a_frac
      REAL(8) HTC_brom,l_beren
      REAL(8) c1,c2,c3,c4,c5,g1       
      REAL(8) qflux_tb_l,qflux_tb_g,qflux_tb,qflux_rad_l
      REAL(8) qflux_fb_l,qflux_fb_g,qflux_fb
      REAL(8) HTC_tb_l,HTC_tb_g,HTC_fb_l,HTC_fb_g
      REAL(8) d_max,d_ave,af_rad,e_f,e_g,e_w,temp,r1,r2,r3,F_wf
!     REAL(8) tdiff
!
      LOGICAL err
!      
      DATA gravity/9.81d0/
      qflux_tb_l=0.0d0
      a_frac=0.0d0
      prop(:)=0.d0     
!
!.....Assign basic properties
!
      prop(2)=cell%p(i)
      cp=cell%cpg(i)
      beta=cell%betag(i)
      IF(cell%quala(i).gt.1.0d-9)THEN
         rho=cell%rhog(i)
         viscos=cell%lviscosg(i)
         cond=cell%lcondg(i)
      ELSE
         prop(1)=DMIN1(973.15d0,(tw+tsat_t)*0.5d0)
         CALL sth2x3_cupid(prop,it,err,                       &
                           st_tbl(ndxstd),                    &
                           st_tbl(ndxstd+nt),                 &
                           st_tbl(ndxstd+nt+np+13*ns+13*ns2))
         IF(err)THEN 
            PRINT *, '#### ERROR: Steam property in tran_film_boiling.f90 at time= ', time, 'at cell= ', jperm(i)
            PAUSE             
         STOP
         ENDIF
         rho=1.0d0/prop(3)                                                   
         IF(prop(9).eq.0.d0)rho=cell%rhog(i) 
         CALL thcond_cupid_rv(prop(1),rho,cond)                                         
         CALL viscos_rv(prop(1),rho,tsat_t,viscos)      
      ENDIF  
      drho=cell%rhol(i)-rho                 
!
!----------------------------------Transition Boiling Modeling----------------------------------------------
!      
!.....Chen Correlation for Transition boiling
!
      IF(tw.le.tsat_t+600.d0)THEN
!
!........Heat franser to liquid phase via wetted area fraction
!
         ag=DMIN1(0.999d0,cell%alphag(i))
         c2=0.05d0/(1.0d0-ag**40.d0)+0.075d0*ag 
         c1=2.4d0*c2 
         c3=0.2d0*c2 
         !g1=mflux_tota *7.373381d-3 
         g1=mflux_tota*1.0d-5
         c4=c1-c2*g1 
         c5=c3*g1 
         a_frac=EXP(-1.34164d0*DMAX1(c4,c5)*DMIN1(15.0d0,SQRT(dt_sat))) 
!
         qflux_tb_l=chf*a_frac
      ENDIF
!
!-------------------------------------Film Boiling Modeling----------------------------------------------
!
!.....Bromley correlation for film boiling: heat transfered to liquid phase(cont.+droplet) via conduction
!
      drho=DMAX1(1.0d-7,drho)
      l_beren=1.58323349d0*(sigma/(gravity*drho))**0.125     !(2*pi)=1.58323349d0**0.25              
      HTC_brom=0.62d0*(gravity*cond**3*rho*drho*(hfg+0.5d0*cp*dt_sat)/(viscos*dt_sat))**0.25d0
      HTC_fb_l=HTC_brom/l_beren**0.25d0
!
!.....Ramp Bromley HTC (0.2<ag<0.999)
!                   
      c1=DMIN1(1.d0,DMAX1(1.2516d0*(0.999-cell%alphag(i)),0.d0))    !  1/(0.999-0.2) = 1.2516                                                
      IF(c1.lt.1.d0)THEN
         c2=(cell%alphag(i)-0.2d0) * 1.2516d0                                                         
         c2=c2*(3.d0*c2-2.d0*c2*c2)     ! Spline fit
         c1=1.d0-c2 
      ENDIF 
      HTC_fb_l=c1*HTC_fb_l 
!
!.....Modify Bromley HTC for subcooling effect (Sudo and Murao correlation)
!         
      HTC_fb_l=HTC_fb_l*(1.d0+0.025d0*DMAX1(tsat_t-cell%tl(i),0.d0))
      HTC_fb_l=HTC_brom     
!
!.....Radiation heat transfer to liquid phase (ignore to gas phase)
!         
      vfg=DMAX1(0.005d0,vfg)
      d_max=SQRT(cell%alphal(i))*dia_rod
      d_ave=7.5d0*sigma/(rho*vfg*vfg)               !We=7.5
      d_ave=DMAX1(1.0d-5,DMIN1(d_max,d_ave))
      af_rad=1.11d0*cell%alphal(i)/d_ave
      e_f=DMIN1(0.75d0,(1.0d0-EXP(-0.9d0*dia_rod*af_rad)))
      !e_f=1.0d0-EXP(-dia_rod*af_rad)                  ! MARS Manual
      e_g=0.02d0
      e_w=0.7d0
      temp=1.0d0-e_f*e_g
      r1=(1-e_g)/(e_g*temp)
      r2=(1-e_f)/(e_f*temp)
      r3=1.0d0/temp+0.111d0                         !0.111=(1-e_w)/e_w
      F_wf=1.0d0/(r2*(1.0d0+r3/r1+r3/r2))
      qflux_rad_l=F_wf*5.67d-8*(tw**4-tsat_t**4)
      HTC_fb_l=HTC_fb_l+qflux_rad_l/dt_sat
!
!.....Gas single phase HTC 
!
      mode=9
      CALL single_phase_HTC(i,mode)
      HTC_fb_g=HTC_tg
!
!.....Interpolation HTC_tb_g for
!      
      HTC_fb_g=HTC_fb_g*cell%alphag(i)
      HTC_fb_g=DMIN1(HTC_fb_g,HTC_tg)
!
!.....Save HTC and heat flux
!
!     Transition boiling
      HTC_tb_l=qflux_tb_l/dt_sat
      HTC_tb_g=HTC_fb_g*(1.0d0-a_frac)
      !qflux_tb_l as itself
      qflux_tb_g=HTC_tb_g*(tw-cell%tg(i))
      qflux_tb=qflux_tb_l+qflux_tb_g
!     Film boiling      
      qflux_fb_l=HTC_fb_l*dt_sat
      qflux_fb_g=HTC_fb_g*(tw-cell%tg(i))
      qflux_fb=qflux_fb_l+qflux_fb_g
!
!.....Selection of HT mechanism via comparing heat fluxes by transition and film boiling
!      
      IF(qflux_tb.lt.qflux_fb)THEN
         mode=7
         HTC_tl=0.0d0
         HTC_tst=HTC_fb_l
         HTC_tg=HTC_fb_g
         HTC_tgp=0.0d0          
         qflux_l=qflux_fb_l
         qflux_g=qflux_fb_g
      ELSE
         mode=5
         HTC_tl=0.0d0
         HTC_tst=HTC_tb_l
         HTC_tg=HTC_tb_g
         HTC_tgp=0.0d0          
         qflux_l=qflux_tb_l
         qflux_g=qflux_tb_g      
      ENDIF
      qflux_t=qflux_l+qflux_g
      HTC_t=HTC_tst+HTC_tg
!
      IF(cell%tl(i).ge.tsat_t-0.05d0)mode=mode+1
!      
!      HTC_tl=0.0d0
!      IF(DABS(dt_sat).gt.1.d-1)THEN
!         HTC_tst=qflux_l/dt_sat
!      ELSE
!         HTC_tst=0.0d0
!         qflux_l=0.0d0
!      ENDIF   
!      tdiff=(tw-cell%tg(i))
!      IF(DABS(tdiff).gt.1.d-1)THEN
!      HTC_tg=qflux_g/(tw-cell%tg(i))
!      ELSE
!         HTC_tg=0.0d0
!         qflux_g=0.0d0
!      ENDIF   
!      HTC_tgp=0.0d0 
!
      END SUBROUTINE trans_film_boiling_HTC
!
