!
      SUBROUTINE trans_film_reflood_HTC_rod(i,mode,quala,hd,ug,p,alphag,alphal,rhog,rhol,tg,tl,vg,vl,condg,condl,visg,visl,cpg,cpl,betag,betal)
!
!     This routine calculates wall heat transfer coefficient and heat flux for transition boiling and film boiling
!     when reflood option is ON
!    
      USE VOL_DATA      
      USE STM_TBL_cupid , ONLY: st_tbl,             &
                                nt,np,ns,ns2,ndxstd
      USE Zmpi          , ONLY: jperm         
      USE Zparam        , ONLY: pi          
      USE Ztimecon      , ONLY: time      
      USE Zwall_HTC     , ONLY: qflux_t,qflux_l,qflux_g,HTC_t,HTC_tl,HTC_tst,HTC_tg,HTC_tgp, &
                                f_direc,c_direc1,c_direc2,dia_rod,                          &
                                mflux_tota,vfg,                                             &
                                hfg,chf,rho,cond,viscos,cp,beta,sigma,                      &
                                tw,tsat_t,dt_sat,zqf,zqf_min,zqf_top,reflood
      USE Zvector       , ONLY: vg_n
!
      IMPLICIT NONE
!      
      INTEGER i,it,mode
!
      REAL(8) gravity
      REAL(8) prop(36),drho
      REAL(8) HTC_brom,l_beren      
      REAL(8) qflux_tb_l,qflux_tb_g,qflux_tb,qflux_rad_l
      REAL(8) qflux_fb_l,qflux_fb_g,qflux_fb
      REAL(8) HTC_tb_l,HTC_tb_g,HTC_fb_l,HTC_fb_g
      REAL(8) HTC_fb_min,HTC_fr,HTC_fbb,HTC_fbb_top,HTC_tb_top
      
      REAL(8) dts_chf,dtw_chf,h_max,d,tw_chf              !for reflood modification
      REAL(8) d_max,d_ave,af_rad,e_f,e_g,e_w,temp,r1,r2,r3,F_wf    !for radiation HTC
      REAL(8) reygas,kby,xmf
!
!.....Additional connectivity input arguments - jrlee
      REAL(8) quala,hd,ug
      REAL(8) p
      REAL(8) alphag,rhog,condg,visg,cpg,betag,tg,vg
      REAL(8) alphal,rhol,condl,visl,cpl,betal,tl,vl
!      
      REAL(8) L1,L2,L3
      REAL(8) ucrit
      PARAMETER (L1=0.1d0)
      PARAMETER (L2=0.2d0)
!     PARAMETER (L3=3.66d0)
      PARAMETER (L3=10.d0)
!
      LOGICAL err
!      
      DATA gravity/9.81d0/
      prop(:)=0.d0     
      qflux_tb_l=0.0D0
!
     !dh=hydraulicd(i)      
!
!.....Assign basic properties
!
      prop(2)=p
      cp=cpg
      beta=betag
      IF(quala.gt.1.0d-9)THEN
         rho=rhog
         viscos=visg
         cond=condg
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
         IF(prop(9).eq.0.d0)rho=rhog
         CALL thcond_cupid_rv(prop(1),rho,cond)                                         
         CALL viscos_rv(prop(1),rho,tsat_t,viscos)      
      ENDIF  
      drho=rhol-rho
!
!-----------------------------------------------------------------------------------------------------------
!
!......Transition boiling for Reflood (Weismann correlation <- Chen correlation)
!       
      IF(zqf.lt.L2)THEN
         tw_chf=tsat_t+40.0d0
         dts_chf=DMAX1(3.0d0,DMIN1(40.0d0,dt_sat))
         dtw_chf=DMAX1(0.0d0,tw-tw_chf)
         h_max=0.5d0*chf/dts_chf
         HTC_tb_l=h_max*DEXP(-0.02d0*dtw_chf)+4500.d0*(mflux_tota/67.8d0)**0.2d0*DEXP(-0.012d0*dtw_chf) 
         HTC_tb_l=DMIN1(HTC_tb_l,h_max)
         IF(zqf.ge.L1)THEN
            HTC_tb_l=((zqf-L1)*1.d-4+(L2-zqf)*HTC_tb_l)*L3
         ENDIF
         qflux_tb_l=HTC_tb_l*dt_sat
      ELSE
         HTC_tb_l=1.d-4         
      ENDIF
!-------------------------------------Film Boiling Modeling----------------------------------------------
!
!.....Minimum film boiling HTC using Bromley correlation
!
!     drho=DMAX1(1.0d-7,drho)
      drho=DMAX1(1.0d-7,rhol-rhog)
      l_beren=1.58323349d0*(sigma/(gravity*drho))**0.125     !(2*pi)=1.58323349d0**0.25              
!     HTC_brom=0.62d0*(gravity*cond**3*cell%rhog(i)*drho*(hfg+0.5d0*cp*dt_sat)/(viscos*dt_sat))**0.25d0
      HTC_brom=0.62d0*(gravity*cond**3*rho*drho*(hfg+0.5d0*cp*dt_sat)/(viscos*dt_sat))**0.25d0
      HTC_fb_l=HTC_brom/l_beren  !**0.25d0
      HTC_fb_min=HTC_fb_l*DMAX1(0.0d0,0.999d0-alphag)
!
!.....Forslund-Rohsenow correlation
!      
      d=3.0d0*sigma/(cell%rhog(i)*DMAX1(0.01,vfg**2))  ! yjm
      d=DMIN1(3.d-3,DMAX1(1.d-4,d))
      HTC_fr=0.1d0*pi*(6.0d0*(0.999d0-alphag)/pi)**0.66666666d0
      HTC_fr=HTC_fr*(gravity*rho*rhol*hfg*cond**3/(dt_sat*viscos*d*(pi/6.d0)**0.33333333d0))**0.25d0   !!!cyj mod. bug fix
!
      reygas=DABS(alphag*rhog*ug)*hd/visg
      kby=(DMAX1((reygas-4000.0d0),0.01d0)/1.0d5)**0.6d0
      HTC_fr=HTC_fr*kby       
!
!.....Modified Bromley correlation for film boiling
!      
      zqf_min=DMAX1(0.005d0,zqf)
      HTC_fbb=DMAX1(0.0d0,(1400.d0-1880.d0*zqf_min))*DMIN1(0.5d0,0.999d0-alphag)+DSQRT(1.0d0-alphag)*HTC_brom/zqf_min**0.25d0  ! yjm
!
      IF(alphag.ge.0.9d0)THEN
         HTC_fbb=HTC_fr
      ELSEIF(alphag.ge.0.6d0)THEN
         xmf=(0.9d0-alphag)/(0.9d0-0.6d0)
         xmf=xmf*(2.0d0-xmf)
         HTC_fbb=xmf*HTC_fbb+(1.0d0-xmf)*HTC_fr
      ELSE
         HTC_fbb=HTC_fbb
      ENDIF      
!
!.....Select the maximum HTC among three HTCs (Minimum by original Bronley, Forslund, and Modified Bromley)
!      
!     HTC_fb_l=DMAX1(HTC_fb_min,HTC_fr,HTC_fbb)      
      HTC_fb_l=DMAX1(HTC_fb_min,HTC_fbb)  ! modified MARS
!
!.....Radiation heat transfer to liquid phase (ignore to gas phase)
!     
      vfg=DMAX1(0.005d0,vfg)
      d_max=DSQRT(alphal)*dia_rod
      d_ave=7.5d0*sigma/(rho*vfg*vfg)               !We=7.5
      d_ave=DMAX1(1.0d-5,DMIN1(d_max,d_ave))
      af_rad=1.11d0*alphal/d_ave
      e_f=DMIN1(0.75d0,(1.0d0-DEXP(-0.9d0*dia_rod*af_rad)))
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
     !CALL single_phase_HTC(i,mode)
      CALL single_phase_HTC_rod(i,mode,alphag,alphal,rhog,rhol,tg,tl,vg,vl,condg,condl,visg,visl,cpg,cpl,betag,betal)
      HTC_fb_g=HTC_tg
!
!.....Interpolation HTC_tb_g for
!      
      HTC_fb_g=HTC_fb_g*cell%alphag(i)
!
!.....Top quenching logic
!
      HTC_tb_g=HTC_fb_g 
      IF(reflood.eq.1.and.zqf_top.lt.L2)THEN
         tw_chf=tsat_t+40.0d0
         dts_chf=DMAX1(1.0d0,DMIN1(40.0d0,dt_sat))
         dtw_chf=DMAX1(0.0d0,tw-tw_chf)
!
!        Top quench Delay logic. (by B.D.Chung) ---------------------------
! 
!        i4 test-reflood bdchung
!        Do not allow top-down quenching if liquid velocity is 
!        upward and vapor velocty is high enough to carry-over the downward liquid film flow 
!        Helmholz instability condition ucrit
!
         ucrit=3.2d0*sqrt(sqrt(9.81d0*cell%sigma(i)*drho)/cell%rhog(i))
!         IF(.TRUE.) then  !Additional modification by DS Choi
!             !critical vg using wallis correlation (slope 1.0)                 
!             !vgcrit=9.81d0*v_da(iv)%DiaEquiv(1)*(v_da(iv)%Rhof-v_da(iv)%Rhog)/v_da(iv)%Rhog 
!             vgcrit=9.81d0*hydraulicd(1)*drho/cell%rhog(i)
!             vgcrit=sqrt(vgcrit)/dmax1(cell%alphag(i),0.01d0) 
!             !limit critical velocity condition for CCFL
!             ucrit=dmin1(ucrit,vgcrit)
!         ENDIF
         IF(vg_n(i,f_direc).gt.0.0d0) then 
            IF(vg_n(i,f_direc).gt.ucrit) chf=0.0d0
         ENDIF
!        Top quench Delay logic. End (by B.D.Chung) ------------------
!
         h_max=0.5d0*chf/dts_chf
         HTC_tb_l=h_max*DEXP(-0.05d0*dtw_chf)+4500.d0*(mflux_tota*0.01474926d0)**0.2d0*DEXP(-0.012d0*dtw_chf)   !1/67.8 = 0.01474926 from Wiesmann
         HTC_tb_l=DMIN1(HTC_tb_l,h_max)
         IF(zqf_top.le.L1)THEN
            HTC_tb_top=HTC_tb_l
         ELSE
            HTC_tb_top=((zqf_top-L1)*DMAX1(1.d-4,HTC_tb_l)+(L2-zqf_top)*HTC_tb_l)*L3
         ENDIF
         qflux_tb_l=HTC_tb_top*dt_sat  
!
         zqf_top=DMAX1(0.005d0,zqf_top)
         HTC_fbb_top=(600.d0-5000.d0*zqf_top)*DMIN1(0.5d0,0.999d0-cell%alphag(i))+  &
                     DSQRT(1.0d0-cell%alphag(i))*HTC_brom/zqf_top**0.25d0
!
!........Select the maximum HTC among three HTCs (Top QF HTC and Bot QF HTC)
!            
         HTC_fb_l=DMAX1(HTC_fbb_top,HTC_fr,HTC_fb_l)
      ENDIF
!
!.....Save HTC and heat flux
!
!     Transition boiling
      HTC_tb_l=qflux_tb_l/dt_sat
      HTC_tb_g=HTC_fb_g
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
      END SUBROUTINE trans_film_reflood_HTC_rod
