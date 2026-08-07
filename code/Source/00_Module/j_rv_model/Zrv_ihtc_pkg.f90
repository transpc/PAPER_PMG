!
MODULE Zrv_ihtc_pkg
!
      USE Zrv_ihtc_models 
      USE Zrv_htc         ! , ONLY: hil_big   
!
      IMPLICIT NONE
!
CONTAINS      
!
!--------------------------------------------------------------
!  HTC: bubbly(SHL)
!--------------------------------------------------------------
!  REAL FUNCTION rv_iht_bbl_shl(i,xag,xal,vg,vl,beta) RESULT(h)
   REAL FUNCTION rv_iht_bbl_shl(i,xag,vg,vl,beta) RESULT(h)
!
      USE Vol_DATA       , ONLY: cell
      USE Zconst2        , ONLY: hydraulicd
!      
      IMPLICIT NONE
!     
      INTEGER i
!     REAL(8) xag,xal,vg,vl,beta
      REAL(8) xag,vg,vl,beta
      REAL(8) vfg,vfg2,alphabub,Reb,DelTsl,dia_max
      REAL(8) cF1,cF2,cF3
      REAL(8) hpz,hlr    
      REAL(8),PARAMETER :: dia_min=0.005d0  
      REAL(8),PARAMETER :: We_crit=5.d0    
!
      DelTsl=cell%Ts(i)-cell%Tl(i)
      IF(xag.le.0.d0.and.DelTsl.ge.0.d0) THEN
         h=0.d0  !single-phase subcooled liquid
         RETURN
      ENDIF
!      
!.....Basic parameters
!
      alphabub=MAX(xag,1.d-5)
      dia_max=MIN(dia_min*alphabub**(1.d0/3.d0),hydraulicd(i))
      vfg=(vg-vl)*beta*beta                                                    
      vfg2=MAX(vfg*vfg,We_crit*cell%sigma(i)/cell%rhol(i)/dia_max)
      vfg=SQRT(vfg2)
!
!.....Reynolds number for bubbly flows            
!
      Reb=(1.d0-alphabub)*cell%rhol(i)*vfg*cell%dbb(i)/cell%lviscosl(i)
!      Reb=We_crit*cell%sigma(i)*(1.d0-alphabub)/cell%lviscosl(i)/vfg  !alternative
!
!.....Plesset-Zwick
!
      hpz=Plesset_Zwick(cell%dbb(i),DelTsl,cell%rhog(i),cell%rhol(i),   &
                        cell%hgsat(i)-cell%hlsat(i),cell%cpl(i),cell%lcondl(i))
!
!.....Lee-Ryley
!
      hlr=Lee_Ryley(cell%dbb(i),cell%lcondl(i),Reb,2.d0,0.74d0)
!
!.....Final HTC for bubbly SHL
!
      cF1=F1(alphabub)
      cF2=F2(alphabub)
      cF3=F3(xag,cell%quala(i),DelTsl)
!      
      h=MAX(hpz*beta,hlr)+0.4d0*ABS(vl)*cell%rhol(i)*cell%cpl(i)*cF1*cF2*cF3
!
   END FUNCTION rv_iht_bbl_shl
!
!--------------------------------------------------------------
!  HTC: bubbly(SCL)
!--------------------------------------------------------------
!  REAL FUNCTION rv_iht_bbl_scl(i,xag,xal,vg,vl,beta) RESULT(h)
   REAL FUNCTION rv_iht_bbl_scl(i,xag,vl,beta) RESULT(h)
!
      USE Vol_DATA       , ONLY: cell
!      
      IMPLICIT NONE
!     
      INTEGER i
!      
!     REAL(8) xag,xal,vg,vl,beta
      REAL(8) xag,vl,beta
      REAL(8) hfg,DelTsl,cF3,alphabub
!
      IF(xag.le.0.d0) THEN
         h=0.d0
         RETURN
      ENDIF
!
      alphabub=MAX(xag,1.d-5)
!
!     alphabub update for 'rc_iht_slg_scl_SB'       
      IF(beta.ne.1.d0) THEN
         alphabub=cell%alpha_bs(i)*beta
      ENDIF   
!
      hfg=cell%hgsat(i)-cell%hlsat(i)
      DelTsl=cell%Ts(i)-cell%Tl(i)
      cF3=F3(xag,cell%quala(i),DelTsl)
!
!     h=Unal(cell%dbb(i),hydraulicd(i),cell%rhog(i),cell%rhol(i),hfg,cell%p(i),vl,alphabub)
      h=Unal(cell%dbb(i),cell%rhog(i),cell%rhol(i),hfg,cell%p(i),vl,alphabub)
      h=h*cF3
!
   END FUNCTION rv_iht_bbl_scl   
!
!--------------------------------------------------------------
!  HTC: bubbly(SHG)
!-------------------------------------------------------------- 
!  REAL FUNCTION rv_iht_bbl_shg(i,xag,xal,vg,vl ) RESULT(h)
   REAL FUNCTION rv_iht_bbl_shg(i,xag) RESULT(h)
!   
      USE Vol_DATA       , ONLY: cell
!      
      IMPLICIT NONE
!     
      INTEGER i
!      
!     REAL(8) xag,xal,vg,vl
      REAL(8) xag
      REAL(8) DelTsg,cF6,cF7
!
      DelTsg=cell%Ts(i)-cell%Tg(i)
      cF6=F6(DelTsg)
      cF7=F7(xag)
!
      h=hig_big*cF6*cF7       
!
   END FUNCTION rv_iht_bbl_shg
!
!--------------------------------------------------------------
!  HTC: bubbly(SCG)
!-------------------------------------------------------------- 
   REAL FUNCTION rv_iht_bbl_scg(i,xag) RESULT(h)
!
      USE Vol_DATA       , ONLY: cell
!      
      IMPLICIT NONE
!     
      INTEGER i
!      
!     REAL(8) xag,xal,vg,vl
      REAL(8) xag
      REAL(8) DelTsg,cF6,cF7
!
      DelTsg=cell%Ts(i)-cell%Tg(i)
      cF6=F6(DelTsg)
      cF7=F7(xag)
!
      h=hig_big*cF6*cF7      
!
   END FUNCTION rv_iht_bbl_scg   
!
!--------------------------------------------------------------
!  HTC: slug SB(SHL)
!-------------------------------------------------------------- 
!  REAL FUNCTION rv_iht_slg_shl_SB(i,xag,xal,vg,vl,beta) RESULT(h)
   REAL FUNCTION rv_iht_slg_shl_SB(i,xag,vg,vl,beta) RESULT(h)
!
      IMPLICIT NONE
!     
      INTEGER i      
!      
!     REAL(8) xag,xal,vg,vl,beta
      REAL(8) xag,vg,vl,beta
!
      h=rv_iht_bbl_shl(i,xag,vg,vl,beta)
!      
   END FUNCTION rv_iht_slg_shl_SB     
!
!--------------------------------------------------------------
! HTC: slug TB(SHL)
!-------------------------------------------------------------- 
!  REAL FUNCTION rv_iht_slg_shl_TB(i,xag,xal,vg,vl,beta) RESULT(h)
   REAL FUNCTION rv_iht_slg_shl_TB(i,xag,vg,vl,beta) RESULT(h)
!
      USE Vol_DATA       , ONLY: cell
      USE Zconst2        , ONLY: hydraulicd 
!
      IMPLICIT NONE
!     
      INTEGER i 
!      
!     REAL(8) xag,xal,vg,vl,beta
      REAL(8) xag,vg,vl,beta
      REAL(8) dtsuph,qwall,qzubr,hbubl,qbubl,hwall,hfg,hshl_sb,hscl_t      
!
!      h=hil_big
!
!....yjm ref. MARS coding
      dtsuph=cell%tl(i)-cell%ts(i)
      qwall=0.0744d0*cell%p(i)**0.79d0*dtsuph**2.d0
!
      qzubr=2.18d4*cell%p(i)**0.34d0
      qzubr=MIN(4.d6,qzubr)
      qwall=MIN(qwall,qzubr)
!
      hfg=MAX(1.d-12,cell%hgsat(i)-cell%hlsat(i))
      hbubl=3.82d0*(cell%lcondl(i)/cell%dsb(i))*(cell%rhol(i)*cell%cpl(i)/(cell%rhog(i)*hfg))*dtsuph
      qbubl=hbubl*dtsuph
!
      qwall=qwall*4.d0/hydraulicd(i)
      qbubl=qbubl*cell%ia_slug_sb(i)
!
      hwall=qwall/MAX(1.d0,dtsuph)
      hbubl=qbubl/MAX(1.d0,dtsuph)
!
!     hshl_sb=rv_iht_slg_shl_SB(i,xag,xal,vg,vl,beta)
      hshl_sb=rv_iht_slg_shl_SB(i,xag,vg,vl,beta)
!     hscl_t=rv_iht_slg_scl_TB(i,xag,xal,vg,vl,beta)
      hscl_t=rv_iht_slg_scl_TB(i,vl)
!
      h=hscl_t*cell%ia_slug_tb(i)+MAX(hshl_sb*cell%ia_slug_sb(i),hbubl)+hwall
!      
   END FUNCTION rv_iht_slg_shl_TB   
!
!--------------------------------------------------------------
!  HTC: slug SB(SCL)
!-------------------------------------------------------------- 
!  REAL FUNCTION rv_iht_slg_scl_SB(i,xag,xal,vg,vl,beta) RESULT(h)
   REAL FUNCTION rv_iht_slg_scl_SB(i,xag,vl,beta) RESULT(h)
!
      IMPLICIT NONE
!     
      INTEGER i
!      
!     REAL(8) xag,xal,vg,vl,beta
      REAL(8) xag,vl,beta
!
!     h=rv_iht_bbl_scl(i,xag,vg,vl,beta)
      h=rv_iht_bbl_scl(i,xag,vl,beta)
!      
   END FUNCTION rv_iht_slg_scl_SB
!
!--------------------------------------------------------------
!  HTC: slug TB(SCL)
!-------------------------------------------------------------- 
!  REAL FUNCTION rv_iht_slg_scl_TB(i,xag,xal,vg,vl,beta) RESULT(h)
   REAL FUNCTION rv_iht_slg_scl_TB(i,vl) RESULT(h)
!   
      USE Vol_DATA  , ONLY: cell
      USE Zconst2   , ONLY: hydraulicd
!      
      IMPLICIT NONE
!     
      INTEGER i
!       
!     REAL(8) xag,xal,vg,vl,beta
      REAL(8) vl
      REAL(8) filmt,refilm,prfilm,nut,nul      
!    
!      h=Sieder_Tate(cell%lcondl(i),cell%rhol(i),cell%lviscosl(i),cell%cpl(i),   &
!                    hydraulicd(i),vg,vl)   
!
!                    cell%dsb(i),hydraulicd(i),vg,vl)   
!
!.....yjm ref. MARS coding /  vg=ug_o, vl=ul_o
      filmt=0.5d0*(1.d0-cell%alphag(i)**0.5d0)*hydraulicd(i)
      filmt=MAX(20.d-6,filmt)
      refilm=ABS(cell%alphal(i)*cell%rhol(i)*vl)*hydraulicd(i)/cell%lviscosl(i) 
      prfilm=cell%cpl(i)*cell%lviscosl(i)/cell%lcondl(i)
      nut=0.023d0*refilm**0.8d0*prfilm**0.4d0
      nul=4.d0*(1.d0+0.022d0*refilm**0.46d0)
      h=(cell%lcondl(i)/filmt)*MAX(nut,nul)
      h=MIN(3.2d5,h)
!         
   END FUNCTION rv_iht_slg_scl_TB   
!
!--------------------------------------------------------------
!  HTC: slug SB(SHG)
!--------------------------------------------------------------    
!  REAL FUNCTION rv_iht_slg_shg_SB(i,xag,xal,vg,vl,beta) RESULT(h)
   REAL FUNCTION rv_iht_slg_shg_SB(i) RESULT(h)
!
      USE Vol_DATA, ONLY: cell
!      
      IMPLICIT NONE
!     
      INTEGER i
!       
!     REAL(8) xag,xal,vg,vl,beta
      REAL(8) DelTs,cF6
!
      DelTs=cell%Ts(i)-cell%Tg(i)
      cF6=F6(DelTs)
!
      h=hig_big*cF6
!      
   END FUNCTION rv_iht_slg_shg_SB   
!
!--------------------------------------------------------------
!  HTC: slug TB(SHG)
!--------------------------------------------------------------    
!  REAL FUNCTION rv_iht_slg_shg_TB(i,xag,xal,vg,vl,beta) RESULT(h)
   REAL FUNCTION rv_iht_slg_shg_TB(i,vg,vl) RESULT(h)
!
      USE Vol_DATA        , ONLY: cell
      USE Zconst2         , ONLY: hydraulicd
!      
      IMPLICIT NONE
!     
      INTEGER i
!       
!     REAL(8) xag,xal,vg,vl,beta
      REAL(8) vg,vl
      REAL(8) Reyg,vfg
!
      vfg=ABS(vg-vl)
      Reyg=cell%rhog(i)*hydraulicd(i)*vfg/cell%lviscosg(i) 
!
      h=Lee_Ryley(hydraulicd(i),cell%lcondg(i),Reyg,2.2,0.82)
!      
   END FUNCTION rv_iht_slg_shg_TB 
!
!--------------------------------------------------------------
!  HTC: slug SB(SCG)
!-------------------------------------------------------------- 
!  REAL FUNCTION rv_iht_slg_scg_SB(i,xag,xal,vg,vl,beta) RESULT(h)
   REAL FUNCTION rv_iht_slg_scg_SB(i) RESULT(h)
!
      IMPLICIT NONE
!     
      INTEGER i
!       
!     REAL(8) xag,xal,vg,vl,beta
!
!     h=rv_iht_slg_shg_SB(i,xag,xal,vg,vl,beta)
      h=rv_iht_slg_shg_SB(i)
!      
   END FUNCTION rv_iht_slg_scg_SB
!
!--------------------------------------------------------------
!  HTC: slug TB(SCG)
!-------------------------------------------------------------- 
!  REAL FUNCTION rv_iht_slg_scg_TB(i,xag,xal,vg,vl,beta) RESULT(h)
   REAL FUNCTION rv_iht_slg_scg_TB(i) RESULT(h)
!
      USE Vol_DATA      , ONLY: cell
!      
      IMPLICIT NONE
!     
      INTEGER i
!       
!     REAL(8) xag,xal,vg,vl,beta
      REAL(8) DelTs,cF6
!
      DelTs=cell%Ts(i)-cell%Tg(i)
      cF6=F6(DelTs)
!
      h=hig_big*cF6
!      
   END FUNCTION rv_iht_slg_scg_TB   
!
!--------------------------------------------------------------
!  HTC: annular liquid film(SHL)
!-------------------------------------------------------------- 
!  REAL FUNCTION rv_iht_anm_shl_lf(i,xag,xal,vg,vl,beta ) RESULT(h)
   REAL FUNCTION rv_iht_anm_shl_lf RESULT(h)
!
      IMPLICIT NONE
!     
!     INTEGER i
!       
!     REAL(8) xag,xal,vg,vl,beta
!
      h=hil_big
!      
   END FUNCTION rv_iht_anm_shl_lf
!
!--------------------------------------------------------------
!  HTC: annular drplet(SHL)
!-------------------------------------------------------------- 
!  REAL FUNCTION rv_iht_anm_shl_drp(i,xag,xal,vg,vl,beta) RESULT(h)
   REAL FUNCTION rv_iht_anm_shl_drp(i) RESULT(h)
! 
      USE Vol_DATA, only: cell
!      
      IMPLICIT NONE
!     
      INTEGER i
!       
!     REAL(8) xag,xal,vg,vl,beta
      REAL(8) DelTs,cF12
      REAL(8) hfg,cF13
!
      DelTs=cell%Ts(i)-cell%Tl(i)
      cF12=F12(DelTs)    
!
      hfg=cell%hgsat(i)-cell%hlsat(i)      
      cF13=F13(DelTs,cell%cpl(i),hfg)
!      
      h=cell%lcondl(i)/MAX(cell%ddrp(i),1.d-8)*cF12*cF13
!      
   END FUNCTION rv_iht_anm_shl_drp      
!
!--------------------------------------------------------------
!  HTC: annular liquid film(SCL)
!-------------------------------------------------------------- 
!  REAL FUNCTION rv_iht_anm_scl_lf(i,xag,xal,vg,vl,beta) RESULT(h)
   REAL FUNCTION rv_iht_anm_scl_lf(i,vl) RESULT(h)
!
      USE Vol_DATA, ONLY: cell
!      
      IMPLICIT NONE
!     
      INTEGER i
!       
!     REAL(8) xag,xal,vg,vl,beta
      REAL(8) vl
!
      h=Theofanous(cell%rhol(i),cell%cpl(i),vl)
!      
   END FUNCTION rv_iht_anm_scl_lf
!
!--------------------------------------------------------------
!  HTC: annular drplet(SCL)
!-------------------------------------------------------------- 
!  REAL FUNCTION rv_iht_anm_scl_drp(i,xag,xal,vg,vl,beta) RESULT(h)
   REAL FUNCTION rv_iht_anm_scl_drp(i) RESULT(h)
!   
      USE Vol_DATA, ONLY: cell
!      
      IMPLICIT NONE
!     
      INTEGER i
!       
!     REAL(8) xag,xal,vg,vl,beta
      REAL(8) DelTs,hfg,cF13
!
      DelTs=cell%Ts(i)-cell%Tl(i)
      hfg=cell%hgsat(i)-cell%hlsat(i)      
      cF13=2.d0+7.d0*MIN(1.d0+cell%cpl(i)*MAX(DelTs,0.d0)/hfg,8.d0)
!
      h=cell%lcondl(i)/cell%ddrp(i)*cF13
!      
   END FUNCTION rv_iht_anm_scl_drp   
!
!--------------------------------------------------------------
!  HTC: annular liquid film(SHG)
!-------------------------------------------------------------- 
!  REAL FUNCTION rv_iht_anm_shg_lf(i,xag,xal,vg,vl,beta) RESULT (h)
   REAL FUNCTION rv_iht_anm_shg_lf(i,xag,vg,vl) RESULT (h)
!   
      USE Vol_DATA  , ONLY: cell
      USE Zconst2   , ONLY: hydraulicd
!      
      IMPLICIT NONE
!     
      INTEGER i
!       
!     REAL(8) xag,xal,vg,vl,beta
      REAL(8) xag,vg,vl
      REAL(8) vfg
!
      vfg=vg-vl
!
      h=Dittus_Boelter(hydraulicd(i),cell%lcondg(i),cell%rhog(i),cell%lviscosg(i),vfg,xag)
!      
   END FUNCTION rv_iht_anm_shg_lf
!
!--------------------------------------------------------------
!  HTC: annular droplet(SHG)
!-------------------------------------------------------------- 
!  REAL FUNCTION rv_iht_anm_shg_drp(i,xag,xal,vg,vl,beta) RESULT(h)
   REAL FUNCTION rv_iht_anm_shg_drp(i,xag,xal,vg,vl) RESULT(h)
!   
      USE Vol_DATA      , ONLY: cell
      USE Zconst2       , ONLY: hydraulicd,ggc
!      
      IMPLICIT NONE
!     
      INTEGER i
!       
!     REAL(8) xag,xal,vg,vl,beta
      REAL(8) xag,xal,vg,vl
      REAL(8) f11,gamma,gamma_s,rey_f,v_crit
      REAL(8) vfg,vfg_s,vfg_hat_2,vfg_hat,vfg_ss,Red
      REAL(8) alp_fd,alp_ad_s,alp_ad,alp_ef,alp_ff
      REAL(8) cF14,a_drp_weight
      REAL(8) dr
!
      xal=MAX(xal,1.d-8)
!
      alp_ad=1.d-4     
      alp_ef=MAX(alp_ad+alp_ad,MIN(2.d-3*cell%rhog(i)/cell%rhol(i),2.d-4))
      gamma=(xal-alp_ad)/(alp_ef-alp_ad)
      IF(xag>cell%alpha_sa(i).AND.xal<alp_ef) THEN
         gamma_s=gamma
      ELSE
         gamma_s=1.d0
      ENDIF
      rey_f=xal*cell%rhol(i)*ABS(vl)*hydraulicd(i)/cell%lviscosl(i) 
      dr=MAX(1.d-5,cell%rhol(i)-cell%rhog(i))
      v_crit=3.2d0*(MAX(cell%sigma(i),1.d-7)*ggc*dr/(cell%rhog(i)*cell%rhog(i)))**0.25       
      f11=gamma_s*MAX(0.d0,1.d0-rey_f**0.25*1.d-4)*EXP(-7.5d-5*(xag*vg/v_crit)**6)
      alp_ff=MAX(0.d0,xal*f11)
      vfg=vg-vl
      IF(xag>cell%alpha_sa(i).AND.xal<alp_ef) THEN
         alp_ad_s=alp_ad*gamma+1.d-5*(1.d0-gamma)
         vfg_s=vfg*(1.d0-f11*gamma)
      ELSE
         alp_ad_s=alp_ad
         vfg_s=vfg*(1.d0-f11)
      ENDIF
      alp_fd=MAX((xal-alp_ff)/(1.d0-alp_ff),alp_ad_s)
      IF(xal<1.d-6) THEN
         vfg_ss=vfg_s*xal*1.d6
      ELSE
         vfg_ss=vfg_s
      ENDIF
!      
      vfg_hat_2=MAX(1.5d0*cell%sigma(i),1.d-10)/(cell%rhog(i)*MIN(0.0025d0*alp_fd**(1.d0/3.d0),hydraulicd(i)))
      vfg_hat_2=MAX(vfg_ss*vfg_ss,vfg_hat_2)      
      vfg_hat=SQRT(MAX(0.d0,vfg_hat_2))
!      
!reflod model by LSJ 180123
!      IF(.not.reflod) THEN
         Red=MAX(1.d-10,1.5d0*cell%sigma(i))*(1.d0-alp_fd)**2.5d0/vfg_hat/cell%lviscosg(i)
!      ELSE
!         Red=rey_reflod 
!      ENDIF 
!      
      h=(2.d0+0.5d0*Red**0.5d0)*cell%lcondg(i)/cell%ddrp(i)
      IF(xal>=alp_ad_s) THEN
         a_drp_weight=1.d0
      ELSE   
         cF14=1.d0-5.d0*MIN(0.2d0,MAX(0.d0,cell%Ts(i)-cell%Tg(i)))
         a_drp_weight=(xal*cF14/alp_ad_s+1.d0-cF14)
      ENDIF
      h=h*a_drp_weight
!      
   END FUNCTION rv_iht_anm_shg_drp     
!
!--------------------------------------------------------------
!  HTC: annular liquid film(SCG)
!-------------------------------------------------------------- 
!  REAL FUNCTION rv_iht_anm_scg_lf(i,xag,xal,vg,vl,beta) RESULT(h)
   REAL FUNCTION rv_iht_anm_scg_lf(i) RESULT(h)
!   
      USE Vol_DATA, ONLY: cell
!      
      IMPLICIT NONE
!     
      INTEGER i
!       
!     REAL(8) xag,xal,vg,vl,beta
      REAL(8) DelTs,cF6
!
      DelTs=cell%Ts(i)-cell%Tg(i)
      cF6=F6(DelTs)
!
      h=hig_big*cF6 
!      
   END FUNCTION rv_iht_anm_scg_lf
!
!--------------------------------------------------------------
!  HTC: annular droplet(SCG)
!--------------------------------------------------------------
   REAL FUNCTION rv_iht_anm_scg_drp(i,xag,xal) RESULT(h)
!   
      USE Vol_DATA       , ONLY: cell
!      
      IMPLICIT NONE
!     
      INTEGER i
!       
      REAL(8) xag,xal
      REAL(8) DelTs,cF6,cF14
      REAL(8) alp_ad,alp_ef,alp_ad_s,a_drp_weight,gamma
!
!.....alp_ad_s for a_drp_weight
!
      alp_ad=1.d-4     
      alp_ef=MAX(alp_ad+alp_ad,MIN(2.d-3*cell%rhog(i)/cell%rhol(i),2.d-4))
      gamma=(xal-alp_ad)/(alp_ef-alp_ad)
      IF(xag>cell%alpha_sa(i).AND.xal<alp_ef) THEN
         alp_ad_s=alp_ad*gamma+1.d-5*(1.d0-gamma)
      ELSE
         alp_ad_s=alp_ad
      ENDIF      
!
!.....cF6 for h      
!
      DelTs=cell%Ts(i)-cell%Tg(i)
      cF6=F6(DelTs)
!
!.....weighted h
!
      h=hig_big*cF6
      IF(xal>=alp_ad_s) THEN
         a_drp_weight=1.d0
      ELSE   
         cF14=1.d0-5.d0*MIN(0.2d0,MAX(0.d0,DelTs))
         a_drp_weight=(xal*cF14/alp_ad_s+1.d0-cF14)
      ENDIF
      h=h*a_drp_weight   
!         
   END FUNCTION rv_iht_anm_scg_drp    
!
!--------------------------------------------------------------
!  HTC: inverted annular gas film(SHL)
!-------------------------------------------------------------- -
!  REAL FUNCTION rv_iht_invann_shl_gf(i,xag,xal,vg,vl,beta) RESULT(h)
   REAL FUNCTION rv_iht_invann_shl_gf RESULT(h)
!   
      IMPLICIT NONE
!     
!     INTEGER i
!       
!     REAL(8) xag,xal,vg,vl,beta
!
      h=hil_big
!         
   END FUNCTION rv_iht_invann_shl_gf
!
!--------------------------------------------------------------
!  HTC: inverted annular small bubble(SHL)
!-------------------------------------------------------------- -
!  REAL FUNCTION rv_iht_invann_shl_sb(i,xag,xal,vg,vl,beta) RESULT(h)
   REAL FUNCTION rv_iht_invann_shl_sb(i,xag,vl,beta) RESULT(h)
!   
      USE Vol_DATA        , ONLY: cell
      USE Zconst2         , ONLY: hydraulicd
      USE Zvector         , ONLY: ul_o,ug_o
!      
      IMPLICIT NONE
!     
      INTEGER i
!       
!     REAL(8) xag,xal,vg,vl,beta
      REAL(8) xag,vl,beta
      REAL(8) DelTsl
      REAL(8) cF17,cF1,cF2,cF3
      REAL(8) We_crit,dia_min,dia_max,vfg,vfg2,Reb,hpz,hlr
      REAL(8) alphaian,alphab,alphabub
!      
!.....Basic Parameters
!
      DelTsl=cell%Ts(i)-cell%Tl(i)
      alphaian=xag
      cF17=F17(i,alphaian)
      alphab=cF17*alphaian
      alphabub=MAX((alphaian-alphab)/(1.d0-alphab),1.d-7) 
      We_crit=5.d0
      dia_min=0.005d0      
      dia_max=MIN(dia_min*alphabub**(1.d0/3.d0),hydraulicd(i))
      vfg=(ug_o(i)-ul_o(i))*beta*beta                                                    
      vfg2=MAX(vfg*vfg,We_crit*cell%sigma(i)/cell%rhol(i)/dia_max)
      vfg=SQRT(vfg2)
!
!.....Reynolds number for bubbly flows            
!
!      IF(.not.reflod) THEN
         Reb=We_crit*cell%sigma(i)*(1.d0-alphabub)*vfg/cell%lviscosl(i)
!      ELSE
!         Reb=rey_reflod
!      ENDIF 
!
!.....Plesset-Zwick
!
      hpz=Plesset_Zwick(cell%dbb(i),DelTsl,cell%rhog(i),cell%rhol(i),    &
                        cell%hgsat(i)-cell%hlsat(i),cell%cpl(i),cell%lcondl(i))
!
!.....Lee-Ryley
!
      hlr=Lee_Ryley(cell%dbb(i),cell%lcondl(i),Reb,2.d0,0.74d0)
!
!.....Final HTC for bubbly SHL
!
      cF1=F1(alphabub)
      cF2=F2(alphabub)
      cF3=F3(xag,cell%quala(i),DelTsl)
      h=MAX(hpz*beta,hlr)+0.4d0*ABS(vl)*cell%rhol(i)*cell%cpl(i)*cF1*cF2*cF3
!         
   END FUNCTION rv_iht_invann_shl_sb
!
!--------------------------------------------------------------
!  HTC: inverted annular gas film(SCL)
!-------------------------------------------------------------- -
!  REAL FUNCTION rv_iht_invann_scl_gf(i,xag,xal,vg,vl,beta) RESULT(h)
   REAL FUNCTION rv_iht_invann_scl_gf(i,xag,vg,vl) RESULT(h)
!   
      USE Vol_DATA        , ONLY: cell
      USE Zconst2         , ONLY: hydraulicd
!      
      IMPLICIT NONE
!     
      INTEGER i
!       
!     REAL(8) xag,xal,vg,vl,beta
      REAL(8) xag,vg,vl
      REAL(8) vfg,alphaian,alpha,DelTsl,cF3    
!
      vfg=vl-vg
      alphaian=xag
      alpha=1.d0-alphaian
      DelTsl=cell%Ts(i)-cell%Tl(i)      
!      
      h=Dittus_Boelter(hydraulicd(i),cell%lcondl(i),cell%rhol(i),cell%lviscosl(i),vfg,alpha)
      cF3=F3(xag,cell%quala(i),DelTsl)
      h=h*cF3
!         
   END FUNCTION rv_iht_invann_scl_gf   
!
!--------------------------------------------------------------
!  HTC: inverted annular small bubble(SCL)
!-------------------------------------------------------------- -
!  REAL FUNCTION rv_iht_invann_scl_sb(i,xag,xal,vg,vl,beta) RESULT(h)
   REAL FUNCTION rv_iht_invann_scl_sb(i,xag,vl) RESULT(h)
!   
      USE Vol_DATA        , ONLY: cell
!      
      IMPLICIT NONE
!     
      INTEGER i
!       
!     REAL(8) xag,xal,vg,vl,beta
      REAL(8) xag,vl
      REAL(8) alphaian,alphab,alphabub
      REAL(8) cF3,cF17,DelTsl,hfg
!      
      alphaian=xag         
      cF17=F17(i,alphaian)
      alphab=cF17*alphaian
      alphabub=MAX((alphaian-alphab)/(1.d0-alphab),1.d-7)   
      DelTsl=cell%Ts(i)-cell%Tl(i)        
!
      hfg=cell%hgsat(i)-cell%hlsat(i)
      DelTsl=cell%Ts(i)-cell%Tl(i)
      cF3=F3(xag,cell%quala(i),DelTsl)
!
!     h=Unal(cell%dbb(i),hydraulicd(i),cell%rhog(i),cell%rhol(i),hfg,cell%p(i),vl,alphabub)
      h=Unal(cell%dbb(i),cell%rhog(i),cell%rhol(i),hfg,cell%p(i),vl,alphabub)
      h=h*cF3  
!         
   END FUNCTION rv_iht_invann_scl_sb     
!
!--------------------------------------------------------------
!  HTC: inverted annular gas film(SHG)
!-------------------------------------------------------------- -
!  REAL FUNCTION rv_iht_invann_shg_gf(i,xag,xal,vg,vl,beta) RESULT(h)
   REAL FUNCTION rv_iht_invann_shg_gf(i,xag) RESULT(h)
!   
      USE Vol_DATA        , ONLY: cell
      USE Zconst2         , ONLY: hydraulicd
!      
      IMPLICIT NONE
!     
      INTEGER i
!       
!     REAL(8) xag,xal,vg,vl,beta
      REAL(8) xag
      REAL(8) cF15,cF17,cF19,cF20,DelTsg,alphab,alphaian
!
      DelTsg=cell%Ts(i)-cell%Tg(i)
      alphaian=xag         
      cF17=F17(i,alphaian)
      alphab=cF17*alphaian
      cF15=SQRT(1.d0-alphab)  
      cF19=2.5d0-DelTsg*(0.2d0-0.1d0*DelTsg)
      cF20=0.5d0*MAX(1.d0-cF15,0.04d0)
!      
      h=cell%lcondg(i)/hydraulicd(i)*cF19/cF20
!         
   END FUNCTION rv_iht_invann_shg_gf   
!
!--------------------------------------------------------------
!  HTC: inverted annular small bubble(SHG)
!-------------------------------------------------------------- -
!  REAL FUNCTION rv_iht_invann_shg_sb(i,xag,xal,vg,vl,beta) RESULT (h)
   REAL FUNCTION rv_iht_invann_shg_sb(i,xag) RESULT (h)
!   
      USE Vol_DATA        , ONLY: cell
!      
      IMPLICIT NONE
!     
      INTEGER i
!       
!     REAL(8) xag,xal,vg,vl,beta
      REAL(8) xag
      REAL(8) DelTsg,cF6,cF7    
!
      DelTsg=cell%Ts(i)-cell%Tg(i)
      cF6=F6(DelTsg) 
      cF7=F7( xag )
!
      h=hig_big*cF6*cF7            
      h=h*cF6
!      
   END FUNCTION rv_iht_invann_shg_sb 
!
!--------------------------------------------------------------
!  HTC: inverted slug annular(SHL)
!-------------------------------------------------------------- -
   REAL FUNCTION rv_iht_invslg_shl_ann(i) RESULT(h)
!   
      USE Vol_DATA        , ONLY: cell
      USE Zconst2         , ONLY: hydraulicd
!      
      IMPLICIT NONE
!     
      INTEGER i
!
      REAL(8) DelTs,cF12,cF13,hfg    
!
      DelTs=cell%Ts(i)-cell%Tl(i)
      hfg=cell%hgsat(i)-cell%hlsat(i)  
       
      cF12=F12(DelTs)
      cF13=F13(DelTs,cell%cpl(i),hfg)
!      
      h=cell%lcondl(i)/hydraulicd(i)*cF12*cF13
!      
   END FUNCTION rv_iht_invslg_shl_ann    
!
!--------------------------------------------------------------
!  HTC: inverted slug droplet(SHL)
!-------------------------------------------------------------- -
   REAL FUNCTION rv_iht_invslg_shl_drp(i) RESULT(h)
!   
      USE Vol_DATA        , ONLY: cell
!      
      IMPLICIT NONE
!     
      INTEGER i
!
      REAL(8) DelTs,cF12,cF13,hfg    
!
      DelTs=cell%Ts(i)-cell%Tl(i)
      hfg=cell%hgsat(i)-cell%hlsat(i)  
      cF12=F12(DelTs)
      cF13=F13(DelTs,cell%cpl(i),hfg)
!      
      h=cell%lcondl(i)/cell%ddrp(i)*cF12*cF13
!      
   END FUNCTION rv_iht_invslg_shl_drp
!
!--------------------------------------------------------------
!  HTC: inverted slug annular(SCL)
!-------------------------------------------------------------- -
   REAL FUNCTION rv_iht_invslg_scl_ann(i) RESULT(h)
!   
      USE Vol_DATA        , ONLY: cell
      USE Zconst2         , ONLY: hydraulicd
!      
      IMPLICIT NONE
!     
      INTEGER i
      REAL(8) DelTs,cF13,hfg    
!
      DelTs=cell%Ts(i)-cell%Tl(i)
      hfg=cell%hgsat(i)-cell%hlsat(i)  
      cF13=F13(DelTs,cell%cpl(i),hfg)
!      
      h=cell%lcondl(i)/hydraulicd(i)*cF13
!      
   END FUNCTION rv_iht_invslg_scl_ann    
!
!--------------------------------------------------------------
!  HTC: inverted slug droplet(SCL)
!-------------------------------------------------------------- -
   REAL FUNCTION rv_iht_invslg_scl_drp(i) RESULT(h)
!   
      USE Vol_DATA        , ONLY: cell
!      
      IMPLICIT NONE
!     
      INTEGER i
      REAL(8) DelTs,cF13,hfg    
!
      DelTs=cell%Ts(i)-cell%Tl(i)
      hfg=cell%hgsat(i)-cell%hlsat(i)  
      cF13=F13(DelTs,cell%cpl(i),hfg)
!      
      h=cell%lcondl(i)/cell%ddrp(i)*cF13
!      
   END FUNCTION rv_iht_invslg_scl_drp
!
!--------------------------------------------------------------
!  HTC: inverted slug annular(SHG)
!--------------------------------------------------------------
!  REAL FUNCTION rv_iht_invslg_shg_ann(i,xag,xal,vg,vl) RESULT(h)
   REAL FUNCTION rv_iht_invslg_shg_ann(i,xag) RESULT(h)
!   
      USE Vol_DATA        , ONLY: cell
      USE Zconst2         , ONLY: hydraulicd
!      
      IMPLICIT NONE
!     
      INTEGER i
!      
!     REAL(8) xag,xal,vg,vl
      REAL(8) xag
      REAL(8) DelTs,cF19,cF22
!
      DelTs=cell%Ts(i)-cell%Tg(i)
      cF19=F19(DelTs)
      cF22=F22(xag)
!      
      h=cell%lcondg(i)/hydraulicd(i)*cF19/cF22
!      
   END FUNCTION rv_iht_invslg_shg_ann    
!
!--------------------------------------------------------------
!  HTC: inverted slug droplet(SHG)
!-------------------------------------------------------------- -
!  REAL FUNCTION rv_iht_invslg_shg_drp(i,xag,xal,vg,vl) RESULT(h)
   REAL FUNCTION rv_iht_invslg_shg_drp(i,xag,vg,vl) RESULT(h)
!   
      USE Vol_DATA        , ONLY: cell
!      
      IMPLICIT NONE
!     
      INTEGER i
!      
!     REAL(8) xag,xal,vg,vl
      REAL(8) xag,vg,vl
      REAL(8) coef,Redrp,vfg,cF21
!
      cF21=F21(i,xag)
      vfg=(vg-vl)*cF21*cF21
      vfg=SQRT(vfg*vfg)
!
!reflod model by LSJ 180123
!      IF(.not.reflod) THEN
         Redrp=cell%rhog(i)*vfg*cell%ddrp(i)/cell%lviscosg(i)
!      ELSE
!         Redrp=rey_reflod
!      ENDIF 
!      
      coef=(2.d0+0.5d0*SQRT(Redrp))
!      
      h=cell%lcondg(i)/cell%ddrp(i)*coef
!      
   END FUNCTION rv_iht_invslg_shg_drp
!
!--------------------------------------------------------------
!  HTC: mist(SHL)
!-------------------------------------------------------------- -
!  REAL FUNCTION rv_iht_mist_shl(i,xag,xal,vg,vl) RESULT(h)
   REAL FUNCTION rv_iht_mist_shl(i,xal) RESULT(h)
!   
      USE Vol_DATA        , ONLY: cell
!      
      IMPLICIT NONE
!     
      INTEGER i
!      
!     REAL(8) xag,xal,vg,vl
      REAL(8) xal
      REAL(8) cF12,cF13,cF23,DelTsl,hfg
!
      DelTsl=cell%Ts(i)-cell%Tl(i)
      hfg=cell%hgsat(i)-cell%hlsat(i)  
!      
      cF12=F12(DelTsl)
      cF13=F13(DelTsl,cell%cpl(i),hfg)
      cF23=F23(xal)
!    
      h=cell%lcondl(i)/cell%ddrp(i)*cF12*cF13*cF23
!      
   END FUNCTION rv_iht_mist_shl
!
!--------------------------------------------------------------
!  HTC: mist(SCL)
!-------------------------------------------------------------- -
!  REAL FUNCTION rv_iht_mist_scl(i,xag,xal,vg,vl) RESULT(h)
   REAL FUNCTION rv_iht_mist_scl(i,xal) RESULT(h)
!   
      USE Vol_DATA        , ONLY: cell
!      
      IMPLICIT NONE
!     
      INTEGER i
!      
!     REAL(8) xag,xal,vg,vl
      REAL(8) xal
      REAL(8) cF13,cF23,DelTsl,hfg
!
      DelTsl=cell%Ts(i)-cell%Tl(i)
      hfg=cell%hgsat(i)-cell%hlsat(i)  
!      
      cF13=F13(DelTsl,cell%cpl(i),hfg)
      cF23=F23(xal)
!          
      h=cell%lcondl(i)/cell%ddrp(i)*cF13*cF23
!      
   END FUNCTION rv_iht_mist_scl   
!
!--------------------------------------------------------------
!  HTC: mist(SHG)
!-------------------------------------------------------------- -
!  REAL FUNCTION rv_iht_mist_shg(i,xag,xal,vg,vl) RESULT(h)
   REAL FUNCTION rv_iht_mist_shg(i,xal,vg,vl) RESULT(h)
!   
      USE Vol_DATA        , ONLY: cell
      USE Zwall_HTC       , ONLY: reflod
!      
      IMPLICIT NONE
!     
      INTEGER i
!      
!     REAL(8) xag,xal,vg,vl
      REAL(8) xal,vg,vl
      REAL(8) cF24,coef,vfg,Redrp,alphadrp,We_crit,DelTsg,blow_f
!
      DelTsg=cell%Ts(i)-cell%Tg(i)
!
!      We_crit=6.d0
      IF(reflod)THEN
         We_crit=4.d0
      ELSE
         We_crit=1.5d0
      ENDIF
      
      cF24=F24(DelTsg,xal)
      alphadrp=MAX(xal,1.d-4)
      vfg=MAX(1.d-6,vg-vl)
!
! yjm  MARS source interphaseHTC Line 1892
      blow_f=(cell%hg(i)-cell%hgsat(i))/(cell%hgsat(i)-cell%hlsat(i))
      blow_f=1.d0/(1.d0+blow_f)**0.7d0      
!
!reflod model by LSJ 180123      
!      IF(.not.reflod) THEN
         Redrp=MAX(We_crit*cell%sigma(i),1.d-10)*(1.d0-alphadrp)**2.5/cell%lviscosg(i)/vfg
!      ELSE
!         Redrp=rey_reflod      
!      ENDIF   
!      
      coef=(2.d0+0.5d0*SQRT(Redrp))
!
      IF(xal.gt.0.d0) THEN
!         h=cell%lcondg(i)/cell%ddrp(i)*coef*cF24      
         h=cell%lcondg(i)/cell%ddrp(i)*coef*cF24*blow_f 
      ELSE
         h=0.d0
      ENDIF   
!      
   END FUNCTION rv_iht_mist_shg
!
!--------------------------------------------------------------
! HTC: mist(SCG)
!-------------------------------------------------------------- -
!  REAL FUNCTION rv_iht_mist_scg(i,xag,xal,vg,vl) RESULT(h)
   REAL FUNCTION rv_iht_mist_scg(i,xal) RESULT(h)
!   
      USE Vol_DATA        , ONLY: cell
!      
      IMPLICIT NONE
!     
      INTEGER i
!      
!     REAL(8) xag,xal,vg,vl
      REAL(8) xal
      REAL(8) DelTsg,cF24,cF6      
      REAL(8), PARAMETER :: p_triple=611.657d0      
!
      DelTsg=cell%Ts(i)-cell%Tg(i)
! 
      cF24=F24(DelTsg,xal)
      cF6=F6(DelTsg)
 !
      IF(xal.eq.0.d0.and.(cell%pps_o(i).lt.p_triple)) THEN
         h=0.d0
      ELSE
         h=hig_big*cF6*cF24
      ENDIF
!
   END FUNCTION rv_iht_mist_scg
!
!--------------------------------------------------------------
!  HTC: vertically stratified(SHL)
!--------------------------------------------------------------
!  REAL FUNCTION rv_iht_vstrat_shl(i,xag,xal,vg,vl) RESULT(h)
   REAL FUNCTION rv_iht_vstrat_shl(i) RESULT(h)
!   
      USE Vol_DATA        , ONLY: cell
      USE Zconst2         , ONLY: hydraulicd,ggc      
!      
      IMPLICIT NONE
!     
      INTEGER i
!      
!     REAL(8) xag,xal,vg,vl
      REAL(8) DelTsl,Nu,Gr,Pr,hyd,bet
!
      DelTsl=MAX(0.1d0,ABS(cell%Ts(i)-cell%Tl(i)))
      hyd=hydraulicd(i)
      bet=MAX( cell%betal(i),1.d-5 )
      Gr=ggc*bet*cell%rhol(i)*cell%rhol(i)*hyd*hyd*hyd*DelTsl/cell%lviscosl(i)/cell%lviscosl(i)
      Pr=cell%lviscosl(i)*cell%cpl(i)/cell%lcondl(i)
      Nu=0.27d0*(Gr*Pr)**0.25
!
      h=Nu*cell%lcondl(i)/hyd
!
   END FUNCTION rv_iht_vstrat_shl
!
!--------------------------------------------------------------
!  HTC: vertically stratified(SHG)
!-------------------------------------------------------------- -
!  REAL FUNCTION rv_iht_vstrat_shg(i,xag,xal,vg,vl) RESULT(h)
   REAL FUNCTION rv_iht_vstrat_shg(i) RESULT(h)
!   
      USE Vol_DATA        , ONLY: cell
      USE Zconst2         , ONLY: hydraulicd,ggc      
!      
      IMPLICIT NONE
!     
      INTEGER i
!      
!     REAL(8) xag,xal,vg,vl
      REAL(8) DelTsg,Nu,Gr,Pr,hyd,bet
!
      DelTsg=MAX(0.1d0,ABS(cell%Ts(i)-cell%Tg(i)))
      hyd=hydraulicd(i)
      bet=MAX( cell%betag(i),1.d-5 )
      Gr=ggc*bet*cell%rhog(i)*cell%rhog(i)*hyd*hyd*hyd*DelTsg/cell%lviscosg(i)/cell%lviscosg(i)
      Pr=cell%lviscosg(i)*cell%cpg(i)/cell%lcondg(i)
      Nu=0.27d0*(Gr*Pr)**0.25
!
      h=Nu*cell%lcondg(i)/hyd
!
   END FUNCTION rv_iht_vstrat_shg
!        
END MODULE Zrv_ihtc_pkg
   
