!------------------------------------------------------------------------
!     Interfacial heat transfer: bubbly 
!------------------------------------------------------------------------
!     SUBROUTINE rv_iht_bbl(i,xag,xal,vg,vl,htc_Hil,htc_Hig)
      SUBROUTINE rv_iht_bbl(i,xag,vg,vl,htc_Hil,htc_Hig)
!
      USE Vol_DATA        , ONLY: cell
      USE Zrv_htc         , ONLY: HTC_MIN
      USE Zrv_ihtc_pkg    , ONLY: rv_iht_bbl_shl,rv_iht_bbl_scl,umbrella_cap, &
                                  pow_interp,rv_iht_bbl_shg,rv_iht_bbl_scg
!      
      IMPLICIT NONE
! 
      INTEGER i          
!      
!     REAL(8) xag,xal,vg,vl      
      REAL(8) xag,vg,vl      
      REAL(8) htc_Hig,htc_Hil
      REAL(8) DelTsg,DelTsl
      REAL(8) Hshl,Hscl,Hshg,Hscg
!
      DelTsl=cell%Ts(i)-cell%Tl(i)
      DelTsg=cell%Ts(i)-cell%Tg(i)
!
!.....SHL HTC: superheated or in 1K interpolation interval near saturation.
!
      IF(DelTsl<1.0) THEN
!        Hshl=rv_iht_bbl_shl(i,xag,xal,vg,vl,1.d0)
         Hshl=rv_iht_bbl_shl(i,xag,vg,vl,1.d0)
         Hshl=Hshl*cell%ia_bubbly(i)
      ELSE
         Hshl=0.0d0
      ENDIF
!
!.....SCL HTC: subcooled or in 1K interpolation interval near saturation.
!
      IF(DelTsl>-1.0) THEN
!        Hscl=rv_iht_bbl_scl(i,xag,xal,vg,vl,1.d0)
         Hscl=rv_iht_bbl_scl(i,xag,vl,1.d0)
         Hscl=Hscl*cell%ia_bubbly(i)
         Hscl=umbrella_cap(Hscl,xag,cell%p(i))
      ELSE
         Hscl=0.0d0
      ENDIF
!      
!.....Power-law interpolation between SHL and SCL withing -1K~+1K
!
      htc_Hil=DMAX1(pow_interp(DelTsl,Hshl,Hscl),HTC_MIN)
!
!.....SHG HTC: superheated or in 1K interpolation interval near saturation.
!
      IF(DelTsg<1.0) THEN
!        Hshg=rv_iht_bbl_shg(i,xag,xal,vg,vl)
         Hshg=rv_iht_bbl_shg(i,xag)
         Hshg=Hshg*cell%ia_bubbly(i)
      ELSE
         Hshg=0.0d0
      ENDIF
!
!.....SCG HTC: subcooled or in 1K interpolation interval near saturation.
!
      IF(DelTsg>-1.d0) THEN
!        Hscg=rv_iht_bbl_scg(i,xag,xal,vg,vl)
         Hscg=rv_iht_bbl_scg(i,xag)
         Hscg=Hscg*cell%ia_bubbly(i)
      ELSE
         Hscg=0.d0
      ENDIF
!      
!.....Power-law interpolation between SHG and SCG withing -1K~+1K
!
      htc_Hig=MAX(pow_interp(DelTsg,Hshg,Hscg),HTC_MIN)
!
      END SUBROUTINE rv_iht_bbl 
!      
!------------------------------------------------------------------------
!     Interfacial heat transfer: slug
!------------------------------------------------------------------------
!     SUBROUTINE rv_iht_slg(i,xag,xal,vg,vl,htc_Hil,htc_Hig)
      SUBROUTINE rv_iht_slg(i,xag,vg,vl,htc_Hil,htc_Hig)
!
      USE Vol_DATA           , ONLY: cell
      USE Zrv_flowmap        , ONLY: alp_tb
      USE Zrv_ihtc_models    , ONLY: F9     
      USE Zrv_ihtc_pkg       , ONLY: rv_iht_slg_scl_SB,rv_iht_slg_scl_TB,umbrella_cap, &
                                     rv_iht_slg_shl_TB,pow_interp,HTC_MIN,rv_iht_slg_shg_SB, &
                                     rv_iht_slg_shg_TB,rv_iht_slg_scg_SB,rv_iht_slg_scg_TB
!      
      IMPLICIT NONE
! 
      INTEGER i
!      
!     REAL(8) xag,xal,vg,vl
      REAL(8) xag,vg,vl
      REAL(8) htc_Hig,htc_Hil      
      REAL(8) DelTsg,DelTsl,cF9
      REAL(8) Hshl,Hscl,Hshg,Hscg
      REAL(8) Hshl_t,Hscl_t,Hshg_t,Hscg_t
!
      DelTsl=cell%Ts(i)-cell%Tl(i)
      DelTsg=cell%Ts(i)-cell%Tg(i)
!
!.....Multiplier to reduce small bubble fraction in non-bubbly flow regimes.
!
      cF9=F9(xag,cell%alpha_bs(i),cell%alpha_sa(i))
!
!.....SHL HTC: superheated or in 1K interpolation interval near saturation.
!
      IF(DelTsl<1.d0) THEN
!         Hshl   = rv_iht_slg_shl_SB(i,xag,xal,vg,vl,cF9)
!         Hshl   = rv_iht_slg_shl_SB(i,xag,vg,vl,cF9)
!        Hshl_t = rv_iht_slg_shl_TB(i,xag,xal,vg,vl,cF9)
         Hshl_t = rv_iht_slg_shl_TB(i,xag,vg,vl,cF9)
!         Hshl   = Hshl*cell%ia_slug_sb(i)+Hshl_t*cell%ia_slug_tb(i)
!
         Hshl=Hshl_t
      ELSE
         Hshl   = 0.d0
      ENDIF
!
!.....SCL HTC: subcooled or in 1K interpolation interval near saturation. 
!
      IF(DelTsl>-1.d0) THEN
!        Hscl   = rv_iht_slg_scl_SB(i,xag,xal,vg,vl,cF9)
         Hscl   = rv_iht_slg_scl_SB(i,xag,vl,cF9)
!        Hscl_t = rv_iht_slg_scl_TB(i,xag,xal,vg,vl,cF9)
         Hscl_t = rv_iht_slg_scl_TB(i,vl)
!
!         IF(itim.eq.24042) print*,'111',itim,hscl,hscl_t
!         
         Hscl   = Hscl*cell%ia_slug_sb(i)+Hscl_t*cell%ia_slug_tb(i)
!
!         IF(itim.eq.24042) print*,'222',itim,cell%ia_slug_sb(i),cell%ia_slug_tb(i)
!         
         Hscl   = umbrella_cap(Hscl,xag,cell%p(i))
!         
!         IF(itim.eq.24042) print*,'333',itim,hscl
      ELSE
         Hscl   = 0.d0
      ENDIF
!      
!.....Power-law interpolation between SHL and SCL withing -1K~+1K
!
      htc_Hil=MAX(pow_interp(DelTsl,Hshl,Hscl),HTC_MIN)
!
!.....SHG HTC: superheated or in 1K interpolation interval near saturation.
!
      IF(DelTsg<1.d0) THEN      
!        Hshg   = rv_iht_slg_shg_SB(i,xag,xal,vg,vl,cF9)
         Hshg   = rv_iht_slg_shg_SB(i)
!        Hshg_t = rv_iht_slg_shg_TB(i,xag,xal,vg,vl,cF9)
         Hshg_t = rv_iht_slg_shg_TB(i,vg,vl)
         Hshg   = Hshg*(1.d0-alp_tb(i))*cell%ia_slug_sb(i)+Hshg_t*cell%ia_slug_tb(i)   
      ELSE
         Hshg   = 0.d0                
      ENDIF
!
!.....SCG HTC: subcooled or in 1K interpolation interval near saturation.
!
      IF(DelTsg>-1.d0) THEN    
!        Hscg   = rv_iht_slg_scg_SB(i,xag,xal,vg,vl,cF9)
         Hscg   = rv_iht_slg_scg_SB(i)
!        Hscg_t = rv_iht_slg_scg_TB(i,xag,xal,vg,vl,cF9)
         Hscg_t = rv_iht_slg_scg_TB(i)
         Hscg   = Hscg*(1.d0-alp_tb(i))*cell%ia_slug_sb(i)+Hscg_t*alp_tb(i)*cell%ia_slug_tb(i)   
      ELSE
         Hscg   = 0.d0               
      ENDIF
!      
!.....Power-law interpolation between SHG and SCG withing -1K~+1K
!
      htc_Hig=MAX(pow_interp(DelTsg,Hshg,Hscg),HTC_MIN)
!
      END SUBROUTINE rv_iht_slg      
!      
!------------------------------------------------------------------------
!     Interfacial heat transfer: annular
!------------------------------------------------------------------------
      SUBROUTINE rv_iht_anm(i,xag,xal,vg,vl,htc_Hil,htc_Hig)
!      
      USE Vol_DATA            , ONLY: cell
      USE Zconst2             , ONLY: ggc
      USE Zrv_ihtc_models     , ONLY: F10 
      USE Zrv_ihtc_pkg        , ONLY: rv_iht_anm_shl_lf,rv_iht_anm_shl_drp,rv_iht_anm_scl_lf, &
                                      umbrella_cap,pow_interp,HTC_MIN,rv_iht_anm_scl_drp,     &
                                      rv_iht_anm_shg_lf,rv_iht_anm_shg_drp,rv_iht_anm_scg_lf, &
                                      rv_iht_anm_scg_drp
!      
      IMPLICIT NONE
!
      INTEGER i
!      
      REAL(8) xag,xal,vg,vl
      REAL(8) htc_Hig,htc_Hil       
      REAL(8) DelTsg,DelTsl
      REAL(8) cF10,lambda
      REAL(8) Hshl,Hscl,Hshg,Hscg
      REAL(8) Hshl_d,Hscl_d,Hshg_d,Hscg_d
!
      DelTsl=cell%Ts(i)-cell%Tl(i)
      DelTsg=cell%Ts(i)-cell%Tg(i)
!
!.....Wavy surface effect
!
!     cF10=F10(xag,cell%rhog(i),cell%rhol(i),cell%sigma(i),vg,vl,ggc,lambda)   
      cF10=F10(xag,cell%rhog(i),cell%rhol(i),cell%sigma(i),vg,ggc,lambda)   
!
!.....SHL HTC: superheated or in 1K interpolation interval near saturation.
!
      IF(DelTsl<1.0) THEN
!        Hshl   = rv_iht_anm_shl_lf(i,xag,xal,vg,vl,1.0)
         Hshl   = rv_iht_anm_shl_lf()
!        Hshl_d = rv_iht_anm_shl_drp(i,xag,xal,vg,vl,1.d0)         
         Hshl_d = rv_iht_anm_shl_drp(i)         
         Hshl   = Hshl*cell%ia_annular_ann(i)*cF10+Hshl_d*cell%ia_annular_drp(i)
      ELSE
         Hshl   = 0.d0
      ENDIF
!
!.....SCL HTC: subcooled or in 1K interpolation interval near saturation. 
!
      IF(DelTsl>-1.d0) THEN
!        Hscl   = rv_iht_anm_scl_lf(i,xag,xal,vg,vl,1.d0)
         Hscl   = rv_iht_anm_scl_lf(i,vl)
!        Hscl_d = rv_iht_anm_scl_drp(i,xag,xal,vg,vl,1.0)         
         Hscl_d = rv_iht_anm_scl_drp(i)         
         Hscl   = Hscl*cell%ia_annular_ann(i)*cF10+Hscl_d*cell%ia_annular_drp(i)
         Hscl   = umbrella_cap(Hscl,xag,cell%p(i))
      ELSE
         Hscl   = 0.d0
      ENDIF
!      
!.....Power-law interpolation
!
      htc_Hil=MAX(pow_interp(DelTsl,Hshl,Hscl),HTC_MIN)
!
!.....SHG HTC: superheated or in 1K interpolation interval near saturation.
!
      IF(DelTsg<1.d0) THEN      
!        Hshg   = rv_iht_anm_shg_lf(i,xag,xal,vg,vl,1.d0)
         Hshg   = rv_iht_anm_shg_lf(i,xag,vg,vl)
!        Hshg_d = rv_iht_anm_shg_drp(i,xag,xal,vg,vl,1.d0)
         Hshg_d = rv_iht_anm_shg_drp(i,xag,xal,vg,vl)
         Hshg   = Hshg*cell%ia_annular_ann(i)*cF10+Hshg_d*cell%ia_annular_drp(i)
      ELSE
         Hshg   = 0.d0
      ENDIF
!
!.....SCG HTC: subcooled or in 1K interpolation interval near saturation.
!
      IF(DelTsg>-1.d0) THEN   
!        Hscg   = rv_iht_anm_scg_lf(i,xag,xal,vg,vl,1.0)
         Hscg   = rv_iht_anm_scg_lf(i)
         Hscg_d = rv_iht_anm_scg_drp(i,xag,xal)
         Hscg   = Hscg*cell%ia_annular_ann(i)*cF10+Hscg_d*cell%ia_annular_drp(i)
      ELSE
         Hscg   = 0.d0
      ENDIF
!
!.....Power-law interpolation
!
      htc_Hig=MAX(pow_interp(DelTsg,Hshg,Hscg),HTC_MIN)
!
   END SUBROUTINE rv_iht_anm  
!      
!------------------------------------------------------------------------
!     Interfacial heat transfer: inverted annular
!------------------------------------------------------------------------
!     SUBROUTINE rv_iht_invann(i,xag,xal,vg,vl,htc_Hil,htc_Hig)
      SUBROUTINE rv_iht_invann(i,xag,vg,vl,htc_Hil,htc_Hig)
!      
      USE Vol_DATA            , ONLY: cell
      USE Zrv_ihtc_pkg        , ONLY: F17,rv_iht_invann_shl_gf,rv_iht_invann_shl_sb, &
                                      rv_iht_invann_scl_gf,rv_iht_invann_scl_sb,     &
                                      umbrella_cap,pow_interp,rv_iht_invann_shg_gf,  &
                                      rv_iht_invann_shg_sb,HTC_MIN
!      
      IMPLICIT NONE
!
      INTEGER i
!      
!     REAL(8) xag,xal,vg,vl
      REAL(8) xag,vg,vl
      REAL(8) htc_Hig,htc_Hil       
      REAL(8) DelTsg,DelTsl
      REAL(8) cF16,cF17,alphaian
      REAL(8) Hshl,Hscl,Hshg,Hscg
      REAL(8) Hshl_b,Hscl_b,Hshg_b,Hscg_b
!
      DelTsl=cell%Ts(i)-cell%Tl(i)
      DelTsg=cell%Ts(i)-cell%Tg(i)
      xag=MAX(0.0d0,MIN(xag,cell%alpha_bs(i)))  ! xag: alphag -> alpha_ian
      alphaian=xag  
      cF17=F17(i,alphaian)
      cF16=MIN(1.d0-cF17,1.d0)      
!
!.....SHL HTC: superheated or in 1K interpolation interval near saturation.
!
      IF(DelTsl<1.d0) THEN
!        Hshl   = rv_iht_invann_shl_gf(i,xag,xal,vg,vl,cF16)
         Hshl   = rv_iht_invann_shl_gf()
!        Hshl_b = rv_iht_invann_shl_sb(i,xag,xal,vg,vl,cF16)         
         Hshl_b = rv_iht_invann_shl_sb(i,xag,vl,cF16)         
         Hshl   = Hshl*cell%ia_invann_ann(i)+Hshl_b*cell%ia_invann_sb(i)
      ELSE
         Hshl   = 0.d0
      ENDIF
!
!.....SCL HTC: subcooled or in 1K interpolation interval near saturation. 
!
      IF(DelTsl>-1.d0) THEN
!        Hscl   = rv_iht_invann_scl_gf(i,xag,xal,vg,vl,cF16)
         Hscl   = rv_iht_invann_scl_gf(i,xag,vg,vl)
!        Hscl_b = rv_iht_invann_scl_sb(i,xag,xal,vg,vl,cF16)         
         Hscl_b = rv_iht_invann_scl_sb(i,xag,vl)         
         Hscl   = Hscl*cell%ia_invann_ann(i)+Hscl_b*cell%ia_invann_sb(i)
         Hscl   = umbrella_cap(Hscl,xag,cell%p(i))
      ELSE
         Hscl   = 0.d0
      ENDIF
!      
!.....Power-law interpolation
!
      htc_Hil=MAX(pow_interp(DelTsl,Hshl,Hscl),HTC_MIN)
!
!.....SHG HTC: superheated or in 1K interpolation interval near saturation.
!
      IF(DelTsg<1.d0) THEN      
!        Hshg   = rv_iht_invann_shg_gf(i,xag,xal,vg,vl,1.d0)
         Hshg   = rv_iht_invann_shg_gf(i,xag)
!        Hshg_b = rv_iht_invann_shg_sb(i,xag,xal,vg,vl,1.0)
         Hshg_b = rv_iht_invann_shg_sb(i,xag)
         Hshg   = Hshg*cell%ia_invann_ann(i)+Hshg_b*cell%ia_invann_sb(i)
      ELSE
         Hshg   = 0.d0
      ENDIF
!
!.....SCG HTC: subcooled or in 1K interpolation interval near saturation.
!
      IF(DelTsg>-1.d0) THEN      
!        Hscg   = rv_iht_invann_shg_gf(i,xag,xal,vg,vl,1.d0)
         Hscg   = rv_iht_invann_shg_gf(i,xag)
!        Hscg_b = rv_iht_invann_shg_sb(i,xag,xal,vg,vl,1.0)
         Hscg_b = rv_iht_invann_shg_sb(i,xag)
         Hscg   = Hscg*cell%ia_invann_ann(i)+Hscg_b*cell%ia_invann_sb(i)
      ELSE
         Hscg   = 0.d0
      ENDIF
!
!.....Power-law interpolation
!
      htc_Hig=MAX(pow_interp(DelTsg,Hshg,Hscg),HTC_MIN)
!
   END SUBROUTINE rv_iht_invann       
!      
!------------------------------------------------------------------------
!     Interfacial heat transfer: inverted slug
!------------------------------------------------------------------------
!     SUBROUTINE rv_iht_invslg(i,xag,xal,vg,vl,htc_Hil,htc_Hig)
      SUBROUTINE rv_iht_invslg(i,xag,vg,vl,htc_Hil,htc_Hig)
!      
      USE Vol_DATA           , ONLY: cell
      USE Zrv_ihtc_pkg       , ONLY: rv_iht_invslg_shl_ann,rv_iht_invslg_shl_drp, &
                                     umbrella_cap,rv_iht_invslg_scl_ann,          &
                                     rv_iht_invslg_scl_drp,pow_interp,HTC_MIN,    &
                                     rv_iht_invslg_shg_ann,rv_iht_invslg_shg_drp
!      
      IMPLICIT NONE
!
      INTEGER i
!      
!     REAL(8) xag,xal,vg,vl
      REAL(8) xag,vg,vl
      REAL(8) htc_Hig,htc_Hil       
      REAL(8) DelTsg,DelTsl
      REAL(8) Hshl,Hscl,Hshg,Hscg
      REAL(8) Hshl_d,Hscl_d,Hshg_d,Hscg_d
!
      DelTsl=cell%Ts(i)-cell%Tl(i)
      DelTsg=cell%Ts(i)-cell%Tg(i)
!
!.....SHL HTC: superheated or in 1K interpolation interval near saturation.
!
      IF(DelTsl<1.d0) THEN
         Hshl   = rv_iht_invslg_shl_ann(i)
         Hshl_d = rv_iht_invslg_shl_drp(i)         
         Hshl   = Hshl*cell%ia_invslg_ann(i)+Hshl_d*cell%ia_invslg_drp(i)
      ELSE
         Hshl   = 0.d0
      ENDIF
!
!.....SCL HTC: subcooled or in 1K interpolation interval near saturation. 
!
      IF(DelTsl>-1.d0) THEN
         Hscl   = rv_iht_invslg_scl_ann(i)
         Hscl_d = rv_iht_invslg_scl_drp(i)
         Hscl   = Hscl*cell%ia_invslg_ann(i)+Hscl_d*cell%ia_invslg_drp(i)
         Hscl   = umbrella_cap(Hscl,xag,cell%p(i))
      ELSE
         Hscl   = 0.d0
      ENDIF
!      
!.....Power-law interpolation
!
      htc_Hil=MAX(pow_interp(DelTsl,Hshl,Hscl),HTC_MIN)
!
!.....SHG HTC: superheated or in 1K interpolation interval near saturation.
!
      IF(DelTsg<1.0) THEN      
!        Hshg   = rv_iht_invslg_shg_ann(i,xag,xal,vg,vl)
         Hshg   = rv_iht_invslg_shg_ann(i,xag)
!        Hshg_d = rv_iht_invslg_shg_drp(i,xag,xal,vg,vl)
         Hshg_d = rv_iht_invslg_shg_drp(i,xag,vg,vl)
         Hshg   = Hshg*cell%ia_invslg_ann(i)+Hshg_d*cell%ia_invslg_drp(i)
      ELSE
         Hshg   = 0.d0
      ENDIF
!
!.....SCG HTC: subcooled or in 1K interpolation interval near saturation.
!
      IF(DelTsg>-1.d0) THEN      
!        Hscg   = rv_iht_invslg_shg_ann(i,xag,xal,vg,vl)
         Hscg   = rv_iht_invslg_shg_ann(i,xag)
!        Hscg_d = rv_iht_invslg_shg_drp(i,xag,xal,vg,vl)
         Hscg_d = rv_iht_invslg_shg_drp(i,xag,vg,vl)
         Hscg   = Hscg*cell%ia_invslg_ann(i)+Hscg_d*cell%ia_invslg_drp(i)
      ELSE
         Hscg   = 0.d0
      ENDIF
!
!.....Power-law interpolation
!
      htc_Hig=MAX(pow_interp(DelTsg,Hshg,Hscg),HTC_MIN)
!
   END SUBROUTINE rv_iht_invslg  
!      
!------------------------------------------------------------------------
!     Interfacial heat transfer: mist
!------------------------------------------------------------------------
      SUBROUTINE rv_iht_mist(i,xag,xal,vg,vl,htc_Hil,htc_Hig)
!      
      USE Vol_DATA           , ONLY: cell
      USE Zrv_ihtc_pkg       , ONLY: rv_iht_mist_shl,rv_iht_mist_scl,pow_interp, &
                                     umbrella_cap,HTC_MIN,rv_iht_mist_shg,       &
                                     rv_iht_mist_scg
!      
      IMPLICIT NONE
!
!     ARgument
      INTEGER, INTENT(IN) :: i      
      REAL(8), INTENT(IN) :: xag, xal
      REAL(8), INTENT(IN) :: vg, vl
      REAL(8) htc_Hig,htc_Hil       
!
!     Local Variables 
      REAL(8) :: DelTsg, DelTsl
      REAL(8) :: Hshl,   Hscl,   Hshg,   Hscg
      
      DelTsl=cell%Ts(i)-cell%Tl(i)
      DelTsg=cell%Ts(i)-cell%Tg(i)
!
!.....SHL HTC: superheated or in 1K interpolation interval near saturation.
!
      IF(DelTsl<1.d0) THEN
!        Hshl=rv_iht_mist_shl(i,xag,xal,vg,vl)
         Hshl=rv_iht_mist_shl(i,xal)
         Hshl=Hshl*cell%ia_mist(i)
      ELSE
         Hshl=0.d0
      ENDIF
!
!.....SCL HTC: subcooled or in 1K interpolation interval near saturation. 
!
      IF(DelTsl>-1.d0) THEN
!        Hscl=rv_iht_mist_scl(i,xag,xal,vg,vl)
         Hscl=rv_iht_mist_scl(i,xal)
         Hscl=Hscl*cell%ia_mist(i)
         Hscl=umbrella_cap(Hscl,xag,cell%p(i))
      ELSE
         Hscl=0.d0
      ENDIF
!      
!.....Power-law interpolation
!
      htc_Hil=DMAX1(pow_interp(DelTsl,Hshl,Hscl),HTC_MIN)
!
!.....SHG HTC: superheated or in 1K interpolation interval near saturation.
!
      IF(DelTsg<1.d0) THEN      
!        Hshg=rv_iht_mist_shg(i,xag,xal,vg,vl)
         Hshg=rv_iht_mist_shg(i,xal,vg,vl)
         Hshg=Hshg*cell%ia_mist(i)
      ELSE
         Hshg=0.d0
      ENDIF
!
!.....SCG HTC: subcooled or in 1K interpolation interval near saturation.
!
      IF(DelTsg>-1.d0) THEN      
!        Hscg=rv_iht_mist_scg(i,xag,xal,vg,vl)
         Hscg=rv_iht_mist_scg(i,xal)
         Hscg=Hscg*cell%ia_mist(i)
      ELSE
         Hscg=0.d0
      ENDIF
!
!.....Power-law interpolation
!
      htc_Hig=DMAX1(pow_interp(DelTsg,Hshg,Hscg),HTC_MIN)
!
   END SUBROUTINE rv_iht_mist
!      
!------------------------------------------------------------------------
!     Interfacial heat transfer: vstrat
!------------------------------------------------------------------------
      SUBROUTINE rv_iht_vstrat(i,xag,xal,vg,vl,htc_Hil,htc_Hig)
!      
      USE Vol_DATA           , ONLY: cell
      USE Zcoord2            , ONLY: cell_leng      
      USE Zrv_ihtc_pkg       , ONLY: F30,F35,rv_iht_vstrat_shl,umbrella_cap, &
                                     HTC_MIN,pow_interp,rv_iht_vstrat_shg
!      
      IMPLICIT NONE
!
      INTEGER i      
!      
      REAL(8) xag,xal,vg,vl,htc_Hig,htc_Hil       
      REAL(8) DelTsg,DelTsl
      REAL(8) Hshl,Hscl,Hshg,Hscg
      REAL(8) cF30,cF35
      REAL(8) ia_vstrat
!      
      DelTsl=cell%Ts(i)-cell%Tl(i)
      DelTsg=cell%Ts(i)-cell%Tg(i)
      cF30=F30(i,xal,vg,vl,DelTsl)
!     cF35=F35(i,xal,vg,vl,DelTsg)
      cF35=F35(i,vg,vl,DelTsg)
!      ia_vstrat=cell%ia_vst_st(i)+cell%ia_vst_sb(i)
      ia_vstrat=2.d0/(cell_leng(i,1)+cell_leng(i,2))
!
!.....SHL HTC: superheated or in 1K interpolation interval near saturation.
!
      IF(DelTsl.lt.1.0) THEN
!        Hshl=rv_iht_vstrat_shl(i,xag,xal,vg,vl)
         Hshl=rv_iht_vstrat_shl(i)
!         Hshl=Hshl*cell%ia_vst_st(i)*(1.d0-cF30)+htc_Hil*cF30
         Hshl=Hshl*ia_vstrat*(1.d0-cF30)+htc_Hil*cF30
      ELSE
         Hshl=0.0d0
      ENDIF
!
!.....SCL HTC: subcooled or in 1K interpolation interval near saturation. 
!
      IF(DelTsl.gt.-1.0) THEN
!        Hscl=rv_iht_vstrat_shl( i, xag, xal, vg, vl )
         Hscl=rv_iht_vstrat_shl( i )
!         Hscl=Hscl*cell%ia_vst_st(i)*(1.d0-cF30)+htc_Hil*cF30
         Hscl=Hscl*ia_vstrat*(1.d0-cF30)+htc_Hil*cF30
         Hscl=umbrella_cap(Hscl,xag,cell%p(i))
      ELSE
         Hscl=0.0d0
      ENDIF
!      
!.....Power-law interpolation: htc_Hil
!
      htc_Hil=DMAX1(pow_interp(DelTsl,Hshl,Hscl),HTC_MIN)
!
!.....SHG HTC: superheated or in 1K interpolation interval near saturation.
!
      IF(DelTsg.lt.1.0) THEN      
!        Hshg=rv_iht_vstrat_shg( i, xag, xal, vg, vl )
         Hshg=rv_iht_vstrat_shg( i )
!         Hshg=Hshg*cell%ia_vst_st(i)*(1.d0-cF35)+htc_Hig*cF35
         Hshg=Hshg*ia_vstrat*(1.d0-cF35)+htc_Hig*cF35
      ELSE
         Hshg=0.0d0
      ENDIF
!
!.....SCG HTC: subcooled or in 1K interpolation interval near saturation.
!
      IF(DelTsg.gt.-1.0) THEN      
!        Hscg=rv_iht_vstrat_shg(i,xag,xal,vg,vl)
         Hscg=rv_iht_vstrat_shg(i)
!         Hscg=Hscg*cell%ia_vst_st(i)*(1.d0-cF35)+htc_Hig*cF35
         Hscg=Hscg*ia_vstrat*(1.d0-cF35)+htc_Hig*cF35
      ELSE
         Hscg=0.0d0
      ENDIF
!
!.....Power-law interpolation: htc_Hig
!
      htc_Hig=DMAX1(pow_interp(DelTsg,Hshg,Hscg),HTC_MIN)
!
   END SUBROUTINE rv_iht_vstrat   
