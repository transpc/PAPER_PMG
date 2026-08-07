!
      SUBROUTINE single_phase_HTC_rod(i,mode,alphag,alphal,rhog,rhol,tg,tl,vg,vl,condg,condl,visg,visl,cpg,cpl,betag,betal)

!
!     This routine calculates wall heat transfer coefficient and heat
!     flux for single phase convection
!     
      USE VOL_DATA
      USE Zwall_HTC    , ONLY: inline_bundle,qflux_t,qflux_l,qflux_g,HTC_t,HTC_tl,HTC_tst,HTC_tg,HTC_tgp,&
                                f_direc,c_direc1,c_direc2,l_plate,dia_rod,pit_dia,&
                                mflux_liqa,mflux_gasa,mflux_tota,&
                                rho,cond,viscos,cp,beta,tw
      USE Zrv_ncell    , ONLY: channel_cell_hts2d,cupid_cell_channel
      USE Zio_unit     , ONLY: unit_log
      USE Zconst2      , ONLY: hydraulicd      
!
      IMPLICIT NONE
!      
      INTEGER i,mode
!
      REAL(8) al_10,ag_10,tc,aratio,al,ag
      REAL(8) HTC_max,HTC_lam,HTC_forced,HTC_free,HTC_freeH,HTC_cross
      REAL(8) C_ditt,N_ditt,Pr,Re,Ra,gravity
      REAL(8) mflux,mflux_cross
!      
      ! Additional connectivity - jrlee
      INTEGER m1,k
      REAL(8) alphag,rhog,condg,visg,cpg,betag,tg,vg
      REAL(8) alphal,rhol,condl,visl,cpl,betal,tl,vl
!      
      DATA gravity/9.81d0/
!      
      al_10=-1.0d0
!
!.....Set fluid property
!
      k=channel_cell_hts2d(i)
      m1=cupid_cell_channel(k)
      IF(mode.eq.9)THEN
         rho=rhog
         tc=tg
         cond=condg
         viscos=visg
         cp=cpg
         beta=betag
      ELSEIF(alphag.gt.0.999d0)then
         al_10=alphal * 1000.0d0
         ag_10=1.0d0-al_10
         rho=rhol*al_10+ag_10*rhog
         tc=tl*al_10+ag_10*tg
         cond=condl*al_10+ag_10*condg
         viscos=visl*al_10+ag_10*visg
         cp=cpl*al_10+ag_10*cpg
         beta=betal*al_10+ag_10*betag
         mflux_liqa=mflux_liqa*al_10+mflux_gasa*ag_10
      ELSE
         rho=rhol
         tc=tl
         cond=condl
         viscos=visl
         cp=cpl
         beta=betal
      ENDIF
!
      IF(mode.eq.0.or.mode.eq.10)THEN  ! Air-water, critical fluid
         mflux=mflux_tota
      ELSEIF(mode.lt.9)THEN           ! Liquid 1-phase, boiling regions
         mflux=mflux_liqa
      ELSE                             ! mode>=9 (gas 1-phase,condensation)
         !stop 'check mode in single_phase_HTC0 wow'
         mflux= mflux_gasa + mflux_liqa*rhog/rhol
      ENDIF
!      
!.....Laminar HTC
!
      HTC_lam=4.36d0*cond/dia_rod
!
!.....Dittus-Boelter (Forced convection)
!
      Pr=Cp*viscos/cond
      Re=mflux*dia_rod/viscos
      Re=mflux*hydraulicd(m1)/viscos
      C_ditt=0.023d0
      n_ditt=0.4d0
      IF(mode.eq.11)n_ditt=0.3d0               ! For Cooling case
      HTC_forced=C_ditt*(Re**0.8d0)*(Pr**n_ditt)*cond/dia_rod
      HTC_forced=C_ditt*(Re**0.8d0)*(Pr**n_ditt)*cond/hydraulicd(m1)      
!
!.....Bundle with in-line rod option
!      
      IF(inline_bundle.eq.1)THEN
!        Multiply G by average area to rod gap area for HTC.   
         aratio=(pit_dia*pit_dia-0.7854d0)/(pit_dia*pit_dia-pit_dia)
!      
         HTC_forced=HTC_forced*pit_dia
!
!........Calculate cross flow mass flux (Default: x and y direction)
!
         al=alphal !DMAX1(cell%alphal(m1),1.0d-15)
         ag=alphag !DMAX1(cell%alphag(m1),1.0d-15)
         mflux_cross=(vg*rhog*ag)**2 + (vg*rhog*ag)**2
         mflux_cross=mflux_cross+ (vl*rhol*al)**2 + (vl*rhol*al)**2
         mflux_cross=DSQRT(mflux_cross)*aratio
         !g=sign(mflux_tota,delgrv*vl_n(i,f_direc))            !!!cyj,
         !save sign of gravity in MARS
!                  
         HTC_cross=0.21d0*Pr**0.4d0*(mflux_cross*dia_rod/viscos)**0.62d0*(cond/dia_rod)
         HTC_forced=DSQRT(HTC_forced*HTC_forced+HTC_cross*HTC_cross)
      ENDIF
!
!.....Churchill-Chu  (Natural convection)
! 
      Ra=Pr*gravity*rho*rho*DABS(beta)*l_plate**3*DMAX1(1.d-5,DABS(tw-tc)) /(viscos*viscos) !apr1400_lbloca_debug
      HTC_free=(0.825d0+0.387d0*Ra**1.66667d-1/(1.0d0+(0.492d0/Pr)**0.5625d0)**0.296296d0)**2*(cond/l_plate)
!
!.....McAdams        (Natural convection-Horizontal pipe)
! 
      HTC_freeH=0.27d0*Ra**0.25d0 *cond/l_plate
!
!.....Use the maximum HTC
!
      HTC_max=DMAX1(HTC_lam,HTC_forced)
      IF(DMAX1(HTC_free,HTC_freeH).gt.HTC_max.and.mode.ne.9)mode=1
      HTC_max=DMAX1(HTC_max,HTC_free,HTC_freeH)
!
!.....Avoid negative HTC
!
     IF(HTC_max.le.0.d0)THEN
         WRITE(*,*)  'Single phase HTC has negative sign. Set to 1.0'
         WRITE(unit_log,*) 'Single phase HTC has negative sign. Set to 1.0'
      ENDIF
      HTC_t=DMAX1(HTC_max,1.0d0)
!
!.....Save HTC and heat flux
!
      qflux_t=HTC_t*(tw-tc)
      IF(mode.eq.9)THEN
         HTC_tl=0.0d0
         HTC_tst=0.0d0
         HTC_tg=HTC_t
         HTC_tgp=0.0d0
         qflux_l=0.0d0
         qflux_g=qflux_t
      ELSE
         IF(alphag.le.0.9d0.or.al_10.lt.0.d0)THEN         !IF ag<0.9,ignore HT to gas phase
            HTC_tl=HTC_t
            HTC_tst=0.0d0
            HTC_tg=0.0d0
            HTC_tgp=0.0d0
            qflux_l=qflux_t
            qflux_g=0.0d0
         ELSE
            HTC_tl=HTC_t*al_10
            HTC_tst=0.0d0
            HTC_tg=HTC_t*ag_10
            HTC_tgp=0.0d0
            qflux_l=HTC_tl*(tw-tl)
            qflux_g=HTC_tg*(tw-tg)
            qflux_t=qflux_l+qflux_g
         ENDIF
      ENDIF
!      
      RETURN
      END SUBROUTINE single_phase_HTC_rod

