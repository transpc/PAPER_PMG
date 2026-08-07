!
      SUBROUTINE single_phase_HTC(i,mode)
!
!     This routine calculates wall heat transfer coefficient and heat flux for single phase convection
!     
      USE VOL_DATA     
      USE Zwall_HTC    , ONLY: inline_bundle,qflux_t,qflux_l,qflux_g,HTC_t,HTC_tl,HTC_tst,HTC_tg,HTC_tgp, &
                                f_direc,c_direc1,c_direc2,dia_rod,pit_dia,            &
                                mflux_liqa,mflux_gasa,mflux_tota,                             &
                                rho,cond,viscos,cp,beta,tw,reflood
      USE Zvector      , ONLY: vl_n,vg_n
      USE Zio_unit     , ONLY: unit_log
!
      IMPLICIT NONE
!      
      INTEGER i,mode
!
      REAL(8) al_10,ag_10,tc,aratio,al,ag  
      REAL(8) HTC_max,HTC_lam,HTC_forced,HTC_free,HTC_freeH,HTC_cross,HTC_max2
      REAL(8) C_ditt,N_ditt,Pr,Re,Ra,gravity,gr,dum
      REAL(8) mflux,mflux_cross
!      
      DATA gravity/9.81d0/
!      
      al_10=-1.0d0       
!
!.....Set fluid property
!
      IF(mode.eq.9)THEN
         rho=cell%rhog(i)
         tc=cell%tg(i)
         cond=cell%lcondg(i)
         viscos=cell%lviscosg(i)
         cp=cell%cpg(i)
         beta=cell%betag(i)
      ELSEIF(cell%alphag(i).gt.0.999d0)then 
!        Only if Liquid fraction<0.001, consider gas phase properties.
         al_10=cell%alphal(i) * 1000.0d0    
         ag_10=1.0d0-al_10 
         rho=cell%rhol(i)*al_10+ag_10*cell%rhog(i) 
         tc=cell%tl(i)*al_10+ag_10*cell%tg(i) 
         cond=cell%lcondl(i)*al_10+ag_10*cell%lcondg(i)
         viscos=cell%lviscosl(i)*al_10+ag_10*cell%lviscosg(i) 
         cp=cell%cpl(i)*al_10+ag_10*cell%cpg(i) 
         beta=cell%betal(i)*al_10+ag_10*cell%betag(i) 
         mflux_liqa=mflux_liqa*al_10+mflux_gasa*ag_10 
      ELSE 
         rho=cell%rhol(i) 
         tc=cell%tl(i) 
         cond=cell%lcondl(i) 
         viscos=cell%lviscosl(i) 
         cp=cell%cpl(i) 
         beta=cell%betal(i) 
      ENDIF  
!
      IF(mode.eq.0.or.mode.eq.10)THEN  ! Air-water, critical fluid
         mflux=mflux_tota
      ELSEIF(mode.lt.9)THEN           ! Liquid 1-phase, boiling regions
         mflux=mflux_liqa
      ELSE                             ! mode>=9 (gas 1-phase, condensation)
         mflux= mflux_gasa + mflux_liqa*cell%rhog(i)/cell%rhol(i)
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
      C_ditt=0.023d0
      n_ditt=0.4d0
      IF(mode.eq.11)n_ditt=0.3d0               ! For Cooling case
      HTC_forced=C_ditt*(Re**0.8d0)*(Pr**n_ditt)*cond/dia_rod
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
         al=DMAX1(cell%alphal(i),1.0d-15)
         ag=DMAX1(cell%alphag(i),1.0d-15)
         mflux_cross=(vg_n(i,c_direc1)*cell%rhog(i)*ag)**2 + (vg_n(i,c_direc2)*cell%rhog(i)*ag)**2 
         mflux_cross=mflux_cross+ (vl_n(i,c_direc1)*cell%rhol(i)*al)**2 + (vl_n(i,c_direc2)*cell%rhol(i)*al)**2 
         mflux_cross=DSQRT(mflux_cross)*aratio 
         !g=sign(mflux_tota,delgrv*vl_n(i,f_direc))            !!!cyj, save sign of gravity in MARS
!                  
         HTC_cross=0.21d0*Pr**0.4d0*(mflux_cross*dia_rod/viscos)**0.62d0 *(cond/dia_rod)
         HTC_forced=DSQRT(HTC_forced*HTC_forced+HTC_cross*HTC_cross)
      ENDIF
!
!.....Churchill-Chu  (Natural convection)
! 
!      Ra=Pr*gravity*rho*rho*DABS(beta)*l_plate**3 *DMAX1(1.d-5,DABS(tw-tc)) /(viscos*viscos) 
!      Ra=Pr*gravity*rho*rho*beta*dia_rod**3 *DMAX1(1.d-5,DABS(tw-tc)) /(viscos*viscos)  ! yjm l_plate -> dia_rod
      Ra=Pr*gravity*rho*rho*DABS(beta)*dia_rod**3 *DMAX1(1.d-5,DABS(tw-tc)) /(viscos*viscos)  ! yjm l_plate -> dia_rod !pik-radiation_component
!      HTC_free=(0.825d0+0.387d0*Ra**1.66667d-1/(1.0d0+(0.492d0/Pr)**0.5625d0)**0.296296d0)**2 *(cond/l_plate)   !     1/6=1.66667, 8/27=0.296296, 9/16=0.5625
      HTC_free=(0.825d0+0.387d0*Ra**1.66667d-1/(1.0d0+(0.492d0/Pr)**0.5625d0)**0.296296d0)**2 *(cond/dia_rod) ! yjm l_plate -> dia_rod
!
!.....McAdams        (Natural convection-Horizontal pipe)
! 
!      HTC_freeH=0.27d0*Ra**0.25d0 *cond/l_plate
      HTC_freeH=0.27d0*Ra**0.25d0 *cond/dia_rod ! yjm l_plate -> dia_rod
!
!.....Use the maximum HTC
!
      HTC_max=DMAX1(HTC_lam,HTC_forced)
      IF(DMAX1(HTC_free,HTC_freeH).gt.HTC_max.and.mode.ne.9)mode=1
      HTC_max=DMAX1(HTC_max,HTC_free,HTC_freeH)
!      
!.....MARS module wall heat transfer single phase
!
      IF(reflood.eq.1.and.cell%alphag(i).gt.0.6d0)THEN
         IF(re.lt.3000.0d0)THEN
            HTC_max2=4.36d0*cond/dia_rod
         ELSEIF(re.lt.1.0d4)THEN
            HTC_max2=10.0d0*cond/dia_rod*(1.0d4-re)/7000.0d0 + 0.023d0*cond/dia_rod*pr**0.4d0*re**0.8d0*pit_dia*(re-3000.0d0)/7000.0d0
         ELSE
            HTC_max2=0.023d0*cond/dia_rod*pr**0.4d0*re**0.8d0*pit_dia
         ENDIF
!         
!        Apply Druker 2f enhancement factor (EPRI-NP3485(1984))
!        two-phase enhancement factor due to two-phase turbulence
!        assume gr/rey**2 = 1.0 and apply Drucker correlation
!
         gr=gravity*DABS(cell%rhol(i)-cell%rhog(i))*cell%rhog(i)*dia_rod**3.0d0/(viscos*viscos) ! l_plate -> dia_rod
         IF(re.lt.0.01d0)THEN
            dum=1.0d0
         ELSE
            dum=DMIN1(1.0d0,gr/(re*re))
         ENDIF
         HTC_max2=HTC_max2*DMIN1(5.0d0,(1.0d0+25.0d0*(1.0d0-cell%alphag(i))*dum)**0.5d0)
!
!        for grid spacer effect
!         IF(floss(1,i).eq.0.0d0.and.floss(2,i).eq.0.0d0.and.floss(3,i).eq.0.0d0)THEN
!            fgrid=1.0d0
!         ELSE
!            ar=0.125d0*1.14d0**0.5d0 ! 1.14 - crossflow form loss (?) 
!            ar=DMIN1(1.0d0,ar)
!            fgrid=1.0d0+5.55d0*ar*DEXP(-0.13d0*gridz(i)/dia_rod)
!         ENDIF
!         HTC_max2=HTC_max2*fgrid
!         HTC_max=HTC_max2*fgrid
!
         HTC_max=HTC_max2
      ELSE
         HTC_max=HTC_max
      ENDIF
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
         IF(cell%alphag(i).le.0.9d0.or.al_10.lt.0.d0)THEN         !IF ag<0.9, ignore HT to gas phase
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
            qflux_l=HTC_tl*(tw-cell%tl(i))
            qflux_g=HTC_tg*(tw-cell%tg(i))  
            qflux_t=qflux_l+qflux_g  
         ENDIF
      ENDIF
!       
!      tdiff=(tw-cell%tl(i))
!      IF(DABS(tdiff).gt.1.d-1)THEN
!         HTC_tl=qflux_l/tdiff
!      ELSE
!         HTC_tl=0.0d0
!         qflux_l=0.0d0
!      ENDIF      
!      HTC_tst=0.0d0
!      tdiff=(tw-cell%tg(i))
!      IF(DABS(tdiff).gt.1.d-1)THEN
!         HTC_tg=qflux_g/(tw-cell%tg(i))
!      ELSE
!         HTC_tg=0.0d0
!         qflux_g=0.0d0
!      ENDIF   
!      HTC_tgp=0.0d0 
!!         
!      IF(mode.eq.9)THEN
!         HTC_t=HTC_tg
!         qflux_t=qflux_g
!      ELSE
!         IF(cell%alphag(i).le.0.9d0.or.al_10.lt.0.d0)THEN         !IF ag<0.9, ignore HT to gas phase
!            HTC_t=HTC_tl
!            qflux_t=qflux_l
!         ELSE
!            qflux_t=qflux_l+qflux_g 
!            IF(ag_10.gt.al_10)THEN
!               HTC_t=HTC_tg/ag_10
!            ELSE
!               HTC_t=HTC_tl/al_10
!            ENDIF   
!         ENDIF
!      ENDIF  
!
      RETURN     
      END SUBROUTINE single_phase_HTC
