!
      SUBROUTINE turb_ke_vis_liq
!
!     This routine calculates liquid and gas turbulent viscosity 
!     using the solutions of k-e equations
!
      USE VOL_DATA              
      USE Zparam    , ONLY: cmu,ke_small,cappa,clog
      USE Zbc_index , ONLY: icell_type
      USE Zconst1   , ONLY: turbubble,turboil,iturb,lowreynolds,turb_phase
      USE Zface     , ONLY: qecell,Kepsilon,Kepsilon_real,liq_only
      USE Ztimecon  , ONLY: alpha_min
      USE Zturb     , ONLY: turb_ke,turb_ke_o,turb_dp,tauw,velt,yplus,walln, &
                             utau,wvis_liq,wvis_gas,ustar_ke,w_real_ke,cmu_real
      USE Zturb        , ONLY: s_macroturb_source
      USE Zvector   , ONLY: vrel_o
      USE Zzone     , ONLY: ncell_fluid
!
      IMPLICIT NONE
!
      INTEGER i
!
      REAL(8) fmu
      REAL(8) vr,kataoka,hfg,alphap,ulprime
      REAL(8) tv_max,tv_mult_max
      REAL(8) pi_Real,a_s 
!
      DATA tv_mult_max/2.0d3/
!
!.....Calculation Cmu for Realizable k-e
!
      IF(iturb.eq.Kepsilon_real)THEN
         DO i=1,ncell_fluid
            pi_Real=1.d0/3.d0*DACOS(DMAX1(DMIN1(DSQRT(6.d0)*w_real_ke(i),1.0d0),-1.0d0))
            a_s=DSQRT(6.d0)*DCOS(pi_real)
            cmu_real(i)=1.0d0/(4.04d0+a_s*ustar_ke(i)*turb_ke(i)/turb_dp(i))          
         ENDDO
      ENDIF      
!
!.....Wall sheer stress and viscosity
!
      DO i=1,ncell_fluid
         IF(icell_type(i).eq.1) THEN
!
!...........Single phase wall sheer and viscosity
!
            IF(iturb.ne.Kepsilon_real)THEN      !std k-e & RNG k-e model
               tauw(i)=cell%rhol(i)*cmu**0.25d0*DSQRT(turb_ke(i))*cappa*velt(i)/DLOG(clog*yplus(i))
               wvis_liq(i)=cell%rhol(i)*cmu**0.25d0*DSQRT(turb_ke(i))*cappa/DLOG(clog*yplus(i))*walln(i)
            ELSE                      !Realizable k-e model
               tauw(i)=cell%rhol(i)*cmu_real(i)**0.25d0*DSQRT(turb_ke(i))*cappa*velt(i)/DLOG(clog*yplus(i))
               wvis_liq(i)=cell%rhol(i)*cmu_real(i)**0.25d0*DSQRT(turb_ke(i))*cappa/DLOG(clog*yplus(i))*walln(i)         
            ENDIF
!
            IF(turboil.eq.1) THEN
!
!...........Wall boiling effect (Kataoka and Serizawa,1997)
!
               Kataoka=0.0d0
               hfg=cell%hgsat(i)-cell%hlsat(i) 
               alphap=cell%alphag(i) 
!               
               IF(turb_ke_o(i).gt.ke_small) THEN
                  ulprime=DSQRT(turb_ke_o(i)) 
                  Kataoka=(cell%rhog(i)*hfg*alphap*ulprime)
               ENDIF
!               
               IF(dabs(Kataoka).gt.1.0d0) THEN
                  KAtaoka=DMIN1(10.d0,DMAX1(0.0D0,6.0d0*qecell(i)/Kataoka))
               ELSE
                  Kataoka=DMIN1(10.d0,DMAX1(0.0D0,6.0d0*qecell(i)/1.0d0))
               ENDIF
!               
               tauw(i)=(1.0D0+Kataoka)*tauw(i)
               wvis_liq(i)=wvis_liq(i)*(1.0d0+Kataoka)
!
!..............Bubble effect (Lahey,2005)
!
               vr=vrel_o(i)
               IF(cell%alphag(i).le.alpha_min) vr=0.d0
               wvis_liq(i)=wvis_liq(i)+DMIN1(10.d0*wvis_liq(i),0.6d0*cell%D1(i)*cell%alphag(i)*vr)
!                
            ENDIF
!
            utau(i)=DSQRT(DABS(tauw(i))/cell%rhol(i))
            wvis_liq(i)=DMAX1(cell%lviscosl(i),wvis_liq(i))
            IF(yplus(i).le.11.63d0) wvis_liq(i)=cell%lviscosl(i)
!
!...........Gas viscosity
!
            IF(turb_phase.eq.liq_only) THEN          
               wvis_gas(i)=wvis_liq(i)*cell%rhog(i)/cell%rhol(i)
               wvis_gas(i)=DMAX1(cell%lviscosg(i),wvis_gas(i))
            ENDIF
!
         ENDIF         
      ENDDO
!
!.....Effective viscosity and conductivity
!
      DO i=1,ncell_fluid
!
!........Turbulence viscosity
!
         IF(turb_dp(i).gt.ke_small) THEN
            IF(iturb.ne.Kepsilon_real)THEN       !std k-e & RNG k-e model
               cell%tviscosl(i)=cell%rhol(i)*cmu*turb_ke(i)*turb_ke(i)/turb_dp(i)
            ELSE                       !Realizable k-e model
               cell%tviscosl(i)=cell%rhol(i)*cmu_real(i)*turb_ke(i)*turb_ke(i)/turb_dp(i)
            ENDIF            
            IF(lowreynolds.eq.1) THEN
               fmu=1.d0-DEXP(-0.0115d0*yplus(i))    
               cell%tviscosl(i)=fmu*cell%tviscosl(i)
            ENDIF              
!
!...........Bubble effect (Lahey,2005)
!
            IF(turbubble.eq.1) THEN
               vr=vrel_o(i)
               IF(cell%alphag(i).le.alpha_min) vr=0.0d0
               cell%tviscosl(i)=cell%tviscosl(i)+DMIN1(10.0d0*cell%tviscosl(i),0.6d0*cell%D1(i)*cell%alphag(i)*vr)
            ENDIF
!            
         ELSE
            cell%tviscosl(i)=0.d0
         ENDIF
!
!........Control the maximum turbulent viscosity by tv_max
!
         IF(s_macroturb_source.ne.'chandesris'.and.s_macroturb_source.ne.'nakayama') THEN 
            tv_max=cell%lviscosl(i)*tv_mult_max
            wvis_liq(i)=DMIN1(wvis_liq(i),tv_max)
            cell%tviscosl(i)=DMIN1(cell%tviscosl(i),tv_max)
         ENDIF   
!         
!........Effective viscosity
!
         cell%eviscosl(i)=cell%lviscosl(i)+cell%tviscosl(i)
         IF(s_macroturb_source.ne.'nakayama'.and.s_macroturb_source.ne.'chandesris') THEN
            IF(icell_type(i).eq.1.and.lowreynolds.eq.0) THEN
               IF(yplus(i).le.11.63d0) THEN
                  cell%tviscosl(i)=0.d0
               ELSE
                  cell%tviscosl(i)=wvis_liq(i)-cell%lviscosl(i)
               ENDIF   
               cell%eviscosl(i)=wvis_liq(i)
            ENDIF          
         ENDIF
!
!........Gas viscosity
!
         IF(turb_phase.eq.liq_only)THEN          
            cell%tviscosg(i)=cell%tviscosl(i)*cell%rhog(i)/cell%rhol(i)
            cell%eviscosg(i)=cell%eviscosl(i)*cell%rhog(i)/cell%rhol(i)
         ENDIF
!
      ENDDO
!
      RETURN
      END SUBROUTINE turb_ke_vis_liq
