!
      SUBROUTINE turb_SST_vis_liq
!
!     This routine calculates liquid and gas turbulent viscosity 
!     using the solutions of k-e equations
!
      USE VOL_DATA     , ONLY: cell
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: cmu,ke_small,cappa,clog
      USE Zbc_index    , ONLY: icell_type
      USE Zconst1      , ONLY: turbubble,turboil,turb_phase
      USE Zface        , ONLY: qecell,liq_only
      USE Ztimecon     , ONLY: alpha_min
      USE Zturb        , ONLY: turb_ke,turb_ke_o,turb_dp,tauw, &
                                wvis_liq,f_b2,strn_ke
      USE Zvector      , ONLY: vrel_o
!
      IMPLICIT NONE
!
      INTEGER i
!
      REAL(8) vr,kataoka,hfg,alphap,ulprime
      REAL(8) tv_max,tv_mult_max
!     REAL(8) por,velo,denom,numer,re,hd,dsqrt_por !LSJ 161122 pporous 
      REAL(8) alpha_star      !CYJ k-w  
      REAL(8),SAVE::alpha_l
!
      DATA tv_mult_max/2.d3/      
!
!.....Additional constants
!
      alpha_l=0.31d0
!
!.....Effective viscosity and conductivity
!
      DO i=1,ncell_fluid
!
!........Turbulence viscosity
!
         IF(turb_dp(i).gt.ke_small) THEN
!
!...........Calculation of invariant measure of the strain rate for k-w model
!
!            re_t=cell%rhol(i)*turb_ke(i)/cell%tviscosl(i)/turb_dp(i)
!            r_k=6.d0
!            alpha_star=(0.072d0/3.d0+re_t/r_k)/(1.d0+re_t/r_k)     ! Fluent 12.0
            alpha_star=1.d0                                         ! Fluent 14.0
            cell%tviscosl(i)=cell%rhol(i)*alpha_l*turb_ke(i)/DMAX1(alpha_l*turb_dp(i)/alpha_star,f_b2(i)*strn_ke(i))              !CYJ k-w
!
!...........Bubble effect (Lahey,2005)
!
            IF(turbubble.eq.1) THEN
               vr=vrel_o(i)
               IF(cell%alphag(i).le.alpha_min) vr=0.d0
               cell%tviscosl(i)=cell%tviscosl(i)+DMIN1(10.d0*cell%tviscosl(i),0.6d0*cell%D1(i)*cell%alphag(i)*vr)
            ENDIF
!
!...........Wall boiling effect (Kataoka and Serizawa,1997)
!           
            IF(icell_type(i).eq.1.and.turboil.eq.1) THEN  
               Kataoka=0.d0
               hfg=cell%hgsat(i)-cell%hlsat(i) 
               alphap=cell%alphag(i) 
!               
               IF(turb_ke_o(i).gt.ke_small) THEN
                  ulprime=DSQRT(turb_ke_o(i)) 
                  Kataoka=(cell%rhog(i)*hfg*alphap*ulprime)
               ENDIF
!               
               IF(dabs(Kataoka).gt.1.d0) THEN
                  Kataoka=DMIN1(10.d0,DMAX1(0.0D0,6.d0*qecell(i)/Kataoka))
               ELSE
                  Kataoka=DMIN1(10.d0,DMAX1(0.0D0,6.d0*qecell(i)/1.d0))
               ENDIF           
               tauw(i)=(1.0D0+Kataoka)*tauw(i) 
               cell%tviscosl(i)=cell%tviscosl(i)*(1.d0+Kataoka)
            ENDIF  
         ELSE
            cell%tviscosl(i)=0.d0
         ENDIF        
!
!........Control the maximum turbulent viscosity by tv_max
!
         tv_max=cell%lviscosl(i)*tv_mult_max
         cell%tviscosl(i)=DMIN1(cell%tviscosl(i),tv_max)
!         
!........Effective viscosity
!
         cell%eviscosl(i)=cell%lviscosl(i)+cell%tviscosl(i)
!         IF(icell_type(i).eq.1.and.lowreynolds.eq.0) THEN
!            IF(yplus(i).le.11.63d0) THEN
!               cell%tviscosl(i)=0.d0
!            ELSE
!               cell%tviscosl(i)=wvis_liq(i)-cell%lviscosl(i)
!              
!            ENDIF   
!            cell%eviscosl(i)=wvis_liq(i)
!         ENDIF                       
!!
!!........Wall sheer stress and viscosity, 
!!
!         IF(icell_type(i).eq.1) THEN
!            tauw(i)=cell%lviscosl(i)*DABS(dvtdn(i))
!            utau(i)=DSQRT(DABS(tauw(i))/cell%rhol(i))
!         ENDIF
         IF(icell_type(i).eq.1)wvis_liq(i)=cell%eviscosl(i)
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
      END SUBROUTINE turb_SST_vis_liq
