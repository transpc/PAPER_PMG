!
      SUBROUTINE init_SST(ke,dp,v)
!
!     This routine initializes turbulence kinetic energy(k) and 
!     dissipation(epsilon). 
!     
      USE VOL_DATA
      USE Zconst1         , ONLY: turb_phase,iVisRatio,vis_ratio
      USE Zparam          , ONLY: ndim,cmu
      USE Zmpi            , ONLY: ncell_fp
      USE Zzone           , ONLY: ncell_fluid
!
      IMPLICIT NONE
!
      INTEGER i,ix
!      
      REAL(8) ke(ncell_fp),dp(ncell_fp),v(ncell_fp,ndim)
      REAL(8) turb_length_scale,turb_intensity,eddy_scale,vmean
!
      turb_intensity=0.05d0
      turb_intensity=turb_intensity*turb_intensity
      turb_length_scale=0.5d0    
      eddy_scale=0.038d0*turb_length_scale       !!! CYJ k-w
!
      DO i=1,ncell_fluid
         vmean=0.0d0
!         
         DO ix=1,ndim
            vmean=vmean+v(i,ix)*v(i,ix)
         ENDDO
!         
         IF(vmean.gt.0.0d0)THEN
            vmean=vmean*turb_intensity
         ELSE
            vmean=1.0d-4*turb_intensity  
         ENDIF
!         
         ke(i)=1.5d0*vmean
         IF(iVisRatio.eq.0)THEN
            dp(i)=cmu**(-0.25d0)*DSQRT(ke(i))/eddy_scale      !!! CYJ k-w
!            utau(i)=cmu**0.25d0*DSQRT(ke(i))
         ELSEIF(iVisRatio.eq.1)THEN
            IF(turb_phase.eq.1)THEN         !gas only
               dp(i)=cell%rhog(i)*ke(i)/(cell%lviscosg(i)*vis_ratio)
            ELSEIF(turb_phase.eq.2)THEN     !liq only
               dp(i)=cell%rhol(i)*ke(i)/(cell%lviscosl(i)*vis_ratio)
            ENDIF
         ENDIF
!         
      ENDDO
!
      RETURN
      END SUBROUTINE init_SST
