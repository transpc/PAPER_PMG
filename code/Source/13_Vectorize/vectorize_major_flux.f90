      SUBROUTINE vectorize_major_flux
!
!     This routine vectorizes volume fluxes.
!
      USE Zvec_param  ,ONLY: nf_flux,nf_tot_nbcon
      USE Zvec_major  ,ONLY: flux_l_nf,flux_g_nf,flux_d_nf,         &
                             flux_l_nf_o,flux_g_nf_o,flux_d_nf_o,   &
                             liq_conv_nf,vap_conv_nf,drp_conv_nf,   &
                             ecnvc_l_nf,ecnvc_g_nf,ecnvc_d_nf,      &
                             al_conv_nf,ad_conv_nf,void_conv_nf,quala_conv_nf, &
                             mflux_l_nf,mflux_g_nf,mflux_d_nf,rad_ir_nbcon
!
      IMPLICIT NONE
!
      LOGICAL,SAVE :: initial
!
      DATA initial/.true./ 
!
      IF(initial)THEN
         ALLOCATE(flux_l_nf(nf_flux),flux_g_nf(nf_flux),flux_d_nf(nf_flux))
         ALLOCATE(flux_l_nf_o(nf_flux),flux_g_nf_o(nf_flux),flux_d_nf_o(nf_flux))
         ALLOCATE(liq_conv_nf(nf_flux),vap_conv_nf(nf_flux),drp_conv_nf(nf_flux))
         ALLOCATE(ecnvc_l_nf(nf_flux),ecnvc_g_nf(nf_flux),ecnvc_d_nf(nf_flux))
         ALLOCATE(al_conv_nf(nf_flux),ad_conv_nf(nf_flux), &
                  void_conv_nf(nf_flux),quala_conv_nf(nf_flux))
         ALLOCATE(rad_ir_nbcon(nf_tot_nbcon))
!choking         
         ALLOCATE(mflux_l_nf(nf_flux),mflux_g_nf(nf_flux),mflux_d_nf(nf_flux))
         mflux_l_nf(:)=0.0d0
         mflux_g_nf(:)=0.0d0
         mflux_d_nf(:)=0.0d0
!                  
         flux_l_nf(:)=0.0d0
         flux_g_nf(:)=0.0d0
         flux_d_nf(:)=0.0d0
         flux_l_nf_o(:)=0.0d0
         flux_g_nf_o(:)=0.0d0
         flux_d_nf_o(:)=0.0d0
!
         liq_conv_nf(:)=0.0d0
         vap_conv_nf(:)=0.0d0
         drp_conv_nf(:)=0.0d0
!
         ecnvc_l_nf(:)=0.0d0 
         ecnvc_g_nf(:)=0.0d0 
         ecnvc_d_nf(:)=0.0d0 
!
         al_conv_nf(:)=0.0d0
         ad_conv_nf(:)=0.0d0
         void_conv_nf(:)=0.0d0
         quala_conv_nf(:)=0.0d0
!
         rad_ir_nbcon(:)=0.0d0
      ENDIF      
!      
      initial=.FALSE.     
!
      RETURN
      END SUBROUTINE vectorize_major_flux
