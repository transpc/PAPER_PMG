       MODULE Zvec_major
! 
!.....velocity, pressure, volume flux
!
      IMPLICIT NONE
      SAVE
!  
      REAL(8), ALLOCATABLE :: flux_l_nf(:),flux_g_nf(:),flux_d_nf(:)
      REAL(8), ALLOCATABLE :: flux_l_nf_o(:),flux_g_nf_o(:),flux_d_nf_o(:)
      REAL(8), ALLOCATABLE :: liq_conv_nf(:),vap_conv_nf(:),drp_conv_nf(:)
      REAL(8), ALLOCATABLE :: ecnvc_l_nf(:),ecnvc_g_nf(:),ecnvc_d_nf(:)
      REAL(8), ALLOCATABLE :: al_conv_nf(:),ad_conv_nf(:),void_conv_nf(:),quala_conv_nf(:)
      REAL(8), ALLOCATABLE :: lbor_conv_nf(:)
      REAL(8), ALLOCATABLE :: alphagj(:),alphalj(:)  !for rv face-upwind interfacial friction model
      REAL(8), ALLOCATABLE :: mflux_l_nf(:),mflux_g_nf(:),mflux_d_nf(:) !choke 
      REAL(8), ALLOCATABLE :: rad_ir_nbcon(:)
!  
      END MODULE Zvec_major
