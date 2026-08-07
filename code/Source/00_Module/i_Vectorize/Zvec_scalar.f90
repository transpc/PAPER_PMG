
      MODULE Zvec_scalar
! 
      IMPLICIT NONE
      SAVE
!
!.....pressure matrix
!  
      REAL(8), ALLOCATABLE::arli_nf(:),argi_nf(:),ardi_nf(:)
!  
!      REAL(8), ALLOCATABLE::rhol_noi(:),rhog_noi(:),rhod_noi(:)
!      REAL(8), ALLOCATABLE::rhol_nok(:),rhog_nok(:),rhod_nok(:)
!      REAL(8), ALLOCATABLE::rhol_inl(:),rhog_inl(:),rhod_inl(:)
!      REAL(8), ALLOCATABLE::rhol_out(:),rhog_out(:),rhod_out(:)
!      REAL(8), ALLOCATABLE::rhol_mcc(:),rhog_mcc(:),rhod_mcc(:)
!
!.....momentum & energy convection
!
!      REAL(8), ALLOCATABLE::ar_liq_noi(:),ar_gas_noi(:),ar_drp_noi(:)
!      REAL(8), ALLOCATABLE::ar_liq_nok(:),ar_gas_nok(:),ar_drp_nok(:)
!      REAL(8), ALLOCATABLE::ar_liq_inl(:),ar_gas_inl(:),ar_drp_inl(:)
!      REAL(8), ALLOCATABLE::ar_liq_out(:),ar_gas_out(:),ar_drp_out(:)                               
!      REAL(8), ALLOCATABLE::ar_liq_mcc(:),ar_gas_mcc(:),ar_drp_mcc(:)
!      REAL(8), ALLOCATABLE::are_liq_noi(:),are_gas_noi(:),are_drp_noi(:)
!      REAL(8), ALLOCATABLE::are_liq_nok(:),are_gas_nok(:),are_drp_nok(:)
!      REAL(8), ALLOCATABLE::are_liq_inl(:),are_gas_inl(:),are_drp_inl(:)
!      REAL(8), ALLOCATABLE::are_liq_out(:),are_gas_out(:),are_drp_out(:)                               
!      REAL(8), ALLOCATABLE::are_liq_mcc(:),are_gas_mcc(:),are_drp_mcc(:)
!
!.....Added by LSJ
!
!     (turb_ke_convection_liq,turb_ke_convection_gas)
      REAL(8), ALLOCATABLE::arli_noi(:),arli_nok(:),argi_noi(:),argi_nok(:)
!
!.....Added by CYJ
!      
      REAL(8), ALLOCATABLE::sfg_non(:,:),sfl_non(:,:),sfd_non(:,:)
      REAL(8), ALLOCATABLE::sfg_inl(:,:),sfl_inl(:,:),sfd_inl(:,:)
      REAL(8), ALLOCATABLE::sfg_out(:,:),sfl_out(:,:),sfd_out(:,:)
      REAL(8), ALLOCATABLE::sfg_mcc(:,:),sfl_mcc(:,:),sfd_mcc(:,:)
!      
      END MODULE Zvec_scalar      
