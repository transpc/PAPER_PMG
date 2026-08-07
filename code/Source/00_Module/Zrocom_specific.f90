      MODULE Zrocom_specific
! 
      IMPLICIT NONE
      SAVE
!
      INTEGER CLOSEloop,injection,size_idx(3)
      INTEGER, ALLOCATABLE:: idx_dc_top(:),idx_dc_bot(:),idx_core(:)
      REAL(8),ALLOCATABLE::tl_space_min(:),tl_space_avg(:),out_mdot(:),in_mdot(:)
      REAL(8),ALLOCATABLE::mdot_bc(:),vtop_liq(:,:),vtop_drp(:,:),          &
                             vtop_gas(:,:),etop_liq(:),etop_drp(:),           &
                             etop_gas(:),ttop_liq(:),ttop_drp(:),ttop_gas(:), &
                             rhotop_liq(:),rhotop_drp(:),rhotop_gas(:),       &
                             alphatop_liq(:),alphatop_drp(:),alphatop_gas(:), &
                             qualatop(:),rad_group(:)
!
      END MODULE Zrocom_specific