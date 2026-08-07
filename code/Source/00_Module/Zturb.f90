      MODULE Zturb
! 
      IMPLICIT NONE
      SAVE
!
      REAL(8),Allocatable::turb_dp(:),turb_dp_o(:),turb_ke(:),turb_ke_o(:),      &
                              turb_dpg(:),turb_dpg_o(:),turb_keg(:),turb_keg_o(:), &
                              pro_ke(:),tauw(:),yplus(:),utau(:),velt(:),          &
                              pro_keg(:),tauwg(:),yplusg(:),utaug(:),veltg(:),     &
                              wvis_liq(:),wvis_gas(:),wcd_liq(:),wcd_gas(:),       & 
                              walln(:),walln2(:),diff_ke(:),diff_dp(:),            &
                              wallnr(:)
      REAL(8),Allocatable::f_b2(:),dvtdn(:)                                 
      REAL(8),Allocatable::yplus_avg(:)  
      REAL(8),Allocatable::strn_ke(:),strn_keg(:)
      REAL(8),Allocatable::ustar_ke(:),ustar_keg(:),w_real_ke(:),w_real_keg(:),cmu_real(:),cmug_real(:)
!
      CHARACTER(30) s_macroturb_source  
!      
      END MODULE Zturb