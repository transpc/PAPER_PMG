      MODULE Zb_condition
!      
      IMPLICIT NONE
      SAVE
!
      REAL(8),Allocatable::pbnd(:),p_fb(:),vb_liq(:,:),vb_drp(:,:),vb_gas(:,:),          &
                             vin_liq(:),vin_drp(:),vin_gas(:),                             &
                             cb_pl(:),cb_pd(:),cb_pg(:),cb_p(:),                           &
                             eb_liq(:),eb_drp(:),eb_gas(:),tb_liq(:),tb_drp(:),tb_gas(:),  &
                             qwall_liq(:),qwall_drp(:),qwall_gas(:),twall(:),              &
                             rhob_liq(:),rhob_drp(:),rhob_gas(:),                          &
                             alphab_liq(:),alphab_drp(:),alphab_gas(:),                    &
                             alpha_liq_nd(:),t_liq_nd(:),rho_liq_nd(:),e_liq_nd(:),        &
                             alpha_gas_nd(:),t_gas_nd(:),rho_gas_nd(:),e_gas_nd(:),        &
                             alpha_drp_nd(:),t_drp_nd(:),rho_drp_nd(:),e_drp_nd(:),        &
                             qualab(:),quala_nd(:),                                        &
                             turb_keb(:),turb_dpb(:),turb_kegb(:),turb_dpgb(:),            &  
                             vb_lold(:,:),vb_gold(:,:),v_wall(:),                          &
                             lvisb_liq(:),lvisb_gas(:)
!
      END MODULE Zb_condition
      
