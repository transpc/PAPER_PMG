!
      SUBROUTINE halden650_lbloca_clean_bc_user
!
      USE Zbc_index    , ONLY: npin
      USE Zb_condition , ONLY: e_gas_nd,e_liq_nd,     &
                               rho_gas_nd,rho_liq_nd, &
                               t_gas_nd,t_liq_nd,     &
                               quala_nd,pbnd
      USE Zncg         , ONLY: tao,cvao_npin,uao_npin,dcva_npin,ra_npin,qn_npin
!
      IMPLICIT NONE                                
!
      tao          =0.d0 
!!!      nvin         =0
!!!      vin_norm(:)  =0.d0
!!!      vin_gas(:)   =0.d0
!!!      vin_liq(:)   =0.d0
!!!      vin_drp(:)   =0.d0
!!!      vb_gas(:,:)  =0.d0
!!!      vb_liq(:,:)  =0.d0
!!!      vb_drp(:,:)  =0.d0
!!!      p_fb(:)      =0.d0
!!!      tb_liq(:)    =0.d0
!!!      tb_gas(:)    =0.d0
!!!      qualab(:)    =0.d0
!!!      eb_liq(:)    =0.d0
!!!      eb_gas(:)    =0.d0
!!!      rhob_liq(:)  =0.d0
!!!      rhob_gas(:)  =0.d0
!!!      cvao_nvin(:) =0.d0
!!!      uao_nvin(:)  =0.d0
!!!      dcva_nvin(:) =0.d0
!!!      ra_nvin(:)   =0.d0
!!!      qn_nvin(:,:) =0.d0
      npin         =0
      pbnd(:)      =0.d0
      t_liq_nd(:)  =0.d0
      t_gas_nd(:)  =0.d0
      quala_nd(:)  =0.d0
      e_liq_nd(:)  =0.d0
      e_gas_nd(:)  =0.d0
      rho_liq_nd(:)=0.d0 
      rho_gas_nd(:)=0.d0 
      cvao_npin(:) =0.d0
      uao_npin(:)  =0.d0
      dcva_npin(:) =0.d0
      ra_npin(:)   =0.d0
      qn_npin(:,:) =0.d0
!      
      END SUBROUTINE halden650_lbloca_clean_bc_user
