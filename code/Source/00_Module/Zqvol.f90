      MODULE Zqvol
! 
      IMPLICIT NONE
      SAVE
!
      REAL(8) q0_solidi
      REAL(8) qwall_origin
      REAL(8) qporous_tot      
      REAL(8),DIMENSION(:),ALLOCATABLE :: q0_gas,q0_drp,q0_liq,            &
                                          qvol_liq,qvol_drp,qvol_gas,      &
                                          H_il,H_ig,gamma,H_gf,gamma_wall, &
                                          hgf_o
      REAL(8),DIMENSION(:),ALLOCATABLE :: q0_ice_solid
      REAL(8),DIMENSION(:),ALLOCATABLE :: qporous_liq,qporous_gas,qvol_ice_solid,qporous_gamma
      REAL(8),DIMENSION(:),ALLOCATABLE :: qrv_gas(:),qrv_liq(:),qrv_gamma   
      REAL(8),DIMENSION(:),ALLOCATABLE :: t_bulk(:),t_plus(:),t_plus_bulk
      REAL(8),DIMENSION(:),ALLOCATABLE :: hil_o(:),hig_o,nsiteden_o,dry_weight
! Qwall_solid
      REAL(8),DIMENSION(:),ALLOCATABLE :: qwall_solid
!
!.....Convective heat transfer term
      REAL(8),DIMENSION(:),ALLOCATABLE :: qconv_sol
      REAL(8),DIMENSION(:),ALLOCATABLE :: htc_convw,tb_convw,ha_convw
!      
      END MODULE Zqvol
