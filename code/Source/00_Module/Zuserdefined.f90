      MODULE Zuserdefined
! 
      IMPLICIT NONE
      SAVE
!.....initial condition and mesh
      LOGICAL udfl_init_variables         !CALL initialize_specific_variables  
      LOGICAL udfl_grid_user              !CALL udfn_grid_user 
      LOGICAL udfl_porous_user            !CALL udfn_porous_user for porosity,permeability
!.....cell properties      
      LOGICAL udfl_set_qvol_porous        !CALL udfn_set_qvol_porous in set_vol_heat_source
      LOGICAL udfl_porous_property        !CALL udfn_porous_property in calc_model for aporous; necessary      
      LOGICAL udfl_mat_prop               !CALL udfn_mat_prop in mat_prop
!.....boundary condition
      LOGICAL udfl_outlet_property        !CALL udfn_outlet_property
      LOGICAL udfl_outlet_press_user      !with udfl_outlet_property
!      LOGICAL l_horizontal_outlet      !CALL udfn_horizontal_outlet
!      LOGICAL l_horizontal_outlet_init !with l_horizontal_outlet
      LOGICAL udfl_tw_profile             !CALL udfn_tw_profile
      LOGICAL udfl_hflux_bc_profile       !입력, set_wall_heat_flux와 함께처리      
      REAL(8),ALLOCATABLE :: vel_bc_profile_inl(:)    !udfl_vel_bc_profile
      REAL(8),ALLOCATABLE :: hflux_bc_profile_chw(:)  !udfl_vel_bc_profile 
      REAL(8),ALLOCATABLE :: hflux_bc_profile_chw_c(:)!udfl_vel_bc_profile      
!.....momentum transfer      
      LOGICAL udfl_mom_drag_i             !udfn_mom_drag_i
      LOGICAL udfl_mom_press_source       !CALL udfn_mom_press_source in calc_momentum
      LOGICAL udfl_mom_wall               !wall frcition 모델, wall_drag외의 모델 사용시, wall_drag에서 vv_prob 제거할 것
      LOGICAL udfl_mom_loss               !CALL udfn_mom_source in calc_momentum
      LOGICAL udfl_mom_film_shear         !CALL udfn_mom_film_shear in calc_momentum
!.....heat transfer      	      
      LOGICAL udfl_calc_HTC_int_i         !CALL udfn_calc_HTC_int_i in int_htc 
      LOGICAL udfl_flashing_hif           !with udfl_calc_HTC_int_i
      LOGICAL udfl_wallHTC_porous         !udfn_heat_wallHTC_porous
      LOGICAL udfl_erg_diff               !CALL udfn_erg_diff in scalar_energy_diffusion
!.....special option
      LOGICAL udfl_model_overwrite        !CALL udfn_model_overwrite in int_swap   
      LOGICAL udfl_update_scalar          !CALL udfn_update_scalar
      LOGICAL udfl_psbt_cfx_model         
!
!.....Restart
      REAL(8) user_rary(100)
      INTEGER user_iary(100)     
!
      LOGICAL :: MG_solver
!
      END MODULE Zuserdefined
