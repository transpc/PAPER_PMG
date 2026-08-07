      MODULE Ztimecon
! 
      IMPLICIT NONE
      SAVE
!
      INTEGER smac,dt_opt,nctrl,ctrl_opt,ictrl,itim_restart,itim_restart_init,itim,itim_last
      INTEGER smac_ctrl(100),dt_opt_ctrl(100),iso_thermal,iter_p
      REAL(8) time,dt_max
      REAL(8) dxmin
      REAL(8) alpha_min
      REAL(8) error_mass,dp_max
      REAL(8) t_end,toutstep,cfl_ratio,cfl_ratio_max,treststep
      REAL(8) t_end_ctrl(0:100),toutstep_ctrl(0:100),cfl_ratio_max_ctrl(100),dt_max_ctrl(100),treststep_ctrl(0:100)
      INTEGER nbline
      LOGICAL repeat_smac
!
      END MODULE Ztimecon