      REAL(8) tot_pwr,sum_rodleng,time_3dk,new_power3D
      REAL(8) t_rx_trip,t_edit_MASTER,dt_edit_MASTER
      REAL(8) ttime_0
      INTEGER(4) n_rm_fuel,n_rm_refl,n_am,nm_botR,nm_topR,c_no_rfl
      INTEGER(4) no_ctrl_rod
      INTEGER(4) njv_hs3d(2000,2)

      COMMON /MASTER_C/tot_pwr,sum_rodleng,time_3dk,new_power3D,    &
             t_rx_trip,t_edit_MASTER,dt_edit_MASTER,                &
             ttime_0,                                               &
             n_rm_fuel,n_rm_refl,n_am,nm_botR,nm_topR,c_no_rfl,     &
             no_ctrl_rod,njv_hs3d

      !tot_pwr: rated core thermal power
      !sum_rodleng: total rod length in meter
      !time_3dk: The time when transient advancement begins
      !t_rx_trip: The time when reactor trip occurs
      !n_rm_fuel: Number of radial meshes for the fuel assembly region
      !n_rm_refl: Number of radial meshes for the reflector regionn
      !n_am: Number of axial meshes
      !nm_botR: Mesh number (index) of the bottom reflector
      !nm_topR: Mesh number (index) of the top reflector
      !c_no_rfl: Channel number of the radial reflector (Not used)
      !no_ctrl_rod: Number of control rod banks
      !njv_hs3d(n,1): n-th inlet junction number
      !njv_hs3d(n,1): n-th exit volume number
