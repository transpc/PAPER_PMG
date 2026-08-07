      MODULE Zsg
! 
      IMPLICIT NONE
      SAVE
!
      INTEGER n_group,max_1d,nr_tube,iboil_model
      REAL(8) tin_1d,pin_1d,pout_1d,ein_1d,visin_1d,rhoin_1d,cpin_1d,eout_1d
      REAL(8) do_tube,di_tube,di_tube2,th_tube,cf,dp_primary,pr_flow_rated,q_rated
      REAL(8) dr_tube,dr_tube2,q_mult_pri,time_sg_heat_tune,pitch,z_econ,dz_fw,relax_hb,tune_hb
!
      INTEGER,ALLOCATABLE::n_1d(:),sd_cell(:,:),mult_cell(:,:),ih(:,:),iavb(:,:),izp(:,:),ihp(:,:)
      INTEGER,ALLOCATABLE::mult_1d_group(:),mult_1d_cell(:),mult_3d_cell1(:),mult_3d_cell2(:)
      REAL(8),ALLOCATABLE::ar_tube(:),h_tube(:,:),vol_1d(:,:),ht_area(:,:)
      REAL(8),ALLOCATABLE::f1_mult(:),f2_mult(:),drde_1d(:,:),q_sd(:,:)
      REAL(8),ALLOCATABLE::vn_1d(:,:),p_1d(:,:),rho_1d(:,:),vis_1d(:,:),tube_length(:),vin_1d(:),hyd_d(:)
      REAL(8),ALLOCATABLE::en_1d(:,:),eo_1d(:,:),t_1d(:,:),cond_1d(:,:),cp_1d(:,:),tl_sg(:,:),tg_sg(:,:),ts_sg(:,:)
      REAL(8),ALLOCATABLE::p_tube(:),vol_tube(:),t_tube(:,:,:),ng_tube(:),pr_flow(:),t_wall(:,:)
!
      INTEGER,ALLOCATABLE::igr(:),j1d(:),idc(:)
      REAL(8),ALLOCATABLE::htcl_sec(:),htcg_sec(:),htcb_sec(:),htc_pr(:),q_pri(:)
!
      END MODULE Zsg