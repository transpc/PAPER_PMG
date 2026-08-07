!      
      SUBROUTINE udfn_sg_allocate
!
      USE Zmpi     ,ONLY:ncell_fp
      USE Zzone    ,ONLY:ncell_fluid
      USE Zsg      ,ONLY:n_group,max_1d,vn_1d,p_1d,rho_1d,vis_1d,vin_1d,en_1d,eo_1d,q_pri, &
                         t_1d,q_sd,drde_1d,cp_1d,cond_1d,pr_flow,htc_pr,                   &
                         htcl_sec,htcg_sec,t_wall,hyd_d,tl_sg,tg_sg,ih,iavb,idc,izp,ihp,   &
                         ht_area,htcb_sec,ts_sg,igr,j1d
!
      IMPLICIT NONE 
!
      ALLOCATE(vn_1d(n_group,max_1d),p_1d(n_group,max_1d),rho_1d(n_group,max_1d))
      ALLOCATE(vis_1d(n_group,max_1d),vin_1d(n_group),ih(n_group,max_1d),iavb(n_group,max_1d))
      ALLOCATE(en_1d(n_group,max_1d),eo_1d(n_group,max_1d),t_1d(n_group,max_1d))
      ALLOCATE(drde_1d(n_group,max_1d),cp_1d(n_group,max_1d),cond_1d(n_group,max_1d))
      ALLOCATE(q_sd(n_group,max_1d),pr_flow(n_group),hyd_d(n_group))
      ALLOCATE(ht_area(n_group,max_1d))
      ALLOCATE(tl_sg(n_group,max_1d),tg_sg(n_group,max_1d),ts_sg(n_group,max_1d),t_wall(n_group,max_1d))
      ALLOCATE(izp(n_group,max_1d),ihp(n_group,max_1d))
!
      ALLOCATE(igr(ncell_fluid),j1d(ncell_fluid),idc(ncell_fp))
      ALLOCATE(htcl_sec(ncell_fp),htcg_sec(ncell_fp),htcb_sec(ncell_fp),htc_pr(ncell_fp),q_pri(ncell_fp))
!
      vn_1d(:,:)=0.0d0
      p_1d(:,:)=0.0d0
      rho_1d(:,:)=0.0d0
      vis_1d(:,:)=0.0d0
      vin_1d(:)=0.0d0
      en_1d(:,:)=0.0d0
      eo_1d(:,:)=0.0d0
      t_1d(:,:)=0.0d0
      drde_1d(:,:)=0.0d0
      cp_1d(:,:)=0.0d0
      cond_1d(:,:)=0.0d0
      q_sd(:,:)=0.0d0
      pr_flow(:)=0.0d0
      t_wall(:,:)=0.0d0
      tl_sg(:,:)=0.0d0
      tg_sg(:,:)=0.0d0
      hyd_d(:)=0.0d0
      ih(:,:)=0.0d0
      iavb(:,:)=0.0d0
      izp(:,:)=0.0d0
      ihp(:,:)=0.0d0
      ht_area(:,:)=0.0d0
      ts_sg(:,:)=0.0d0
!
      idc(:)=0
      igr(:)=0
      j1d(:)=0
      htcl_sec(:)=0.0d0
      htcg_sec(:)=0.0d0
      htcb_sec(:)=0.0d0
      htc_pr(:)=0.0d0
      q_pri(:)=0.0d0
!
      RETURN
      END SUBROUTINE udfn_sg_allocate
