!
      SUBROUTINE allocate_model_var
!
!     This routine allocates variables containing the information on models.
!
      USE Zmpi         , ONLY: ncell_fp
      USE Zparam       , ONLY: ndim,nb_max
      USE Zbc_index    , ONLY: num_wall_group
      USE Zboron       , ONLY: cboronb_liq,cboronb
      USE Zconst1      , ONLY: iat,iturb,nlift,nwlf,ntdf,iheatpart,fric_face,rv_htmodel_forCFD
      USE Zface        , ONLY: q1cell,qqcell,qecell,qclcell,qcgcell,ndensitycell,twall_model
      USE Ziat         , ONLY: ia_conv,ia_old,ia,ia_b,dsm_b,dbubble_init, &
                                iat_size,iat_coal,iat_break,iat_nucl
      USE Zmodel       , ONLY: drift_c0,drift_c1,cb_bubble,h_il_cfd,h_ig_cfd,h_gf_cfd,vfgl_cfd,vfgd_cfd,resist
      USE Znormal      , ONLY: sa_walll      
      USE Zndforce     , ONLY: c_bface,c_bface_indx,dvdxl,d_bfc,f_wl,cwlf,clift,ctd,            &
                               vfgl_o,relax_cd,                                                 &
                               face_wall_group,cell_closewall,cell_closewall_indx,dis_closewall
      USE Zqvol        , ONLY: t_bulk,t_plus,t_plus_bulk,nsiteden_o,dry_weight,                 &
                               htc_convw,tb_convw,ha_convw
      USE Zrv_model    , ONLY: rv_model
      USE Zturb        , ONLY: turb_dp,turb_dp_o,turb_ke,turb_ke_o,turb_dpg,turb_dpg_o,         &
                               turb_keg,turb_keg_o,diff_ke,diff_dp,                             &
                               pro_ke,tauw,yplus,utau,velt,pro_keg,pro_keg,tauwg,yplusg,        &
                               utaug,veltg,strn_ke,strn_keg,                                    &
                               wvis_liq,wvis_gas,wcd_liq,wcd_gas,walln,walln2,wallnr,f_b2,dvtdn
      USE Zturbzeq     , ONLY: cell_hindex,tleng,tlengs,vorticity,chheight
      USE Zvector      , ONLY: face_fr_l,face_fr_g
      USE Zzone        , ONLY: ncell_fluid,ncell_cond_all,ncell_cond
      USE Zrad_comp    , ONLY: qrad_rod,qrad_sol,qrad_flu     
!
      IMPLICIT NONE
!
      INTEGER n, n1
!
      n=ncell_fp
      n1=ncell_fluid
!
!.....Ziat
!
      IF(iat.gt.0)THEN      
         ALLOCATE(ia_conv(n),ia_b(nb_max),dsm_b(nb_max))
         ALLOCATE(iat_size(n),iat_coal(n),iat_break(n),iat_nucl(n))
         ia_conv(:)=0.0d0
         ia_b(:)=0.0d0
         dsm_b(:)=0.0d0
         iat_size(:)=0.0d0
         iat_coal(:)=0.0d0
         iat_break(:)=0.0d0
         iat_nucl(:)=0.0d0
         dbubble_init=0.0d0
      ENDIF
      ALLOCATE(ia_old(n),ia(n))
      ia_old(:)=0.0d0
      ia(:)=0.0d0
!
!.....Zturb
!
      IF(iturb.gt.0)THEN
         ALLOCATE(turb_dpg(n),turb_dpg_o(n),turb_keg(n),turb_keg_o(n),    &
                  diff_ke(n),diff_dp(n), pro_ke(n),tauw(n),velt(n),       &
                  pro_keg(n),tauwg(n),yplusg(n),utaug(n),veltg(n),        &
                  f_b2(n),dvtdn(n),strn_ke(n),strn_keg(n))
         turb_dpg(:)=0.0d0
         turb_dpg_o(:)=0.0d0
         turb_keg(:)=0.0d0
         turb_keg_o(:)=0.0d0
         diff_ke(:)=0.0d0
         diff_dp(:)=0.0d0
         pro_ke(:)=0.0d0
         tauw(:)=0.0d0
         velt(:)=0.0d0
         pro_keg(:)=0.0d0
         tauwg(:)=0.0d0
         yplusg(:)=0.0d0
         utaug(:)=0.0d0
         veltg(:)=0.0d0
         f_b2(:)=0.0d0
         dvtdn(:)=0.0d0
         strn_ke(:)=0.0d0
         strn_keg(:)=0.0d0     
      ENDIF
      IF(iturb.ge.0)THEN
         ALLOCATE(turb_dp(n),turb_dp_o(n),turb_ke(n),turb_ke_o(n),    &
                  wvis_liq(n),wvis_gas(n),wcd_liq(n),wcd_gas(n))
         turb_dp(:)=0.0d0
         turb_dp_o(:)=0.0d0
         turb_ke(:)=0.0d0
         turb_ke_o(:)=0.0d0         
         wvis_liq(:)=0.0d0
         wvis_gas(:)=0.0d0
         wcd_liq(:)=0.0d0
         wcd_gas(:)=0.0d0  
      ENDIF      
      IF(iturb.ge.0.or.iheatpart.gt.0)THEN
         ALLOCATE(yplus(n),utau(n))
         yplus(:)=0.0d0
         utau(:)=0.0d0  
      ENDIF      
      ALLOCATE(walln(n),walln2(n),wallnr(n))
      walln(:)=0.0d0
      walln2(:)=0.0d0
!
!.....Zturbzeq
!
      IF(iturb.eq.0)THEN
         ALLOCATE(cell_Hindex(n))
         ALLOCATE(tleng(n),vorticity(n),ChHeight(5,0:1000))
         ALLOCATE(dvdxl(n,ndim,ndim))
         cell_Hindex(:)=0
         tlengs=0.0d0
         tleng(:)=0.0d0
         vorticity(:)=0.0d0
         ChHeight(:,:)=0.0d0
         dvdxl(:,:,:)=0.0d0
      ENDIF
!
!.....Relaxation
!
      IF(relax_cd.ne.0.0d0)THEN
         ALLOCATE(vfgl_o(n))
         vfgl_o(:)=0.0d0
      ENDIF      
!
!.....Zndforce
!
      IF(nlift.ne.-1.0d0)THEN    !Nlift can have a negative value.
         ALLOCATE(Clift(n))
         Clift(:)=0.0d0
      ENDIF
      IF(ntdf.ne.-1.0d0)THEN
         ALLOCATE(Ctd(n))
         Ctd(:)=0.0d0
      ENDIF
      IF(nwlf.ne.-1.0d0)THEN
         ALLOCATE(Cwlf(n),F_wl(n,ndim))
         F_wl(:,:)=0.0d0
         Cwlf(:)=0.0d0
      ENDIF      
! 
!.....Find walls
!
      IF(iturb.ge.0.or.nwlf.ne.-1.0d0)THEN
         ALLOCATE(c_bface(n),c_bface_indx(n))
         !ALLOCATE(nwb(n))
         ALLOCATE(d_bfc(n))
         ALLOCATE(face_wall_group(n,num_wall_group),cell_closewall(n,num_wall_group))
         ALLOCATE(cell_closewall_indx(n,num_wall_group),dis_closewall(n,num_wall_group))
         c_bface(:)=0
         c_bface_indx(:)=0
         !nwb(:)=0
         d_bfc(:)=0.0d0
         face_wall_group(:,:)=0
         cell_closewall(:,:)=0
         cell_closewall_indx(:,:)=0
         dis_closewall(:,:)=0.0d0   
      ENDIF 
!
!.....RV int_swap
!      
      IF(rv_model.gt.0)THEN
         ALLOCATE(h_il_cfd(n),h_ig_cfd(n),h_gf_cfd(n),vfgl_cfd(n),vfgd_cfd(n))
         h_il_cfd(:)=0.0d0
         h_ig_cfd(:)=0.0d0
         h_gf_cfd(:)=0.0d0
         vfgl_cfd(:)=0.0d0
         vfgd_cfd(:)=0.0d0     
      ENDIF
!
!.....Heat Partitioning Model
!    
      IF(iheatpart.gt.0)THEN
         ALLOCATE(q1cell(n1),qqcell(n1),qclcell(n1),qcgcell(n1),ndensitycell(n1))
         ALLOCATE(t_bulk(n1),t_plus(n1),t_plus_bulk(n1),nsiteden_o(n1),dry_weight(n1))
         q1cell(:)=0.0d0
         qqcell(:)=0.0d0
         qclcell(:)=0.0d0
         qcgcell(:)=0.0d0
         ndensitycell(:)=0.0d0
         t_bulk(:)=0.0d0
         t_plus(:)=0.0d0
         t_plus_bulk(:)=0.0d0   
         nsiteden_o(:)=0.0d0      
         dry_weight(:)=0.0d0          
      ENDIF
      IF(iheatpart.gt.0.or.ncell_cond_all.gt.0)THEN
         ALLOCATE(qecell(n1))
         qecell(:)=0.0d0
      ENDIF
      IF(iheatpart.gt.0.or.rv_htmodel_forCFD.gt.0)THEN
         ALLOCATE(sa_walll(n1))
         sa_walll(:)=0.d0
      ENDIF
!
!.....Fric_face
!     
      IF(fric_face)THEN
         ALLOCATE(face_fr_l(n,ndim),face_fr_g(n,ndim))
         face_fr_l(:,:)=0.0d0
         face_fr_g(:,:)=0.0d0      
      ENDIF
!
!.....Zboron: use in somaFlow.in
!
      ALLOCATE(cboronb_liq(nb_max),cboronb(nb_max))
      cboronb_liq(:)=0.0d0
      cboronb(:)=0.0d0
!
!.....Zmodel: use in calc_momentum
!
      ALLOCATE(drift_c0(n),drift_c1(n))
      drift_c0(:)=1.0d0
      drift_c1(:)=1.0d0
      ALLOCATE(Cb_bubble(n))
      Cb_bubble(:)=0.0d0
      Twall_Model=0
!   
!.....Radiation_Component
!      
      ALLOCATE(qrad_flu(ncell_fluid),qrad_rod(ncell_fluid),qrad_sol(ncell_cond))       
      qrad_flu=0.0d0
      qrad_rod=0.0d0
      qrad_sol=0.0d0
!
!.....resist in SMR
!      
      ALLOCATE(resist(ncell_fluid))
      resist=0.0d0
!      
!.....Convective heat transfer model on the solid wall
!      
      ALLOCATE(htc_convw(10),tb_convw(10),ha_convw(10))
      htc_convw(:)=0.0d0
      tb_convw(:)=0.0d0
      ha_convw(:)=0.0d0 
!      
      RETURN
      END SUBROUTINE allocate_model_var
