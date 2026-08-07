!
      SUBROUTINE allocate_preliminary_var(ncell_fluid,ncell_fp)
!
!     This routine allocates variables containing geomeric information.
!
      USE Zmpi         , ONLY: maxmt_fluid,maxmt_fp
      USE Zparam       , ONLY: ndim,nb_max,mesh_openfoam
      USE Zbc_index    , ONLY: npb,index_flux,index_property
      USE Zb_condition , ONLY: pbnd,p_fb,vb_liq,vb_drp,vb_gas,vin_liq,vin_drp,vin_gas,  &
                                cb_pl,cb_pd,cb_pg,cb_p,                                 &
                                eb_liq,eb_drp,eb_gas,tb_liq,tb_drp,tb_gas,              &
                                qwall_liq,qwall_drp,qwall_gas,twall,                    &
                                rhob_liq,rhob_drp,rhob_gas,                             &
                                alphab_liq,alphab_drp,alphab_gas,                       &
                                alpha_liq_nd,t_liq_nd,rho_liq_nd,e_liq_nd,              &
                                alpha_gas_nd,t_gas_nd,rho_gas_nd,e_gas_nd,              &
                                alpha_drp_nd,t_drp_nd,rho_drp_nd,e_drp_nd,              &
                                qualab,quala_nd,                                        &
                                turb_keb,turb_dpb,turb_kegb,turb_dpgb,                  & 
                                vb_lold,vb_gold,v_wall,lvisb_liq,lvisb_gas
      USE Zgradoption  , ONLY: ifrink
      USE Zconst2      , ONLY: hydraulicd,sl       !LSJ 161123 porous
      USE Zcoord1      , ONLY: xloc
      USE Zcoord2      , ONLY: xfc,xfc_min,xfc_max,xloc_xfc_min,xloc_xfc_max, &
                               xloc_xfc_radius_min,xloc_xfc_radius_max
      USE Zcoord3      , ONLY: sv,vol,porosity,permeability
      USE Zncg         , ONLY: ncg_species,qn_cell0,qn_nvin,qn_npin
      USE Znode        , ONLY: nd_max,nmax_vertex,num_cell_node,cell_node,node_face_cell, &
                               num_nd,rwcn,dxr
      USE Zporous      , ONLY: sgap
      USE Zzone        , ONLY: nzone,nmaterial,icore
      USE Zuserdefined , ONLY: user_iary,user_rary
!
      IMPLICIT NONE
!
!     input
      INTEGER ncell_fluid,ncell_fp
!     local variable 
      INTEGER n,n1
!
      n=ncell_fp
      n1=ncell_fluid
!
!.....Zconst2
!
!    bug  read_grid need sl to be ndim line 401 copied from sl_tmp(nn,ndim)
      ALLOCATE(hydraulicd(n),sl(n1,ndim),sgap(n,ndim))
      hydraulicd(:)=0.0d0 
      sl(:,:)=0.0d0       
      sgap(:,:)=0.0d0
!
!.....Zcoord1
!
      ALLOCATE(xloc(n,ndim))
      xloc(:,:)=0.0d0
!
!.....Zcoord2
!
      IF(maxmt_fluid.gt.0) THEN
         ALLOCATE(xfc(maxmt_fluid,ndim))
         ALLOCATE(xfc_min(maxmt_fluid,ndim),xfc_max(maxmt_fluid,ndim),xloc_xfc_min(maxmt_fluid,ndim),xloc_xfc_max(maxmt_fluid,ndim))
      ELSE
!     when maxmt_fluid=0 on some procs avoid troubles with check bounds
         ALLOCATE(xfc(1,ndim))
         ALLOCATE(xfc_min(1,ndim),xfc_max(1,ndim),xloc_xfc_min(1,ndim),xloc_xfc_max(1,ndim))
      ENDIF
      ALLOCATE(xloc_xfc_radius_min(ncell_fluid),xloc_xfc_radius_max(ncell_fluid))
      xfc(:,:)=0.0d0
!
!.....Zcoord3
!
      IF(maxmt_fluid.gt.0) THEN
         ALLOCATE(sv(maxmt_fluid,ndim),permeability(maxmt_fluid))
      ELSE
!     when maxmt_fluid=0 on some procs avoid troubles with check bounds
         ALLOCATE(sv(1,ndim),permeability(1))
      ENDIF
      ALLOCATE(vol(n),porosity(n))
      sv(:,:)=0.0d0
      vol(:)=0.0d0
      porosity(:)=0.0d0
      permeability(:)=0.0d0
!
!.....Znormal
!
!
!.....Zzone
!
      ALLOCATE(nzone(n),nmaterial(n))
      nzone(:)=0
      nmaterial(:)=0.0d0
!
!.....Zbc_index
!
      ALLOCATE(npb(n))
      npb(:)=0
!
!.....Zb_condition
!
      ALLOCATE(pbnd(nb_max),p_fb(nb_max),                                                 &
                vb_liq(nb_max,ndim),vb_drp(nb_max,ndim),vb_gas(nb_max,ndim), &
                vin_liq(nb_max),vin_drp(nb_max),vin_gas(nb_max),                           &
                cb_pl(nb_max),cb_pd(nb_max),cb_pg(nb_max),cb_p(nb_max),                    &
                eb_liq(nb_max),eb_drp(nb_max),eb_gas(nb_max),                              &
                tb_liq(nb_max),tb_drp(nb_max),tb_gas(nb_max),                              &
                qwall_liq(nb_max),qwall_drp(nb_max),qwall_gas(nb_max),twall(nb_max),       &
                rhob_liq(nb_max),rhob_drp(nb_max),rhob_gas(nb_max),                        &
                alphab_liq(nb_max),alphab_drp(nb_max),alphab_gas(nb_max),                  &
                alpha_liq_nd(nb_max),t_liq_nd(nb_max),rho_liq_nd(nb_max),e_liq_nd(nb_max), &
                alpha_gas_nd(nb_max),t_gas_nd(nb_max),rho_gas_nd(nb_max),e_gas_nd(nb_max), &
                alpha_drp_nd(nb_max),t_drp_nd(nb_max),rho_drp_nd(nb_max),e_drp_nd(nb_max), &
                qualab(nb_max),quala_nd(nb_max),                                           &
                turb_keb(nb_max),turb_dpb(nb_max),turb_kegb(nb_max),turb_dpgb(nb_max),     &
                vb_lold(nb_max,ndim),vb_gold(nb_max,ndim),v_wall(ndim),                    &
                lvisb_liq(nb_max),lvisb_gas(nb_max))
      pbnd(:)=0.0d0
      p_fb(:)=0.0d0
      vb_liq(:,:)=0.0d0
      vb_drp(:,:)=0.0d0
      vb_gas(:,:)=0.0d0
      vin_liq(:)=0.0d0
      vin_drp(:)=0.0d0
      vin_gas(:)=0.0d0
      cb_pl(:)=0.0d0
      cb_pd(:)=0.0d0
      cb_pg(:)=0.0d0
      cb_p(:)=0.0d0
      eb_liq(:)=0.0d0
      eb_drp(:)=0.0d0
      eb_gas(:)=0.0d0
      tb_liq(:)=0.0d0
      tb_drp(:)=0.0d0
      tb_gas(:)=0.0d0
      qwall_liq(:)=0.0d0
      qwall_drp(:)=0.0d0
      qwall_gas(:)=0.0d0
      twall(:)=0.0d0
      rhob_liq(:)=0.0d0
      rhob_drp(:)=0.0d0
      rhob_gas(:)=0.0d0
      alphab_liq(:)=0.0d0
      alphab_drp(:)=0.0d0
      alphab_gas(:)=0.0d0
      alpha_liq_nd(:)=0.0d0
      t_liq_nd(:)=0.0d0
      rho_liq_nd(:)=0.0d0
      e_liq_nd(:)=0.0d0
      alpha_gas_nd(:)=0.0d0
      t_gas_nd(:)=0.0d0
      rho_gas_nd(:)=0.0d0
      e_gas_nd(:)=0.0d0
      alpha_drp_nd(:)=0.0d0
      t_drp_nd(:)=0.0d0
      rho_drp_nd(:)=0.0d0
      e_drp_nd(:)=0.0d0
      qualab(:)=0.0d0
      quala_nd(:)=0.0d0
      turb_keb(:)=0.0d0
      turb_dpb(:)=0.0d0
      turb_kegb(:)=0.0d0
      turb_dpgb(:)=0.0d0
      vb_lold(:,:)=0.0d0
      vb_gold(:,:)=0.0d0
      v_wall(:)=0.0d0
      lvisb_liq(:)=0.0d0
      lvisb_gas(:)=0.0d0
!
!.....Zncg
!
      ALLOCATE(ncg_species(8),qn_cell0(8),qn_nvin(nb_max,8),qn_npin(nb_max,8)) 
      ncg_species(:)=0.0d0
      qn_cell0(:)=0.0d0
      qn_nvin(:,:)=0.0d0
      qn_npin(:,:)=0.0d0
!
!.....Znode
!
      IF(mesh_openfoam.eq.1.and.ifrink.ge.1)THEN
         ALLOCATE(num_cell_node(n),cell_node(nd_max,n),node_face_cell(nmax_vertex,maxmt_fp),dxr(nd_max,n,ndim))
         ALLOCATE(num_nd(maxmt_fp),rwcn(nd_max,n))
         num_cell_node(:)=0
         cell_node(:,:)=0
         node_face_cell(:,:)=0
         num_nd(:)=0
         rwcn(:,:)=0.0d0
         dxr(:,:,:)=0.0d0
      ENDIF
!
!.....rv
!
      ALLOCATE(icore(n))
      icore(:)=0
!    
!.....user array for restart
!  
      user_rary(:)=-1.0d0
      user_iary(:)=-1
!            
!.....index for vector flux
!         
      ALLOCATE(index_flux(nb_max),index_property(nb_max))
      index_flux(:)=0   
      index_property(:)=0
!      
      RETURN
      END SUBROUTINE allocate_preliminary_var
