!
      MODULE Zrestart
!      
!.....Typed data
!
      IMPLICIT NONE
      SAVE
!
      TYPE CELL_DATA_rst
         REAL(8) alphag,alphal,alphad,quala
         REAL(8) rhog,rhol,rhod,rhom
         REAL(8) eg,el,ed,p,pps
         REAL(8) hg,hl,hgsat,hlsat,egsat,elsat,quals,ha
         REAL(8) tg,tl,td,ts
         REAL(8) drhogdp,drhogde,drhogdx
         REAL(8) drholdp,drholde
         REAL(8) dtgdp,dtgde,dtgdx
         REAL(8) dtldp,dtlde,dtsdp,dtsde,dtsdx
         REAL(8) lviscosg,lviscosl,lviscosd  
         REAL(8) tviscosg,tviscosl,tviscosd 
         REAL(8) eviscosg,eviscosl,eviscosd  
         REAL(8) vFgl,vFgd
         REAL(8) entr,dentr,yeta
         REAL(8) lcondg,lcondl,sigma,betag,betal,cpg,cpl
         REAL(8) condg,condl  
         REAL(8) alphag_o,alphal_o,alphad_o,quala_o
         REAL(8) eg_o,el_o,ed_o,p_o
         REAL(8) tg_o,tl_o,td_o,ts_o
         REAL(8) aint1,aint2,aint3,D1,D2,Ddepart,Dlift
         INTEGER regime 
         REAL(8) limiter
         REAL(8) vfwg,vfwl
         REAL(8) vfwg_x,vfwg_y,vfwg_z,vfwl_x,vfwl_y,vfwl_z         
         REAL(8) rhog_o
         REAL(8) cboron,cboron_o,rhol_o
         REAL(8) twall             
         REAL(8) fwkl,fwkg 
         REAL(8) estm,estm_o,pps_o
      ENDTYPE
!      
      TYPE(CELL_DATA_rst),ALLOCATABLE::cell_rst(:)
!
      TYPE WAL_DATA_rst 
         REAL(8) twall_partition
         REAL(8) wall_fluxl_diff,wall_fluxg_diff,wall_fluxd_diff
         REAL(8) ddepartw,ratio_evap
      ENDTYPE
!      
      TYPE(WAL_DATA_rst),ALLOCATABLE::face_rst(:)
!
      TYPE SOL_DATA_rst 
         REAL(8) tsol
         REAL(8) tsol_o
         REAL(8) tsol_max,tpellet_surf
      ENDTYPE
!      
      TYPE(SOL_DATA_rst),ALLOCATABLE::solid_rst(:)
!
!.....Normal cell data
!
      REAL(8),ALLOCATABLE :: p_rst(:)
      REAL(8),ALLOCATABLE :: gamma_rst(:)
      REAL(8),ALLOCATABLE :: turb_ke_o_rst(:),turb_dp_o_rst(:),turb_keg_o_rst(:),turb_dpg_o_rst(:)
      REAL(8),ALLOCATABLE :: turb_ke_rst(:),turb_dp_rst(:),turb_keg_rst(:),turb_dpg_rst(:)
      REAL(8),ALLOCATABLE :: Cwlf_rst(:),Clift_rst(:),Ctd_rst(:)
      REAL(8),ALLOCATABLE :: H_il_rst(:),H_ig_rst(:),H_gf_rst(:)
      REAL(8),ALLOCATABLE :: del_eg_rst(:),del_el_rst(:)
      REAL(8),ALLOCATABLE :: del_x_rst(:)
      REAL(8),ALLOCATABLE :: del_ag_rst(:),del_ad_rst(:)
      REAL(8),ALLOCATABLE :: del_rhog_rst(:),del_rhol_rst(:)
      REAL(8),ALLOCATABLE :: iat_nucl_rst(:)
      REAL(8),ALLOCATABLE :: solid_tmp1(:),solid_tmp2(:)
      REAL(8),ALLOCATABLE :: solid_tmp3(:),solid_tmp4(:)
      REAL(8),ALLOCATABLE :: qqcell_rst(:),qecell_rst(:),qclcell_rst(:),qcgcell_rst(:)     
!
!.....i,j array data
!
      REAL(8),ALLOCATABLE :: fluxvol_g_rst(:,:)
      REAL(8),ALLOCATABLE :: fluxvol_l_rst(:,:)
      REAL(8),ALLOCATABLE :: fluxvol_d_rst(:,:)
!
!.....i,ix array data
!
      REAL(8),ALLOCATABLE :: vg_n_rst(:,:)
      REAL(8),ALLOCATABLE :: vl_n_rst(:,:)
      REAL(8),ALLOCATABLE :: vd_n_rst(:,:)
      REAL(8),ALLOCATABLE :: vg_o_rst(:,:)
      REAL(8),ALLOCATABLE :: vl_o_rst(:,:)
      REAL(8),ALLOCATABLE :: vd_o_rst(:,:)
!
      ENDMODULE Zrestart