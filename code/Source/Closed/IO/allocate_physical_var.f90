!
      SUBROUTINE allocate_physical_var
!
!     This routine allocates variables containing physical information.
!
      USE Zmpi          , ONLY: ncell_fp
      USE Zzone         , ONLY: ncell_fluid
      USE Zparam        , ONLY: ndim,nb_max
      USE Zare          , ONLY: ar_liq,ar_drp,ar_gas,are_liq,are_drp,are_gas
      USE Zconst2       , ONLY: gfactor,i_repeat,iprn,dt,dt_old, stime_hup,stime_hflat, &
                                 imomwall
      USE Zdaint_ag     , ONLY: daint1_ag,daint2_ag,daint1_ag_bc,daint2_ag_cm,          &
                                 aint_01b,aint_bc,aint_09b,aint_01d,aint_cm,aint_09d,   &
                                 D1_01,D1_bc,D1_09,D2_cm,D2_09
      USE Zdel_scalar   , ONLY: del_eg,del_el,del_x,del_ag,del_ad,del_rhog,del_rhol,err_stm
      USE Zdhda         , ONLY: dHldag,dHgdag,dHfgdag,dHldtl,dHgdtg,dHfgdtg
      USE Zenergy_diff  , ONLY: ediff_liq,ediff_gas,ediff_drp
      USE Zm_src        , ONLY: src_liq,src_gas,src_drp
      USE Zmass_diff    , ONLY: mdiff_gas,mediff_gas       
      USE Zpress        , ONLY: p,pp,pp_o,dpdx,dpdx_o
      USE Zpress_coeff  , ONLY: coefp_l,coefp_d,coefp_g,coefm_l,coefm_g
      USE Zqvol         , ONLY: q0_gas,q0_drp,q0_liq,qvol_liq,qvol_drp,qvol_gas,h_il,h_ig,h_gf,hgf_o, &
                                 gamma,gamma_wall,qporous_liq,qporous_gas,qporous_gamma,q0_ice_solid
      USE Zscalar_coeff , ONLY: sb
      USE Zscalar_src   , ONLY: ha_dp_c,ha_dt_c
      USE Zvector       , ONLY: vrel_o,ul_o,ug_o,                            &
                                vl_n,vg_n,vd_n,vl_o,vg_o,vd_o,vl_t,vd_t,vg_t
      USE Zvoid         , ONLY: dagdx,gamma_void
      USE Zzone         , ONLY: num_max_zone
      USE Zqvol         , ONLY: qwall_solid,qrv_liq,qrv_gas,qrv_gamma
!
      IMPLICIT NONE
!
      INTEGER n,n1
!
      n=ncell_fp
      n1=ncell_fluid
!
!.....Zvector
!
      ALLOCATE(     &
               vrel_o(n),ul_o(n),ug_o(n))
      ALLOCATE(vl_n(n,ndim),vd_n(n,ndim),vg_n(n,ndim),     &
               vl_o(n,ndim),vd_o(n,ndim),vg_o(n,ndim),     &
               vl_t(n,ndim),vd_t(n,ndim),vg_t(n,ndim))
      vl_n(:,:)=0.0d0
      vd_n(:,:)=0.0d0
      vg_n(:,:)=0.0d0
      vl_o(:,:)=0.0d0
      vd_o(:,:)=0.0d0
      vg_o(:,:)=0.0d0
      vl_t(:,:)=0.0d0
      vd_t(:,:)=0.0d0
      vg_t(:,:)=0.0d0
      vrel_o(:)=0.0d0
      ul_o(:)=0.0d0
      ug_o(:)=0.0d0
!
!.....Zpress
!
      ALLOCATE(p(n),dpdx(n,ndim),pp(n),pp_o(n))
      ALLOCATE(dpdx_o(n,ndim))
      p(:)=0.0d0
      dpdx(:,:)=0.0d0
      pp(:)=0.0d0
      pp_o(:)=0.0d0
      dpdx_o(:,:)=0.0d0      
!
!.....Zpress_coeff
!
      ALLOCATE(coefp_l(n),coefp_d(n),coefp_g(n))
      coefp_l(:)=0.0d0
      coefp_d(:)=0.0d0
      coefp_g(:)=0.0d0
      ALLOCATE(coefm_l(n),coefm_g(n))
      coefm_l(:)=0.0d0
      coefm_g(:)=0.0d0
!
!  Zscalar_coeff
!
      ALLOCATE(sb(n,6))
      sb(:,:)=0.0d0
!
!.....Zscalar_src
!
      ALLOCATE(Ha_dp_c(n),Ha_dt_c(n))
      Ha_dp_c(:)=0.0d0
      Ha_dt_c(:)=0.0d0
!
!.....Zm_src
!
      ALLOCATE(src_liq(n1,ndim),src_gas(n1,ndim),src_drp(n1,ndim))
      src_liq(:,:)=0.0d0
      src_gas(:,:)=0.0d0
      src_drp(:,:)=0.0d0      
!
!.....Zenergy_diff
!
      ALLOCATE(ediff_liq(n),ediff_gas(n),ediff_drp(n))
      ediff_liq(:)=0.0d0
      ediff_gas(:)=0.0d0
      ediff_drp(:)=0.0d0
!
!.....Zenergy_conv
!
!
!.....Zdel_scalar
!
      ALLOCATE(del_eg(n),del_el(n),del_x(n),del_ag(n),del_ad(n),del_rhog(n),del_rhol(n))
      del_eg(:)=0.0d0
      del_el(:)=0.0d0
      del_x(:)=0.0d0
      del_ag(:)=0.0d0
      del_ad(:)=0.0d0
      del_rhog(:)=0.0d0
      del_rhol(:)=0.0d0
      ALLOCATE(err_stm(n))
      err_stm(:)=0
!
!.....Zdhda
!
      ALLOCATE(dHldag(n),dHgdag(n),dHfgdag(n),dHldtl(n),dHgdtg(n),dHfgdtg(n))
      dHldag(:)=0.0d0
      dHgdag(:)=0.0d0
      dHfgdag(:)=0.0d0
      dHldtl(:)=0.0d0
      dHgdtg(:)=0.0d0
      dHfgdtg(:)=0.0d0
!
!.....Zare
!
      ALLOCATE(ar_liq(n),ar_drp(n),ar_gas(n),are_liq(n),are_drp(n),are_gas(n))
      ar_liq(:)=0.0d0
      ar_drp(:)=0.0d0
      ar_gas(:)=0.0d0
      are_liq(:)=0.0d0
      are_drp(:)=0.0d0
      are_gas(:)=0.0d0
!
!.....Zmass_conv
!
!
!.....Zvoid
!
      ALLOCATE(dagdx(n1,ndim),gamma_void(n))
      dagdx(:,:)=0.0d0
      gamma_void(:)=0.0d0
!
!.....Zdaint_ag
!
      ALLOCATE(daint1_ag(n),daint2_ag(n),daint1_ag_bc(n),daint2_ag_cm(n))
      ALLOCATE(aint_01b(n),aint_bc(n),aint_09b(n),aint_01d(n),aint_cm(n),aint_09d(n))
      ALLOCATE(D1_01(n),D1_bc(n),D1_09(n),D2_cm(n),D2_09(n))
      daint1_ag(:)=0.0d0
      daint1_ag_bc(:)=0.0d0
      daint2_ag(:)=0.0d0
      daint2_ag_cm(:)=0.0d0
      aint_01b(:)=0.0d0
      aint_bc(:)=0.0d0
      aint_09b(:)=0.0d0
      aint_01d(:)=0.0d0
      aint_cm(:)=0.0d0
      aint_09d(:)=0.0d0
      D1_01(:)=0.0d0
      D1_bc(:)=0.0d0
      D1_09(:)=0.0d0
      D2_cm(:)=0.0d0
      D2_09(:)=0.0d0
!
!.....Zqvol
!
      ALLOCATE(q0_gas(num_max_zone),q0_drp(num_max_zone),q0_liq(num_max_zone), &
                qvol_liq(n),qvol_drp(n),qvol_gas(n),                           &
                H_il(n),H_ig(n),gamma(n),H_gf(n),gamma_wall(n),hgf_o(n))
      ALLOCATE(qporous_liq(n),qporous_gas(n),qporous_gamma(n))
      ALLOCATE(qrv_liq(n),qrv_gas(n),qrv_gamma(n))
      ALLOCATE(q0_ice_solid(num_max_zone))
      ALLOCATE(qwall_solid(nb_max))
      q0_gas(:)=0.0d0
      q0_drp(:)=0.0d0
      q0_liq(:)=0.0d0
      qvol_liq(:)=0.0d0
      qvol_drp(:)=0.0d0
      qvol_gas(:)=0.0d0
      H_il(:)=0.0d0
      H_ig(:)=0.0d0
      gamma(:)=0.0d0
      H_gf(:)=0.0d0
      gamma_wall(:)=0.0d0
      qporous_liq(:)=0.0d0
      qporous_gas(:)=0.0d0
      qporous_gamma(:)=0.0d0
      qrv_liq(:)=0.0d0
      qrv_gas(:)=0.0d0
      qrv_gamma(:)=0.0d0
      q0_ice_solid(:)=0.0d0
      qwall_solid(:)=0.d0
      hgf_o(:)=0.0d0
!
!.....Zconst2
!
      ALLOCATE(gfactor(n))
      i_repeat=0
      iprn=0
      dt=0.0d0
      dt_old=0.0d0
      stime_hup=0.0d0
      stime_hflat=0.0d0
      imomwall=0.0d0
      gfactor(:)=0.0d0
!
!     Implement Zmass_diff 2015.07.29 JHLee (SNU)
!
      Allocate(mdiff_gas(n),mediff_gas(n))
      mdiff_gas(:)=0.0d0
      mediff_gas(:)=0.0d0
!
      RETURN
      END SUBROUTINE allocate_physical_var
!
