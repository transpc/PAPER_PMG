      MODULE Vol_DATA
!
      TYPE CELL_DATA
         SEQUENCE
         REAL(8),DIMENSION(:),ALLOCATABLE :: alphag,alphal,alphad,quala,ha,                &
                                             rhog,rhol,rhod,rhom,rhomr,rhoa,               &
                                             eg,el,ed,p,pps,                               &
                                             hg,hl,hgsat,hlsat,egsat,elsat,quals,          &
                                             tg,tl,td,ts,tst,                              &
                                             drhogdp,drhogde,drhogdx,                      &
                                             drholdp,drholde,                              &
                                             dtgdp,dtgde,dtgdx,                            &
                                             dtldp,dtlde,dtsdp,dtsde,dtsdx,                &
                                             lviscosg,lviscosl,lviscosd,                   &
                                             tviscosg,tviscosl,tviscosd,                   &
                                             eviscosg,eviscosl,eviscosd,                   &
                                             vFgl,vFgd,                                    &
                                             entr,dentr,yeta,                              &
                                             lcondg,lcondl,sigma,betag,betal,cpg,cpl,      &
                                             condg,condl,                                  &
                                             alphag_o,alphal_o,alphad_o,quala_o,           &
                                             eg_o,el_o,ed_o,p_o,                           &
                                             tg_o,tl_o,td_o,ts_o,                          &
                                             aint1,aint2,aint3,D1,D2,Ddepart,Dlift,        &
                                             limiter,                                      &
                                             vfwg,vfwl,                                    &
                                             rhog_o,                                       &
                                             cboron,cboron_o,rhol_o,                       &
                                             twall,                                        &
                                             fwkl,fwkg,                                    &
                                             estm,estm_o,pps_o,                            &
!         
                                             vfwg_x,vfwg_y,vfwg_z,vfwl_x,vfwl_y,vfwl_z,    &
!
!........Natural Convection: Twall-bc, Buoyancy Coefficient(Rodi.)
!
                                             T_top,T_bot,                                  &
                                             ced33,mdiff,                                  &
!
!........Film model variables
!
                                             film_thickness,film_shear,                    &
!
!======================================================================
! CUPID-RV variable
!======================================================================
!
!........Regime criteria
!
                                             wf_dry,wf_VST,                                &
                                             alpha_bs,alpha_de,alpha_sa,alpha_cd,alpha_gs, &
!
!........Interfacial area
!
                                             ia_bubbly,                                    &
                                             ia_slug_tb,ia_slug_sb,                        &
                                             ia_churn,                                     &
                                             ia_annular_drp,ia_annular_ann,                &
                                             ia_mpr,                                       &
                                             ia_invann_ann,ia_invann_sb,                   &
                                             ia_invchn,                                    &
                                             ia_invslg_drp,ia_invslg_ann,                  &
                                             ia_mist,                                      &
                                             ia_mpo,                                       &
                                             ia_VST,                                       &
                                             ia_vst_st,ia_vst_sb,                          &
!
!........Bubble diameter
!
                                             dbb,                                          &      ! bubble dia
                                             dsb,dtb,                                      &      ! slug dia
                                             ddrp,                                         &      ! drop dia
!         dbmin,dbmax
!         ddmin,ddmax
!         
!........cell length, hyrdaulic diameter -> hydraulicD(:) in Zconst2.f90
!
                                             length,                                       &
!
!........flow direction
!         1: upward
!         2: downward
!         3: count-current
!
                                             fdir,                                         &
!
!........Drift flux model
!
                                             c1,c0,                                        &
!         
!........face-upwind void fraction for RV interfacial friction model
!
                                             alphagf,alphalf         
!
         INTEGER,DIMENSION(:),ALLOCATABLE :: regime,idummyV, &
                                             vst
!               
!........reflood flag
!
         LOGICAL,DIMENSION(:),ALLOCATABLE :: ireflod
!
      ENDTYPE
!      
      TYPE(CELL_DATA) :: cell
!
      END MODULE
