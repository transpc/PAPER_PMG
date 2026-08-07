      SUBROUTINE read_flow(nout)
!      
      USE VOL_DATA        , ONLY: cell
      USE WALL_DATA       , ONLY: face
      USE SOLID_DATA      , ONLY: solid
      USE Zmpi            , ONLY: ncell_fp,ncell_ps,i3perm,metis
      USE Zzone           , ONLY: num_fzone,num_pzone,num_szone,ncell_fluid,ncell_cond,num_max_zone,ncell_fluid_all
      USE Zcore           , ONLY: np,myrank
      USE Zparam          , ONLY: nin_max,nb_max,ndim,cmu,nb_max,nin_max,nb_sym,nb_mars
      USE STM_TBL_cupid   , ONLY: nfluid
      USE Z2nd_order      , ONLY: mass_conv_2nd,mom_conv_2nd,eng_conv_2nd,qula_conv_2nd,      &
                                  ncg_conv_2nd,boron_conv_2nd,turb_conv_2nd,faclim,vkt,       &
                                  limiter_mass,limiter_mom,limiter_eng,limiter_xn,limiter_qn, &
                                  limiter_boron,limiter_turb
      USE Zb_condition    , ONLY: vb_liq,vb_gas,alphab_gas,qualab,p_fb,tb_liq,tb_gas,vb_drp,     &
                                  vin_liq,vin_gas,vin_drp,pbnd,t_liq_nd,t_gas_nd,alpha_gas_nd,   &
                                  quala_nd,twall,alphab_liq,rhob_liq,rhob_gas,cb_pl,cb_pg,cb_pd, &
                                  turb_keb,turb_dpb,turb_kegb,turb_dpgb,lvisb_liq,lvisb_gas
      USE Zbc_index       , ONLY: nvin,vin_norm,npin,num_wall_group,select_wall_group,vin_mfr, &
                                  i_horizontal_outlet
      USE Zbicg           , ONLY: psolve,eps_bicg,max_bicg,min_bicg,relax_u,relax_p, &
                                  lev_type,levmpi_type,lev,levmpi,                   &
                                  lev_type_c,levmpi_type_c,lev_c,levmpi_c
      USE Zboron          , ONLY: cboronb,cboronb_liq
      USE Zconst1         , ONLY: vv_prob,restart,save_option,noutput,restart,                          &
                                  iprofbc,parallel,mtopol,mdrag,mhtc,iat,iturb,ntdf,nlift,turb_phase,   &
                                  nwlf,iheatpart,mboron,mdiffoff,cplmars,cplmaster,lsquareoff,          &
                                  turbubble,turboil,topolsurface,buoyancy_turb,lowreynolds,             &
                                  mdiffscheme,fric_face,nd_face,p_work_face,wconden,iVisRatio,          &
                                  vis_ratio,lwconden_alphal0,lrestart_changed_nbcon,lrestart_overwrite, &
                                  i_macroturb_source,i_turb_zero,i_wall_fric,i_turb_disp,i_wall_lub,    &
                                  i_bubble_diameter,i_drop_diameter,i_rv_int_fric,i_critical_flow,      &
                                  i_subchannel_fric,i_subchannel_fric_axial,i_subchannel_fric_cross,    &
                                  i_subchannel_mixing,i_2p_multiplier,rv_htmodel_forCFD
      USE Zconst2         , ONLY: iprn,grav,dt,lgravity,                      &
                                  stime_hup,stime_hflat,stime_vup,stime_vflat
      USE Zcoord1         , ONLY: xloc
      USE Zdecoupled      , ONLY: al_min_c,ag_min_c,al_min_e,ag_min_e
      USE Zface           , ONLY: Kepsilon_RNG, Kepsilon_real
      USE Zflowregime     , ONLY: vFgl_1_min,vFgl_2_min,vFgl_3_min
      USE Zgradoption     , ONLY: iavgtype,irc_damp,ifrink,iter_grad,non_orth_diff,non_orth,grav_grad, &
                                  non_orth_eng,non_orth_turb,non_orth_iter
      USE Zheat_partition , ONLY: kfactor     
      USE Ziat            , ONLY: dbubble_init,                         &
                                  r_db_min,r_db_max,r_ddrop,r_dh_hibiki
      USE Zncg            , ONLY: n_ncg_sp,imp_ncg,ncg_species,qn_cell0,qn_nvin,qn_npin,ncg_diff, &
                                  i_ncg_vis
      USE Zmodel          , ONLY: drag_coeff_a,drag_coeff_b,h_il_coeff_a,h_il_coeff_b,h_il_min,     &
                                  h_fg_min,h_ig_min,h_ig_coeff_a,h_ig_coeff_b,h_fg_coeff_b,         &
                                  h_fg_coeff_a,cd_min_user,cd_min_ag99,h_il_min_user,h_il_min_ag99, &
                                  dtl,dtg
      USE Zimplicit       , ONLY: imp_mom_conv,imp_mom_diff,imp_scalar_conv,imp_scalar_diff,imp_ke_diff, &
                                  max_iter_mom,ag_min_m,al_min_m,eps_imp_mom,imp_alpha,eps_imp_alpha,    &
                                  max_iter_alpha,eps_imp_scalar,max_iter_scalar,iter_scalar,             &
                                  skip_imp_scalar,iter_mom,imp_ke_conv,                                  &
                                  eps_imp_ke,eps_imp_dp,max_iter_ke,max_iter_dp,                         &
                                  imp_boron_trans,max_iter_boron,eps_imp_boron
      USE Zndforce        , ONLY: relax_hik,relax_cd
      USE Zqvol           , ONLY: q0_gas,q0_liq,q0_ice_solid,qwall_solid,&
                                  htc_convw, tb_convw, ha_convw
      USE Zrv_model       , ONLY: rv_model,rv_fw_reg,rv_ht_str,rv_ht_w,rv_fric_i,rv_ht_i,rv_fric_w, &
                                  rv_choke,rv_gapcond,lfric_swap,lhtc_swap,lfricw_swap,rv_mcp,rv_valve,ia_option
      USE Zscalar_coeff   , ONLY: l_th_equil,l_min_hik,alphag_min,alphal_min
      USE Ztimecon        , ONLY: ctrl_opt,cfl_ratio,t_end,toutstep,cfl_ratio_max,smac,dt_opt, &
                                  dt_max,nctrl,t_end_ctrl,toutstep_ctrl,cfl_ratio_max_ctrl,    &
                                  smac_ctrl,dt_opt_ctrl,dt_max_ctrl,alpha_min,treststep,       &
                                  treststep_ctrl,itim,iso_thermal
      USE Zturbzeq        , ONLY: tlengs
      USE Zmodel          , ONLY: rad_model,abs_coeff,wall_emiss,max_iter_rad,eps_imp_rad, &
                                  i_droplet,i_weight,i_fs_temp_intpol
      USE Zwall_HTC       , ONLY: inline_bundle,reflood,                                            &
                                  hyd_core,l_plate,dia_rod,pit_dia,h_bundle,base_bundle,dia_rod_cfd
      USE Zporous         , ONLY: vfporous
!.....added by pik      
      USE Zrv_wall_fric   , ONLY: rough_wall   
      USE Zporous         , ONLY: l_subchannel,l_mixing_vane,l_spacer_grid,l_2p_multiplier_evvd, &
                                  udfi_subchannel_flowdir,&
                                  kloss_grid,kloss_cross,ftm,ka,beta
      USE viewData_common , ONLY: viewField
      USE Ztplot          , ONLY: tplot_num,tplot_cell,tplot_dt,tplot_prop,tplot_x,tplot_y,tplot_z, &
                                  tplot_cell_loc,tplot_cell_rank     
      USE Zio_unit        , ONLY: unit_somaflow,unit_tplotv,unit_tplots,unit_log,unit_rv,unit_chn
      USE unitManager     , ONLY: createUnit
      USE Zdel_scalar     , ONLY: limit_eng_src_opt,limit_iht_opt,suspend_iht_opt,suspend_erg_opt,smac3_pres_eng, &
                                  prn_div_eng,ag_min_hig,al_min_hil,qu_min_hgf,dsrc, &
                                  stmtbl_repeat_for_nc,max_ihtc_opt,relax_interface_dtemp,max_ihtc_opt_coeff
      USE Zrv_choke       , ONLY: choke_throat_area,fzone_throat,choke_pout,env_press_option,relax_choke,time_cflow_on
      USE Zmcp            
      USE Zvalve
      USE Zpdrop          , ONLY: npid,time_pid_on,time_pid_off,num_dp_region       
!
      IMPLICIT NONE
!      
!.....Input
      INTEGER :: nout   
!.....Local variables
      INTEGER :: i,ix,ii,err
      INTEGER :: n_constTwall,n_constqwall      
      INTEGER :: izone
      INTEGER :: num_fzone_count,num_szone_count
      REAL(8) :: vb_liq_size,vb_gas_size
      REAL(8) :: cvm,c_vm,arhob_gas,arhob_liq,arhob_drp
      REAL(8) :: rhob_m,rhob_vm,denom
      REAL(8) :: tsol0(num_max_zone),sumGravity      
      REAL(8) :: tsol0_temp(num_max_zone),q0_ice_solid_temp(num_max_zone)
      CHARACTER*50 :: f_tplotv,f_tplots
      CHARACTER*20 :: fn
!.....Local allocatable arrays      
      REAL(8),DIMENSION(:),ALLOCATABLE :: p0,tl0,                    &
                                          tg0,quala0, a0_g, cboron0, &
                                          distance,distance_min
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: v0,xloc_all
      
      NAMELIST /problem_description/ vv_prob, ndim, num_fzone_count, num_szone_count,nfluid
      NAMELIST /initial_condition/ grav, v0, p0, tl0, tg0, a0_g, quala0, cboron0, q0_gas, q0_liq, &
         tsol0, q0_ice_solid,&      
         stime_hup,stime_hflat,stime_vup,stime_vflat 
      NAMELIST /boundary_condition/ nvin, vin_norm, vb_liq, vb_gas, alphab_gas, qualab, cboronb, p_fb, tb_liq, tb_gas, vin_mfr, &
         npin, pbnd, t_liq_nd, t_gas_nd, alpha_gas_nd, quala_nd,i_horizontal_outlet,                                            &
         n_constTwall, twall, n_constqwall, qwall_solid,&
         htc_convw, tb_convw, ha_convw
      !NAMELIST /boundary_condition/ nvin, vin_norm, vb_liq, vb_gas, alphab_gas, qualab, cboronb, p_fb, tb_liq, tb_gas, vin_mfr, &
      !   npin, pbnd, t_liq_nd, t_gas_nd, alpha_gas_nd, quala_nd,i_horizontal_outlet,                                            &
      !   n_constTwall, n_constqwall !ÀÌ½ÂÁØ      
      NAMELIST /n_scheme/ parallel, psolve, metis, &
         lev,lev_type,levmpi,levmpi_type,          &
         lev_c,lev_type_c,levmpi_c,levmpi_type_c,  &
         eps_bicg, max_bicg, min_bicg,             &
         iavgtype, irc_damp, grav_grad, iter_grad, &
         non_orth, non_orth_iter, non_orth_diff, non_orth_eng, non_orth_turb, &
         p_work_face, fric_face, nd_face,                                     &
         relax_u, relax_p,                                                    &
         alpha_min, lsquareoff
      NAMELIST /imp_scheme/ imp_mom_conv, imp_mom_diff, imp_alpha, imp_scalar_conv, imp_scalar_diff, imp_ke_diff,          &
         iter_scalar,skip_imp_scalar,iter_mom,                                                                             &
         stmtbl_repeat_for_nc,max_ihtc_opt,max_ihtc_opt_coeff,relax_interface_dtemp,                                       & 
         limit_eng_src_opt,limit_iht_opt,suspend_iht_opt,suspend_erg_opt,smac3_pres_eng,prn_div_eng,ag_min_hig,al_min_hil, &
         qu_min_hgf,dsrc,al_min_c,ag_min_c,al_min_e,ag_min_e,al_min_m,ag_min_m,                                            &
         eps_imp_mom,eps_imp_alpha,eps_imp_scalar,eps_imp_ke,eps_imp_dp,eps_imp_rad,                                       &
         max_iter_mom,max_iter_alpha,max_iter_scalar,max_iter_ke,max_iter_dp,max_iter_rad,                                 &
         imp_boron_trans,max_iter_boron,eps_imp_boron
      NAMELIST /second_order/ mass_conv_2nd,mom_conv_2nd,eng_conv_2nd,qula_conv_2nd,ncg_conv_2nd,boron_conv_2nd,turb_conv_2nd, &
         limiter_mass,limiter_mom,limiter_eng,limiter_xn,limiter_qn,limiter_boron,limiter_turb,                                &
         faclim,vkt      
      NAMELIST /turb/ iturb,turb_phase,lowreynolds,buoyancy_turb,tlengs,turbubble,turboil,iVisRatio,vis_ratio, &
         i_macroturb_source,i_turb_zero,i_wall_fric
      NAMELIST /int_ht/ mtopol,mHTC,&
         H_il_coeff_a,H_il_coeff_b,H_il_min,&
         H_ig_coeff_a,H_ig_coeff_b,H_ig_min,&
         H_fg_coeff_a,H_fg_coeff_b,H_fg_min,&
         dtl,dtg,&
         H_il_min_user,H_il_min_ag99,&
         relax_hik,l_th_equil,l_min_hik,alphag_min,alphal_min,&
         i_weight
      NAMELIST /int_mom/ mdrag,Ntdf,Nlift,Nwlf,&
         Drag_coeff_a,Drag_coeff_b,&
         vFgl_1_min,vFgl_2_min,vFgl_3_min,&
         Cd_min_user,Cd_min_ag99,relax_cd,&
         i_turb_disp,i_wall_lub         
      NAMELIST /iat_heatpart/ iat,dbubble_init,iheatpart,kfactor,wconden,rv_htmodel_forCFD,dia_rod_cfd,rad_model,abs_coeff,wall_emiss,&
         r_db_min,r_db_max,r_ddrop,i_bubble_diameter,r_dh_hibiki,i_drop_diameter
      NAMELIST /ncg/ n_ncg_sp,ncg_diff,imp_ncg,i_ncg_vis,&
         ncg_species,qn_cell0,qn_nvin,qn_npin                    
      NAMELIST /rv_models/ rv_model,rv_fw_reg,rv_ht_str,rv_ht_w,rv_fric_i,rv_ht_i,rv_fric_w,rv_choke,rv_gapcond,rv_mcp,rv_valve
      NAMELIST /rv_parameters/ inline_bundle,reflood,i_rv_int_fric,i_critical_flow,rough_wall, &
         hyd_core,l_plate,dia_rod,pit_dia,h_bundle,base_bundle,ia_option
      NAMELIST /rv_valvemodel/ num_valve,fzone_valve,time_valve_closed
      NAMELIST /rv_chokemodel/ fzone_throat,choke_throat_area,env_press_option,choke_pout,relax_choke,time_cflow_on
      NAMELIST /rv_mcpmodel/ num_mcp,fzone_mcp,frac_tabl,han,hvn,had,hvd,hat,hvt,har,hvr, &
         rated_pump_speed,rated_pump_hd,rwinit,relax_flow, &
         num_mcp_transient,mcp_transient_start,speed_pump, &
         time_mcp_ramping,vflow_direct,mcp_vflow,time_mcp_on,time_mcp_off
      NAMELIST /rv_pdropmodel/ npid,num_dp_region,time_pid_on,time_pid_off      
      NAMELIST/subchannel/l_subchannel,l_mixing_vane,l_spacer_grid,l_2p_multiplier_evvd,    &
                           i_subchannel_fric,i_subchannel_mixing,i_subchannel_fric_axial,   &
                           i_subchannel_fric_cross,i_2p_multiplier,udfi_subchannel_flowdir, &
                           kloss_grid,kloss_cross,ftm,ka,beta
      NAMELIST /misc_option/ mboron,mdiffoff,cplmars,cplmaster,topolsurface,mdiffscheme,iso_thermal,ifrink, & 
!.....add      
      vfporous,i_droplet,i_fs_temp_intpol
      NAMELIST /time_control/ ctrl_opt,dt,cfl_ratio,t_end,toutstep,treststep,cfl_ratio_max,smac,dt_opt,dt_max, &
         nctrl,t_end_ctrl,toutstep_ctrl,treststep_ctrl,cfl_ratio_max_ctrl,smac_ctrl,dt_opt_ctrl,dt_max_ctrl
      NAMELIST /post/ save_option,noutput,iprn,restart,itim,&
         tplot_dt,tplot_num,tplot_prop,viewField
      NAMELIST /ipara/ nb_max,nin_max,nb_sym,nb_mars
!
      CALL default_flow
!
      read(unit_somaflow, nml=problem_description)
      read(unit_somaflow, nml=n_scheme)
      read(unit_somaflow, nml=turb)
      read(unit_somaflow, nml=int_ht)
      read(unit_somaflow, nml=int_mom)
      read(unit_somaflow, nml=iat_heatpart)
      read(unit_somaflow, nml=rv_models)
      read(unit_somaflow, nml=misc_option)
      read(unit_somaflow, nml=post)
      read(unit_somaflow, nml=ipara)
!
!.....Subchannel parameters for Subchannel-scale calculations
!
      unit_chn=901
      OPEN(unit_chn,file='subchannel_parameters.in',status='old',iostat=err)
      IF(err.eq.0)then
         IF(myrank.eq.0)WRITE(*       ,"(11x,a)")'Reading subchannel_parameters.in...'
         IF(myrank.eq.0)WRITE(unit_log,"(11x,a)")'Reading subchannel_parameters.in...'
         READ(unit_chn, nml=subchannel)
      ENDIF
      CLOSE(unit_chn)
!     
      
      num_wall_group=0
      select_wall_group=1            ! 1: surface vector, 2: coordinate, 3: nbcon     
!     
      IF(nout.eq.-1) THEN
         ALLOCATE(grav(ndim))
         CALL user_def_inp(0) 
         rewind(unit_somaflow)
         RETURN
      ENDIF  
!     
      ALLOCATE(v0(num_max_zone,ndim),p0(num_max_zone),tl0(num_max_zone))
      ALLOCATE(tg0(num_max_zone),quala0(num_max_zone),a0_g(num_max_zone), cboron0(num_max_zone))
!     
      rewind(unit_somaflow)
      read(unit_somaflow, nml=initial_condition)
      DO i=1,num_fzone_count
         v0(i,:)   =v0(1,:)
         p0(i)     =p0(1)
         tl0(i)    =tl0(1)
         tg0(i)    =tg0(1)
         a0_g(i)   =a0_g(1)
         quala0(i) =quala0(1)
         cboron0(i)=cboron0(1)
         q0_gas(i) =q0_gas(1)
         q0_liq(i) =q0_liq(1)
      ENDDO
      DO i=1,num_szone_count
         tsol0(i)=tsol0(1)
         q0_ice_solid(i)=q0_ice_solid(1)
      ENDDO      
      
      rewind(unit_somaflow)
      read(unit_somaflow, nml=initial_condition)
!      
      sumGravity=dabs(grav(1))+dabs(grav(2))
      IF(ndim.gt.2)sumGravity=sumGravity+dabs(grav(ndim))
      lgravity = .false.
      IF(sumGravity .ge. 1.d-3) lgravity = .true.
      IF(num_fzone_count.ne.num_fzone) THEN 
         PRINT*,'ERROR! Initial fluid + porous zone number in somaFlow.in is different from that in somaGrid.in!'
         PRINT*,num_fzone_count,num_fzone
         PAUSE
         STOP
      ENDIF 
!     
      IF(num_szone_count.ne.num_pzone+num_szone) THEN 
         PRINT*,'ERROR! Initial solid + porous zone number in somaFlow.in is different from that in somaGrid.in!'
         PRINT*,num_szone_count,num_pzone+num_szone
         PAUSE
         STOP
      ENDIF  
!     
      ix=1
      tsol0_temp(:)=tsol0(:)
      q0_ice_solid_temp(:)=q0_ice_solid(:)
      tsol0(:)=0.d0
      q0_ice_solid(:)=0.d0
      DO izone=num_fzone-num_pzone+1,num_fzone+num_szone
         tsol0(izone)=tsol0_temp(ix)
         q0_ice_solid(izone)=q0_ice_solid_temp(ix)
         ix=ix+1
      ENDDO
!     
      iprofbc=0
!      
!.....Parallel (0:BICGSTAB,1:CSR-MPI)
!.....Psolve (1:CG solver+diag.precond,2:CG solver+ILU precond,3:BICG solver+Diag. precond,4:BICG solver+ILU precond.)      
!
!      
      read(unit_somaflow, nml=boundary_condition)
      read(unit_somaflow, nml=n_scheme)
      read(unit_somaflow, nml=imp_scheme)
      read(unit_somaflow, nml=second_order)
      read(unit_somaflow, nml=turb)
      read(unit_somaflow, nml=int_ht)
      read(unit_somaflow, nml=int_mom)
      read(unit_somaflow, nml=iat_heatpart)
      read(unit_somaflow, nml=ncg)
      read(unit_somaflow, nml=rv_models)
      read(unit_somaflow, nml=misc_option)
      read(unit_somaflow, nml=time_control)
      read(unit_somaflow, nml=post)
      IF(rv_htmodel_forCFD.ne.0.and.dia_rod.le.0) STOP '### dia_rod shoud be larger than zero when rv_htmodel_forCFD=1.'
!
!.....Additional parameters for CUPID-RV calculations
! 
      IF(rv_model.eq.1) THEN
         unit_rv=900
         OPEN(unit_rv,file='rv_parameters.in',status='old',iostat=err)
         IF(err.gt.0) STOP '### rv_parameters.in is required when rv_model=1.'
         READ(unit_rv, nml=rv_parameters)
         READ(unit_rv, nml=rv_valvemodel)
         READ(unit_rv, nml=rv_chokemodel)
         READ(unit_rv, nml=rv_mcpmodel)   
         READ(unit_rv, nml=rv_pdropmodel)           
         IF(dia_rod.le.0) STOP '### dia_rod shoud be larger than zero when rv_model=1.'
      ENDIF      
!
      CALL default_flow_after_allocation
      CALL default_flow_convert_to_string
!
      IF(np.gt.1.and.parallel.eq.0)parallel=1
!     
      DO i=1,nvin
         DO ix=1,ndim
            vb_drp(i,ix)=vb_gas(i,ix)
         ENDDO
         IF(vin_norm(i).gt.0)THEN
           vin_liq(i)=-vb_liq(i,1)
           vin_gas(i)=-vb_gas(i,1)
           vin_drp(i)=-vb_drp(i,1)
         ENDIF
      ENDDO
!     
      DO i=n_constTwall+2,3,-1
         Twall(i) = Twall(i-2)
      ENDDO      
      DO i=n_constqwall+nin_max,nin_max+1,-1
         qwall_solid(i) = qwall_solid(i-nin_max)
      ENDDO   
!     
      Cd_min_user=MAX(0.1d0,Cd_min_user)
      Cd_min_ag99=MAX(0.1d0,Cd_min_ag99)
      
      iter_grad=max(1,iter_grad)
      non_orth_iter=max(1,non_orth_iter)
      IF((grav_grad.eq.2.or.non_orth_diff.eq.2).and.(ifrink.ne.1.and.ifrink.ne.2))THEN
         IF(myrank.eq.0)THEN
            WRITE(*,"(5x,a)")'At IO, set frink method to on when, at NS,' 
            WRITE(*,"(5x,a)")'pressure interpolation of gravity weighted is on' 
            WRITE(*,"(5x,a)")'cross diffusion of frink gradient is on!' 
            WRITE(*,"(5x,a)")'Set ifrink to 1 or 2 when grav_grad==2 or non_orth_diff==2 !' 
         ENDIF
         PAUSE
         STOP
      ENDIF
!     
      IF(imp_ke_diff.eq.0)THEN
         imp_ke_diff=0
         imp_ke_conv=0
      ELSEIF(imp_ke_diff.eq.2)THEN
         imp_ke_diff=1
         imp_ke_conv=0
      ELSEIF(imp_ke_diff.eq.3)THEN
         imp_ke_diff=0
         imp_ke_conv=1 
      ELSE
         imp_ke_diff=1
         imp_ke_conv=1         
      ENDIF       
!   
      IF(n_ncg_sp.eq.1)THEN
         qn_cell0(1) =1.d0
         qn_nvin(:,1)=1.d0
         qn_npin(:,1)=1.d0
      ENDIF
!     
      IF(mdiffoff.eq.1) imp_mom_diff=0
!
!.....treat restart
!     
      lrestart_changed_nbcon=.false.
      lrestart_overwrite=.false.
      IF(restart.ge.100)THEN 
         lrestart_changed_nbcon=.true.
         restart=mod(restart,100)
      ENDIF
      IF(restart.ge.10)THEN
         lrestart_overwrite=.true.
         restart=mod(restart,10)
      ENDIF 
!      
      IF(tplot_num.gt.100) THEN
         print*,'>> Warning: Positions must be less than 10'
         print*,'>> Warning: Change tplot_num in read_flow.f90'
         STOP
      ENDIF  
!     
      IF(tplot_num.gt.0) THEN
!
!........Open 'timeplot.dat' in CPU0
!
         IF(myrank.eq.0) THEN
            DO i=1,tplot_num
               WRITE(fn,*)i               
               fn=ADJUSTL(fn)
               f_tplotv='tplot_vector'//TRIM(fn)//'.dat'              
               f_tplots='tplot_scalar'//TRIM(fn)//'.dat'
               unit_tplotv(i)=createUnit("tplotv")
               unit_tplots(i)=createUnit("tplots")
               unit_tplotv(i)=550+i
               unit_tplots(i)=660+i
               IF(restart.eq.0)THEN
                  OPEN(unit_tplotv(i),file=f_tplotv)
                  OPEN(unit_tplots(i),file=f_tplots)
                  ! OPEN(550+i,file=f_tplotv)
                  ! OPEN(660+i,file=f_tplots)
               ELSE
                  OPEN(unit_tplotv(i),file=f_tplotv,position='append')
                  OPEN(unit_tplots(i),file=f_tplots,position='append')
                  ! OPEN(550+i,file=f_tplotv,position='append')
                  ! OPEN(660+i,file=f_tplots,position='append')
               ENDIF   
            ENDDO   
         ENDIF   
!
!........Read global location (or cell number) from somaFlow.in
!               
         ALLOCATE(tplot_cell(tplot_num)) 
         IF(tplot_prop.eq.0) THEN
            IF(ndim.eq.2) ALLOCATE(tplot_x(tplot_num),tplot_y(tplot_num)) 
            IF(ndim.eq.3) ALLOCATE(tplot_x(tplot_num),tplot_y(tplot_num),tplot_z(tplot_num)) 
            DO i=1,tplot_num
               IF(ndim.eq.2) READ(unit_somaflow,*) tplot_x(i),tplot_y(i)
               IF(ndim.eq.3) READ(unit_somaflow,*) tplot_x(i),tplot_y(i),tplot_z(i) 
            ENDDO    
         ELSEIF(tplot_prop.eq.1)then
            DO i=1,tplot_num
               READ(unit_somaflow,*) tplot_cell(i)
               IF(tplot_cell(i).gt.ncell_fluid_all) THEN
                  IF(myrank.eq.0) THEN
                     print*,'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
                     print*,'>>>>> For tplot_prop=1, cell index should be less than ncell_fluid. <<<<<'
                     print*,'>>>>>                 Change the [cell index].                      <<<<<'
                     print*,'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
                  ENDIF   
                  print*,''
                  print*,'>>>>> Myrank=',myrank,', Max. cell index=',ncell_fluid,'  <<<<<'
                  PAUSE
                  STOP
               ENDIF
            ENDDO           
         ENDIF
!
!........Search the nearest point (global viewpoint) when tplot_prop=0
!
         IF(tplot_prop.eq.0) THEN      
            ALLOCATE(xloc_all(ncell_fluid_all,ndim),distance(tplot_num),distance_min(tplot_num)) 
            distance_min=10.d10
            DO ix=1,ndim                        
               CALL allgatherv_r(xloc(1,ix),xloc_all(1,ix),ncell_fluid,ncell_fluid_all,0) 
            ENDDO         
            DO i=1,ncell_fluid_all
            DO ii=1,tplot_num
               distance(ii)=(tplot_x(ii)-xloc_all(i,1))*(tplot_x(ii)-xloc_all(i,1))+(tplot_y(ii)-xloc_all(i,2))*(tplot_y(ii)-xloc_all(i,2))
               IF(ndim.eq.3) distance(ii)=distance(ii)+(tplot_z(ii)-xloc_all(i,3))*(tplot_z(ii)-xloc_all(i,3))
               distance(ii)=dsqrt(distance(ii))
               IF(distance(ii).le.distance_min(ii)) THEN 
                  tplot_cell(ii)=i
                  distance_min(ii)=distance(ii)
               ENDIF   
            ENDDO
            ENDDO
         ENDIF   
!
         ALLOCATE(tplot_cell_loc(tplot_num),tplot_cell_rank(tplot_num)) 
         tplot_cell_loc(:)=0
         tplot_cell_rank(:)=0
         DO i=1,tplot_num
            tplot_cell_loc(i)=i3perm(tplot_cell(i))
            if(tplot_cell_loc(i).lt.1.or.tplot_cell_loc(i).gt.ncell_fluid)tplot_cell_loc(i)=0
            IF(tplot_cell_loc(i).gt.0)THEN
               tplot_cell_rank(i)=myrank
               WRITE(*,"(11x,a,1i8,1i6,1i3)")'tplot_cell,tplot_cell_loc,rank=',tplot_cell(i),tplot_cell_loc(i),tplot_cell_rank(i)
            ENDIF   
         ENDDO  
         CALL allreducei_i(tplot_cell_rank,tplot_num)
!     
      ENDIF  
!
!.....treat wconden
!
      lwconden_alphal0=.false.
      IF(IABS(wconden).ge.10)THEN
          lwconden_alphal0=.true.
          IF(wconden.ge.10)THEN
             wconden=wconden-10
          ELSEIF(wconden.le.10)THEN
             wconden=wconden+10
          ENDIF   
      ENDIF 
!
!.....treat rv_model
!
      IF(rv_model.lt.1)THEN
         rv_fric_i=0
         rv_ht_i=0
         rv_fric_w=0
      ENDIF
!
      lfric_swap=.false.
      lhtc_swap=.false.
      lfricw_swap=.false.
      IF(rv_fric_i.lt.0)lfric_swap=.true.
      IF(rv_ht_i.lt.0)lhtc_swap=.true.
      IF(rv_fric_w.lt.0)lfricw_swap=.true.
!      
      rv_fric_i=IABS(rv_fric_i)
      rv_ht_i=IABS(rv_ht_i)      
      rv_fric_w=IABS(rv_fric_w)
!      
      IF(myrank.eq.0)THEN
         IF(lfric_swap)WRITE(*,"(11x,a)")'Int. drag model for open and porous media is used.'
         IF(lhtc_swap)WRITE(*,"(11x,a)")'Int. HTC model for open and porous media is used .'
         IF(lfricw_swap)WRITE(*,"(11x,a)")'Wall friction model for open and porous media is used.'
         IF(lfric_swap)WRITE(unit_log,"(11x,a)")'Int. drag model for open and porous media is used.'
         IF(lhtc_swap)WRITE(unit_log,"(11x,a)")'Int. HTC model for open and porous media is used .'
         IF(lfricw_swap)WRITE(unit_log,"(11x,a)")'Wall friction model for open and porous media is used.'
      ENDIF
!
!----------------------------------------------------------------------
!-----------------------End of Read somaFlow.in------------------------
!----------------------------------------------------------------------
!
!.....Allocate global variables
!
      IF(nout.eq.0)THEN !next-20170512
         IF(parallel.ge.1) THEN
            include '../00_Module/c_Solver/allocate_cell_par.h'
            include '../00_Module/c_Solver/allocate_solid_par.h'
!           IF(ncell_ps.gt.0) ALLOCATE(solid(ncell_ps))
         ELSE
            include '../00_Module/c_Solver/allocate_cell_ser.h'
            include '../00_Module/c_Solver/allocate_solid_ser.h'
!           IF(ncell_cond.gt.0) ALLOCATE(solid(ncell_cond))
         ENDIF   
         include '../00_Module/c_Solver/allocate_face.h'
      ENDIF
      
!     ALLOCATE(face(1:ncell_fp))
!
!.....Initialize thermal property for initial run or restart run
!     restart=0; initial run
!     restart=1; restart run from user inputted itim
!     restart=2; restart run from the last itim
!      
      IF(restart.eq.0) THEN
         REWIND(unit_somaflow)
         READ(unit_somaflow,NML = post)
         restart=mod(restart,10)
         CALL user_def_inp(0) 
         CALL initialize_variables(v0,p0,tl0,tg0,a0_g,quala0,tsol0,cboron0)           
         IF(nout.eq.0)CALL read_tpfh2o          
         CALL ncg_cell
         CALL property_calc(0)
         CALL init_steamtable        
      ELSE
         REWIND(unit_somaflow)
         READ(unit_somaflow,NML = post)
         restart=mod(restart,10)
         CALL read_tpfh2o
         CALL ncg_cell
         CALL restart_read(nout)
         CALL init_steamtable          
      ENDIF
!
      IF(restart.eq.0) itim=0
!
!.....Transform Mass flow rate into Velocity
!
      DO i=1,nvin
         IF(vin_mfr(i).gt.0) CALL inlet_velocity_from_mass_flow_rate(i)
      ENDDO
!
!.....Miscellaneous coefficients
!
      DO i=1,nvin
         cvm=c_vm(alphab_gas(i))
         arhob_liq=alphab_liq(i)*rhob_liq(i)
         arhob_drp=0.d0
         arhob_gas=alphab_gas(i)*rhob_gas(i)
         rhob_m=arhob_liq+arhob_gas
         rhob_vm=cvm*(arhob_liq+arhob_gas)
         denom=rhob_gas(i)*rhob_liq(i)+rhob_m*rhob_vm
         cb_pl(i)=(rhob_gas(i)+rhob_vm)/denom
         cb_pg(i)=(rhob_liq(i)+rhob_vm)/denom
         cb_pd(i)=cb_pg(i)
      ENDDO
!
!.....Calculate coefficients for turbulence model
!
      DO i=1,nvin
         vb_liq_size=0.d0
         vb_gas_size=0.d0
         DO ix=1,ndim
            vb_liq_size=vb_liq_size+(vb_liq(i,ix)**2)
            vb_gas_size=vb_gas_size+(vb_gas(i,ix)**2)
         ENDDO
!         IF(vin_norm(i).gt.0)THEN
!            vb_liq_size=vin_liq(i)**2
!            vb_gas_size=vin_gas(i)**2
!         ENDIF
         turb_keb(i) =1.5d0 * 0.05d0**2.d0*vb_liq_size 
         turb_kegb(i)=1.5d0 * 0.05d0**2.d0*vb_gas_size 
!
         IF(iVisRatio.eq.0)THEN       !Turbulence Intensity and Length Scale
            IF(iturb.eq.1)THEN        
               turb_dpb(i) =cmu**(-0.25)*turb_keb(i)**0.5d0/(0.07d0 * 0.5d0)     !SST k-w
               turb_dpgb(i)=cmu**(-0.25)*turb_kegb(i)**0.5d0/(0.07d0 * 0.5d0)
!               turb_dpb(i) =cmu**(-0.25)*turb_keb(i)**0.5d0/(0.42d0*d_bfc(i))
!               turb_dpgb(i)=cmu**(-0.25)*turb_kegb(i)**0.5d0/(0.42d0*d_bfc(i))            
            ELSEIF(iturb.ge.2.and.iturb.le.4)THEN
               turb_dpb(i) =0.09d0**0.75*turb_keb(i)**1.5d0/(0.07d0 *0.5d0)      ! k-e
               turb_dpgb(i)=0.09d0**0.75*turb_kegb(i)**1.5d0/(0.07d0 *0.5d0) 
            ENDIF
         ELSEIF(iVisRatio.eq.1)THEN   !Turbulene Intensity and Viscosity Ratio
            IF(iturb.eq.1)THEN
               turb_dpgb(i)=rhob_gas(i)*turb_kegb(i)/(lvisb_gas(i)*vis_ratio)           !SST k-omega
               turb_dpb(i) =rhob_liq(i)*turb_keb(i)/ (lvisb_liq(i)*vis_ratio)
            ELSEIF(iturb.ge.2.and.iturb.le.4)THEN
               turb_dpgb(i)=rhob_gas(i)*cmu*turb_kegb(i)**2/(lvisb_gas(i)*vis_ratio)    ! k-epsilon
               turb_dpb(i) =rhob_liq(i)*cmu*turb_keb(i)**2/ (lvisb_liq(i)*vis_ratio)
            ENDIF
         ENDIF         
      ENDDO
!         
      CALL initialize_topology_criteria
!     
      DEALLOCATE(v0,p0,tl0)
      DEALLOCATE(tg0,quala0,a0_g)
!
      IF(cplmaster.gt.0) CALL read_cupid_master
!
!.....Save boron concentration
!       
      DO i=1,nvin
         cboronb_liq(i)=cboronb(i)
      ENDDO
!      
      END SUBROUTINE read_flow
!
!--------------------------------------------------------------------------------
!
      SUBROUTINE inlet_velocity_from_mass_flow_rate(ii)
!
!     Transform Mass flow rate into Velocity
!
      USE Zcore           , ONLY: np
      USE Znum_cell       , ONLY: istart_nf,istart_nbcon_nf
      USE Zvec_index      , ONLY: nbcon_nf
      USE Zb_condition    , ONLY: vb_liq,vb_gas,rhob_liq,rhob_gas
      USE Zbc_index       , ONLY: vin_mfr
      USE Zvec_geo        , ONLY: saa_nf
!
      IMPLICIT NONE
!
!.....Input
      INTEGER ii
!.....Local variables
      INTEGER i,k
      INTEGER nf_number,istart,isize,i1,istart2,i2
      REAL(8) area_in
!
      area_in=0.d0
      nf_number=2
      istart=istart_nf(1,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      isize =istart_nf(2,nf_number)
      DO i=1,isize
         i1=istart+i
         i2=istart2+i
         k=nbcon_nf(i2)
         IF(k.eq.ii) area_in=area_in+saa_nf(i1)
      ENDDO
      IF(np.gt.1) CALL allreducei_r1(area_in)
!
      IF(vin_mfr(ii).eq.1)THEN        !g/s  --> m/s
         vb_gas(ii,:)=vb_gas(ii,:)/(rhob_gas(ii)*area_in*1.d3)
         vb_liq(ii,:)=vb_liq(ii,:)/(rhob_liq(ii)*area_in*1.d3)
      ELSEIF(vin_mfr(ii).eq.2)THEN    !kg/s --> m/s
         vb_gas(ii,:)=vb_gas(ii,:)/(rhob_gas(ii)*area_in)
         vb_liq(ii,:)=vb_liq(ii,:)/(rhob_liq(ii)*area_in)
      ENDIF
!
      END SUBROUTINE inlet_velocity_from_mass_flow_rate
!------------------------------------------------------------------------------------
      SUBROUTINE default_flow
!
      USE Zcore            ,ONLY: myrank,np
      USE Zparam           , ONLY: nin_max,nb_max,ndim,cmu,nb_max,nin_max,nb_sym,nb_mars
      USE STM_TBL_cupid    , ONLY: nfluid
      USE Zturbzeq         , ONLY: tleng,tlengs
      USE Zgrad_ls_c3d     , ONLY: lsindex
      USE Zimplicit        , ONLY: imp_boron_trans,max_iter_boron,eps_imp_boron
      USE Zmodel           ,ONLY: i_droplet,i_weight,s_wall_fric,i_fs_temp_intpol
      USE Zndforce         ,ONLY: s_turb_disp,s_wall_lub  
      USE Zncg             ,ONLY: i_ncg_vis                      
      USE Ziat             ,ONLY: r_db_min,r_db_max,r_ddrop,s_bubble_diameter,r_dh_hibiki,s_drop_diameter 
      USE Zbicg            ,ONLY: lev_type,levmpi_type,lev,levmpi,        &
                                  lev_type_c,levmpi_type_c,lev_c,levmpi_c
      USE ZturbZeq         ,ONLY: s_turb_zero
      USE Zporous          ,ONLY: l_subchannel,l_mixing_vane,l_spacer_grid,l_2p_multiplier_evvd,    &
                                   s_subchannel_fric,s_subchannel_mixing,s_subchannel_fric_axial,   &
                                   s_subchannel_fric_cross,s_2p_multiplier,udfi_subchannel_flowdir, &
                                   vfporous
      USE Zturb            ,ONLY: s_macroturb_source
      USE Zrv_wall_fric    , ONLY: rough_wall      
      USE Zrv_int_friction ,ONLY: s_rv_int_fric
      USE Zrv_choke        ,ONLY: choke,s_critical_flow,fzone_throat, &
                                  choke_pout,relax_choke,choke_throat_area,env_press_option, &
                                  time_cflow_on,choke,choke_update
      USE Zpdrop          , ONLY: time_pid_on,time_pid_off,num_dp_region  
      USE Zbc_index        ,ONLY: l_horizontal_outlet_init
      USE Zwall_HTC        ,ONLY: hyd_core,l_plate,dia_rod,pit_dia,h_bundle,base_bundle,inline_bundle,reflood
      USE Zconst1          ,ONLY: i_macroturb_source,i_turb_zero,i_wall_fric,i_turb_disp,i_wall_lub,                     &
                                  i_bubble_diameter,i_drop_diameter,i_rv_int_fric,i_critical_flow,                       &
                                  i_subchannel_fric,i_subchannel_fric_axial,i_subchannel_fric_cross,i_subchannel_mixing, &
                                  i_2p_multiplier,rv_htmodel_forCFD,                                                     &
                                  iturb
      USE Zconst2          , ONLY: stime_hup,stime_hflat,hydraulicd,gfactor, &
                                   stime_vup,stime_vflat,                    &
                                   hydraulicd_init
      USE Zdel_scalar      , ONLY: stmtbl_repeat_for_nc,max_ihtc_opt,relax_interface_dtemp,max_ihtc_opt_coeff
      USE Zrv_model        , ONLY: ia_option       
      USE Zmcp             
      USE Zvalve
!
      IMPLICIT NONE      
!
!      CHARACTER(30) s_macroturb_source,s_turb_zero,s_wall_fric,s_turb_disp,s_wall_lub,                    &
!                     s_bubble_diameter,s_drop_diameter,s_rv_int_fric,s_critical_flow,s_subchannel_fric,   &
!                     s_subchannel_fric_axial, s_subchannel_fric_cross,s_subchannel_mixing,s_2p_multiplier
!
!.....Initialzation_1       
      choke=.FALSE.              !initially, no crtical flow
      l_horizontal_outlet_init=.FALSE. !input-find pressure boundary cells only when l_horizontal_outlet=.true.  
      
!      lsindex(:)=1               !least square multiplier
!      gfactor(:)=1.d0           !gravity multiplier to momentum equation         
!      tleng(:)=tlengs            !length scale for zero equation turbulence model 
!      if(hydraulicd_init.eq.0)THEN
!         hydraulicd(:)=hyd_core
!      endif
!      
!     NAMELIST /problem_description/
      nfluid=1
!.....NAMELIST/ipara/  
      nb_max=10
      nin_max=4                   !4=general, 5=UPTF   
      nb_sym=101
      nb_mars=201
!.....NAMELIST /n_scheme/ 
      lev_type=0
      levmpi_type=0
      lev=0
      levmpi=0
      lev_type_c=0
      levmpi_type_c=0
      lev_c=0
      levmpi_c=0
!
!..... NAMELIST/imp_scheme/ 
      imp_boron_trans=0      
      eps_imp_boron=0.1D-07 
      max_iter_boron=100
      stmtbl_repeat_for_nc=0
      max_ihtc_opt=0
      relax_interface_dtemp=0
      max_ihtc_opt_coeff=0.5d0
!      
!.....NAMELIST/intial_condition/      
      stime_vup=0.d0         
      stime_vflat=0.d0                            
      stime_hup=0.d0            
      stime_hflat=0.d0
!      
!.....NAMELIST/turb/          
      s_macroturb_source='none' !'none','chandesris','nakayama' !<--turbulence macroscopic
      s_turb_zero='constant'    !'constant','telluride','cfx','vanDriest','baldwin_Lomax','noto','user' 
      s_wall_fric='none'        !'none','kakac','darcy','udf'  
!      
!.....NAMELIST/int_ht/      
       i_weight=0
!      
!.....NAMELIST/int_mom/  
      s_turb_disp='cfx'         !'cfx','gosman','burns','udf'
      s_wall_lub='antal&cfx'    !'antal&cfx','antal&antal','tomiyama','udf' 
!      
!.....NAMELIST/iat_heatpart/      
      rv_htmodel_forCFD=0
      s_bubble_diameter='yoneda' !'yoneda','trac','modified_yoneda','viswanathan','hibiki','r_dh_hibiki'    
      r_dh_hibiki=25.52d-3
      r_db_min=0.0001d0               
      r_db_max=0.05d0                 
      s_drop_diameter='constant' 
      r_ddrop=0.001d0
!      
!.....NAMELIST/ncg/  
      i_ncg_vis=0           !slow and exact viscosity
!      
!.....NAMELIST/misc_option/
      vfporous=1
      i_droplet=0           !droplet field
      i_fs_temp_intpol=0  !temperature interpolation of fluid and solid interface
!      
!.....NAMELIST/rv_parameters/ 
      inline_bundle=0
      reflood=0
      s_rv_int_fric='drift_flux'         !drift_flux,constant,drag_coeff       
      s_critical_flow='Henry-Fauske-mod' !Henry-Fauske-mod,Henry-Fauske,Murdock-Bauman  
      rough_wall=2.d-6                  !wall roughness for single phase flow
      hyd_core=0.009731d0    ! Hydraulic diameter
      l_plate=3.6576d0       ! Length of rectangle to calculate HTC for natural convection by Chulchill-Chu 
      dia_rod=0.013046d0     ! Diameter of fuel rod
      pit_dia=1.3263d0       ! Ratio of pitch and Diameter of fuel rod                         
      h_bundle=3.6576d0      ! Axial length of bundle for Reflood calculation                              
      base_bundle=0.d0       ! height of the botton of bundle; base heigh   
      ia_option=1            ! RV intetfacial area model (1=org, 2=mod)      
!      
!.....NAMELIST/rv_valvemodel/       
      IF(.not.ALLOCATED(fzone_valve)) ALLOCATE(fzone_valve(5,2))               !max number of valves is 5.
      IF(.not.ALLOCATED(time_valve_closed)) ALLOCATE(time_valve_closed(5,2))     
      IF(.not.ALLOCATED(num_valveface)) ALLOCATE(num_valveface(5)) 
      IF(.not.ALLOCATED(valve_closed)) ALLOCATE(valve_closed(5))
      fzone_valve(:,:)=0
      time_valve_closed(:,:)=0.d0
      num_valveface(:)=0
      valve_closed(:)=0
!
!.....NAMELIST/rv_chokemodel/       
      IF(.not.ALLOCATED(fzone_throat)) ALLOCATE(fzone_throat(1,2))   
      fzone_throat(1,:)=0    !number of fluid zone1 (inside the broken face)  
      choke_throat_area=0.1  !choked area (m2)
      env_press_option=0           
      choke_pout=1.d5        !environment pressure (out of the broken face)
      relax_choke=0.99       !relaxation of choking velocity  
      time_cflow_on=0.d0
      choke=.FALSE.
      choke_update=0
!      
!.....NAMELIST/rv_mcpmodel/ 
      init_mcp=1             !default steady state condition (not an input parameter)
      IF(.not.ALLOCATED(fzone_mcp)) ALLOCATE(fzone_mcp(10,2))  !max num_mcp=10 is assumed.
      fzone_mcp(:,:)=1       ! number of fluid zone (1:num_mcp, 1:2)
      IF(.not.ALLOCATED(mcp_on)) ALLOCATE(mcp_on(10))  !max num_mcp=10 is assumed.
      mcp_on(:)=0         ! mcp is in operation (True)   
      IF(.not.ALLOCATED(num_mcpface)) ALLOCATE(num_mcpface(20))  !max num_mcp=10 is assumed.
      num_mcpface(:)=0
! 
      frac_tabl(:)=0.d0
      han(:)=0.d0
      hvn(:)=0.d0
      had(:)=0.d0
      hvd(:)=0.d0      
      hat(:)=0.d0
      hvt(:)=0.d0      
      har(:)=0.d0
      hvr(:)=0.d0   
!      
      IF(.not.ALLOCATED(rated_pump_speed)) ALLOCATE(rated_pump_speed(10)) !max num_mcp=10 is assumed.
      rated_pump_speed(:)=0.d0
      IF(.not.ALLOCATED(rated_pump_hd)) ALLOCATE(rated_pump_hd(10))  !max num_mcp=10 is assumed.
      rated_pump_hd(:)=0.d0
      IF(.not.ALLOCATED(rwinit)) ALLOCATE(rwinit(10)) !max num_mcp=10 is assumed. 
      rwinit(:)=0.d0
      relax_flow=0.003d0
      IF(.not.ALLOCATED(mcp_transient_start)) ALLOCATE(mcp_transient_start(20)) !max num_mcp_transient=20 is assummed.
      mcp_transient_start(:)=0.d0
      IF(.not.ALLOCATED(speed_pump)) ALLOCATE(speed_pump(10,20)) !max num_mcp=10, max num_mcp_transient=20 is assummed.
      speed_pump=0.d0    
!
!.....NAMELIST/rv_pdropmodel/
      time_pid_on=0.d0
      time_pid_off=0.d0
      IF(.not.ALLOCATED(num_dp_region)) ALLOCATE(num_dp_region(20,4)) !num_dp_region(max npid=20,4)
      num_dp_region(:,:)=0       
!      
!.....NAMELIST/subchannel/     
      l_subchannel=.false.
      l_spacer_grid=.FALSE.            
      l_mixing_vane=.false. 
      s_subchannel_fric='none'        !'none','aniso_fric_sem','aniso_fric_imp','aniso_fric_exp'
      s_subchannel_fric_axial='none'  !'none','matra','chandesris','takeda'
      s_subchannel_fric_cross='none'  !'none','simple_formloss'
      s_subchannel_mixing='none'      !'none','em','evvd' 
      s_2p_multiplier='none'          !'none','default','armand'
      l_2p_multiplier_evvd=.false. 
      udfi_subchannel_flowdir=ndim    !(1,2,3)
!
      RETURN
      ENTRY default_flow_after_allocation
!
!.....Initialzation_2      
!      MG_solver = .false.        !true when parallel==2
!      choke=.FALSE.              !initially, no crtical flow
      lsindex(:)=1               !least square multiplier
      gfactor(:)=1.d0           !gravity multiplier to momentum equation         
      IF(iturb.eq.0) tleng(:)=tlengs            !length scale for zero equation turbulence model
      IF(np.gt.1)CALL broadcast_i1(hydraulicd_init)
      IF(hydraulicd_init.eq.0)THEN
         hydraulicd(:)=hyd_core
      ELSE   
         IF(myrank.eq.0)WRITE(*,"(11x,a)")'hydraulicd(:) was already initialized in udfn_porous_user!'
      ENDIF   
!      
      RETURN
      ENTRY default_flow_convert_to_string      
      
      SELECT CASE(i_macroturb_source)
      CASE(0)
         s_macroturb_source='none' !'none','chandesris','nakayama' !<--turbulence macroscopic
      CASE(1)
         s_macroturb_source='chandesris'
      CASE(2)
         s_macroturb_source='nakayama'
      ENDSELECT
!      
      SELECT CASE(i_turb_zero)
      CASE(0)
         s_turb_zero='none'    !'none','constant','telluride','cfx','vanDriest','baldwin_Lomax','noto','user' 
      CASE(1)
         s_turb_zero='constant'
      CASE(2)
         s_turb_zero='telluride'
      CASE(3)
         s_turb_zero='cfx'
      CASE(4)
         s_turb_zero='vanDriest'
      CASE(5)
         s_turb_zero='baldwin_Lomax'
      CASE(6)
         s_turb_zero='noto'
      CASE(7)
         s_turb_zero='user'
      ENDSELECT
!      
      SELECT CASE(i_wall_fric)
      CASE(0)
         s_wall_fric='none'        !'none','kakac','darcy','udf'  
      CASE(1)
         s_wall_fric='kakac'
      CASE(2)
         s_wall_fric='darcy'
      CASE(3)
         s_wall_fric='udf'
      ENDSELECT
!      
      SELECT CASE(i_turb_disp)
      CASE(0)
         s_turb_disp='none'         !'none','cfx','gosman','burns','udf'
      CASE(1)
         s_turb_disp='cfx'
      CASE(2)
         s_turb_disp='gosman'
      CASE(3)
         s_turb_disp='burns'
      CASE(4)
         s_turb_disp='udf'
      ENDSELECT
!      
      SELECT CASE(i_wall_lub)
      CASE(0)
         s_wall_lub='none'    !'none','antal&cfx','antal&antal','tomiyama','udf' 
      CASE(1)
         s_wall_lub='antal&cfx'
      CASE(2)
         s_wall_lub='antal&antal'
      CASE(3)
         s_wall_lub='tomiyama'
      CASE(4) 
         s_wall_lub='udf'
      ENDSELECT
!      
      SELECT CASE(i_bubble_diameter)
      CASE(0)
         s_bubble_diameter='none' !'none','yoneda','trac','modified_yoneda','viswanathan','hibiki','r_dh_hibiki'    
      CASE(1)
         s_bubble_diameter='yoneda'  
      CASE(2)
         s_bubble_diameter='trac'
      CASE(3)
         s_bubble_diameter='modified_yoneda'    
      CASE(4)
         s_bubble_diameter='viswanathan'
      CASE(5)
         s_bubble_diameter='hibiki'
      CASE(6)
         s_bubble_diameter='yun'
      ENDSELECT
!      
      SELECT CASE(i_drop_diameter)
      CASE(0)
         s_drop_diameter='none' 
      CASE(1)
         s_drop_diameter='constant' 
      ENDSELECT
!      
      SELECT CASE(i_rv_int_fric)
      CASE(0)
         s_rv_int_fric='none'         !'none',drift_flux,constant,drag_coeff       
      CASE(1)
         s_rv_int_fric='drift_flux'  
      CASE(2)
         s_rv_int_fric='constant'
      CASE(3)
         s_rv_int_fric='drag_coeff'
      ENDSELECT
!      
      SELECT CASE(i_critical_flow)
      CASE(0)
         s_critical_flow='none' !'none',Henry-Fauske-mod,Henry-Fauske,Murdock-Bauman  
      CASE(1)
         s_critical_flow='Henry-Fauske-mod'
      CASE(2)
         s_critical_flow='Henry-Fauske'
      CASE(3)
         s_critical_flow='Murdock-Bauman'
      ENDSELECT
!      
      SELECT CASE(i_subchannel_fric)
      CASE(0)
         s_subchannel_fric='none'        !'none','aniso_fric_sem','aniso_fric_imp','aniso_fric_exp'
      CASE(1)
         s_subchannel_fric='aniso_fric_sem'
      CASE(2)
         s_subchannel_fric='aniso_fric_imp'
      CASE(3)
         s_subchannel_fric='aniso_fric_exp'
      ENDSELECT
!      
      SELECT CASE(i_subchannel_fric_axial)
      CASE(0)
         s_subchannel_fric_axial='none'  !'none','matra','chandesris','takeda'
      CASE(1)
         s_subchannel_fric_axial='matra'
      CASE(2)
         s_subchannel_fric_axial='chandesris'
      CASE(3)
         s_subchannel_fric_axial='takeda'
      CASE(4)
         s_subchannel_fric_axial='ctf'
      ENDSELECT
!      
      SELECT CASE(i_subchannel_fric_cross)
      CASE(0)
         s_subchannel_fric_cross='none'  !'none','simple_formloss'
      CASE(1)
         s_subchannel_fric_cross='simple_formloss'
      ENDSELECT
!      
      SELECT CASE(i_subchannel_mixing)
      CASE(0)
         s_subchannel_mixing='none'      !'none','em','evvd' 
      CASE(1)
         s_subchannel_mixing='em'
      CASE(2)
         s_subchannel_mixing='evvd'
      ENDSELECT
!      
      SELECT CASE(i_2p_multiplier)
      CASE(0)
         s_2p_multiplier='none'          !'none','default','armand'      
      CASE(1)
         s_2p_multiplier='default'
      CASE(2)
         s_2p_multiplier='armand'
      ENDSELECT
!                                                            
      END SUBROUTINE default_flow
      

