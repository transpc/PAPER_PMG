!
      SUBROUTINE user_def_inp(iflag)
!
!     This routine defines user/problem specific input parameters & options
!
      USE VOL_DATA        , ONLY: cell     
      USE Solid_DATA      , ONLY: solid    
      USE Wall_DATA       , ONLY: face
      USE Zmpi            , ONLY: jperm     
      USE Zzone           , ONLY: ncell_fluid,nzone
      USE Zcore           , ONLY: myrank,np
      USE Znum_cell       , ONLY: istart_nf !kuhn
      USE Zvec_index      , ONLY: left_nf   !kuhn
      USE Zvec_param      , ONLY: nf_ctw,nf_fsw
      USE Znum_cell       , ONLY: i_neigh,i_neigh
      USE Zparam          , ONLY: pi,ndim,nin_max
      USE Zb_condition    , ONLY: pbnd,rhob_liq,rhob_gas,vb_liq,vb_gas,                       &
                                  p_fb,tb_gas,tb_liq,alphab_gas,alphab_liq,alphab_drp,        &
                                  alpha_gas_nd,alpha_liq_nd,alpha_drp_nd,qualab,eb_liq,eb_gas
      USE Zbc_index       , ONLY: nvin,nbcon,npb
      USE Zconst1         , ONLY: vv_prob,iturb,restart,fric_face
      USE Zconst2         , ONLY: stime_hflat,hydraulicd,gfactor,dt,dtr, &
                                  stime_vup,stime_vflat,grav,ggc
      USE Zcoord1         , ONLY: xloc
      USE Zgrad_ls_c3d    , ONLY: lsindex
!      USE Zmodel          , ONLY: vfgl_large 
      USE Zncg            , ONLY: tao,cvao_nvin,uao_nvin,dcva_nvin,ra_nvin,qn_nvin,qn_cell,n_ncg_sp
      USE Ztimecon        , ONLY: time,itim
      USE Zturbzeq        , ONLY: tleng,tlengs 
      USE Zwall_HTC       , ONLY: l_plate,dia_rod,pit_dia
      USE Zuserdefined    , ONLY: udfl_porous_user,                                          &
                                  udfl_grid_user,                                            &
                                  udfl_init_variables,                                       &
                                  udfl_calc_HTC_int_i,udfl_erg_diff,udfl_flashing_hif,       &
                                  udfl_porous_property,udfl_mat_prop,                        &
                                  udfl_mom_wall,                                             &
                                  udfl_mom_drag_i,udfl_mom_drag_i,                           &
                                  udfl_outlet_property,udfl_mom_press_source,                &
                                  udfl_mom_loss,udfl_outlet_press_user,udfl_set_qvol_porous, &
                                  udfl_update_scalar,udfl_wallHTC_porous,                    &
                                  udfl_mom_film_shear,                                       &
                                  udfl_hflux_bc_profile,                                     &
                                  udfl_tw_profile,                                           &
                                  user_rary,user_iary,                                       &
                                  udfl_model_overwrite
      USE Zmodel          , ONLY: wVertical
      USE viewData_common , ONLY: viwunit 
      USE Zsiphon         , ONLY: topen
      USE Zvec_index_solid, ONLY: dfilm_fsw,vfilm_fsw,dfilm_ctw,vfilm_ctw !kuhn
      !OPR1000
      USE Zconst1         , ONLY: cplmars
      USE Zrv_model       , ONLY: rv_model,rv_ht_str !SMR
      USE Zuserdefined    , ONLY: udfl_psbt_cfx_model
      
!
!.....will be removed
      USE Zmodel           ,ONLY: s_wall_fric
      USE Ziat             ,ONLY: s_bubble_diameter
      USE ZturbZeq         ,ONLY: s_turb_zero
      USE Zporous          ,ONLY: l_mixing_vane
      USE Zcheck_scalar   , ONLY: eps_rho,eps_p,eps_eng,eps_vol
      USE Zmcp
      USE Zvector         , ONLY: ul_o
      USE Zncg            , ONLY: tao,cvao_cell,uao_cell,dcva_cell,ra_cell,qn_cell,n_ncg_sp    
      USE KSMR
      USE NuScale        
      use zdel_scalar,only:max_ihtc_opt_coeff
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: iflag
!.....Local variables
      INTEGER :: i,j,ii,j0
      INTEGER :: err=0    
      INTEGER :: nf_number,istart,isize,jj,i1 !kuhn
      INTEGER,SAVE :: time_bc,time_bc_next,iprn
      LOGICAL,SAVE :: INITIAL=.TRUE.
!      LOGICAL,SAVE :: INITIAL_h=.TRUE.
      LOGICAL,SAVE :: initial_sgp=.TRUE.

      LOGICAL,SAVE :: initialKuhn=.TRUE.
      REAL(8) :: xx,yy,rr,rr0,pcnv
      REAL(8) :: Vlin,Pin,Tgin,Tlin,Tdin,Agin,Adin,Pout,Agout,tmp1,tmp2
!.....Local arrays
      REAL(8),SAVE :: zkuhn(13),zTwall(12)
      REAL(8) qn_cell0(n_ncg_sp)
      REAL(8) :: dpdz,vel,rcs_flow,rcs_rho,rcs_area,rcs_leng,rcs_head,resist
!
      IF(iflag.eq.0)THEN
!      
         udfl_init_variables=.FALSE. !initial condition
         udfl_grid_user=.FALSE.      !mesh 
         udfl_porous_user=.FALSE.    !porosity,permeability
!         
!--------VnV input head----------------------------------------------------------------------------
         IF(vv_prob.eq.'manometric'         )udfl_init_variables=.TRUE.
         IF(vv_prob.eq.'PAFS-POOL'          )THEN
            udfl_init_variables=.TRUE.
            udfl_porous_user=.TRUE.            
         ENDIF   
         IF(vv_prob.eq.'moving_wall'        )udfl_init_variables=.TRUE.
         IF(vv_prob.eq.'fluidic_device'     )udfl_init_variables=.TRUE.
         IF(vv_prob.eq.'fluidic_device'     )udfl_grid_user=.TRUE.
         IF(vv_prob.eq.'T_blowdown'     .or.&
            vv_prob.eq.'therm_diff'         )udfl_init_variables=.TRUE.
         IF(vv_prob.eq.'sgp_separator'      )THEN
            udfl_porous_user=.TRUE.
            udfl_init_variables=.TRUE.
         ENDIF
         IF(vv_prob.eq.'3D_boiling'         ) udfl_porous_user=.TRUE.  
         IF(vv_prob.eq.'fs_31203_3D'    .or.&
            vv_prob.eq.'fs_31302_3D'    .or.&            
            vv_prob.eq.'fs_31701_3D'    .or.&            
            vv_prob.eq.'fs_31805_3D'    .or.&            
            vv_prob.eq.'rbht1196_3d'        )THEN  
            udfl_porous_user   =.TRUE.  
            udfl_init_variables=.TRUE.
         ENDIF    
         IF(vv_prob.eq.'fs_31203'    .or.&
            vv_prob.eq.'fs_31302'    .or.&            
            vv_prob.eq.'fs_31701'    .or.&            
            vv_prob.eq.'fs_31805'    .or.&  
            vv_prob.eq.'rbht1196_1d' .or.&  
            vv_prob.eq.'ismr_rad'    .or.&  
            vv_prob.eq.'ismr_2d')THEN
            udfl_init_variables=.TRUE.            
         ENDIF   
         IF(vv_prob.eq.'UTPF-RV'             )udfl_init_variables=.TRUE.
         IF(vv_prob.eq.'kuhn_111' .or. vv_prob.eq.'Nuscale-03Pool') udfl_init_variables=.TRUE.
         IF(vv_prob.eq.'apr1400_mc_rv' .or. &
            vv_prob.eq.'opr1000_mc_rv')THEN
            udfl_porous_user=.TRUE.
            udfl_grid_user=.TRUE.
         ENDIF   
         IF(vv_prob.eq.'atlas_mc_porous')THEN
            udfl_init_variables=.TRUE.  
            udfl_porous_user=   .TRUE.
         ENDIF     
         IF(vv_prob.eq.'manosteam_mcc'      )udfl_init_variables=.TRUE. 
         IF(vv_prob.eq.'apr1400_lbloca'     .or.&
            vv_prob.eq.'opr1000_mc_rv_lbloca')THEN
            udfl_porous_user=.TRUE.
            udfl_grid_user=.TRUE.            
         ENDIF
         IF(vv_prob.eq.'siphon'              )udfl_init_variables=.TRUE.          
         IF(vv_prob.eq.'opr1000_rv' .or. &
            vv_prob.eq.'apr1400_rv' )THEN
            udfl_porous_user=.TRUE.
            udfl_grid_user=.TRUE.
         ENDIF     
         IF(vv_prob.eq.'halden650_5')udfl_porous_user=.TRUE.
!                                  
!--------User input head----------------------------------------------------------------------------
!
!........User1_YHY
!
!........User2_PIK
         IF(vv_prob.eq.'opr1000_rv_lbloca')THEN
            udfl_porous_user=.TRUE.
            udfl_grid_user=.TRUE.
         ENDIF 
         IF(vv_prob.eq.'pwr_mc_poro'    .or.&
            vv_prob.eq.'apr1400_mc_poro'    .or.&
            vv_prob.eq.'opr1000_mc_poro'    )THEN
            udfl_init_variables=.TRUE.  
            udfl_porous_user=   .TRUE.
         ENDIF  
         IF(vv_prob.eq.'icarus2002') THEN  
            udfl_porous_user=.TRUE. 
            udfl_init_variables=.TRUE.
         ENDIF
!
!........User3_LJR
!
         IF(vv_prob.eq.'OPR1000_fullvessel_1x1'       .or. &
            vv_prob.eq.'OPR1000_fullcore_modmesh02_rv'      ) then
            udfl_porous_user=.TRUE.
            udfl_grid_user=.TRUE.
         ENDIF
!
         IF(vv_prob.eq.'OPR1000_fullcore_modmesh02_rv_vessel' .or. &
            vv_prob.eq.'OPR1000_single_assem'                       ) then
            udfl_porous_user   =.true.
            udfl_grid_user     =.true.
           !udfl_init_variables=.TRUE.
         ENDIF
!
!........User5_CYJ
!
         IF(vv_prob.eq.'h2p1_0'.or.vv_prob.eq.'h2p1_0x'.or. &
            vv_prob.eq.'h2p1_1'.or.vv_prob.eq.'h2p1_1x'.or. &
            vv_prob.eq.'h2p1_2'.or.vv_prob.eq.'h2p1_2x'.or. &
            vv_prob.eq.'h2p1_3'.or.vv_prob.eq.'h2p1_3x'.or. &
            vv_prob.eq.'h2p1_4'.or.vv_prob.eq.'h2p1_4x'.or. &
            vv_prob.eq.'VD_h2p1_0')THEN
            udfl_init_variables=.TRUE.
         ENDIF
!
         IF(vv_prob.eq.'ST2-CT-01' .or. &
            vv_prob.eq.'ST2-CT-02' .or. &
            vv_prob.eq.'ST2-CT-03'       )THEN
            udfl_init_variables=.TRUE.
            udfl_mom_press_source=.TRUE.
            udfl_mom_drag_i      =.TRUE.
            udfl_erg_diff        =.TRUE.
            udfl_update_scalar   =.TRUE.
            udfl_tw_profile      =.TRUE.
         ENDIF
!
!........User6_LSJ
!
         IF(vv_prob.eq.'KSMR') THEN
            !l_subchannel=.TRUE.            
         END IF  
         IF(vv_prob.eq.'KSMR-SG-porous') THEN
            udfl_porous_user=.true. !porous for SG
         END IF   
         IF(vv_prob.eq.'KSMR-PZR') THEN
            udfl_init_variables=.TRUE.
         END IF          
         IF(vv_prob.eq.'Nuscale-RRV'.or.vv_prob.eq.'Nuscale-02')THEN
            udfl_init_variables=.TRUE.
         ENDIF      
         IF(vv_prob.eq.'Nuscale-RVV'.or.vv_prob.eq.'Nuscale-PZR')THEN
            udfl_init_variables=.TRUE.
            udfl_porous_user=.true. !for icore        
         ENDIF   
         IF(vv_prob.eq.'Nuscale-05')THEN
            udfl_init_variables=.TRUE.
            udfl_porous_user=.true. !for icore        
         ENDIF          
!
!........User6_DSJ
!
!------------------------------------------------------------------------                     
         RETURN
      ENDIF                    
!
!.....Initialize user-defined SUBROUTINE
!
!
!.....cell property
      udfl_set_qvol_porous=.FALSE.     !volumetric heat source; udfn_set_qvol_porous & set_vol_heat_source
      udfl_porous_property=.FALSE.     !heat transfer area in porous media; udfn_porous_property & calc_models
      udfl_mat_prop=.FALSE.            !solid property; udfn_mat_prop & mat_prop
!
!.....boundary condition 
      udfl_outlet_property=.FALSE.     !change pressure boundary properties; udfn_outlet_property & set_outlet_property
      udfl_outlet_press_user=.FALSE.   !indicator of setting outlet pressure by user; pressure_solve
      udfl_tw_profile=.FALSE.          !profiled twall; imp_eng_diffusion, radiation_model, wall_condensation_steam &udfn_tw_profile
      udfl_hflux_bc_profile=.FALSE.    !udfn_hflux_bc_profile_chw,udfn_hflux_bc_profile_chw_c
!
!.....momentum transfer      
      udfl_mom_wall=.FALSE.            !wall friction; udfn_mom_wall & wall_drag
      udfl_mom_drag_i=.FALSE.          !interfacial drag; udfn_mom_drag_i & int_drag
      udfl_mom_press_source=.FALSE.    !momentum pressure source; udfn_mom_press_source & calc_momentum
      udfl_mom_loss=.FALSE.            !momentum source;  udfn_mom_source & calc_momentum
      udfl_mom_film_shear=.FALSE.      !momentum source of wall shear; udfn_mom_film_shear & calc_momentum --> input
!
!.....heat transfer      
      udfl_calc_HTC_int_i=.FALSE.      !interfacial heat transfer; udfn_calc_HTC_int_i & int_htc
      udfl_flashing_hif=.FALSE.        !indicator of flashing only when udfl_calc_HTC_int_i=.true.; udfn_calc_HTC_int_i
      udfl_wallHTC_porous=.FALSE.      !heat transfer in porous media; udfn_heat_wallHTC_porous & heat_partition
      udfl_erg_diff=.FALSE.            !energy diffusion; udfn_erg_diff & scalar_energy_diffusion 
!      
!.....overwrite heat and momentum transfer model           
      udfl_model_overwrite=.FALSE.     !overwrite all model; udfn_model_overwrite & int_swap(11)
!
!.....change scalar_update
      udfl_update_scalar=.FALSE.       !change update_scalar; udfn_update_scalar & scalar_update 
!
!.....change heat partition model
      udfl_psbt_cfx_model=.FALSE.      !specific model for heat partition; heat_partition_2_bi_lo, heat_partition_2_bi_lo_err
!
!-----VnV input tail-----------------------------------------------------------------------------
!
!
!.....UPTF-RV problem has basically to be set as  nin_max=5, reflood=0, domain_decomposition=1.
      IF (vv_prob.eq.'UPTF-RV') THEN  
         IF(nin_max.ne.5) THEN
            PRINT*,'UPTF problem should set nin_max=5 and check reflood,domain_decomposition as well'
            PAUSE
            STOP
         ENDIF
         CALL UPTF_vel_bc(time)
         hydraulicd(:)=0.5d0   
     ELSEIF (vv_prob.eq.'rbht1196_3d'.or.vv_prob.eq.'rbht1196_1d') THEN
!                  
         IF(initial) THEN 
            OPEN(50,file='ht_bc_transient_1d.in',status='old',iostat=err)               
            initial=.false.     
            time_bc_next=0.d0 
         ENDIF
         IF(err.eq.0)THEN         
            IF(time.ge.time_bc_next) THEN
               Print*,'--Change the BC transient @ t=',time
               READ(50,*) time_bc,time_bc_next,Vlin,Pin,Tgin,Tlin,Tdin,Agin,Adin,Pout,Agout
!
!              inlet                       
               p_fb(1)=Pin
               tb_gas(1)=Tgin
               tb_liq(1)=Tlin
               alphab_gas(1)=Agin
               alphab_liq(1)=1.d0-alphab_gas(1)
               alphab_drp(1)=0.d0
               qualab(1)=0.d0 !0.99d0 !0.d0
               qn_cell0(:)=qn_cell(1,:)
               CALL convert_temp2erg(p_fb(1),tb_liq(1),tb_gas(1),qualab(1),eb_liq(1), &
                                      eb_gas(1),rhob_liq(1),rhob_gas(1),tmp1,tmp2,    &
                                      tao,cvao_nvin(1),uao_nvin(1),dcva_nvin(1),ra_nvin(1),qn_cell0)
               vb_liq(1,ndim)=Vlin/0.004864d0/rhob_liq(1)   !m/s
               vb_gas(1,ndim)=vb_liq(1,ndim)               
!                  
!              outlet                   
               pbnd(1)=Pout
               alpha_gas_nd(1)=Agout
               alpha_liq_nd(1)=0.d0 !1.d0-alpha_gas_nd(1)
               alpha_drp_nd(1)=0.d0
!
            ENDIF  
         ENDIF
!                    
     ELSEIF (vv_prob.eq.'fs_31203'.or.vv_prob.eq.'fs_31203_3D'.or.&
             vv_prob.eq.'fs_31302'.or.vv_prob.eq.'fs_31302_3D'.or.&
             vv_prob.eq.'fs_31701'.or.vv_prob.eq.'fs_31701_3D'.or.&
             vv_prob.eq.'fs_31805'.or.vv_prob.eq.'fs_31805_3D') THEN
         udfl_mom_loss=.TRUE. 
         eps_rho=8.d-3
         eps_p=30.d3
         eps_eng=1.d-1
         eps_vol=1.d-6
!         
!.....Gap Conductance Test: LSJ
!
      ELSEIF (vv_prob.eq.'gap_conductance') THEN
!
      ELSEIF (vv_prob.eq.'2D_boiling'.or.vv_prob.eq.'3D_boiling') THEN
         udfl_set_qvol_porous=.TRUE.
         hydraulicd(:)=0.2d0
!
      ELSEIF(vv_prob.eq.'nat_conv_krane') THEN
!
      ELSEIF (vv_prob.eq.'flashing') THEN
         udfl_calc_HTC_int_i=.TRUE.
         udfl_flashing_hif=.TRUE.         
         IF (time.le.10.d0) THEN
            pbnd(1)=1.d0-(1.0d0-0.854d0)*time/10.d0  ! before : 0.5 -> 0.854
            pbnd(1)=pbnd(1)*1.0d6
         ENDIF
!         
      ELSEIF (vv_prob.eq.'cavitation') THEN
!         udfl_outlet_property=.TRUE.
!     
      ELSEIF (vv_prob.eq.'plume') THEN
         udfl_outlet_property=.TRUE.
         DO i=1,ncell_fluid
            j0=i_neigh(i)-1
            DO j=i_neigh(i),i_neigh(i+1)-1
               IF(nbcon(j).eq.1)THEN
                  IF(xloc(i,1).lt.0.038d0.or.xloc(i,1).gt.0.042d0) nbcon(j)=-1
               ENDIF
            ENDDO
         ENDDO            
!      
      ELSEIF (vv_prob.eq.'2D_loca') THEN
!     
         udfl_set_qvol_porous=.TRUE.
         IF(time.gt.300.) CALL sbloca_user
!
      ELSEIF (vv_prob.eq.'2D_laminar') THEN
         DO i=1,nvin
            rhob_liq(i)=1000.0d0
         ENDDO
         DO i=1,ncell_fluid
            cell%lviscosl(i)=0.1d0
            cell%rhol(i)=1000.0d0
         ENDDO
         DO i=1,ncell_fluid
            j0=i_neigh(i)-1
            DO j=i_neigh(i),i_neigh(i+1)-1
               IF(nbcon(j).eq.5) npb(i)=1
            ENDDO
         ENDDO
!
      ELSEIF (vv_prob.eq.'2D_conduction') THEN
         udfl_mat_prop=.true.
         DO i=1,nvin
            rhob_liq(i)=1000.0d0
         ENDDO
         DO i=1,ncell_fluid
            cell%lviscosl(i)=0.1d0
            cell%lcondl(i)=3000.0d0
            cell%rhol(i)=1000.0d0
         ENDDO
!        
      ELSEIF (vv_prob.eq.'Horizontal_flow') THEN  
         udfl_mat_prop=.true.
         DO i=1,nvin
            rhob_liq(i)=1000.0d0
         ENDDO
         DO i=1,ncell_fluid
            cell%lviscosl(i)=0.1d0
            cell%lcondl(i)=3000.0d0
            cell%rhol(i)=1000.0d0
         ENDDO
!         l_horizontal_outlet=.TRUE.
!
      ELSEIF (vv_prob.eq.'separation') THEN
         DO i=1,ncell_fluid
            cell%rhol(i)=1000.0d0
         ENDDO
!
      ELSEIF (vv_prob.eq.'dam_break') THEN
         IF(itim.eq.1) THEN
            DO i=1,ncell_fluid
              IF(xloc(i,1).gt.0.05d0.or.xloc(i,2).gt.0.1d0) THEN
                  cell%alphal(i)=0.001d0
                  cell%alphag(i)=0.999d0-1.d-8
                  cell%alphad(i)=1.d-8
               ELSE
                  cell%alphal(i)=0.999d0-1.d-8
                  cell%alphag(i)=0.001d0
                  cell%alphad(i)=1.d-8
             ENDIF
            ENDDO
            udfl_mom_drag_i=.TRUE.
         ENDIF         

      ELSEIF(vv_prob.eq.'3D_air/water') THEN
!
      ELSEIF(vv_prob.eq.'block_porous') THEN
         udfl_porous_property=.true.
!         
      ELSEIF(vv_prob.eq.'annul_porous') THEN
         udfl_set_qvol_porous=.true. 
         udfl_porous_property=.true.         
!
      ELSEIF (vv_prob.eq.'boron_trans') THEN
         CALL initialize_time_variant_inlet
!
      ELSEIF (vv_prob.eq.'rocom' .or. vv_prob.eq.'rocom_mc') THEN
         time=time-dt
         IF(vv_prob.eq.'rocom'.and.iturb.eq.2.and.time.le.1.d0) THEN
            dt=MIN(dt,1.d-3) 
            dtr=1.d0/dt
         ENDIF
         time=time+dt
         udfl_mom_loss=.TRUE.
         IF(vv_prob.eq.'rocom')THEN
             udfl_outlet_property=.true.
             udfl_outlet_press_user=.true. 
         ENDIF
         IF(initial)THEN
            initial=.FALSE.
            DO i=1,ncell_fluid 
               ii=0
               j0=i_neigh(i)-1
               DO j=i_neigh(i),i_neigh(i+1)-1
                  IF(nbcon(j).lt.0)ii=ii+1
               ENDDO
               IF(ii.ge.4)then !it is a leg 
                  lsindex(i)=0          
                  gfactor(i)=0.0d0 
                  WRITE(*,*)'          1-D cell=',jperm(i)
               ENDIF    
            ENDDO             
         ENDIF  
!
      ELSEIF (vv_prob.eq.'fluidic_device') then
         udfl_mom_press_source=.TRUE.
         udfl_erg_diff=.TRUE. 
         udfl_init_variables=.TRUE. 
         udfl_mom_drag_i=.TRUE.
         udfl_update_scalar=.TRUE.
         udfl_flashing_hif=.TRUE.
         udfl_mom_wall=.TRUE.                           
!
      ELSEIF (vv_prob.eq.'PAFS-POOL') THEN
         IF(INITIAL)THEN
            IF(myrank.eq.0) OPEN(41,file='project.flavia.res')          
            CALL GID_out_cell_lsj(time,0)         
            CALL initialize_specific_variables_solid
            CALL initialize_least_square_option  !just for free surface
            INITIAL=.false.
         ENDIF  
         udfl_set_qvol_porous=.TRUE.
         udfl_mat_prop=.TRUE.
         udfl_mom_press_source=.TRUE.
         udfl_erg_diff=.TRUE. 
         udfl_porous_property=.TRUE.
         udfl_init_variables=.TRUE. 
         udfl_mom_drag_i=.TRUE.
         udfl_update_scalar=.TRUE.
         udfl_flashing_hif=.TRUE.
         udfl_mom_wall=.TRUE.                           
         udfl_calc_HTC_int_i=.TRUE.
!
      ELSEIF(vv_prob.eq.'stern') THEN
         udfl_mom_wall=.true.
         udfl_wallHTC_porous=.true.
         udfl_porous_property=.true.
!
      ELSEIF(vv_prob.eq.'PSBT_sngl')THEN
         udfl_psbt_cfx_model=.true.      
         udfl_calc_HTC_int_i=.true.   
         IF(initial)THEN
            initial=.false.
         ENDIF                 
!
      ELSEIF (vv_prob.eq.'sgp_separator') THEN
!
         udfl_outlet_property=.true.
         udfl_mom_loss=.TRUE.
         udfl_set_qvol_porous=.true.
         CALL udfn_sg_pre
         IF(initial_sgp)THEN
            CALL udfn_sg_input
            initial_sgp=.false.
         ENDIF
!
      ELSEIF (vv_prob.eq.'multi_ncg') THEN
!
         IF(time.gt.100.d0)THEN
            qn_nvin(:,1)=0.0d0 
            qn_nvin(:,2)=1.0d0 
            IF(INITIAL) THEN
                CALL ncg_cell
                CALL init_steamtable
                INITIAL=.FALSE.
            ENDIF
         ENDIF
!         
      ELSEIF (vv_prob.eq.'mult_ncg_2nd_conv') THEN
!
         IF(time.gt.500.d0)THEN
            qn_nvin(:,1)=0.0d0 
            qn_nvin(:,3)=1.0d0
            vb_gas(1,3)=0.5d0
            IF((time-500.d0).ge.0.d0)vb_gas(1,3)=DMIN1(0.5d0,0.5d0*(time-500.d0)) !New SMAC3
            IF(INITIAL) THEN
                CALL ncg_cell
                CALL init_steamtable
                INITIAL=.FALSE.
            ENDIF
         ENDIF
!         
      ELSEIF (vv_prob.eq.'moving_wall')THEN
!
         alpha_drp_nd(1)=0.0d0
!
      ELSEIF (vv_prob.eq.'choking_edward') THEN
!
         hydraulicd(:)=0.0675d0         
!        
      ELSEIF (vv_prob.eq.'Nuscale-03Pool') THEN
!
         ii=0
         pcnv=0.0d0
         DO i=1,ncell_fluid
           IF(nzone(i).eq.1)THEN
              pcnv=pcnv+cell%p(i)
              ii=ii+1
           ENDIF
         ENDDO
         
         IF(np.gt.1)THEN
            CALL allreducei_i1(ii)
            CALL allreducei_r1(pcnv)
         ENDIF
         pcnv=pcnv/DBLE(ii)
!
         IF(pcnv.lt.10.0D6)THEN
            vb_gas(1,2)=1.0d-1 
         ELSE
            vb_gas(1,2)=0.0d0
         ENDIF
!
         IF(nf_ctw.gt.0)THEN
            nf_number=6
            istart=istart_nf(1,nf_number)
            isize =istart_nf(2,nf_number)
         ELSEIF(nf_fsw.gt.0)THEN
            nf_number=5
            istart=istart_nf(1,nf_number)
            isize =istart_nf(2,nf_number)         
         ENDIF          
         DO i=1,isize
            IF(nf_ctw.gt.0)THEN
               dfilm_ctw(i)=0.001d0
               vfilm_ctw(i)=0.1d0
            ELSEIF(nf_fsw.gt.0)THEN
               dfilm_fsw(i)=0.0001d0 !0.3d0
               vfilm_fsw(i)=0.0001d0 !0.4d0
            ENDIF
         ENDDO
!
      ELSEIF(vv_prob.eq.'kuhn_111')THEN
         udfl_tw_profile=.TRUE.
!........only for VFT14_kuhn_fs
         IF(initialKuhn)THEN
            zkuhn=(/2.418d0, 2.315d0, 2.181d0, 2.047d0, 1.897d0, 1.709d0, &
                    1.531d0, 1.313d0, 1.097d0, 0.837d0, 0.569d0, 0.237d0, 0.0d0/)
            zTwall=(/366.15d0, 366.15d0, 365.45d0, 364.25d0, 364.45d0, 363.75d0, &
                     362.55d0, 362.05d0, 361.85d0, 361.85d0, 361.85d0, 361.85d0/)
            initialKuhn=.false.
         ENDIF
!
         DO i=1,ncell_fluid
            IF(wVertical(i).eq.1)THEN
               DO j=1,12
                  IF(xloc(i,3).lt.zkuhn(j).and.xloc(i,3).gt.zkuhn(j+1))THEN
                     face%twall_partition(i)=zTwall(j)
               ENDIF
               ENDDO
            ENDIF
         ENDDO
!
         IF(INITIAL)THEN
            zkuhn=(/2.418d0, 2.315d0, 2.181d0, 2.047d0, 1.897d0, 1.709d0, &
                    1.531d0, 1.313d0, 1.097d0, 0.837d0, 0.569d0, 0.237d0, 0.0d0/)
            INITIAL=.false.
         ENDIF
!
         IF(iflag.ne.1)RETURN
!        
         IF(nf_ctw.gt.0)THEN
            nf_number=6
            istart=istart_nf(1,nf_number)
            isize =istart_nf(2,nf_number)
         ELSEIF(nf_fsw.gt.0)THEN
            nf_number=5
            istart=istart_nf(1,nf_number)
            isize =istart_nf(2,nf_number)         
         ENDIF
!         
         DO i=1,isize
            i1=istart+i
            ii=left_nf(i1)
            DO jj=1,12
               IF(xloc(ii,3).lt.zkuhn(jj).and.xloc(ii,3).gt.zkuhn(jj+1))THEN
                  tmp1=0.0001d0+0.0001d0*(jj-1)
                  tmp1=MIN(tmp1,0.001)
               ENDIF
            ENDDO  
            IF(nf_ctw.gt.0)THEN
               dfilm_ctw(i)=tmp1
               vfilm_ctw(i)=0.5d0
            ELSEIF(nf_fsw.gt.0)THEN
               dfilm_fsw(i)=tmp1
               vfilm_fsw(i)=0.5d0
            ENDIF
         ENDDO    
!
      ELSEIF (vv_prob.eq.'apr1400_mc_rv'.or.&
              vv_prob.eq.'opr1000_mc_rv') THEN
        IF(initial)THEN
            iprn=0
            !r1,2,3,4=4.05,2.25,6.05,8.55
            IF(myrank.eq.0)WRITE(*,"(11x,a)")'lsindex and gfactor=0 for 1D or LEG cells!'
            initial=.FALSE.
            DO i=1,ncell_fluid 
               ii=0
               j0=i_neigh(i)-1
               DO j=i_neigh(i),i_neigh(i+1)-1
                  IF(nbcon(j).lt.0)ii=ii+1
               ENDDO
               IF(ii.ge.4)then 
                  lsindex(i)=0          
                  gfactor(i)=0.0d0 
               ENDIF    
            ENDDO 
            IF(vv_prob.eq.'opr1000_mc_rv')THEN !1.75=downcommer, 2.05=vessel
               rr0=2.05d0
            ELSEIF(vv_prob.eq.'apr1400_mc_rv')THEN !
               rr0=2.31d0
            ENDIF
            DO i=1,ncell_fluid 
               xx=xloc(i,1)
               yy=xloc(i,2)
               rr=xx*xx+yy*yy
               rr=DSQRT(rr)
               IF(rr.gt.rr0)then
                  lsindex(i)=0          
                  gfactor(i)=0.0d0 
               ENDIF    
            ENDDO
         ENDIF   
         udfl_mom_wall=.true. 
!
      ELSEIF(vv_prob.eq.'copain_porous')THEN
         udfl_porous_property=.true.
         solid%tsol(:)=300.06D0
!
      ELSEIF (vv_prob.eq.'apr1400_lbloca') THEN
         eps_rho=8.d-3 !8.d-3 !apr1400_lbloca_debug
         eps_p=20.d3   !20.d3
         eps_eng=1.d-1 !1.0d-1
         eps_vol=1.d-6 !1.0d-6
         stime_vup=0.0d0       
         stime_vflat=1.0d0         
         CALL apr1400_lbloca_ctrl(iflag)
!
      ELSEIF (vv_prob.eq.'halden650_5') THEN
         udfl_set_qvol_porous=.true.         
         DO i=1,ncell_fluid 
            IF(i.ge.40.and.i.le.43)THEN
              lsindex(i)=0          
              gfactor(i)=0.0d0 
            ENDIF
         ENDDO          
         CALL halden650_lbloca_ctrl(iflag)
!         
      ELSEIF (vv_prob.eq.'atlas_mc_porous') THEN
        IF(initial)THEN
            initial=.FALSE.
            CALL udfn_mat_prop
            HydraulicD(:)=0.01263 
            DO i=1,ncell_fluid 
               ii=0
               j0=i_neigh(i)-1
               DO j=i_neigh(i),i_neigh(i+1)-1
                  IF(nbcon(j).lt.0)ii=ii+1
               ENDDO
               IF(ii.ge.4)then 
                  lsindex(i)=0          
                  gfactor(i)=0.0d0 
                  WRITE(*,*)'          1-D cell=',jperm(i)
               ENDIF    
            ENDDO 
        ENDIF 
        udfl_set_qvol_porous=.true.
        udfl_porous_property=.true.
        udfl_mat_prop=.TRUE. 
        udfl_wallHTC_porous=.TRUE.
        udfl_model_overwrite=.TRUE.        
!
      ELSEIF (vv_prob.eq.'mass_check_mcc'.or.&
              vv_prob.eq.'single_channel_merged') THEN
!
     ELSEIF (vv_prob.eq.'check_couple_mcc'.or.&
             vv_prob.eq.'single_channel_dif_3D' ) THEN
!    
      ELSEIF(vv_prob.eq.'siphon')THEN
!
         udfl_outlet_property=.true. 
         udfl_model_overwrite=.TRUE. !CALL average_spatially(cell%vfgl)
         IF(initial)THEN
            initial=.false.
            topen=1.0d0
         ENDIF
         DO i=1,ncell_fluid 
            IF(nzone(i).eq.1)THEN !small pipe
               tleng(i)=0.03d0
            ELSEIF(nzone(i).eq.2)THEN !large pipe
               tleng(i)=tlengs
            ELSEIF(nzone(i).eq.3)THEN !tank
               tleng(i)=1.775d0
               tleng(i)=0.2d0   
            ENDIF   
          ENDDO 
!           
      ELSEIF (vv_prob.eq.'opr1000_rv') THEN
        IF(initial)THEN
            iprn=0
            !r1,2,3,4=4.05,2.25,6.05,8.55
            IF(myrank.eq.0)WRITE(*,"(11x,a)")'lsindex and gfactor=0 for 1D or LEG cells!'
            initial=.FALSE.
            HydraulicD(:)=0.01263 
            DO i=1,ncell_fluid 
               ii=0
               j0=i_neigh(i)-1
               DO j=i_neigh(i),i_neigh(i+1)-1
                  IF(nbcon(j).lt.0)ii=ii+1
               ENDDO
               IF(ii.ge.4)then 
                  lsindex(i)=0          
                  gfactor(i)=0.0d0 
               ENDIF    
            ENDDO 
            DO i=1,ncell_fluid 
               xx=xloc(i,1)
               yy=xloc(i,2)
               rr=xx*xx+yy*yy
               rr=DSQRT(rr)
               IF(rr.gt.4.05d0)then 
                  lsindex(i)=0          
                  gfactor(i)=0.0d0 
               ENDIF    
            ENDDO
         ENDIF   
         l_plate=3.6576d0         
         dia_rod=0.013046d0 
         pit_dia=0.012852d0
         udfl_mom_wall=.true. 
         udfl_model_overwrite=.TRUE.         
         IF(restart.ne.0)CALL opr1000_coolant_transient_user
!         
      ELSEIF (vv_prob.eq.'apr1400_rv') THEN
        IF(initial)THEN
            iprn=0
            !r1,2,3,4=4.05,2.25,6.05,8.55
            IF(myrank.eq.0)WRITE(*,"(11x,a)")'lsindex and gfactor=0 for 1D or LEG cells!'
            initial=.FALSE.
            HydraulicD(:)=0.01263 
            DO i=1,ncell_fluid 
               ii=0
               j0=i_neigh(i)-1
               DO j=i_neigh(i),i_neigh(i+1)-1
                  IF(nbcon(j).lt.0)ii=ii+1
               ENDDO
               IF(ii.ge.4)then 
                  lsindex(i)=0          
                  gfactor(i)=0.0d0 
               ENDIF    
            ENDDO 
            DO i=1,ncell_fluid 
               xx=xloc(i,1)
               yy=xloc(i,2)
               rr=xx*xx+yy*yy
               rr=DSQRT(rr)
               IF(rr.gt.7.3d0)then 
                  lsindex(i)=0          
                  gfactor(i)=0.0d0 
               ENDIF    
            ENDDO
         ENDIF   
         l_plate=3.6576d0         
         dia_rod=0.013046d0 
         pit_dia=0.012852d0
         udfl_mom_wall=.true. 
      ELSEIF (vv_prob.eq.'apr1400core_rv') THEN
         !l_plate=3.6576d0         
         !dia_rod=0.013046d0 
         !pit_dia=0.012852d0
        !udfl_mom_wall=.true. 
         udfl_model_overwrite=.TRUE.         
!         
      ELSEIF(vv_prob.eq.'OPR1000_fullcore_modmesh02_rv') THEN
         udfl_mom_wall=.true.

      ELSEIF(vv_prob.eq.'OPR1000_fullcore_modmesh02_rv_vessel') THEN
!
         IF(cplmars.ne.0)then
            IF(initial)THEN
               iprn=0
               !r1,2,3,4=4.05,2.25,6.05,8.55
               IF(myrank.eq.0)then
                  WRITE(*,"(11x,a)")'lsindex and gfactor=0 for 1D or LEG cells!'
               ENDIF !'
               initial=.FALSE.
               DO i=1,ncell_fluid
                  ii=0
                  j0=i_neigh(i)-1
                  DO j=i_neigh(i),i_neigh(i+1)-1
                     IF(nbcon(j).lt.0)ii=ii+1
                  ENDDO
                  IF(ii.ge.4)then
                     lsindex(i)=0
                     gfactor(i)=0.0d0
                  ENDIF
               ENDDO
               DO i=1,ncell_fluid
                  xx=xloc(i,1)
                  yy=xloc(i,2)
                  rr=xx*xx+yy*yy
                  rr=DSQRT(rr)
                  IF(rr.gt.4.05d0)then
                     lsindex(i)=0
                     gfactor(i)=0.0d0
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         udfl_mom_wall=.true.

         IF(rv_model.eq.0)then
            udfl_set_qvol_porous=.true.
         ENDIF
!
         IF(cplmars.eq.0 .and. restart.ne.0)then
            CALL opr1000_coolant_transient_user
         ENDIF
!
      ELSEIF(vv_prob.eq.'OPR1000_single_assem') THEN
         udfl_erg_diff=.TRUE. 

         IF(l_mixing_vane)then
            fric_face=1
         ELSE
            fric_face=0
         ENDIF

!  
!-----User input tail---------------------------------------------------------------------------------------
!
!........User1_YHY
      ELSEIF(vv_prob.eq.'check_Hik')THEN
         udfl_calc_HTC_int_i=.true.
         s_bubble_diameter='hibiki'
         CALL udfn_check_Hik
!
!........User2_PIK
      ELSEIF (vv_prob.eq.'opr1000_rv_lbloca') THEN
        IF(initial)THEN
            iprn=0
            !r1,2,3,4=4.05,2.25,6.05,8.55
            IF(myrank.eq.0)WRITE(*,"(11x,a)")'lsindex and gfactor=0 for 1D or LEG cells!'
            initial=.FALSE.
            DO i=1,ncell_fluid 
               ii=0
               j0=i_neigh(i)-1
               DO j=i_neigh(i),i_neigh(i+1)-1
                  IF(nbcon(j).lt.0)ii=ii+1
               ENDDO
               IF(ii.ge.4)then 
                  lsindex(i)=0          
                  gfactor(i)=0.0d0 
               ENDIF    
            ENDDO 
            DO i=1,ncell_fluid 
               xx=xloc(i,1)
               yy=xloc(i,2)
               rr=xx*xx+yy*yy
               rr=DSQRT(rr)
               IF(rr.gt.4.05d0)then 
                  lsindex(i)=0          
                  gfactor(i)=0.0d0 
               ENDIF    
            ENDDO
         ENDIF 
         user_iary(31)=1
         user_rary(31)=100.0d0           
         user_iary(32)=2         
         udfl_mom_wall=.true. 
         udfl_model_overwrite=.TRUE.         
         CALL apr1400_lbloca_ctrl(iflag) 
!         
      ELSEIF (vv_prob.eq.'pwr_mc_poro'.or.vv_prob.eq.'apr1400_mc_poro'.or.vv_prob.eq.'opr1000_mc_poro') THEN
        IF(initial)THEN
            initial=.FALSE.
            CALL udfn_mat_prop
            HydraulicD(:)=0.01263 
            DO i=1,ncell_fluid 
               ii=0
               j0=i_neigh(i)-1
               DO j=i_neigh(i),i_neigh(i+1)-1
                  IF(nbcon(j).lt.0)ii=ii+1
               ENDDO
               IF(ii.ge.4)then 
                  lsindex(i)=0          
                  gfactor(i)=0.0d0 
                  WRITE(*,*)'          1-D cell=',jperm(i)
               ENDIF    
            ENDDO 
        ENDIF 
        s_turb_zero='telluride'
        udfl_set_qvol_porous=.true.
        udfl_porous_property=.true.
        s_wall_fric='darcy'           
        udfl_mat_prop=.TRUE. 
        udfl_wallHTC_porous=.TRUE.
        udfl_model_overwrite=.TRUE.        
!      
      ELSEIF(vv_prob.eq.'opr1000_mc')THEN
!            
         stime_vup=0.0d0       
         stime_vflat=1.0d0
         udfl_mom_wall=.true. !2019_07_17_debug-critical        
         udfl_set_qvol_porous=.true.
        IF(initial)THEN
            initial=.FALSE.
            DO i=1,ncell_fluid 
               ii=0
               j0=i_neigh(i)-1
               DO j=i_neigh(i),i_neigh(i+1)-1
                  IF(nbcon(j).lt.0)ii=ii+1
               ENDDO
               IF(ii.ge.4)then 
                  lsindex(i)=0          
                  gfactor(i)=0.0d0
               ENDIF    
            ENDDO 
            DO i=1,ncell_fluid 
               xx=xloc(i,1)
               yy=xloc(i,2)
               rr=xx*xx+yy*yy
               rr=DSQRT(rr)
               IF(rr.gt.4.05d0)then 
                  lsindex(i)=0          
                  gfactor(i)=0.0d0 
               ENDIF    
            ENDDO 
        ENDIF
!
      ELSEIF(vv_prob.eq.'icarus2002') THEN
!         
        udfl_mom_loss=.TRUE.      !icarus2002-next, mpi problem 
!
!........User3_LJR
!
!........User4_LSJ
!
      ELSE IF(vv_prob.eq.'KSMR') THEN
         udfl_mom_wall=.true. !KSMR pdrop model (calc_model.f90)
         udfl_set_qvol_porous=.true.
         IF(INITIAL) THEN
            INITIAL=.false.
!            
            ggc = DSQRT(DOT_PRODUCT(grav,grav ) )            
!            
            CALL ReadZone()
            CALL InitializeTemperature()
         END IF
         IF(rv_ht_str.ne.1) CALL AssignCoreHeat()
         call AssignSGHeat()
         IF(.false.) CALL ForceCoreInletTemperatrue()
         CALL PrintRcsFlowRate()
         CALL PrintTemperatureCoreOut()    
!         
      ELSE IF(vv_prob.eq.'KSMR-SG-porous') THEN
         udfl_mom_loss=.true. !SG anisotropic friction model (calc_momentum.f90)          
         udfl_mom_wall=.true. !KSMR pdrop model (calc_model.f90)
         udfl_set_qvol_porous=.true.
         IF(INITIAL) THEN
! sg volume             
            IF(np.gt.1) CALL allreducei_r1(vol_sg)  
            INITIAL=.false.
!            
            ggc = DSQRT(DOT_PRODUCT(grav,grav ) )            
!            
            CALL ReadZone()
            CALL InitializeTemperature()
         END IF
         
         CALL AssignCoreHeat()
         call AssignSGHeat()
!         IF(.true.) CALL ForceCoreInletTemperatrue()
         CALL PrintRcsFlowRate()
         CALL PrintTemperatureCoreOut()   
         CALL PrintPressure()

     ELSE IF(vv_prob.eq.'KSMR-SG-pid') THEN
!         udfl_mom_loss=.true. !SG anisotropic friction model (calc_momentum.f90)          
         udfl_mom_wall=.true. !KSMR pdrop model (calc_model.f90)
         udfl_set_qvol_porous=.true.
         IF(INITIAL) THEN
            INITIAL=.false.
            ggc = DSQRT(DOT_PRODUCT(grav,grav ) )            
            CALL ReadZone()
            CALL InitializeTemperature()
         END IF
         
         CALL AssignCoreHeat()
         call AssignSGHeat()
!         IF(.true.) CALL ForceCoreInletTemperatrue()
         CALL PrintRcsFlowRate()
         CALL PrintTemperatureCoreOut()            
         CALL PrintPressure()         
         
      ELSE IF(vv_prob.eq.'KSMR-PZR') THEN
         udfl_mom_wall=.true. !KSMR pdrop model (calc_model.f90)
         udfl_set_qvol_porous=.true.
         IF(INITIAL) THEN
            INITIAL=.false.
!            
            DO i=1,ncell_fluid
               qn_cell0(:)=qn_cell(i,:)
               CALL convert_temp2erg(cell%p(i),cell%tl(i),cell%tg(i),cell%quala(i),cell%el(i),cell%eg(i), &
                                     cell%rhol(i),cell%rhog(i),cell%pps_o(i),cell%estm_o(i),           &
                                     tao,cvao_cell(i),uao_cell(i),dcva_cell(i),ra_cell(i),qn_cell0)            
            ENDDO              
!            
            ggc = DSQRT(DOT_PRODUCT(grav,grav ) )            
!            
            CALL ReadZone()
            CALL InitializeTemperature()
         END IF
         CALL AssignCoreHeat()
         call AssignSGHeat()
!         IF(.true.) CALL ForceCoreInletTemperatrue()
         CALL PrintRcsFlowRate()
         CALL PrintTemperatureCoreOut()   
!
      ELSEIF(vv_prob.eq.'Nuscale-RRV'.or.vv_prob.eq.'Nuscale-02') THEN
         IF(initial) THEN          
            DO i=1,ncell_fluid
               qn_cell0(:)=qn_cell(i,:)
               CALL convert_temp2erg(cell%p(i),cell%tl(i),cell%tg(i),cell%quala(i),cell%el(i),cell%eg(i), &
                                     cell%rhol(i),cell%rhog(i),cell%pps_o(i),cell%estm_o(i),           &
                                     tao,cvao_cell(i),uao_cell(i),dcva_cell(i),ra_cell(i),qn_cell0)            
            ENDDO  
            initial=.false.
         ENDIF         
!         
      ELSE IF(vv_prob.eq.'Nuscale-RVV'.or.vv_prob.eq.'Nuscale-PZR') THEN
         udfl_mom_wall=.true.
         udfl_set_qvol_porous=.true.
!
         udfl_calc_HTC_int_i=.true.
         udfl_flashing_hif=.TRUE.         
!         
         IF(INITIAL) THEN
            INITIAL=.false.
!
            DO i=1,ncell_fluid
               qn_cell0(:)=qn_cell(i,:)
               CALL convert_temp2erg(cell%p(i),cell%tl(i),cell%tg(i),cell%quala(i),cell%el(i),cell%eg(i), &
                                     cell%rhol(i),cell%rhog(i),cell%pps_o(i),cell%estm_o(i),           &
                                     tao,cvao_cell(i),uao_cell(i),dcva_cell(i),ra_cell(i),qn_cell0)            
            ENDDO             
!            
            ggc = DSQRT(DOT_PRODUCT(grav,grav ) )            
!            
            CALL ReadZone_NuScale()
            CALL InitializeTemperature_NuScale()
         END IF
         CALL AssignCoreHeat_NuScale()
         call AssignSGHeat_NuScale()
!         IF(.true.) CALL ForceCoreInletTemperatrue_NuScale()
         CALL PrintRcsFlowRate_NuScale()
         CALL PrintTemperatureCoreOut_NuScale()    
!         
      ELSE IF(vv_prob.eq.'Nuscale-05') THEN
         udfl_mom_loss=.true. !SG anisotropic friction model (calc_momentum.f90)
         udfl_mom_wall=.true. !NuScale pdrop model (calc_model.f90)
         udfl_set_qvol_porous=.true.
!
         udfl_calc_HTC_int_i=.true.
         udfl_flashing_hif=.TRUE.         
!         
         IF(INITIAL) THEN
            INITIAL=.false.
!
            DO i=1,ncell_fluid
               qn_cell0(:)=qn_cell(i,:)
               CALL convert_temp2erg(cell%p(i),cell%tl(i),cell%tg(i),cell%quala(i),cell%el(i),cell%eg(i), &
                                     cell%rhol(i),cell%rhog(i),cell%pps_o(i),cell%estm_o(i),           &
                                     tao,cvao_cell(i),uao_cell(i),dcva_cell(i),ra_cell(i),qn_cell0)            
            ENDDO             
!            
            ggc = DSQRT(DOT_PRODUCT(grav,grav ) )            
!            
            CALL ReadZone_NuScale()
            CALL InitializeTemperature_NuScale()
         END IF
         CALL AssignCoreHeat_NuScale()
         call AssignSGHeat_NuScale()
!         IF(.true.) CALL ForceCoreInletTemperatrue_NuScale()
         CALL PrintRcsFlowRate_NuScale()
         CALL PrintTemperatureCoreOut_NuScale()            
         
      ELSEIF (vv_prob.eq.'mcp_ftype2') THEN     
         rcs_flow = rwinit(1)
         rcs_rho  = cell%rhol(1)  
         rcs_head = rcs_rho*9.8d0*24.d0 !rated_pump_hd(1)*1.d0 !rho*g*h
         rcs_area = 0.25d0 !m2 
         rcs_leng = 12.d0 !m
         dpdz     = rcs_head / rcs_leng
         vel      = rcs_flow / rcs_rho / rcs_area              
         
         DO i=1,ncell_fluid
            resist = dpdz/vel/vel
            cell%vfwl(i) = resist*ul_o(i)
         ENDDO              
!                 
!
!........User5_CYJ
!   
      ELSEIF(vv_prob.eq.'h2p1_0'.or.vv_prob.eq.'h2p1_0x'.or. &
             vv_prob.eq.'h2p1_1'.or.vv_prob.eq.'h2p1_1x'.or. &
             vv_prob.eq.'h2p1_2'.or.vv_prob.eq.'h2p1_2x'.or. &
             vv_prob.eq.'h2p1_3'.or.vv_prob.eq.'h2p1_3x'.or. &
             vv_prob.eq.'h2p1_4'.or.vv_prob.eq.'h2p1_4x'.or. &
            vv_prob.eq.'VD_h2p1_0')THEN
         udfl_tw_profile=.TRUE.
         CALL udfn_H2P1_inputTg
         CALL udfn_H2P1_inputTw
!
      ELSEIF(vv_prob.eq.'ST2-CT-01' .or. &
             vv_prob.eq.'ST2-CT-02' .or. &
             vv_prob.eq.'ST2-CT-03'       )THEN
         udfl_tw_profile=.TRUE.
         udfl_mom_press_source=.TRUE.
         udfl_mom_drag_i      =.TRUE.
         udfl_erg_diff        =.TRUE.
         udfl_update_scalar   =.TRUE.
         CALL udfn_CUBE_inputBC
         CALL udfn_CUBE_inputTW
         
!
!........User6_DSJ
!---------------------------------------------------------------------------------------------------
!
!     
!-----End of user input tail---------------------------------------------------------------------------------------
   
      ENDIF !IF(vv_prob.~
!    
!.....Ramping the wall heat flux B.C in time
!
      IF(time.le.stime_hflat+dt)THEN
         CALL wall_heat_flux
      ENDIF
!
!.....Ramping the velocity B.C. in time
!
      IF(time.le.stime_vflat+dt)THEN
         CALL udfn_vel_bc_ramp
      ENDIF
!
      END SUBROUTINE user_def_inp
