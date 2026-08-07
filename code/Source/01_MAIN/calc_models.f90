!
      SUBROUTINE calc_models
!
!     This routine includes model-involved calls;
!     1) IAT model
!     2) Turbulence model
!     3) Heat partitioning model
!     4) Mass diffusivity model 2015.07.27 JHLee (SNU)
!
      USE Zconst1      , ONLY: iat,wconden
      USE Zcore        , ONLY: myrank
      USE Zncg         , ONLY: ncg_diff
      USE Zuserdefined , ONLY: udfl_porous_property
      USE Zmodel       , only: rad_model
      USE Zuserdefined , ONLY: udfl_mom_wall           
      USE Zrv_model    , ONLY: rv_model,rv_ht_i,rv_fric_i,free_model,lfric_swap,lhtc_swap,rv_mcp,rv_choke,rv_valve
      USE Zqvol        , ONLY: qporous_liq,qporous_gas,qporous_gamma,gamma_wall       
      USE Zio_unit     , ONLY: unit_log
      USE Zrad_comp    , ONLY: rad_comp_mod
      
!
      IMPLICIT NONE     
!
      LOGICAL initial
      DATA initial/.TRUE./   
!
      CALL vectorize_scalar_upwind
!
!.....fluxBC model: choke model, mcp model, valve model
!      
      IF(rv_valve.eq.1.or.rv_choke.eq.1.or.rv_mcp.eq.1) CALL fluxBC_init        
!
      IF(rv_model.lt.1.or.(rv_fric_i.lt.1.or.rv_ht_i.lt.1).or.(lfric_swap.or.lhtc_swap))THEN
         free_model=.TRUE.      
      ELSE
         free_model=.FALSE.
      ENDIF 
!
      qporous_liq(:)=0.0d0
      qporous_gas(:)=0.0d0 
      qporous_gamma(:)=0.0d0
      gamma_wall(:)=0.0d0       
!
!.....Assign heat transfer area of porous body in user-defined function
!
      IF(udfl_porous_property) CALL udfn_porous_property
!      
      IF(free_model)THEN      
!      
!........Interfacial area transport
!
         IF(iat.gt.0) CALL IAT_calc
!
!........Flow regime, Interfacial Area, Entrainment/De-entrainment
!
         CALL int_area
      
      ENDIF
!
!.....Turbulence models
!
      CALL turbulence_mod        
      
!      
      IF(free_model)THEN   
!
!........Interface drag coefficients
!
         CALL int_drag
!
!........Calculate interface heat transfer coefficients
!
         CALL int_htc
!
!........Wall drag coefficients
!
         IF(udfl_mom_wall)THEN
            CALL udfn_mom_wall
         ELSE   
            CALL wall_drag
         ENDIF
!         
      ELSE
!
!........Set variables which are not initialized due to skipping above open media model.
!      
         CALL int_swap(0)
!
!........SMR
!
         IF(udfl_mom_wall) CALL udfn_mom_wall  
!           
      ENDIF
!    
!.....Store interfacial drag and heat transfer coefficients of open media.
!
      IF(rv_model.eq.1) CALL int_swap(1)
!      
!.....Wall condensation model      
!      
      IF(wconden.ne.0) CALL wall_condensation_model 
!
!.....Radiation model
!
      IF(rad_model.eq.1) CALL radiation_model
!
!.....Radation model between zones
!      
      IF(rad_comp_mod.eq.1)CALL radiation_component
!
!.....Convective heat transfer from solid wall
!
      CALL convective_heat_solid      
!
!.....Heat Partitioning Model
!
      CALL heat_partition
!
!.....Set the volumetric heat source for fluid or porous from input
!
      CALL set_vol_heat_source       
!      
!.....Mass diffusivity model 2015.07.27 JHLee (SNU) 
!
      IF(ncg_diff.gt.0) CALL mass_diff_mod      
!
!.....Reactor vessel models
!
      IF(rv_model.eq.1) CALL reactor_vessel_models
!
!.....Swap interfacial drag and heat transfer coefficients of open media and RV media.
!
      CALL int_swap(11)
!
      IF(initial)THEN
         initial=.FALSE.
         IF(.not.free_model.and.myrank.eq.0)WRITE(*,"(11x,a)")'Open media model is off.'
         IF(.not.free_model.and.myrank.eq.0)WRITE(unit_log,"(11x,a)")'Open media model is off.'
      ENDIF      
!
      RETURN
      END SUBROUTINE calc_models
