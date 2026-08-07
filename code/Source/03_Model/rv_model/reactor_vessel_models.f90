!
      SUBROUTINE reactor_vessel_models
!
      USE Zqvol       , ONLY: qporous_gas,qporous_liq,gamma_wall   
      USE Zqvol       , ONLY: qporous_gamma,qrv_gas,qrv_liq,qrv_gamma    
      USE Zrv_model   , ONLY: rv_fw_reg,rv_ht_str,rv_ht_w,rv_fric_i,rv_ht_i,rv_fric_w
      USE Zrv_hts_1d  , ONLY: ncell_hts_1d    
      USE Zrv_hts_2d  , ONLY: nrod_2d  
      USE Zmpi        , ONLY: ncell_fp
!
      IMPLICIT NONE
!
      INTEGER i
      REAL(8) qporous_liq_save(ncell_fp),qporous_gas_save(ncell_fp),gamma_wall_save(ncell_fp),qporous_gamma_save(ncell_fp)
!
      CALL rv_power_transient 
!
      IF(rv_fw_reg.eq.1) CALL rv_flow_regime
!
      DO i=1,ncell_fp
         qporous_liq_save(i)=qporous_liq(i) 
         qporous_gas_save(i)=qporous_gas(i)
         gamma_wall_save(i) =gamma_wall(i)
         qporous_gamma_save(i)=qporous_gamma(i)
         qporous_liq(i)=0.0d0
         qporous_gas(i)=0.0d0
         gamma_wall(i) =0.0d0
      ENDDO
!      
      IF(rv_ht_w.eq.1)THEN
         IF(ncell_hts_1d.gt.0) CALL rv_wall_ht_1d(1) !left boundary
         IF(ncell_hts_1d.gt.0) CALL rv_wall_ht_1d(2) !right boundary
         IF(nrod_2d     .gt.0) CALL rv_wall_ht_2d
         qrv_liq(:)=qporous_liq(:)
         qrv_gas(:)=qporous_gas(:)
         qrv_gamma(:)=qporous_gamma(:)  
      ELSE
         qporous_liq(:)=qrv_liq(:)
         qporous_gas(:)=qrv_gas(:)
         qporous_gamma(:)=qrv_gamma(:)                
      ENDIF
!      
      DO i=1,ncell_fp
         qporous_liq(i)=qporous_liq(i)+qporous_liq_save(i) 
         qporous_gas(i)=qporous_gas(i)+qporous_gas_save(i)
         gamma_wall(i) =gamma_wall(i) +gamma_wall_save(i)
         qporous_gamma(i)=qporous_gamma(i)+qporous_gamma_save(i)
      ENDDO
!
      IF(rv_ht_str.eq.1) CALL rv_heat_structure
!
      IF(rv_fric_i.ne.0) CALL rv_int_friction   !1=MARS manual model, 2=MARS code model
!
      IF(rv_ht_i.eq.1)   CALL rv_int_ht
!
      IF(rv_fric_w.eq.1) CALL rv_wall_friction
!
      RETURN
      END SUBROUTINE reactor_vessel_models
