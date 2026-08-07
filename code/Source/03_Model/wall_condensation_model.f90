!
      SUBROUTINE wall_condensation_model
!
!     The latent heat transfer due to the condensation is reflected only.
!     Liquid film modeling is not considered.
!     Condensation wall condition : twall
!
      USE Zmpi            , ONLY: ncell_fp
      USE Zconst1         , ONLY: wconden,iturb
      USE Zqvol           , ONLY: qporous_liq,qporous_gas,qporous_gamma, &
                                  gamma_wall
      USE Zmodel          , ONLY: qconden
      
!
      IMPLICIT NONE
!
      REAL(8) qporous_liq_save(ncell_fp),qporous_gas_save(ncell_fp),qporous_gamma_save(ncell_fp)
      REAL(8) gamma_wall_save(ncell_fp)
!      
      qporous_liq_save(:)  =qporous_liq(:) 
      qporous_gas_save(:)  =qporous_gas(:)
      qporous_gamma_save(:)=qporous_gamma(:)
!      
      qporous_liq(:)  =0.0d0
      qporous_gas(:)  =0.0d0      
      qporous_gamma(:)=0.0d0 
!
      gamma_wall_save(:)=gamma_wall(:)
      gamma_wall(:)     =0.0d0 
      qconden(:)        =0.0d0           
!
!.....Calculate wall condensation
!
      IF(wconden.eq.1) THEN                  ! Resolved Boundary Layer Approach
         CALL wall_condensation_multi_diff
         CALL wall_condensation_RBLA    
      ELSEIF(wconden.eq.2) THEN              ! Heat and Mass Transfer Analogy : Wall Function Approach
         CALL wall_condensation_multi_diff
         IF(iturb.le.-1)STOP '### HMTA condensation model sould be used with a proper turbulence model (iturb>=1) ###'
         CALL wall_condensation_HMTA
      ELSEIF(wconden.eq.3) THEN              ! Uchida correlation
         CALL wall_condensation_uchida
      ELSEIF(wconden.eq.4) THEN
         CALL wall_condensation_steam
      ELSEIF(wconden.eq.-1) THEN
         CALL wall_condensation_multi_diff
         CALL wall_condensation_RBLA_porous
      ELSEIF(wconden.eq.-2) THEN
         CALL wall_condensation_multi_diff
         CALL wall_condensation_HMTA_porous
      ELSEIF(wconden.eq.-3) THEN
         CALL wall_condensation_uchida_porous
      ENDIF
!
      qporous_liq(:)  =qporous_liq(:)+qporous_liq_save(:) 
      qporous_gas(:)  =qporous_gas(:)+qporous_gas_save(:)
      qporous_gamma(:)=qporous_gamma(:)+qporous_gamma_save(:)
!
      gamma_wall(:)   =gamma_wall(:)+gamma_wall_save(:)
!      
      END SUBROUTINE wall_condensation_model
