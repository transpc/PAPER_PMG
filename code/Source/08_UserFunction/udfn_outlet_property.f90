!
      SUBROUTINE udfn_outlet_property
!   
!     User-defined outlet properties
!
      USE Zconst1 , ONLY: vv_prob
!
!.....cavitation
!      
      IF(vv_prob.eq.'cavitation') CALL cavitation_outlet_property_user
!
!.....plume
!            
      IF(vv_prob.eq.'plume') CALL plum_outlet_property_user
!
!.....rocom
!      
      IF(vv_prob.eq.'rocom')THEN
         CALL rocom_outlet_bc_user
         CALL rocom_control_bc_user
         CALL rocom_periodic_threee_inout_user
         CALL outlet_horizontal_user
      ENDIF
!
!.....atlas
!      
      IF(vv_prob.eq.'atlas')CALL outlet_horizontal_user
!
!.....sgp_separator
!
      IF(vv_prob.eq.'sgp_separator') CALL udfn_sg_outlet_property
!
!.....siphon
!
      IF(vv_prob.eq.'siphon')CALL siphon_control_pbnd
!      
      RETURN 
      END SUBROUTINE udfn_outlet_property
