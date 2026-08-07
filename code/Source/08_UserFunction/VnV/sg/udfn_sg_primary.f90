!
      SUBROUTINE udfn_sg_primary 
!
      IMPLICIT NONE
!
!.....Calculate flow velocity inside utube
!
      CALL udfn_1d_momentum
!
!.....Calculate the primary and secodary side heat transfer coefficients
!
      CALL udfn_sg_htc
!
!.....Calculate heat conduction of utube
!
      CALL udfn_heat_tube
!
!.....Calculate energy of the primary coolant inside utube
!
      CALL udfn_1d_energy
!
!.....Update the primary coolant properties
!
      CALL udfn_prop_liquid
!
      RETURN
      END SUBROUTINE udfn_sg_primary
!