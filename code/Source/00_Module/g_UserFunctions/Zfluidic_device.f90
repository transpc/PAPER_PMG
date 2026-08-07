      Module Zfluidic_device
! 
      IMPLICIT NONE
      SAVE
!      
      REAL(8) flow_si,k_sp,k_dp
!   bug flow_sp,flow_fd  never computer creates a problem in udfn_mom_wall.f90
!   do not know if zero value is correct ???
      REAL(8) :: flow_sp=0.d0
      REAL(8) :: flow_fd=0.d0
!
      ENDMODULE Zfluidic_device