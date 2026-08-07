module KSMR
!   
   INTEGER, PARAMETER :: n_zone=12
   INTEGER, ALLOCATABLE :: zone_comp(:),zone_comp_all(:)
!   REAL(8), ALLOCATABLE :: resist(:)
   REAL(8):: zone_comp_area(n_zone), zone_comp_height(n_zone)
   REAL(8), SAVE :: sg_dp,core_dp,dc_dp
   REAL(8), SAVE :: vol_sg
   REAL(8), PARAMETER :: T_in =530.5d0 !538.15 !568.5d0 ! 560.d0 !538.15
   REAL(8), PARAMETER :: T_out=560.5d0 !594.15 !568.5d0 ! 560.d0 !594.15
   REAL(8), PARAMETER :: core_power = 259.d6
   REAL(8), PARAMETER :: time_power_ramping = 1.d0 !100.d0 !100.d0
   REAL(8), PARAMETER :: t_trip=2000.d0
   REAL(8):: power_ans73
!
end module KSMR
