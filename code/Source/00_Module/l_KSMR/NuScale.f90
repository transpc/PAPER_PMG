module NuScale
!   
   INTEGER, PARAMETER :: n_zone=10
   INTEGER, ALLOCATABLE :: zone_comp(:),zone_comp_all(:)
!   REAL(8), ALLOCATABLE :: resist(:)
   REAL(8):: zone_comp_area(n_zone), zone_comp_height(n_zone)
   REAL(8), PARAMETER :: T_in =568.d0 !.15d0 !568.d0 !538.15d0
   REAL(8), PARAMETER :: T_out=568.d0 !594.15d0 !568.d0 !594.15d0
   REAL(8), PARAMETER :: core_power =192.d6
   REAL(8), PARAMETER :: time_power_ramping = 300.d0
   REAL(8), PARAMETER :: time_power_off = 300.d100
   REAL(8), PARAMETER :: t_trip=300.d0
   REAL(8):: power_ans73
!
end module NuScale
