      MODULE Zrv_mpi
      
      IMPLICIT NONE
      SAVE
!      
!.....Serial decomposition
!
      INTEGER niut_fuel_rod,ncell_fuel_rod_p
      INTEGER,DIMENSION(:),ALLOCATABLE :: celem_fuel_rod,celem_fluid_core,celem_hts_1d
      INTEGER,DIMENSION(:),ALLOCATABLE :: ri_fuel_rod,si_fuel_rod,iut_fuel_rod,rintf_fuel_rod,sintf_fuel_rod,jperm_fuel_rod, &
                                          jperm_fluid_core,jperm_hts_1d
!
!.....For ALLGATHERV
!
      INTEGER,DIMENSION(:),ALLOCATABLE :: ncell_fluid1_core,ncell_fluid1_core_dsp,jjperm_fluid_core
      INTEGER,DIMENSION(:),ALLOCATABLE :: ncell_fuel1_rod,ncell_fuel1_rod_dsp,jjperm_fuel_rod
      INTEGER,DIMENSION(:),ALLOCATABLE :: ncell_fuel1_rod_2d,ncell_fuel1_rod_2d_dsp
!
      END MODULE Zrv_mpi
