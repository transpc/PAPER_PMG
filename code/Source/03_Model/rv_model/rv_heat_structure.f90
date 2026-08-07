!
      SUBROUTINE rv_heat_structure
!
      USE Zrv_hts_1d,    ONLY:ncell_hts_1d
      USE Zrv_hts_2d,    ONLY:nrod_2d
      USE Zrv_ncell,     ONLY:ncell_fuel_rod
!
      IMPLICIT NONE
!
!     LOGICAL,SAVE :: initial=.true.
!
!     IF(initial)THEN
!
!........Calculate geometric parameters
!
!         IF(ncell_hts_1d.gt.0) CALL rv_geo_hts_1d
!         IF(nrod_2d.gt.0) CALL rv_geo_hts_2d
!
!........rezone only if relood option is on.
!
!         IF(reflood.ge.1)CALL rv_rezone
!         write(*,*)reflood
!
!         initial=.false.
!     ENDIF
!
!.....Calculate heat conduction in the heat structure
!
      IF(ncell_hts_1d.gt.0) CALL rv_hts_1d
!      IF(nqvol.gt.0) CALL rv_power_transient 
      IF(nrod_2d.gt.0) THEN
         IF(ncell_fuel_rod.gt.0) CALL rv_hts_2d
      ENDIF
!
      END SUBROUTINE rv_heat_structure
