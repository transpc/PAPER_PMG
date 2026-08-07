!
      SUBROUTINE set_vol_heat_source
!
!     This routine allocate volumetric heat source
!
      USE Zconst2      , ONLY: stime_hup,stime_hflat
      USE Zqvol        , ONLY: q0_liq,q0_gas,q0_ice_solid,qvol_gas,qvol_liq,qvol_ice_solid
      USE Ztimecon     , ONLY: time
      USE Zuserdefined , ONLY: udfl_set_qvol_porous
      USE Zzone        , ONLY: nzone,ncell_fluid,ncell_cond,nzone_c
!
      IMPLICIT NONE
!      
      INTEGER(4) i
!
!.....User defined volumetric heat source
!
      IF(udfl_set_qvol_porous)THEN
         CALL udfn_set_qvol_porous
         RETURN
      ENDIF
!
!.....Volumetric heat source for the solid zone
!
      DO i=1,ncell_cond
         qvol_ice_solid(i)=0.d0
         IF(time.ge.stime_hup.and.time.lt.stime_hflat)THEN
            qvol_ice_solid(i)=q0_ice_solid(nzone_c(i))*(time-stime_hup)/(stime_hflat-stime_hup)
         ELSEIF(time.lt.stime_hup)THEN
            qvol_ice_solid(i)=0.0d0
         ELSE
            qvol_ice_solid(i)=q0_ice_solid(nzone_c(i))
         ENDIF
      ENDDO      
!
!.....Volumetric heat source for the fluid zone
!
      DO i=1,ncell_fluid
         qvol_gas(i)=q0_gas(nzone(i))
         qvol_liq(i)=q0_liq(nzone(i))
      ENDDO      
!
      RETURN
      END SUBROUTINE set_vol_heat_source