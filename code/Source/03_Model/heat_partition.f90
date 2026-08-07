!
      SUBROUTINE heat_partition
!
!     This routine calculates diffusive heat fluxes between cell and cell face 
!     whose property is solid-fluid interface, constant temperature, constant 
!     heatflux(nbcon=-2,-3&-4,-5). Do not calculate heat partition for 
!     adiabatic wall and fluid cell w/o wall.
!
      USE VOL_DATA                 
      USE Zconst1       , ONLY: iheatpart,rv_htmodel_forCFD
      USE Zmodel        , ONLY: use_porous
      USE Zuserdefined  , ONLY: udfl_wallHTC_porous
!
      IMPLICIT NONE
!
!
!.....Use of heat flux partitioning model when iheatpart > 0.
!
      IF(iheatpart.ne.0)THEN        
         !CALL heat_partition_fluid   ! Newton-Raphson method
         CALL heat_partition_2_bi_lo  ! Bi-section method
      ENDIF
!
!.....Use of RV heat transfer model package for Fluid-solid interface heat transfer (not RV heat structure)
!     
      IF(rv_htmodel_forCFD.gt.0)THEN
         CALL non_rv_wall_HT_2d
      ENDIF      
!
!.....Heat partition in porous cells, nmaterial(i)<0
!
      IF(udfl_wallHTC_porous)THEN
!
!........Heat partition for liquid only
!           
         CALL udfn_heat_wallHTC_porous
      ELSEIF(use_porous.ne.0)THEN
!
!........Heat partition for wall boiling 
!   
         CALL heat_partition_porous
      ENDIF
!
      END SUBROUTINE heat_partition
