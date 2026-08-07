MODULE Zrv_htc
!
   IMPLICIT NONE
   SAVE
!
   REAL, PARAMETER :: hig_big = 1.0D4   !< big htc for SCG, W/m^2K
   REAL, PARAMETER :: hil_big = 3.0D6   !< big htc for SHL TB, W/m^2K
   REAL, PARAMETER :: HTC_MIN = 1.0D-30 !< Min. htc allowed, W/m^2K
!  
END MODULE Zrv_htc
