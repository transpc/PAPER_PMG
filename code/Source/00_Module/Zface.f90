      MODULE Zface
! 
      IMPLICIT NONE
      SAVE
!
      INTEGER,PARAMETER :: Free_slip=-2,Laminar=-1, Zequation=0, SST=1, Kepsilon=2, Kepsilon_RNG=3, Kepsilon_real=4
      INTEGER,PARAMETER :: LES_WALE=100
      INTEGER,PARAMETER :: gas_only=1, liq_only=2, both=3
      INTEGER Twall_Model
      REAL(8),Allocatable::q1cell(:),qqcell(:),qecell(:),qclcell(:),qcgcell(:),ndensitycell(:)
!
      END MODULE Zface
