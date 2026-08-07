      MODULE Zinlet_bc
! 
      IMPLICIT NONE
      SAVE
!
      INTEGER,PARAMETER::inary=100
      INTEGER,Allocatable::index_trace(:),incell_index(:)
      REAL(8),Allocatable::vgy_in(:),vly_in(:),dsm_in(:),   &
                              alphag_in(:),ia_in(:)
!
      END MODULE Zinlet_bc