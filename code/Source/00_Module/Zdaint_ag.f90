      MODULE Zdaint_ag
! 
      IMPLICIT NONE
      SAVE
!
      REAL(8),Allocatable::daint1_ag(:), daint2_ag(:),daint1_ag_bc(:),daint2_ag_cm(:)
      REAL(8),Allocatable::aint_01b(:),aint_bc(:),aint_09b(:)
      REAL(8),Allocatable::aint_01d(:),aint_cm(:),aint_09d(:)
      REAL(8),Allocatable::D1_01(:),D1_bc(:),D1_09(:),D2_cm(:),D2_09(:)
!
      END MODULE Zdaint_ag