      MODULE Zpress
! 
      IMPLICIT NONE
      SAVE
!
      INTEGER flag
      REAL(8),Allocatable::p(:),dpdx(:,:),pp(:),pp_o(:)
      REAL(8),Allocatable::dpdx_o(:,:) 
!
      END MODULE Zpress
