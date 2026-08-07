      MODULE Zscalar_coeff
! 
      IMPLICIT NONE
      SAVE
!
      INTEGER l_th_equil,l_min_hik
      REAL(8) alphag_min,alphal_min
      REAL(8),Allocatable :: sb(:,:)
      REAL(8),Allocatable :: sfg_nf(:,:),sfl_nf(:,:),sfd_nf(:,:)
      REAL(8),Allocatable :: sfg_non_k(:,:),sfl_non_k(:,:),sfd_non_k(:,:)
      REAL(8),Allocatable :: sfg6_nf(:),sfl6_nf(:),sfd6_nf(:)
      REAL(8),Allocatable :: sfg6_non_k(:),sfl6_non_k(:),sfd6_non_k(:)
!
      END MODULE Zscalar_coeff
!
      MODULE Zcheck_scalar
! 
      IMPLICIT NONE
      SAVE
!
      REAL(8) :: eps_rho=8.d-3
      REAL(8) :: eps_p  =30.d3
      REAL(8) :: eps_eng=1.d-1
      REAL(8) :: eps_vol=1.d-6
!
      END MODULE Zcheck_scalar
