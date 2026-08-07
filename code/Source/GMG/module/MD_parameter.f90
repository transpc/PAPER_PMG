! ====================================================================!
!  module for CSR 
!---------------------------------------------------------------------+
      module MD_parameter
      implicit none
! ---
      integer :: ndim,nvert_max,nface_max,nvert,nface,maxit,mdf_matrix,isol,     &
                 nsemi,ICG,nv_max
      INTEGER(4) :: ip_inter,nvpe
      integer :: nf_max,ndom,ipar,iVcy,isemi
      real(8):: crit,teta,teta_p,alpha
      INTEGER(4) :: ele_type
      
!---
      save
    end module