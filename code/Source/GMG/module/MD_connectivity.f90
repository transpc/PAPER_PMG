 
!---------------------------------------------------------------------+
      module MD_connectivity
      implicit none
!---
      INTEGER(4)  nnz_neigh, nnz_neighc, nnz_tmp
      INTEGER(4),dimension(:),allocatable::live
      INTEGER(4),DIMENSION(:),ALLOCATABLE:: ia_neigh, ja_neigh
      INTEGER(4),DIMENSION(:),ALLOCATABLE::ia_neighc, ja_neighc     
      INTEGER(4),DIMENSION(:),ALLOCATABLE::ia_tmp, ja_tmp       

!---
      save
      end module
