! ====================================================================!
!  module for CSR 
!---------------------------------------------------------------------+
      module MD_matrix
      implicit none
! ---
      integer nnz
      integer,dimension(:),allocatable :: ia,ja,ju
      real*8,dimension(:),allocatable :: au,alu, diagr
      
      real*8,dimension(:),allocatable:: u, b,ut

!... For PBCG...
      INTEGER :: ngroup
      INTEGER, allocatable :: ia_a1(:), ju_a1(:), ja_a1(:), iend1(:), indx_a1(:)
      INTEGER, allocatable :: ia_a2(:), ju_a2(:), ja_a2(:), iend2(:), indx_a2(:)
      REAL(8), ALLOCATABLE :: alu0(:), alu1(:)
      INTEGER, ALLOCATABLE :: iap(:,:), iaa(:), jaa(:), nbgroup(:,:), jap(:)
      
! for pbcg-ali
      INTEGER ::  nnz_l,nnz_u
      INTEGER, allocatable :: ia_l(:),ia_u(:), ja_l(:), ja_u(:)
      INTEGER, allocatable :: iend(:)
      REAL(8), ALLOCATABLE :: diag_l(:), diag_u(:)
      REAL(8), ALLOCATABLE :: alu_l(:),alu_u(:)
      
      
!---
      save
    end module
     
