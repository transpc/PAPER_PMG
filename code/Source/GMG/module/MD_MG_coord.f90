      MODULE MD_MG_coord
! 
      IMPLICIT NONE
      SAVE
!
      INTEGER(4):: nmax1,nmax2
      INTEGER(4) :: nnode1,nnode2
      INTEGER(4) :: nelem1,nelem2
      INTEGER(4) :: nnode1gl,nnode2gl
      INTEGER(4),DIMENSION (:), ALLOCATABLE :: nnei1
      INTEGER(4),DIMENSION (:,:), ALLOCATABLE :: inei1
      INTEGER(4),DIMENSION (:), ALLOCATABLE :: icoarse,icoarse1
!      INTEGER(4),DIMENSION (:), ALLOCATABLE :: imapc1,imapc2
!     
      
      REAL(8),DIMENSION(:,:), ALLOCATABLE :: coord1,coord2
      
! NEW
      INTEGER(4),DIMENSION (:), ALLOCATABLE ::  ialv
      INTEGER(4) nnods,ncolc,ncolf
      REAL(8),DIMENSION(:,:), ALLOCATABLE :: coordc
      INTEGER(4),DIMENSION (:), ALLOCATABLE :: imapc,icoarsef,inmax,imapt
!
      END MODULE MD_MG_coord
