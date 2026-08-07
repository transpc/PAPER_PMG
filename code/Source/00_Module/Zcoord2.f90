      MODULE Zcoord2
!      
      IMPLICIT NONE
      SAVE
!
      REAL(8),DIMENSION(:),ALLOCATABLE :: fac,fac1,fac_c,fac1_c
      REAL(8),DIMENSION(:),ALLOCATABLE :: xloc_xfc_radius_min,xloc_xfc_radius_max
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: xfc,cell_leng,xfc_wallcell
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: xfc_min,xfc_max,xloc_xfc_min,xloc_xfc_max
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: prnvar
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: xfc_c
!
      END MODULE Zcoord2
      
