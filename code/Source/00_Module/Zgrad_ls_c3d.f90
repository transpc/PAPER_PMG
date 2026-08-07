      MODULE Zgrad_ls_c3d
!      
      IMPLICIT NONE
      SAVE
!
      INTEGER,Allocatable::lsindex(:)
      REAL(8),Allocatable::a11_3(:),a12_3(:),a13_3(:),a22_3(:),a23_3(:),a33_3(:),a_3(:),det_3(:) 
!
      END MODULE Zgrad_ls_c3d