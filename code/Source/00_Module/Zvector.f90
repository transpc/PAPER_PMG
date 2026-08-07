      MODULE Zvector
! 
      IMPLICIT NONE
      SAVE
!
      REAL(8),Allocatable:: vrel_o(:),ul_o(:),ug_o(:) 
      REAL(8),Allocatable::vl_n(:,:),vd_n(:,:),vg_n(:,:),vl_o(:,:),vd_o(:,:),vg_o(:,:), &
                           vl_t(:,:),vd_t(:,:),vg_t(:,:)
      REAL(8),Allocatable::vl_f_non(:,:),vg_f_non(:,:)
      REAL(8),ALLOCATABLE::face_fr_l(:,:),face_fr_g(:,:)
!
      END MODULE Zvector
