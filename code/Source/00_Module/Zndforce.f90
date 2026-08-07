      MODULE Zndforce
! 
      IMPLICIT NONE
      SAVE
!
      INTEGER,Allocatable::c_bface(:),c_bface_indx(:),nwb(:)
      INTEGER,Allocatable::face_wall_group(:,:),cell_closewall(:,:),cell_closewall_indx(:,:)
      REAL(8) relax_hik,relax_cd
      REAL(8),Allocatable::dvdxl(:,:,:)
      REAL(8),Allocatable::d_bfc(:),F_wl(:,:)
      REAL(8),Allocatable::Cwlf(:),Clift(:),Ctd(:),Vfgl_o(:)
      REAL(8),Allocatable::dis_closewall(:,:)      
!
      CHARACTER(30) s_turb_disp
      CHARACTER(30) s_wall_lub
!
      END MODULE Zndforce
