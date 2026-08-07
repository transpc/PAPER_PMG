      Module Zrad_comp !pik-radiation_component
! 
      IMPLICIT NONE
      SAVE
!
!.....Constant
!     
      REAL(8),PARAMETER::sigma_sb=5.670d-8
!
!.....Control
!     
      INTEGER:: rad_comp_mod=0
!
!.....input
!      
      INTEGER :: input_opt      
      INTEGER :: nset,nsize_max
      INTEGER,ALLOCATABLE:: nproperty(:,:)
      INTEGER,ALLOCATABLE:: nsize(:,:)
      INTEGER,ALLOCATABLE:: cell_idx_tmp(:,:,:),cell_idx(:,:,:),core_idx(:,:,:)       
      REAL(8),ALLOCATABLE:: epsil(:,:),rho(:,:),viewf(:,:,:)
      REAL(8),ALLOCATABLE:: area(:,:)
!      
!.....output
!      
      REAL(8),ALLOCATABLE::heat_rad(:,:),heat_rad1(:,:),heat_rad2(:,:)  
      REAL(8),ALLOCATABLE::qrad_flu(:),qrad_rod(:),qrad_sol(:)
!
!.....intermediate
!      
      REAL(8),ALLOCATABLE::t1(:,:),t2(:,:)
      REAL(8),ALLOCATABLE::a(:),b(:),c(:),d(:),dterm(:)
!
!.....Special output
      INTEGER:: ncell_out
      INTEGER,ALLOCATABLE:: cell_out_tmp(:,:),cell_out(:,:),nproperty_out_tmp(:)
!      
      INTEGER:: nloc_out
      REAL(8),ALLOCATABLE:: xloc_out_tmp(:,:),frac_out_tmp(:,:)
!        
      INTEGER,ALLOCATABLE:: core_cell_cupid(:)
!
      ENDMODULE Zrad_comp
