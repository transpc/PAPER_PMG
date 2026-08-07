      MODULE Zncg
! 
      IMPLICIT NONE
      SAVE
!
      INTEGER n_ncg_sp,imp_ncg,ncg_diff,i_ncg_vis
      INTEGER,ALLOCATABLE::ncg_species(:)
      REAL(8) tao
      REAL(8),ALLOCATABLE::cvao_cell(:),uao_cell(:),dcva_cell(:),ra_cell(:),qn_cell(:,:),qn_cell_o(:,:),qn_cell0(:),wmole_gas(:),advn_cell(:)
      REAL(8),ALLOCATABLE::cvao_nvin(:),uao_nvin(:),dcva_nvin(:),ra_nvin(:),qn_nvin(:,:)
      REAL(8),ALLOCATABLE::cvao_npin(:),uao_npin(:),dcva_npin(:),ra_npin(:),qn_npin(:,:)
!
      END MODULE Zncg
