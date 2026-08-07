      MODULE Zrv_ncell
      
      IMPLICIT NONE
      SAVE
!      
!
      INTEGER ncell_fuel_rod,ncell_fluid_core,num_ch,nz_fine,num_ch_1d,ncell_fuel_rod_all,ncell_fluid_core_all
      INTEGER,ALLOCATABLE::neigh_fuel_rod(:,:),cupid_cell_hts2d(:),n_channel_fluid(:),nz_fluid(:)
      INTEGER,ALLOCATABLE::cupid_cell_channel(:),channel_cell_hts2d(:)
      INTEGER,ALLOCATABLE::nz0_fuel_rod(:),nz_fuel_rod(:),nrod_fuel_rod(:)
      REAL(8),ALLOCATABLE::num_fuel_rod(:)
!
      !Additional connectivity - jrlee
      !MASTER mapping
      INTEGER,ALLOCATABLE::assem_nx0(:),assem_ny0(:),assem_nz0(:),assem_nxy0(:)
      INTEGER,ALLOCATABLE::assem_nx1(:),assem_ny1(:),assem_nz1(:),assem_nxy1(:)
      INTEGER,ALLOCATABLE::assem_nx2(:),assem_ny2(:),assem_nz2(:),assem_nxy2(:)

      INTEGER,ALLOCATABLE::master_to_rod(:,:,:,:),nf_input_cell(:)
      REAL(8),ALLOCATABLE::p3d_cupid(:),qvol_mas(:)
      REAL(8),ALLOCATABLE::dnbr_ce1(:),dnbr_cupid1(:)
      REAL(8),ALLOCATABLE::dnbr_ce2(:),dnbr_cupid2(:)
      INTEGER,ALLOCATABLE::master_to_assem0(:,:),master_to_assem1(:,:),master_to_assem1_rod(:,:),master_to_assem1_cell(:,:)

      !MASTER mapping - yhy modification
      INTEGER,ALLOCATABLE::asm_nx(:),asm_ny(:),asm_nz(:),asm_ni(:),asm_ni2(:)
      INTEGER,ALLOCATABLE::chn_nx(:),chn_ny(:),chn_nz(:)
      

      ENDMODULE Zrv_ncell
