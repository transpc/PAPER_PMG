      MODULE Zrv_hts_1d
! 
      IMPLICIT NONE
      SAVE
!
      INTEGER :: ncell_hts_1d,ncell_hts_1d_all,nr_1d,ng_hts_1d,imp_cond_1d
      INTEGER,ALLOCATABLE :: n_ch_1d(:),ht_geo_1d(:),nmat_1d(:),bcl_1d(:),bcr_1d(:),ig_hts_1d(:)
      INTEGER,ALLOCATABLE :: cupid_cell_1d(:)
      REAL(8) :: dz_1d
      REAL(8),ALLOCATABLE :: ri_1d(:,:),ht_area_left_1d(:),ht_area_right_1d(:)
      REAL(8),ALLOCATABLE :: t_hts_1d(:,:),vl_1d(:,:),vr_1d(:,:),sl_1d(:,:),sr_1d(:,:),slw_1d(:),srw_1d(:)
      REAL(8),ALLOCATABLE :: hfluxl_1d(:),hfluxr_1d(:),twl_1d(:),twr_1d(:),                   &
                             hll_1d(:),hlr_1d(:),hstl_1d(:),hstr_1d(:),hspl_1d(:),hspr_1d(:),hgl_1d(:),hgr_1d(:), &
                             tll_1d(:),tlr_1d(:),tstl_1d(:),tstr_1d(:),tspl_1d(:),tspr_1d(:),tgl_1d(:),tgr_1d(:)
      REAL(8),ALLOCATABLE :: num_rod_1d(:)                              
!
      END MODULE Zrv_hts_1d   
