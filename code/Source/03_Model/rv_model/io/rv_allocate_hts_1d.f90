!
      SUBROUTINE rv_allocate_hts_1d
!
      USE Zrv_hts_1d,    ONLY:ncell_hts_1d,nr_1d,vl_1d,vr_1d,sl_1d,sr_1d,                    &
                               hfluxl_1d,hfluxr_1d,twl_1d,twr_1d,ng_hts_1d,t_hts_1d,         &
                               cupid_cell_1d,ig_hts_1d,ht_area_left_1d,ht_area_right_1d,     &
                               hll_1d,hstl_1d,hspl_1d,hgl_1d,tll_1d,tstl_1d,tspl_1d,tgl_1d,  &
                               hlr_1d,hstr_1d,hspr_1d,hgr_1d,tlr_1d,tstr_1d,tspr_1d,tgr_1d
      USE Zwall_HTC,     ONLY:mode_1d
!
      IMPLICIT NONE
!
      INTEGER n,nr,ng
!
      nr=nr_1d
      ng=ng_hts_1d
!
      ALLOCATE(vl_1d(ng,nr),vr_1d(ng,nr),sl_1d(ng,nr),sr_1d(ng,nr))
      ALLOCATE(hfluxl_1d(ng),hfluxr_1d(ng),twl_1d(ng),twr_1d(ng))
!
      vl_1d(:,:)=0.0d0
      vr_1d(:,:)=0.0d0
      sl_1d(:,:)=0.0d0
      sr_1d(:,:)=0.0d0
      hfluxl_1d(:)=0.0d0
      hfluxr_1d(:)=0.0d0
      twl_1d(:)=0.0d0
      twr_1d(:)=0.0d0
!
      n=ncell_hts_1d
      ALLOCATE(t_hts_1d(n,nr),cupid_cell_1d(n),ig_hts_1d(n),ht_area_left_1d(n),ht_area_right_1d(n))

      t_hts_1d(:,:)=0
      cupid_cell_1d(:)=0
      ig_hts_1d(:)=0
      ht_area_right_1d(:)=0.0d0
      ht_area_right_1d(:)=0.0d0

!
      ALLOCATE(hll_1d(n),hstl_1d(n),hspl_1d(n),hgl_1d(n),tll_1d(n),tstl_1d(n),tspl_1d(n),tgl_1d(n))
      hll_1d(:)=0.0d0
      hstl_1d(:)=0.0d0
      hspl_1d(:)=0.0d0
      hgl_1d(:)=0.0d0
      tll_1d(:)=0.0d0
      tstl_1d(:)=0.0d0
      tspl_1d(:)=0.0d0
      tgl_1d(:)=0.0d0
!
      ALLOCATE(hlr_1d(n),hstr_1d(n),hspr_1d(n),hgr_1d(n),tlr_1d(n),tstr_1d(n),tspr_1d(n),tgr_1d(n))
      hstr_1d(:)=0.0d0
      hspr_1d(:)=0.0d0 
      hgr_1d(:)=0.0d0
      tlr_1d(:)=0.0d0
      tstr_1d(:)=0.0d0
      tspr_1d(:)=0.0d0  
      tgr_1d(:)=0.0d0        
      hlr_1d(:)=0.0d0
!
!.....Zwall_HEAT
!
      ALLOCATE(mode_1d(n))
      mode_1d(:)=0
!
      RETURN
      END SUBROUTINE rv_allocate_hts_1d
