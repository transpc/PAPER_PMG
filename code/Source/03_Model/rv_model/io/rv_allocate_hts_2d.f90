!
      SUBROUTINE rv_allocate_hts_2d
!
      USE Zrv_hts_2d,     ONLY:nrod_2d,dxl_2d,dxr_2d
!
      USE Zrv_hts_2d,     ONLY:nz0_2d,nr_2d,z_fuel,azl_2f,azr_2f,arl_2f,arr_2f,vl_2f,vr_2f,ht_area_fuel,v_2f
      USE Zrv_hts_2d,     ONLY:arw_2f
      USE Zrv_hts_2d,     ONLY:dz_fuel0
      USE Zrv_ncell,      ONLY:num_ch,nz_fine
!
      IMPLICIT NONE
!
      INTEGER :: nc,nz,nr,nd
!
      nc=num_ch
      nd=nrod_2d
      nz=nz0_2d*nz_fine
      nr=nr_2d
!
      ALLOCATE(z_fuel(nz),dz_fuel0(nz))
      ALLOCATE(dxl_2d(nr),dxr_2d(nr),ht_area_fuel(nc))
      z_fuel(:)=0.d0
      dz_fuel0(:)=0.d0
      dxl_2d(:)=0.d0
      dxr_2d(:)=0.d0
!
      !modified for various height
      ALLOCATE(azl_2f(nz,nr),azr_2f(nz,nr),arl_2f(nz,nr),arr_2f(nz,nr),vl_2f(nz,nr),vr_2f(nz,nr))
      azl_2f=0.d0
      azr_2f=0.d0
      arl_2f=0.d0
      arr_2f=0.d0
      vl_2f=0.d0
      vr_2f=0.d0
      ALLOCATE(v_2f(nz,nr))
      v_2f=0.d0
      
      ALLOCATE(arw_2f(nz))
      arw_2f=0.d0
      
!
      END SUBROUTINE rv_allocate_hts_2d
