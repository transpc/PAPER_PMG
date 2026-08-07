      MODULE Zrv_hts_2d
! 
      IMPLICIT NONE
!
      LOGICAL :: l_ht_str_2d_qcell      
!
      REAL(8) :: power_2d
!
!............................................................................
!
      INTEGER :: nrod_2d,nz0_2d,nr_2d,ht_geo_2d
      INTEGER :: nqvol
      INTEGER,DIMENSION(:),ALLOCATABLE :: wet_bi,wet_ti,n_ch_frap,n_ch_frap_input
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: nmat_2d
!     
      REAL(8),DIMENSION(:),ALLOCATABLE :: ri_2d,dxl_2d,dxr_2d,z_fuel,ht_area_fuel,                     &
                                          qvol_time,                                                   &
                                          hlr_2f,hgr_2f,hstr_2f,hspr_2f,tlr_2f,tgr_2f,tstr_2f,tspr_2f, &
                                          wet_b,wet_t,ztop_q,                                          &
                                          rweight_power
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: t_fuel,qvol_2f0,qvol_2f,qvol_norm_2d
!
      !Modfied 
      !arw_2f/dz_fuel/dz_half_fuel
      REAL(8)::height_fuel,qvol_tot
      REAL(8),DIMENSION(:),ALLOCATABLE :: arw_2f,                        &
                                          dz_fuel,dz_fuel0,dz_half_fuel, &
                                          ht_area_fuel0
      REAL(8),DIMENSION(:,:),ALLOCATABLE:: azl_2f,arl_2f,vl_2f, &
                                           azr_2f,arr_2f,vr_2f, &
                                           v_2f
      
      ! OPR1000
      REAL(8),DIMENSION(:),ALLOCATABLE :: dz_fuel00
!
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: qcell_to_qvol_mul
!
!.....rod bundle having different radial dimension 
      INTEGER::ri_2d_opt
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: ri_2d_input,rod_dxr_2d,rod_alw_2f,rod_arw_2f
      REAL(8),DIMENSION(:,:,:),ALLOCATABLE :: rod_azr_2,rod_arr_2f,rod_azr_2f,rod_vr_2f,rod_azl_2f,rod_arl_2f,rod_vl_2f,rod_v_2f
!
      REAL(8),DIMENSION(:),ALLOCATABLE :: qvol_tot_channel, &  !individual power for each rod      
                                          twall_fuel

      END MODULE Zrv_hts_2d
