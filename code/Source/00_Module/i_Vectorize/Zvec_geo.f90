      MODULE Zvec_geo
! 
      IMPLICIT NONE
      SAVE
!     geo_var
      REAL(8),DIMENSION(:),ALLOCATABLE :: sap_nf,sa_nf,saa_nf
      REAL(8),DIMENSION(:),ALLOCATABLE :: dji_nf,djia_nf,djir_non
      REAL(8),DIMENSION(:),ALLOCATABLE :: fac_nf,fac1_nf,fac_non,fac_fsw, &
                                          fac1_non,fac1_fsw,f0,f1,        &
                                          perm_non,perm_out,perm_nf,      &
                                          sad_non
      REAL(8),DIMENSION(:),ALLOCATABLE :: djia_non_k,fac1_non_k
!
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: xn_nf,sv_nf,svp_nf,xfc_nf, &
                                            dxfc_nf,dxfc_non_k,dji_x_nf
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: xloc_m_non_i,xloc_m_non_k,dnj_non
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: xn_non_k,dji_x_non_k,sv_non_k
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: xloc_m_non_ki,xloc_m_non_kk
!     Solid
      REAL(8),DIMENSION(:),ALLOCATABLE :: fac_c_nf,fac1_c_nf, &
                                          sap_c_nf,dji_a_c_nf
!
      END MODULE Zvec_geo
