      MODULE Zporous
! 
      IMPLICIT NONE
      SAVE
!
      INTEGER udfi_subchannel_flowdir
      INTEGER vfporous !,diffporous         !LSJ 161121 porous (modified to iavgtype=2)
      !INTEGER iter2max
      !INTEGER ipor_homogen,nporous,lporous       !somaporous.in (1=homogeneous poro, 2=inhomogeneous poro)
      !INTEGER iper_homogen,npermea,lpermea       !somaporous.in (1=homogeneous perm, 2=inhomogeneous poro)
      !INTEGER nhyd,lhyd                          !somaporous.in (hydraulic diameter)
      !INTEGER nsubf,lsubf                        !somaporous.in (sl)
      !INTEGER iso_fric,imp_fric,axial_fric,cross_fric     !somaporous.in (1=isotropic fric, 2=anisotropic fric)  
      !INTEGER nspacer,lspacer
      !INTEGER iturbmix_onoff,iturbmix                     !somaporous.in (0=not use, 1=use)
      !INTEGER, ALLOCATABLE:: indxporous(:),indxpermea(:),indxhyd(:),indxsubf(:),indxspacer(:)
!
      REAL(8) kloss_simple,kloss_cross,kloss_grid
      REAL(8), ALLOCATABLE:: fric_model_gas(:,:),fric_model_liq(:,:),fric_model_drp(:,:)
      REAL(8), ALLOCATABLE:: sgap(:,:)
      REAL(8), ALLOCATABLE:: s_ij_non_i(:),s_ij_non_k(:)
      REAL(8), ALLOCATABLE:: s_gapij_non_i(:),s_gapij_non_k(:)
      REAL(8), ALLOCATABLE:: poro(:),xmin_por(:,:),xmax_por(:,:)
      REAL(8), ALLOCATABLE:: permea(:,:),xmin_per(:,:),xmax_per(:,:)
      REAL(8), ALLOCATABLE:: hyd(:),xmin_hyd(:,:),xmax_hyd(:,:)
      REAL(8), ALLOCATABLE:: subf(:,:),xmin_subf(:,:),xmax_subf(:,:)
      REAL(8), ALLOCATABLE:: kspacer(:),xmin_spacer(:,:),xmax_spacer(:,:)
      
!
!     EVVD for scalar_matrix.f90
      REAL(8) ftm,ka,beta
      REAL(8), ALLOCATABLE:: tm_mas_l(:),tm_mas_g(:),vd_mas_l(:),vd_mas_g(:)
      REAL(8), ALLOCATABLE:: tm_eng_l(:),tm_eng_g(:),vd_eng_l(:),vd_eng_g(:)
      REAL(8), ALLOCATABLE:: cell_area(:)
!     Mixing vane in scalar convection
!     1: Momentum(not used), 2: Mass, 3: Energy
      REAL(8), ALLOCATABLE:: mixing_vane_l(:,:)       
!
!.....Subchannel
      LOGICAL l_subchannel                       
      LOGICAL l_mixing_vane    
      LOGICAL l_spacer_grid      
      LOGICAL l_2p_multiplier_evvd 
      CHARACTER(30) s_subchannel_fric
      CHARACTER(30) s_subchannel_mixing     
      CHARACTER(30) s_subchannel_fric_axial 
      CHARACTER(30) s_subchannel_fric_cross 
      CHARACTER(30) s_2p_multiplier       
!      
      INTEGER,Allocatable:: chn_type_tmp(:),chn_type(:)
      INTEGER,Allocatable:: mv_fac(:),mv_non(:)
      REAL(8),Allocatable:: q_CHF_Biasi(:),APR_DNBR(:),q_local(:),qvol_liq_DNBR(:,:),qvol_DNBR(:),qvol_DNBR_dummy(:),rod_power(:,:)

      INTEGER n_sg
      REAL(8),Allocatable:: h_sg(:),sg_loc(:),mv_loc(:)

      INTEGER n_mv
      REAL(8) h_mv

      INTEGER nz_nk,nz_th0
      REAL(8),ALLOCATABLE:: dz_nk(:),hz_nk(:)
      REAL(8) qsum00
      
      END MODULE Zporous
