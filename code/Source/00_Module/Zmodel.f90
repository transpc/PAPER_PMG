      MODULE Zmodel
! 
      IMPLICIT NONE
      SAVE
!
      REAL(8) Drag_coeff_a,Drag_coeff_b,Drag_coeff_Cb
      REAL(8) H_il_coeff_a,H_il_coeff_b,H_il_min
      REAL(8) H_ig_coeff_a,H_ig_coeff_b,H_ig_min
      REAL(8) H_fg_coeff_a,H_fg_coeff_b,H_fg_min
      REAL(8) Cd_min_user,Cd_min_ag99
      REAL(8) H_il_min_user,H_il_min_ag99
      REAL(8) dtl,dtg
      REAL(8),ALLOCATABLE:: drift_c0(:),drift_c1(:),cb_bubble(:)
      REAL(8),ALLOCATABLE:: h_il_cfd(:),h_ig_cfd(:),h_gf_cfd(:),vfgl_cfd(:),vfgd_cfd(:)
!.....radiation      
      INTEGER rad_model
      INTEGER max_iter_rad                                
      REAL(8) abs_coeff
      REAL(8) wall_emiss
      REAL(8) eps_imp_rad   
      REAL(8),ALLOCATABLE :: rad_source(:)
!.....condensation                         
      INTEGER,ALLOCATABLE :: wVertical(:)                                                
      REAL(8),ALLOCATABLE :: coef_diff(:)       
      REAL(8),ALLOCATABLE :: qconden(:) ,qrad(:)
!
!.....wall frction
      REAL(8),ALLOCATABLE :: vfwl_k(:)
      REAL(8),ALLOCATABLE :: fsar(:,:)
!      
!.....SMR      
      REAL(8),ALLOCATABLE :: resist(:)       

!.....Hymeres2
      REAL(8),ALLOCATABLE :: molefr(:),cube_tw(:)
      REAL(8),ALLOCATABLE :: h2p1_tw(:)
      REAL(8) h2p1_tp1
      REAL(8) h2p1_tp2
!.....input
      INTEGER i_weight                 !volume faction weighted drag and ihtc
      INTEGER i_droplet                !drop field
      CHARACTER(30) s_wall_fric        !wall friction model
      INTEGER i_fs_temp_intpol         !interpolate fluid-solid inteface temperature
      INTEGER use_porous               !index to check whether the porous volumes are included
!             
      END MODULE Zmodel
