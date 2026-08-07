
      MODULE Zvec_turb
! 
      IMPLICIT NONE
      SAVE
!  
      REAL(8),DIMENSION(:),ALLOCATABLE :: ke_non_i,ke_non_k,dp_non_i,dp_non_k
      REAL(8),DIMENSION(:),ALLOCATABLE :: ke_inl_i,ke_inl_k,dp_inl_i,dp_inl_k
      REAL(8),DIMENSION(:),ALLOCATABLE :: ke_out_i,ke_out_k,dp_out_i,dp_out_k
      REAL(8),DIMENSION(:),ALLOCATABLE :: ke_mcc_i,ke_mcc_k,dp_mcc_i,dp_mcc_k
!
      END MODULE Zvec_turb
