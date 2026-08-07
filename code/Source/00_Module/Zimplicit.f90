      MODULE Zimplicit
! 
      IMPLICIT NONE
      SAVE
!
      INTEGER imp_mom_conv,imp_mom_diff,imp_scalar_conv,imp_scalar_diff,imp_ke_diff,imp_ke_conv
      INTEGER imp_alpha,iter_scalar,skip_imp_scalar,iter_mom
      INTEGER max_iter_mom,max_iter_alpha,max_iter_scalar,max_iter_ke,max_iter_dp
      REAL(8) ag_min_m,al_min_m,eps_imp_mom,eps_imp_alpha,eps_imp_scalar,eps_imp_ke,eps_imp_dp,dp_eps
      INTEGER imp_boron_trans
      INTEGER max_iter_boron
      REAL(8) eps_imp_boron
!
      END MODULE Zimplicit