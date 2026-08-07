!
      SUBROUTINE initialize_topology_criteria
!
!     Define topology map criteria for alpha, gamma
!      
      USE Zconst1     , ONLY: vv_prob
      USE Zflowregime , ONLY: alphag_bc,alphag_cm,gamma_1,gamma_2
!
      IMPLICIT NONE 
!
      alphag_bc= 0.7d0
      alphag_cm= 0.90d0
      gamma_1=   0.2d0
      gamma_2=   0.45d0
!      
      IF(vv_prob.eq.'PAFS-POOL'.or.vv_prob.eq.'fluidic_device'.or.  &
         vv_prob.eq.'check_Hik'.or.vv_prob.eq.'DIVA-NEW'.or.     &
         vv_prob.eq.'UPTF-RV'.or.vv_prob.eq.'Nuscale-RVV'.or.vv_prob.eq.'Nuscale-PZR') THEN
         alphag_bc=0.3d0
      ENDIF
!      
      RETURN
      END SUBROUTINE initialize_topology_criteria