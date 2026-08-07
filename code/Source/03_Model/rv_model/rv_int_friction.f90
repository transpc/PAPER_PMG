!
      SUBROUTINE rv_int_friction
!
      USE VOL_DATA 
      USE Zzone         , ONLY: ncell_fluid
!
      IMPLICIT NONE
!      
      REAL(8),ALLOCATABLE:: ia_sum(:)      
!
      ALLOCATE(ia_sum(ncell_fluid))
      ia_sum=0.0d0
!      
      CALL rv_flow_regime       
      CALL rv_int_fric_model(ncell_fluid,cell%vfgl,cell%vfgd,ia_sum)
!      
      DEALLOCATE(ia_sum)
!
      RETURN
      END SUBROUTINE rv_int_friction
