!
      SUBROUTINE fluxBC_main
!
      USE Zrv_model    , ONLY: rv_mcp,rv_choke
!
      IMPLICIT NONE
!
!......choke model
!
      IF(rv_choke.eq.1) CALL choke_main   
!
!......mcp model
!      
      IF(rv_mcp.eq.1) CALL mcp_main
!
!......Update u* with Neumann condition
!      
      CALL barrier_mpi
!      
      CALL fluxBC_ustar   
!                                                                       
      RETURN
      END SUBROUTINE fluxBC_main