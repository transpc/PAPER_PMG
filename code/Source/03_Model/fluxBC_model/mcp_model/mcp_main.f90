!
      SUBROUTINE mcp_main
!
      USE Zmcp ,ONLY: time_mcp_off,vflow_direct
      USE Ztimecon ,ONLY: time 
!
      IMPLICIT NONE
!      
      IF(vflow_direct.eq.1.and.time.ge.time_mcp_off) RETURN      
!
!.....Define MCP flow rate
!
      CALL mcp_flow
!
      RETURN
      END SUBROUTINE mcp_main