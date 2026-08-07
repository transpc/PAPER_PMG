!
      SUBROUTINE mcp_flow
!
!     Henry-Fauske Critical Flow Model       
!
      USE Ztimecon        , ONLY: time      
      USE Zmcp            
!
      IMPLICIT NONE
     
      INTEGER :: tt
      REAL(8) :: direction 
      REAL(8),ALLOCATABLE,SAVE :: mcp_mflowrate(:)
      LOGICAL,SAVE :: initial
      DATA initial /.true./
!
!......'mcp_on' is always 1 when rv_mcp=1 and becomes 0 only speed_pump=0 rpm
!       However, it is directly controlled when vflow_direct=1 is at mcp_onoff.f90.
!      
      IF(vflow_direct.eq.0) mcp_on(:)=1                    !MCP on/off in global domain
!
!......initializaton mcp model      
!
      IF(.not.ALLOCATED(mcp_mflowrate)) ALLOCATE(mcp_mflowrate(num_mcp))
      IF(initial) THEN
         mcp_mflowrate(:)=0.d0
         initial=.false.
      ENDIF   
!
!......mcp model: homologous curve or direct input
!      
      IF(vflow_direct.eq.0) THEN
         CALL mcp_model(mcp_mflowrate) !homologous curve
      ELSE
         mcp_mflowrate(:)=mcp_vflow
      ENDIF   
!
!......Ramping of MCP flow rate
!      
      IF(time.lt.time_mcp_ramping) THEN
         mcp_mflowrate(:)=mcp_mflowrate(:)*time/time_mcp_ramping          !MCP total mass flow rate (kg/s, absolute value) 
      ELSE
         mcp_mflowrate(:)=mcp_mflowrate(:)                             !MCP total mass flow rate (kg/s, absolute value)  
      ENDIF      
!
!......Define fluxl_mcpface, fluxg_mcpface, fluxd_mcpface
!       
      fluxl_mcpface(:)=0.d0
      fluxg_mcpface(:)=0.d0
      fluxd_mcpface(:)=0.d0      
      DO tt=1,num_mcp          !flux of MCPs define in global domain
         IF(mcp_on(tt).eq.0) CYCLE
!  
         IF(fzone_mcp(tt,1)-fzone_mcp(tt,2).gt.0) THEN
            direction=-1.d0        !MCP flow direction (high to low fluid zone number)
         ELSE
            direction=1.d0         !MCP flow direction (low to high fluid zone number)
         ENDIF
         fluxl_mcpface(tt)=direction*mcp_mflowrate(tt)  !volume flow rate (m3/sec)
!         IF(myrank.eq.7) print*,myrank,tt,fluxl_mcpface(tt)
      ENDDO
!            
      END SUBROUTINE mcp_flow
