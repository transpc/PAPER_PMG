MODULE Zmcp
!
   IMPLICIT NONE
   SAVE
!
   INTEGER :: num_mcp         !total number of mcps
   INTEGER :: num_mcploc      !number of mcps in each subdomain
   INTEGER :: cellnum_mcpface
   INTEGER :: icell_mcpface(1000)
!
!..Finding MCP faces
!  
   INTEGER, ALLOCATABLE :: fzone_mcp(:,:)  !fluid zones of each throat (1)=(left zone), (2)=(right zone) (1:num_mcp,1:2)
   INTEGER, ALLOCATABLE :: icell_mcp(:,:)  !inner cell number of throat faces (1:num_mcp,left zone)
   INTEGER, ALLOCATABLE :: ocell_mcp(:,:)  !outter cell number of throat faces (1:num_mcp,right zone)
   INTEGER, ALLOCATABLE :: num_mcpface(:)
   INTEGER, ALLOCATABLE :: n_face_mcp(:,:)   !non-face index of tt-th mcp (1:num_mcp,1:num_mcpface(tt))
   INTEGER, ALLOCATABLE :: dir_face_mcp(:,:) !non face direction of tt-th mcp (1:num_mcp,1:num_mcpface(tt)
   INTEGER, ALLOCATABLE :: nzone_tmp(:)     !global nzone for MPI calculation
   REAL(8), ALLOCATABLE :: num_mcpface_global(:),num_mcpface_tmp(:)  !number of mcp face (1:num_mcp)
!
!..MCP mapping: local subdomain to global domain
!      
   INTEGER, ALLOCATABLE :: mapping_mcp(:)   
!
!..nonk array for nf_nonk array
!      
   INTEGER, ALLOCATABLE :: nonk_mcp(:,:)
!
!..Flux (or velocity) on MCP
!   
   REAL(8) vl_mcp,vg_mcp,vl_mcp_o,vg_mcp_o
   REAL(8),ALLOCATABLE:: fluxl_mcpface(:),fluxg_mcpface(:),fluxd_mcpface(:)  !global mcp flux condition
!
!..Area of each MCP
!      
   REAL(8),ALLOCATABLE:: mcp_area(:),mcp_vol(:) !mcp area (1:num_mcp)
   REAL(8),ALLOCATABLE:: mcp_area_global(:),mcp_vol_global(:) !mcp area (1:num_mcp)
!
!..MCP model on/off
!       
   INTEGER, ALLOCATABLE :: mcp_on(:)   !mcp in operation (1:num_mcp)
   INTEGER :: init_mcp                 !mcp is steady state calculate (1), transient calculation (0)  
   INTEGER :: vflow_direct      
!
!..MCP model
!  
   REAL(8), dimension(11) :: frac_tabl
   REAL(8), dimension(11) :: han
   REAL(8), dimension(11) :: hvn
   REAL(8), dimension(11) :: had
   REAL(8), dimension(11) :: hvd
   REAL(8), dimension(11) :: hat
   REAL(8), dimension(11) :: hvt
   REAL(8), dimension(11) :: har
   REAL(8), dimension(11) :: hvr 
   
   REAL(8), ALLOCATABLE :: rwinit(:)           !< initial flow rate per MCP, kg/s.
   REAL(8), ALLOCATABLE :: head_mult(:)
   REAL(8), ALLOCATABLE :: dp_pump(:)          !< delta p across the pumps.   
   REAL(8), ALLOCATABLE :: rated_vol_flow(:)   !< m3/s.
   REAL(8), ALLOCATABLE :: rated_pump_speed(:) !< RPM.
   REAL(8), ALLOCATABLE :: rated_pump_hd(:)    !< head, m.
   REAL(8) :: rated_rho                        !< liquid density, kg/m3.
   REAL(8) :: relax_flow
   INTEGER :: num_mcp_transient                       !number of MCP transient intervals
   REAL(8), ALLOCATABLE :: mcp_transient_start(:)    !mcp transient calculation start time (sec)
   REAL(8), ALLOCATABLE :: speed_pump(:,:)           !mcp speed control (rpm) (1:num_mcp, 1:mcp_transient_strat)
!
   REAL(8) :: mcp_vflow,time_mcp_ramping   
   REAL(8) :: time_mcp_on,time_mcp_off
!
END MODULE Zmcp
