!
      SUBROUTINE fluxBC_init
!
      USE Zinterface
      USE Zzone         , ONLY: ncell_fluid
      USE Zmcp          , ONLY: mcp_on,vflow_direct,time_mcp_on,time_mcp_off
      USE Zvalve        , ONLY: valve_closed,time_valve_closed
      USE Ztimecon      , ONLY: time
      USE Zrv_model     , ONLY: rv_choke,rv_mcp,rv_valve
      USE Zbc_index     , ONLY: ngrad
      USE Zgrad_ls_c3d  , ONLY: lsindex
      USE Zvalve        , ONLY: num_valveloc,num_valveface,n_face_valve,mapping_valve, &
                             sap_nf_o,sad_non_o,sa_nf_o
      USE Zvec_index    , ONLY: left_nf,right_non
      USE Znum_cell     , ONLY: istart_nf      
      USE Zcore 
!
      IMPLICIT NONE
!
      INTEGER :: tt,zz,ii,i1,i,istart,ir,kk,nf_number
      LOGICAL, SAVE::initial
      DATA initial /.true./     
!
!......Define fluxBC faces
!      
      IF(initial) THEN
         IF(rv_choke.eq.1) CALL choke_define_fluxBC
         IF(rv_mcp.eq.1) THEN
             CALL mcp_define_fluxBC
             CALL mcp_mapping
         ENDIF     
         IF(rv_valve.eq.1) CALL valve_define_fluxBC
         initial=.false.
      ENDIF
!
!      mcp model is used when rv_mcp=1. For vflow_direct=1, mcp model operates time_mcp_on<time<time_mcp_off
      IF(rv_mcp.eq.1) THEN
         IF(vflow_direct.eq.1) THEN
            IF(time.ge.time_mcp_on.and.time.lt.time_mcp_off) THEN     !Local MCP on/off is controlled by 0 RPM
               mcp_on(:)=1                                       !MCP on/off in global domain
            ELSE
               mcp_on(:)=0
            ENDIF
         ENDIF  
      ENDIF   
      
!      valve model is used when time_valve_closed(:,1)<time<time_valve_closed(:,2)
      IF(rv_valve.eq.1) THEN
         nf_number=0
         istart=istart_nf(1,nf_number) 
         sap_nf_o=0.d0
         sa_nf_o=0.d0
         sad_non_o=0.d0
         sad_non_o=0.d0
         DO zz=1,num_valveloc
            tt=mapping_valve(zz) 
            IF(time.ge.time_valve_closed(tt,1).and.time.le.time_valve_closed(tt,2)) THEN
               valve_closed(tt)=1 !closed
               DO i=1,num_valveface(zz)
                  i1=n_face_valve(zz,i) 
                  ii=left_nf(i1)
                  ir=i1-istart
                  kk=right_non(ir)
                  ngrad(ii)=1
                  IF(kk.le.ncell_fluid) THEN
                     ngrad(kk)=1                      
                  ENDIF
!                  
!temporary                  
                  IF(tt.eq.2) THEN
                     lsindex(ii)=0 !for lsquareoff=0
                     ngrad(ii)=0
                     lsindex(kk)=0 !for lsquareoff=0
                     ngrad(kk)=0
                  ENDIF
!                  
               ENDDO                 
            ELSE
               valve_closed(tt)=0 !open
               DO i=1,num_valveface(zz)
                  i1=n_face_valve(zz,i) 
                  ii=left_nf(i1)
                  lsindex(ii)=0
                  ngrad(ii)=0
                  ir=i1-istart
                  kk=right_non(ir)
                  IF(kk.le.ncell_fluid) THEN
                     lsindex(kk)=0 
                     ngrad(kk)=0
                  ENDIF  
               ENDDO                
            ENDIF   
            
         ENDDO      
      ENDIF   
      
      CALL barrier_mpi
!      
      RETURN
      END SUBROUTINE fluxBC_init
