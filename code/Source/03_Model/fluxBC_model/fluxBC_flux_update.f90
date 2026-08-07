!
      SUBROUTINE fluxBC_flux_update_ice(flux_l_nf,flux_g_nf,flux_d_nf,fluxt_l_nf,fluxt_g_nf,fluxt_d_nf)
!      
      USE Zvec_param    , ONLY: nf_flux
      USE Znum_cell     , ONLY: istart_nf
      USE Zvec_index    , ONLY: left_nf,right_non
      USE Zrv_choke     , ONLY: num_throatface,n_face_throat,num_throatface,choke,dir_face_throat,vl_choke,vg_choke
      USE Zvec_geo      , ONLY: saa_nf,xn_nf
      USE Zvector       , ONLY: vl_n,vg_n      
      USE Zmcp          , ONLY: num_mcploc,num_mcpface,n_face_mcp,mcp_on,mapping_mcp,mcp_area_global, &
                                fluxl_mcpface,fluxg_mcpface,fluxd_mcpface,dir_face_mcp
      USE Zvalve        , ONLY: valve_closed,num_valveloc,num_valveface,n_face_valve,mapping_valve      
!
      IMPLICIT NONE
!
      INTEGER :: nf_number,istart      
      INTEGER :: i,i1,zz,tt,ii,kk,ir
      REAL(8) :: fluxl_throatface,fluxg_throatface,fluxd_throatface
      REAL(8) :: flux_l_nf(nf_flux),flux_g_nf(nf_flux),flux_d_nf(nf_flux)
      REAL(8) :: fluxt_l_nf(nf_flux),fluxt_g_nf(nf_flux),fluxt_d_nf(nf_flux)  
!
      nf_number=0
      istart=istart_nf(1,nf_number)         
!      
!......choke model
!          
      IF(.not.choke) GOTO 100
!      
!1) update volume flow rate          
      DO i=1,num_throatface
         i1=n_face_throat(i)
         fluxl_throatface=flux_l_nf(i1)          
         fluxg_throatface=flux_g_nf(i1)    
         fluxd_throatface=flux_d_nf(i1)
         flux_l_nf(i1)=vl_choke*saa_nf(i1)*dir_face_throat(i)
         flux_g_nf(i1)=vg_choke*saa_nf(i1)*dir_face_throat(i)
         flux_d_nf(i1)=vl_choke*saa_nf(i1)*dir_face_throat(i)
         IF(flux_l_nf(i1).gt.0) flux_l_nf(i1)=MIN(fluxl_throatface,flux_l_nf(i1))
         IF(flux_g_nf(i1).gt.0) flux_g_nf(i1)=MIN(fluxg_throatface,flux_g_nf(i1))
         IF(flux_d_nf(i1).gt.0) flux_d_nf(i1)=MIN(fluxd_throatface,flux_d_nf(i1))
         IF(flux_l_nf(i1).lt.0) flux_l_nf(i1)=MAX(fluxl_throatface,flux_l_nf(i1))
         IF(flux_g_nf(i1).lt.0) flux_g_nf(i1)=MAX(fluxg_throatface,flux_g_nf(i1))
         IF(flux_d_nf(i1).lt.0) flux_d_nf(i1)=MAX(fluxd_throatface,flux_d_nf(i1))
         IF(flux_l_nf(i1).eq.0.d0) flux_l_nf(i1)=0.d0
         IF(flux_g_nf(i1).eq.0.d0) flux_g_nf(i1)=0.d0
         IF(flux_d_nf(i1).eq.0.d0) flux_d_nf(i1)=0.d0              
!
         fluxt_l_nf(i1)=flux_l_nf(i1)
         fluxt_g_nf(i1)=flux_g_nf(i1)
         fluxt_d_nf(i1)=flux_d_nf(i1)
!            
      ENDDO
!
!2) update u*
      DO i=1,num_throatface
         i1=n_face_throat(i)
         ii=left_nf(i1)
         ir=i1-istart
         kk=right_non(ir)
         vl_n(ii,:)=vl_choke*dir_face_throat(i)*xn_nf(i1,:)
         vg_n(ii,:)=vg_choke*dir_face_throat(i)*xn_nf(i1,:)
         vl_n(kk,:)=vl_n(ii,:)
         vg_n(kk,:)=vg_n(ii,:)
      ENDDO    
!      
!
100 CONTINUE       
!         
!......MCP model
!      
      DO zz=1,num_mcploc                !local loop
         tt=mapping_mcp(zz)             !mapping local MCP to global MCP number
         IF(mcp_on(tt).eq.0) CYCLE      !mcp_on: global MCP number  
          DO i=1,num_mcpface(zz)
             i1=n_face_mcp(zz,i) 
!             
             flux_l_nf(i1)=fluxl_mcpface(tt)*dir_face_mcp(zz,i)/mcp_area_global(tt)*saa_nf(i1)  !global MCP flux to local MCP flux on face
             flux_g_nf(i1)=fluxg_mcpface(tt)*dir_face_mcp(zz,i)/mcp_area_global(tt)*saa_nf(i1)
             flux_d_nf(i1)=fluxd_mcpface(tt)*dir_face_mcp(zz,i)/mcp_area_global(tt)*saa_nf(i1)   
!             
             fluxt_l_nf(i1)=flux_l_nf(i1)
             fluxt_g_nf(i1)=flux_g_nf(i1)
             fluxt_d_nf(i1)=flux_d_nf(i1)   
          ENDDO
      ENDDO  
!
!......valve model
!         
      DO zz=1,num_valveloc                
         tt=mapping_valve(zz)             
         IF(valve_closed(tt).eq.0) CYCLE   
         DO i=1,num_valveface(zz)
            i1=n_face_valve(zz,i) 
! 
            flux_l_nf(i1)=0.d0
            flux_g_nf(i1)=0.d0
            flux_d_nf(i1)=0.d0
!            
            fluxt_l_nf(i1)=0.d0
            fluxt_g_nf(i1)=0.d0
            fluxt_d_nf(i1)=0.d0
         ENDDO
      ENDDO         
!                                                                       
      RETURN
      END SUBROUTINE fluxBC_flux_update_ice
!
!!!    
!
      SUBROUTINE fluxBC_flux_update(flux_l_nf,flux_g_nf,flux_d_nf)
!
      USE Zvec_param    , ONLY: nf_flux
      USE Zvec_index    , ONLY: left_nf,right_non
      USE Znum_cell     , ONLY: istart_nf      
      USE Zrv_choke     , ONLY: num_throatface,n_face_throat,num_throatface,choke,dir_face_throat,vl_choke,vg_choke
      USE Zvec_geo      , ONLY: saa_nf,xn_nf
      USE Zvector       , ONLY: vl_n,vg_n      
      USE Zmcp          , ONLY: num_mcploc,num_mcpface,n_face_mcp,mcp_on,mapping_mcp,mcp_area_global, &
                                fluxl_mcpface,fluxg_mcpface,fluxd_mcpface,dir_face_mcp
      USE Zvalve        , ONLY: valve_closed,num_valveloc,num_valveface,n_face_valve,mapping_valve    
!
      IMPLICIT NONE
!
      INTEGER :: nf_number,istart
      INTEGER :: i,i1,zz,tt,ii,kk,ir
      REAL(8) :: fluxl_throatface,fluxg_throatface,fluxd_throatface
      REAL(8) :: flux_l_nf(nf_flux),flux_g_nf(nf_flux),flux_d_nf(nf_flux)
!
      nf_number=0
      istart=istart_nf(1,nf_number)         
!      
!......choke model
!          
      IF(.not.choke) GOTO 100
!      
!1) update volume flow rate       
      DO i=1,num_throatface
         i1=n_face_throat(i)
         fluxl_throatface=flux_l_nf(i1)          
         fluxg_throatface=flux_g_nf(i1)    
         fluxd_throatface=flux_d_nf(i1)
         flux_l_nf(i1)=vl_choke*saa_nf(i1)*dir_face_throat(i)
         flux_g_nf(i1)=vg_choke*saa_nf(i1)*dir_face_throat(i)
         flux_d_nf(i1)=vl_choke*saa_nf(i1)*dir_face_throat(i)
         IF(flux_l_nf(i1).gt.0) flux_l_nf(i1)=MIN(fluxl_throatface,flux_l_nf(i1))
         IF(flux_g_nf(i1).gt.0) flux_g_nf(i1)=MIN(fluxg_throatface,flux_g_nf(i1))
         IF(flux_d_nf(i1).gt.0) flux_d_nf(i1)=MIN(fluxd_throatface,flux_d_nf(i1))
         IF(flux_l_nf(i1).lt.0) flux_l_nf(i1)=MAX(fluxl_throatface,flux_l_nf(i1))
         IF(flux_g_nf(i1).lt.0) flux_g_nf(i1)=MAX(fluxg_throatface,flux_g_nf(i1))
         IF(flux_d_nf(i1).lt.0) flux_d_nf(i1)=MAX(fluxd_throatface,flux_d_nf(i1))
         IF(flux_l_nf(i1).eq.0.d0) flux_l_nf(i1)=0.d0
         IF(flux_g_nf(i1).eq.0.d0) flux_g_nf(i1)=0.d0
         IF(flux_d_nf(i1).eq.0.d0) flux_d_nf(i1)=0.d0   
      ENDDO
!
!2) update u*
      DO i=1,num_throatface
         i1=n_face_throat(i)
         ii=left_nf(i1)
         ir=i1-istart
         kk=right_non(ir)
         vl_n(ii,:)=vl_choke*dir_face_throat(i)*xn_nf(i1,:)
         vg_n(ii,:)=vg_choke*dir_face_throat(i)*xn_nf(i1,:)
         vl_n(kk,:)=vl_n(ii,:)
         vg_n(kk,:)=vg_n(ii,:)
      ENDDO          
!
100 CONTINUE     
!      
!.....MCP model
!      
      DO zz=1,num_mcploc                !local loop
         tt=mapping_mcp(zz)             !mapping local MCP to global MCP number
         IF(mcp_on(tt).eq.0) CYCLE      !mcp_on: global MCP number   
          DO i=1,num_mcpface(zz)
             i1=n_face_mcp(zz,i) 
             flux_l_nf(i1)=fluxl_mcpface(tt)*dir_face_mcp(zz,i)/mcp_area_global(tt)*saa_nf(i1)  !global MCP flux to local MCP flux on face
             flux_g_nf(i1)=fluxg_mcpface(tt)*dir_face_mcp(zz,i)/mcp_area_global(tt)*saa_nf(i1)
             flux_d_nf(i1)=fluxd_mcpface(tt)*dir_face_mcp(zz,i)/mcp_area_global(tt)*saa_nf(i1) 
          ENDDO
      ENDDO 
!
!......valve model
!         
      DO zz=1,num_valveloc                
         tt=mapping_valve(zz)             
         IF(valve_closed(tt).eq.0) CYCLE   
         DO i=1,num_valveface(zz)
            i1=n_face_valve(zz,i) 
            flux_l_nf(i1)=0.d0
            flux_g_nf(i1)=0.d0
            flux_d_nf(i1)=0.d0
         ENDDO
      ENDDO  
!      
      RETURN
      END SUBROUTINE fluxBC_flux_update
