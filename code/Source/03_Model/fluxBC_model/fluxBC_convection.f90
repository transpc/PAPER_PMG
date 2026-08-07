!
      SUBROUTINE fluxBC_convection_ice(fluxl_c_nf,fluxg_c_nf,fluxd_c_nf,mfluxl_cup_nf,mfluxg_cup_nf,mfluxd_cup_nf)
!
      USE Zparam       , ONLY: ndim     
      USE Zvec_param   , ONLY: nf_flux  
      USE Znum_cell    , ONLY: istart_nf
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Zvec_major   , ONLY: flux_l_nf,flux_g_nf,flux_d_nf      
      USE Zvec_geo     , ONLY: saa_nf,xn_nf
      USE Zare         , ONLY: ar_gas,ar_liq,ar_drp
      USE Zrv_choke    , ONLY: num_throatface,n_face_throat,choke
      USE Zmcp         , ONLY: num_mcploc,num_mcpface,n_face_mcp,mcp_on,mapping_mcp
      USE Zvalve       , ONLY: valve_closed,num_valveloc,num_valveface,n_face_valve,mapping_valve
!
      IMPLICIT NONE
!
      INTEGER :: nf_number,istart
      INTEGER :: i,i1,ii,ir,kk,ix,zz,tt
      REAL(8) :: a_l,a_g,a_d,b_l,b_g,b_d
      REAL(8) :: agrk,alrk,adrk,agri,alri,adri
      REAL(8) :: vl(ndim),vg(ndim),vd(ndim)
      REAL(8) :: fluxl_c_nf(nf_flux,ndim),fluxg_c_nf(nf_flux,ndim),fluxd_c_nf(nf_flux,ndim)      
      REAL(8) :: mfluxl_cup_nf(nf_flux),mfluxg_cup_nf(nf_flux),mfluxd_cup_nf(nf_flux)
!
      nf_number=0
      istart=istart_nf(1,nf_number)
!
      IF(.not.choke) GOTO 100      
!
!......choke model
!
      DO i=1,num_throatface
         i1=n_face_throat(i)
         ii=left_nf(i1)                
         ir=i1-istart                !i1=istart+i 
         kk=right_non(ir)
               
         a_l=min(flux_l_nf(i1),0.d0) !flux_l_nf for MCP face is already updated in pressure_matrix subroutines
         b_l=max(flux_l_nf(i1),0.d0)
         a_g=min(flux_g_nf(i1),0.d0)
         b_g=max(flux_g_nf(i1),0.d0)
         a_d=min(flux_d_nf(i1),0.d0)
         b_d=max(flux_d_nf(i1),0.d0)               
!
         agrk=ar_gas(kk)
         alrk=ar_liq(kk)
         adrk=ar_drp(kk)
         agri=ar_gas(ii)
         alri=ar_liq(ii)
         adri=ar_drp(ii)          
!
         mfluxl_cup_nf(i1)=a_l*alrk+b_l*alri
         mfluxg_cup_nf(i1)=a_g*agrk+b_g*agri
         mfluxd_cup_nf(i1)=a_d*adrk+b_d*adri
!
         DO ix=1,ndim
            vl(ix)=flux_l_nf(i1)*xn_nf(i1,ix)/saa_nf(i1)
            vg(ix)=flux_g_nf(i1)*xn_nf(i1,ix)/saa_nf(i1)
            vd(ix)=flux_d_nf(i1)*xn_nf(i1,ix)/saa_nf(i1)
!                
            fluxl_c_nf(i1,ix)=(a_l*alrk)*vl(ix)+(b_l*alri)*vl(ix)
            fluxg_c_nf(i1,ix)=(a_g*agrk)*vg(ix)+(b_g*agri)*vg(ix)
            fluxd_c_nf(i1,ix)=(a_d*adrk)*vd(ix)+(b_d*adri)*vd(ix)            
         ENDDO
      ENDDO  
!
100   CONTINUE
!
!......mcp model
!      
      DO zz=1,num_mcploc                !local loop
         tt=mapping_mcp(zz)             !mapping local MCP to global MCP number
         IF(mcp_on(tt).eq.0) CYCLE       !mcp_on: global MCP number  
         DO i=1,num_mcpface(zz)
            i1=n_face_mcp(zz,i) 
            ii=left_nf(i1)                
            ir=i1-istart                !i1=istart+i 
            kk=right_non(ir)
!                
            a_l=min(flux_l_nf(i1),0.d0)  !flux_l_nf for MCP face is already updated in pressure_matrix subroutines
            b_l=max(flux_l_nf(i1),0.d0)
            a_g=min(flux_g_nf(i1),0.d0)
            b_g=max(flux_g_nf(i1),0.d0)
            a_d=min(flux_d_nf(i1),0.d0)
            b_d=max(flux_d_nf(i1),0.d0)
!
            agrk=ar_gas(kk)
            alrk=ar_liq(kk)
            adrk=ar_drp(kk)
            agri=ar_gas(ii)
            alri=ar_liq(ii)
            adri=ar_drp(ii)          
!
            mfluxl_cup_nf(i1)=a_l*alrk+b_l*alri
            mfluxg_cup_nf(i1)=a_g*agrk+b_g*agri
            mfluxd_cup_nf(i1)=a_d*adrk+b_d*adri
!
            DO ix=1,ndim
               vl(ix)=flux_l_nf(i1)*xn_nf(i1,ix)/saa_nf(i1)
               vg(ix)=flux_g_nf(i1)*xn_nf(i1,ix)/saa_nf(i1)
               vd(ix)=flux_d_nf(i1)*xn_nf(i1,ix)/saa_nf(i1)
!               
               fluxl_c_nf(i1,ix)=(a_l*alrk)*vl(ix)+(b_l*alri)*vl(ix)
               fluxg_c_nf(i1,ix)=(a_g*agrk)*vg(ix)+(b_g*agri)*vg(ix)
               fluxd_c_nf(i1,ix)=(a_d*adrk)*vd(ix)+(b_d*adri)*vd(ix)
            ENDDO   
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
            mfluxl_cup_nf(i1)=0.d0
            mfluxg_cup_nf(i1)=0.d0
            mfluxd_cup_nf(i1)=0.d0
            fluxl_c_nf(i1,:)=0.d0
            fluxg_c_nf(i1,:)=0.d0
            fluxd_c_nf(i1,:)=0.d0
         ENDDO
      ENDDO                 
!            
      END SUBROUTINE fluxBC_convection_ice
    
    
!
      SUBROUTINE fluxBC_convection_smac3(diag_g_nf,diag_l_nf,off_diag_l_non_i,off_diag_g_non_i,off_diag_l_non_k,off_diag_g_non_k,src_l,src_g)
!
      USE Zzone        , ONLY: ncell_fluid,ncell_fluid_pad
      USE Zparam       , ONLY: ndim     
      USE Zvec_param   , ONLY: nf_non,nf_fluxk1,nf_nonk
      USE Znum_cell    , ONLY: istart_nf,istart_nfs,nf_number_id
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Zvec_major   , ONLY: flux_l_nf,flux_g_nf
      USE Zvec_geo     , ONLY: saa_nf,xn_nf
      USE Zare         , ONLY: ar_gas,ar_liq
      USE Zrv_choke    , ONLY: num_throatface,n_face_throat,choke,nonk_throat
      USE Zmcp         , ONLY: num_mcploc,num_mcpface,n_face_mcp,mcp_on,mapping_mcp,nonk_mcp
!
      IMPLICIT NONE
!           
      INTEGER :: nf_number,istart
      INTEGER :: i,i1,ii,ir,kk,ix,zz,tt,nv,istart0,i0,ik
      REAL(8) :: vl,vg
      REAL(8) :: diagl,diagg
      REAL(8),DIMENSION(ncell_fluid_pad,ndim) :: src_g,src_l
      REAL(8),DIMENSION(nf_fluxk1) :: diag_g_nf,diag_l_nf      
      REAL(8),DIMENSION(nf_non) :: off_diag_g_non_i,off_diag_l_non_i
      REAL(8),DIMENSION(nf_nonk) :: off_diag_g_non_k,off_diag_l_non_k      
!
      IF(.not.choke) GOTO 100
!
!......choke model
!
      nv=0
      nf_number=nf_number_id(nv)
      istart=istart_nf(1,nf_number)      
      istart0=istart_nfs(0)          
      DO i=1,num_throatface
         i1=n_face_throat(i)
         i0=istart0-istart+i1  !i0=istart0+i=istart0+(i1-start) 
         ii=left_nf(i1) 
         ir=i1-istart
         kk=right_non(ir)
         ik=nonk_throat(i)
!         
         off_diag_l_non_i(ir)=off_diag_l_non_i(ir)+diag_l_nf(i0)
         off_diag_g_non_i(ir)=off_diag_g_non_i(ir)+diag_g_nf(i0)
         
         kk=right_non(i1-istart) !i1=istart+i
         IF(kk.le.ncell_fluid) THEN  
            off_diag_l_non_k(ik)=off_diag_l_non_k(ik)+diag_l_nf(ik)
            off_diag_g_non_k(ik)=off_diag_g_non_k(ik)+diag_g_nf(ik)
         ENDIF 
!              
         IF(flux_l_nf(i1).lt.0.d0) THEN !!left cell
            DO ix=1,ndim
               vl=flux_l_nf(i1)*xn_nf(i1,ix)/saa_nf(i1)  !xn_nf=vector mode (+ when ii<kk), svp(i1,:)=outward when inlet, 
               vg=flux_g_nf(i1)*xn_nf(i1,ix)/saa_nf(i1)
!                 
               diagl=flux_l_nf(i1)*ar_liq(kk)
               diagg=flux_g_nf(i1)*ar_gas(kk)                    
!
               src_l(ii,ix)=src_l(ii,ix)-diagl*vl  !left cell receives. alprhouf is moved to RHS
               src_g(ii,ix)=src_g(ii,ix)-diagg*vg                 
            ENDDO  
!            
!         ELSE
!            DO ix=1,ndim
!               vl=flux_l_nf(i1)*xn_nf(i1,ix)/saa_nf(i1)  !xn_nf=vector mode (+ when ii<kk), svp(i1,:)=outward when inlet, 
!               vg=flux_g_nf(i1)*xn_nf(i1,ix)/saa_nf(i1)
!!                  
!               diagl=flux_l_nf(i1)*ar_liq(ii) !flux on kk-cell is inverse direction of ii-cell
!               diagg=flux_g_nf(i1)*ar_gas(ii)                   
!!
!               src_l(kk,ix)=src_l(kk,ix)+diagl*vl
!               src_g(kk,ix)=src_g(kk,ix)+diagg*vg               
!            ENDDO                 
!            
         ENDIF   
      ENDDO   
!
100 CONTINUE      
!
!......mcp model
!         
      nv=0
      nf_number=nf_number_id(nv)
      istart=istart_nf(1,nf_number)      
      istart0=istart_nfs(0)    
      DO zz=1,num_mcploc                !local loop
         tt=mapping_mcp(zz)             !mapping local MCP to global MCP number
         IF(mcp_on(tt).eq.0) CYCLE      !mcp_on: global MCP number  
         DO i=1,num_mcpface(zz)
            i1=n_face_mcp(zz,i) 
            i0=istart0-istart+i1  !i0=istart0+i=istart0+(i1-start) 
            ii=left_nf(i1)  
            ir=i1-istart
            kk=right_non(ir)     !i1=istart+i
            ik=nonk_mcp(zz,i)
            
            off_diag_l_non_i(ir)=off_diag_l_non_i(ir)+diag_l_nf(i0)
            off_diag_g_non_i(ir)=off_diag_g_non_i(ir)+diag_g_nf(i0) 
            
            kk=right_non(i1-istart) !i1=istart+i
            IF(kk.le.ncell_fluid) THEN  
               off_diag_l_non_k(ik)=off_diag_l_non_k(ik)+diag_l_nf(ik)
               off_diag_g_non_k(ik)=off_diag_g_non_k(ik)+diag_g_nf(ik)
            ENDIF            
!
            IF(flux_l_nf(i1).lt.0.d0) THEN
               DO ix=1,ndim
                  vl=flux_l_nf(i1)*xn_nf(i1,ix)/saa_nf(i1)  !xn_nf=vector mode (+ when ii<kk), svp(i1,:)=outward when inlet, 
                  vg=flux_g_nf(i1)*xn_nf(i1,ix)/saa_nf(i1)
!                  
                  diagl=flux_l_nf(i1)*ar_liq(kk)
                  diagg=flux_g_nf(i1)*ar_gas(kk)                    
!                  
                  src_l(ii,ix)=src_l(ii,ix)-diagl*vl
                  src_g(ii,ix)=src_g(ii,ix)-diagg*vg                     
               ENDDO 
!               
!            ELSE
!               DO ix=1,ndim
!                  vl=flux_l_nf(i1)*xn_nf(i1,ix)/saa_nf(i1)  !xn_nf=vector mode (+ when ii<kk), svp(i1,:)=outward when inlet, 
!                  vg=flux_g_nf(i1)*xn_nf(i1,ix)/saa_nf(i1)
!!                  
!                  diagl=flux_l_nf(i1)*ar_liq(ii) 
!                  diagg=flux_g_nf(i1)*ar_gas(ii)   
!!
!                  src_l(kk,ix)=src_l(kk,ix)+diagl*vl
!                  src_g(kk,ix)=src_g(kk,ix)+diagg*vg
!               ENDDO   
!               
            ENDIF            
!                   
         ENDDO   
      ENDDO   
!      
      END SUBROUTINE fluxBC_convection_smac3
    
    
