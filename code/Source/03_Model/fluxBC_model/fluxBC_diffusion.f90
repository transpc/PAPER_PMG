
!
      SUBROUTINE fluxBC_diffusion_ice(fluxl_diff_nf,fluxg_diff_nf,fluxd_diff_nf)
!
      USE Vol_DATA     , ONLY: cell
      USE Zzone        , ONLY: ncell_fluid      
      USE Zparam       , ONLY: ndim     
      USE Zvec_param   , ONLY: nf_non,nf_inl,nf_adw,nf_fsw,nf_ctw,nf_chw,nf_nonk      
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Znum_cell    , ONLY: istart_nf,istart_nfs,nf_number_id,right_nb_k
      USE Zvec_major   , ONLY: flux_l_nf,flux_g_nf,flux_d_nf      
      USE Zvec_geo     , ONLY: saa_nf,xn_nf,sad_non
      USE Zvector      , ONLY: vg_o,vl_o,vd_o     
      USE Zrv_choke    , ONLY: num_throatface,n_face_throat,choke,nonk_throat
      USE Zmcp         , ONLY: num_mcploc,num_mcpface,n_face_mcp,mcp_on,mapping_mcp,nonk_mcp
      USE Zvalve       , ONLY: valve_closed,num_valveloc,num_valveface,n_face_valve,mapping_valve,nonk_valve
!
      IMPLICIT NONE
!
      INTEGER :: nv,nf_number,istart,istart0
      INTEGER :: i,i1,ii,ir,ix,zz,tt,i0,kk,ik,k,kk0
      REAL(8) :: avisli,avisgi,avisdi
      REAL(8) :: vl,vg,vd,vlf,vgf,vdf
      REAL(8) :: fluxl_diff_nf(nf_nonk+nf_non+nf_inl+nf_adw+nf_fsw+nf_ctw+nf_chw,ndim) 
      REAL(8) :: fluxg_diff_nf(nf_nonk+nf_non+nf_inl+nf_adw+nf_fsw+nf_ctw+nf_chw,ndim) 
      REAL(8) :: fluxd_diff_nf(nf_nonk+nf_non+nf_inl+nf_adw+nf_fsw+nf_ctw+nf_chw,ndim) 
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
         DO ix=1,ndim
            vlf=flux_l_nf(i1)*xn_nf(i1,ix)/saa_nf(i1)  !xn_nf=vector mode (+ when ii<kk), svp(i1,:)=outward when inlet, 
            vgf=flux_g_nf(i1)*xn_nf(i1,ix)/saa_nf(i1)            
            vdf=flux_d_nf(i1)*xn_nf(i1,ix)/saa_nf(i1)  !xn_nf=vector mode (+ when ii<kk), svp(i1,:)=outward when inlet, 
!            
            vl=vlf-vl_o(ii,ix)
            vg=vgf-vg_o(ii,ix)
            vd=vdf-vd_o(ii,ix)
            avisli=cell%alphal(ii)*cell%eviscosl(ii)
            avisgi=cell%alphag(ii)*cell%eviscosg(ii)
            avisdi=cell%alphad(ii)*cell%eviscosd(ii)
            fluxl_diff_nf(i0,ix)=avisli*vl*sad_non(ir)
            fluxg_diff_nf(i0,ix)=avisgi*vg*sad_non(ir)
            fluxd_diff_nf(i0,ix)=avisdi*vd*sad_non(ir)
!            
            IF(kk.le.ncell_fluid) THEN  
               ik=nonk_throat(i)      ! k=right_nb_k(ik) --> left_nf(k), right_non(k)        
               k=right_nb_k(ik)
               kk0=right_non(k)  
!               
               vl=vlf-vl_o(kk0,ix)
               vg=vgf-vg_o(kk0,ix)
               vd=vdf-vd_o(kk0,ix)            
               avisli=cell%alphal(kk0)*cell%eviscosl(kk0)
               avisgi=cell%alphag(kk0)*cell%eviscosg(kk0)
               avisdi=cell%alphad(kk0)*cell%eviscosd(kk0)
!               
               fluxl_diff_nf(ik,ix)=avisli*vl*sad_non(ik)
               fluxg_diff_nf(ik,ix)=avisgi*vg*sad_non(ik)
               fluxd_diff_nf(ik,ix)=avisdi*vd*sad_non(ik)               
            ENDIF
!            
         ENDDO
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
            kk=right_non(ir) !i1=istart+i
            DO ix=1,ndim
               vlf=flux_l_nf(i1)*xn_nf(i1,ix)/saa_nf(i1)  !xn_nf=vector mode (+ when ii<kk), svp(i1,:)=outward when inlet, 
               vgf=flux_g_nf(i1)*xn_nf(i1,ix)/saa_nf(i1)            
               vdf=flux_d_nf(i1)*xn_nf(i1,ix)/saa_nf(i1)  !xn_nf=vector mode (+ when ii<kk), svp(i1,:)=outward when inlet, 
!            
               vl=vlf-vl_o(ii,ix)
               vg=vgf-vg_o(ii,ix)
               vd=vdf-vd_o(ii,ix)
               avisli=cell%alphal(ii)*cell%eviscosl(ii)
               avisgi=cell%alphag(ii)*cell%eviscosg(ii)
               avisdi=cell%alphad(ii)*cell%eviscosd(ii)
               fluxl_diff_nf(i0,ix)=avisli*vl*sad_non(ir)
               fluxg_diff_nf(i0,ix)=avisgi*vg*sad_non(ir)
               fluxd_diff_nf(i0,ix)=avisdi*vd*sad_non(ir)
!            
               IF(kk.le.ncell_fluid) THEN               
                  ik=nonk_mcp(zz,i)      ! k=right_nb_k(ik) --> left_nf(k), right_non(k)        
                  k=right_nb_k(ik)
                  kk0=right_non(k)  
!                  
                  vl=vlf-vl_o(kk0,ix)
                  vg=vgf-vg_o(kk0,ix)
                  vd=vdf-vd_o(kk0,ix)            
                  avisli=cell%alphal(kk0)*cell%eviscosl(kk0)
                  avisgi=cell%alphag(kk0)*cell%eviscosg(kk0)
                  avisdi=cell%alphad(kk0)*cell%eviscosd(kk0)
!                  
                  fluxl_diff_nf(ik,ix)=avisli*vl*sad_non(ik)
                  fluxg_diff_nf(ik,ix)=avisgi*vg*sad_non(ik)
                  fluxd_diff_nf(ik,ix)=avisdi*vd*sad_non(ik)                  
               ENDIF   
            ENDDO            
         ENDDO
      ENDDO  
!
!......valve model
!         
      nv=0
      nf_number=nf_number_id(nv)
      istart=istart_nf(1,nf_number)      
      istart0=istart_nfs(0)      
      DO zz=1,num_valveloc                
         tt=mapping_valve(zz)             
         IF(valve_closed(tt).eq.0) CYCLE   
         DO i=1,num_valveface(zz)
            i1=n_face_valve(zz,i) 
            i0=istart0-istart+i1  !i0=istart0+i=istart0+(i1-start) 
            ii=left_nf(i1) 
            ir=i1-istart                !i1=istart+i 
            kk=right_non(ir)
            DO ix=1,ndim
               vlf=0.d0
               vgf=0.d0
               vdf=0.d0
!            
               vl=vlf-vl_o(ii,ix)
               vg=vgf-vg_o(ii,ix)
               vd=vdf-vd_o(ii,ix)
               avisli=cell%alphal(ii)*cell%eviscosl(ii)
               avisgi=cell%alphag(ii)*cell%eviscosg(ii)
               avisdi=cell%alphad(ii)*cell%eviscosd(ii)
               fluxl_diff_nf(i0,ix)=avisli*vl*sad_non(ir)
               fluxg_diff_nf(i0,ix)=avisgi*vg*sad_non(ir)
               fluxd_diff_nf(i0,ix)=avisdi*vd*sad_non(ir)
!            
               IF(kk.le.ncell_fluid) THEN    
                  ik=nonk_valve(zz,i)      ! k=right_nb_k(ik) --> left_nf(k), right_non(k) 
                  k=right_nb_k(ik)
                  kk0=right_non(k)
!                  
                  vl=vlf-vl_o(kk0,ix)
                  vg=vgf-vg_o(kk0,ix)
                  vd=vdf-vd_o(kk0,ix)
                  avisli=cell%alphal(kk0)*cell%eviscosl(kk0)
                  avisgi=cell%alphag(kk0)*cell%eviscosg(kk0)
                  avisdi=cell%alphad(kk0)*cell%eviscosd(kk0)
                  fluxl_diff_nf(ik,ix)=avisli*vl*sad_non(ik)
                  fluxg_diff_nf(ik,ix)=avisgi*vg*sad_non(ik)
                  fluxd_diff_nf(ik,ix)=avisdi*vd*sad_non(ik)                  
               ENDIF
!               
            ENDDO 
         ENDDO
      ENDDO                 
!            
    END SUBROUTINE fluxBC_diffusion_ice
    
    
!
      SUBROUTINE fluxBC_diffusion_smac3(diag_l_nf,diag_g_nf,off_diag_l_non_i,off_diag_g_non_i,off_diag_l_non_k,off_diag_g_non_k,diag_l,diag_g,src_l,src_g,iter)
!
      USE Vol_DATA     , ONLY: cell
      USE Zzone        , ONLY: ncell_fluid,ncell_fluid_pad
      USE Zparam       , ONLY: ndim     
      USE Zvec_param   , ONLY: nf_non,nf_inl,nf_adw,nf_fsw,nf_ctw,nf_chw ,nf_nonk     
      USE Znum_cell    , ONLY: istart_nf,istart_nfs,nf_number_id,right_nb_k
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Zvec_major   , ONLY: flux_l_nf,flux_g_nf    
      USE Zvec_geo     , ONLY: saa_nf,xn_nf,djir_non
      USE Zrv_choke    , ONLY: num_throatface,n_face_throat,choke,nonk_throat
      USE Zmcp         , ONLY: num_mcploc,num_mcpface,n_face_mcp,mcp_on,mapping_mcp,nonk_mcp
      USE Zvalve       , ONLY: valve_closed,num_valveloc,num_valveface,n_face_valve,mapping_valve,nonk_valve
!
      IMPLICIT NONE
!
      INTEGER :: nv,nf_number,istart,istart0
      INTEGER :: i,i1,ii,ir,ix,zz,tt,i0,kk,iter,ik,k,kk0
      REAL(8) :: cf_g,cf_l
      REAL(8) :: vlf,vgf
      REAL(8),DIMENSION(ncell_fluid) :: diag_g,diag_l
      REAL(8),DIMENSION(ncell_fluid_pad,ndim) :: src_g,src_l
      REAL(8),DIMENSION(nf_non) :: off_diag_g_non_i,off_diag_l_non_i
      REAL(8),DIMENSION(nf_nonk) :: off_diag_g_non_k,off_diag_l_non_k
      REAL(8),DIMENSION(nf_nonk+nf_non+nf_inl+nf_adw+nf_fsw+nf_ctw+nf_chw) :: diag_g_nf,diag_l_nf
!      
      IF(iter.eq.1) THEN 
          
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
         ik=nonk_throat(i)      ! k=right_nb_k(ik) --> left_nf(k), right_non(k) 
!         
         diag_g_nf(i0)=0.d0
         diag_l_nf(i0)=0.d0
         off_diag_g_non_i(ir)=0.d0
         off_diag_l_non_i(ir)=0.d0      

         IF(kk.le.ncell_fluid) THEN   
            off_diag_g_non_k(ik)=0.d0
            off_diag_l_non_k(ik)=0.d0 
         ENDIF           
!         
         DO ix=1,ndim
            vlf=flux_l_nf(i1)*xn_nf(i1,ix)/saa_nf(i1)  !xn_nf=vector mode (+ when ii<kk), svp(i1,:)=outward when inlet, 
            vgf=flux_g_nf(i1)*xn_nf(i1,ix)/saa_nf(i1)            
!            
            cf_g=cell%alphag(ii)*cell%eviscosg(ii)*saa_nf(i1)
            cf_l=cell%alphal(ii)*cell%eviscosl(ii)*saa_nf(i1)
            cf_g=cf_g*djir_non(ir)
            cf_l=cf_l*djir_non(ir)
            diag_g(ii)=diag_g(ii)+cf_g
            diag_l(ii)=diag_l(ii)+cf_l
            src_g(ii,ix)=src_g(ii,ix)+cf_g*vgf 
            src_l(ii,ix)=src_l(ii,ix)+cf_l*vlf
            
            IF(kk.le.ncell_fluid) THEN
               ik=nonk_throat(i)      ! k=right_nb_k(ik) --> left_nf(k), right_non(k)        
               k=right_nb_k(ik)
               kk0=right_non(k)
               cf_g=cell%alphag(kk0)*cell%eviscosg(kk0)*saa_nf(ik)
               cf_l=cell%alphal(kk0)*cell%eviscosl(kk0)*saa_nf(ik)
               cf_g=cf_g*djir_non(ik)
               cf_l=cf_l*djir_non(ik)
               diag_g(kk0)=diag_g(kk0)+cf_g
               diag_l(kk0)=diag_l(kk0)+cf_l
               src_g(kk0,ix)=src_g(kk0,ix)+cf_g*vgf 
               src_l(kk0,ix)=src_l(kk0,ix)+cf_l*vlf            
            ENDIF   

         ENDDO
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
         IF(mcp_on(tt).eq.0) CYCLE       !mcp_on: global MCP number  
         DO i=1,num_mcpface(zz)
            i1=n_face_mcp(zz,i) 
            i0=istart0-istart+i1  !i0=istart0+i=istart0+(i1-start) 
            ii=left_nf(i1)  
            ir=i1-istart
            kk=right_non(ir) !i1=istart+i
            ik=nonk_mcp(zz,i)      ! k=right_nb_k(ik) --> left_nf(k), right_non(k) 
!            
            diag_g_nf(i0)=0.d0
            diag_l_nf(i0)=0.d0
            off_diag_g_non_i(ir)=0.d0
            off_diag_l_non_i(ir)=0.d0
            
            IF(kk.le.ncell_fluid) THEN   
               off_diag_g_non_k(ik)=0.d0
               off_diag_l_non_k(ik)=0.d0
            ENDIF   
!            
            DO ix=1,ndim
               vlf=flux_l_nf(i1)*xn_nf(i1,ix)/saa_nf(i1)  !xn_nf=vector mode (+ when ii<kk), svp(i1,:)=outward when inlet, 
               vgf=flux_g_nf(i1)*xn_nf(i1,ix)/saa_nf(i1)            
!            
               cf_g=cell%alphag(ii)*cell%eviscosg(ii)*saa_nf(i1)
               cf_l=cell%alphal(ii)*cell%eviscosl(ii)*saa_nf(i1)
               cf_g=cf_g*djir_non(ir)
               cf_l=cf_l*djir_non(ir)
               diag_g(ii)=diag_g(ii)+cf_g
               diag_l(ii)=diag_l(ii)+cf_l
               src_g(ii,ix)=src_g(ii,ix)+cf_g*vgf 
               src_l(ii,ix)=src_l(ii,ix)+cf_l*vlf
            
               IF(kk.le.ncell_fluid) THEN
                  ik=nonk_mcp(zz,i)      ! k=right_nb_k(ik) --> left_nf(k), right_non(k)        
                  k=right_nb_k(ik)
                  kk0=right_non(k)                     
                  cf_g=cell%alphag(kk0)*cell%eviscosg(kk0)*saa_nf(ik)
                  cf_l=cell%alphal(kk0)*cell%eviscosl(kk0)*saa_nf(ik)
                  cf_g=cf_g*djir_non(ik)
                  cf_l=cf_l*djir_non(ik)
                  diag_g(kk0)=diag_g(kk0)+cf_g
                  diag_l(kk0)=diag_l(kk0)+cf_l
                  src_g(kk0,ix)=src_g(kk0,ix)+cf_g*vgf 
                  src_l(kk0,ix)=src_l(kk0,ix)+cf_l*vlf            
               ENDIF   
!   
            ENDDO            
         ENDDO
      ENDDO  
!
!......valve model
!         
      nv=0
      nf_number=nf_number_id(nv)
      istart=istart_nf(1,nf_number)      
      istart0=istart_nfs(0)         
      DO zz=1,num_valveloc                
         tt=mapping_valve(zz)             
         IF(valve_closed(tt).eq.0) CYCLE   
         DO i=1,num_valveface(zz)
            i1=n_face_valve(zz,i) 
            i0=istart0-istart+i1  !i0=istart0+i=istart0+(i1-start)  
            ii=left_nf(i1) 
            ir=i1-istart                !i1=istart+i 
            kk=right_non(ir)
!
            diag_g_nf(i0)=0.d0
            diag_l_nf(i0)=0.d0
            off_diag_g_non_i(ir)=0.d0
            off_diag_l_non_i(ir)=0.d0  
            IF(kk.le.ncell_fluid) THEN  
               ik=nonk_valve(zz,i) 
               off_diag_g_non_k(ik)=0.d0
               off_diag_l_non_k(ik)=0.d0
            ENDIF               
!            
            DO ix=1,ndim
               vlf=0.d0 
               vgf=0.d0 
!            
               cf_g=0.d0 
               cf_l=0.d0 
               cf_g=cf_g*djir_non(ir)
               cf_l=cf_l*djir_non(ir)
               diag_g(ii)=diag_g(ii)+cf_g
               diag_l(ii)=diag_l(ii)+cf_l
               src_g(ii,ix)=src_g(ii,ix)+cf_g*vgf 
               src_l(ii,ix)=src_l(ii,ix)+cf_l*vlf
            
               IF(kk.le.ncell_fluid) THEN
                  ik=nonk_valve(zz,i)      ! k=right_nb_k(ik) --> left_nf(k), right_non(k)        
                  k=right_nb_k(ik)
                  kk0=right_non(k)                      
                  cf_g=0.d0
                  cf_l=0.d0
                  cf_g=cf_g*djir_non(ik)
                  cf_l=cf_l*djir_non(ik)
                  diag_g(kk0)=diag_g(kk0)+cf_g
                  diag_l(kk0)=diag_l(kk0)+cf_l
                  src_g(kk0,ix)=src_g(kk0,ix)+cf_g*vgf 
                  src_l(kk0,ix)=src_l(kk0,ix)+cf_l*vlf            
               ENDIF   
!   
            ENDDO 
         ENDDO
      ENDDO 
      
      ELSE !iter > 1
       
      IF(.not.choke) GOTO 200         
!
!......choke model
!
!      DO i=1,num_throatface
!         i1=n_face_throat(i)
!         ii=left_nf(i1)   
!         ir=i1-istart               
!         kk=right_non(ir)
!         i0=istart0+ir
!         diag_g_nf(i0)=0.d0
!         diag_l_nf(i0)=0.d0
!         off_diag_g_non_i(ir)=0.d0
!         off_diag_l_non_i(ir)=0.d0         
!         off_diag_g_non_k(ir)=0.d0
!         off_diag_l_non_k(ir)=0.d0
!!         
!         DO ix=1,ndim
!            vlf=flux_l_nf(i1)*xn_nf(i1,ix)/saa_nf(i1)  !xn_nf=vector mode (+ when ii<kk), svp(i1,:)=outward when inlet, 
!            vgf=flux_g_nf(i1)*xn_nf(i1,ix)/saa_nf(i1)            
!!            
!            cf_g=cell%alphag(ii)*cell%eviscosg(ii)*saa_nf(i1)
!            cf_l=cell%alphal(ii)*cell%eviscosl(ii)*saa_nf(i1)
!            cf_g=cf_g*djir_non(ir)
!            cf_l=cf_l*djir_non(ir)
!            diag_g(ii)=diag_g(ii)+cf_g
!            diag_l(ii)=diag_l(ii)+cf_l
!            src_g(ii,ix)=src_g(ii,ix)+cf_g*(vg_n(ii,ix)-vgf)
!            src_l(ii,ix)=src_l(ii,ix)+cf_l*(vl_n(ii,ix)-vlf)
!            
!            IF(kk.le.ncell_fluid) THEN
!               cf_g=cell%alphag(kk)*cell%eviscosg(kk)*saa_nf(i1)
!               cf_l=cell%alphal(kk)*cell%eviscosl(kk)*saa_nf(i1)
!               cf_g=cf_g*djir_non(ir)
!               cf_l=cf_l*djir_non(ir)
!               diag_g(kk)=diag_g(kk)+cf_g
!               diag_l(kk)=diag_l(kk)+cf_l
!               off_diag_g_non_i(ir)=0.d0
!               off_diag_l_non_i(ir)=0.d0
!               src_g(kk,ix)=src_g(kk,ix)+cf_g*(vg_n(ii,ix)-vgf)
!               src_l(kk,ix)=src_l(kk,ix)+cf_l*(vl_n(ii,ix)-vlf)            
!            ENDIF   
!
!         ENDDO
!      ENDDO
!
200   CONTINUE
!
!......mcp model
!      
!      nv=0
!      nf_number=nf_number_id(nv)
!      istart=istart_nf(1,nf_number)      
!      istart0=istart_nfs(0)     
!      DO zz=1,num_mcploc                !local loop
!         tt=mapping_mcp(zz)             !mapping local MCP to global MCP number
!         IF(mcp_on(tt).eq.0) CYCLE       !mcp_on: global MCP number  
!         DO i=1,num_mcpface(zz)
!            i1=n_face_mcp(zz,i) 
!            i0=istart0-istart+i1  !i0=istart0+i=istart0+(i1-start) 
!            ii=left_nf(i1)  
!            ir=i1-istart
!            ik=nonk_mcp(zz,i)      ! k=right_nb_k(ik) --> left_nf(k), right_non(k)        
!            k=right_nb_k(ik)
!            kk=right_non(k) 
!            
!            diag_g_nf(i0)=0.d0
!            diag_l_nf(i0)=0.d0
!            off_diag_g_non_i(ir)=0.d0
!            off_diag_l_non_i(ir)=0.d0         
!            off_diag_g_non_k(ik)=0.d0
!            off_diag_l_non_k(ik)=0.d0            
!!            
!            DO ix=1,ndim
!               vlf=flux_l_nf(i1)*xn_nf(i1,ix)/saa_nf(i1)  !xn_nf=vector mode (+ when ii<kk), svp(i1,:)=outward when inlet, 
!               vgf=flux_g_nf(i1)*xn_nf(i1,ix)/saa_nf(i1)            
!!            
!               diag_g_nf(i0)=0.d0
!               diag_l_nf(i0)=0.d0
!               off_diag_g_non_i(ir)=0.d0
!               off_diag_l_non_i(ir)=0.d0         
!               off_diag_g_non_k(ik)=0.d0
!               off_diag_l_non_k(ik)=0.d0
!!         
!               cf_g=cell%alphag(ii)*cell%eviscosg(ii)*saa_nf(i1)
!               cf_l=cell%alphal(ii)*cell%eviscosl(ii)*saa_nf(i1)
!               cf_g=cf_g*djir_non(ir)
!               cf_l=cf_l*djir_non(ir)
!               diag_g(ii)=diag_g(ii)+cf_g
!               diag_l(ii)=diag_l(ii)+cf_l
!               src_g(ii,ix)=src_g(ii,ix)+cf_g*(vg_n(ii,ix)-vgf) 
!               src_l(ii,ix)=src_l(ii,ix)+cf_l*(vl_n(ii,ix)-vlf)
!            
!               IF(kk.le.ncell_fluid) THEN
!                  cf_g=cell%alphag(kk)*cell%eviscosg(kk)*saa_nf(i1)
!                  cf_l=cell%alphal(kk)*cell%eviscosl(kk)*saa_nf(i1)
!                  cf_g=cf_g*djir_non(ir)
!                  cf_l=cf_l*djir_non(ir)
!                  diag_g(kk)=diag_g(kk)+cf_g
!                  diag_l(kk)=diag_l(kk)+cf_l
!                  src_g(kk,ix)=src_g(kk,ix)+cf_g*(vg_n(ii,ix)-vgf)
!                  src_l(kk,ix)=src_l(kk,ix)+cf_l*(vl_n(ii,ix)-vlf)
!               ENDIF   
!!   
!            ENDDO            
!         ENDDO
!      ENDDO  
!
!......valve model
!         
!      DO zz=1,num_valveloc                
!         tt=mapping_valve(zz)             
!         IF(valve_closed(tt).eq.0) CYCLE   
!         DO i=1,num_valveface(zz)
!            i1=n_face_valve(zz,i) 
!            ii=left_nf(i1) 
!            ir=i1-istart                !i1=istart+i 
!            i0=istart0+ir
!            diag_g_nf(i0)=0.d0
!            diag_l_nf(i0)=0.d0
!            off_diag_g_non_i(ir)=0.d0
!            off_diag_l_non_i(ir)=0.d0         
!            off_diag_g_non_k(ir)=0.d0
!            off_diag_l_non_k(ir)=0.d0            
!!            
!            DO ix=1,ndim
!               vlf=0.d0 
!               vgf=0.d0 
!!            
!               diag_g_nf(i0)=0.d0
!               diag_l_nf(i0)=0.d0
!               off_diag_g_non_i(ir)=0.d0
!               off_diag_l_non_i(ir)=0.d0         
!               off_diag_g_non_k(ir)=0.d0
!               off_diag_l_non_k(ir)=0.d0
!!         
!               cf_g=0.d0 
!               cf_l=0.d0 
!               cf_g=cf_g*djir_non(ir)
!               cf_l=cf_l*djir_non(ir)
!               diag_g(ii)=diag_g(ii)+cf_g
!               diag_l(ii)=diag_l(ii)+cf_l
!               src_g(ii,ix)=src_g(ii,ix)+cf_g*vgf 
!               src_l(ii,ix)=src_l(ii,ix)+cf_l*vlf
!            
!               IF(kk.le.ncell_fluid) THEN
!                  cf_g=0.d0
!                  cf_l=0.d0
!                  cf_g=cf_g*djir_non(ir)
!                  cf_l=cf_l*djir_non(ir)
!                  diag_g(kk)=diag_g(kk)+cf_g
!                  diag_l(kk)=diag_l(kk)+cf_l
!                  src_g(kk,ix)=src_g(kk,ix)+cf_g*(vg_n(ii,ix)-vgf) 
!                  src_l(kk,ix)=src_l(kk,ix)+cf_l*(vl_n(ii,ix)-vlf)            
!               ENDIF   
!!   
!            ENDDO 
!         ENDDO
!      ENDDO           
          
      ENDIF    
!            
      END SUBROUTINE fluxBC_diffusion_smac3
