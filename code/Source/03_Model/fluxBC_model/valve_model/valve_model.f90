!
      SUBROUTINE valve_model_vst(i,vst)
!      
      USE Zvalve       , ONLY: num_valveloc,valve_closed,num_valveface,mapping_valve,icell_valve,ocell_valve         
!      
      IMPLICIT NONE      
!
      INTEGER i,vst
      INTEGER zz,tt,ii,kk,i0
!    
      DO zz=1,num_valveloc                
         tt=mapping_valve(zz)             
         IF(valve_closed(tt).eq.0) CYCLE   
         DO i0=1,num_valveface(zz)
            ii=icell_valve(zz,i0)
            kk=ocell_valve(zz,i0)
            IF(i.eq.ii.or.i.eq.kk) vst=0
         ENDDO
      ENDDO          
!         
      RETURN
      END SUBROUTINE valve_model_vst
!
!-------------------------------------------------------------------
!    
      SUBROUTINE valve_model_imp_end_diffusion1(cf_g_non,cf_l_non)
!      
      USE Znum_cell    , ONLY: istart_nf,nf_number_id
      USE Zvec_param   , ONLY: nf_non
      USE Zvalve       , ONLY: num_valveloc,valve_closed,num_valveface,n_face_valve,mapping_valve
!      
      IMPLICIT NONE      
!
      INTEGER zz,tt,i1,ir,i,nf_number,nv,istart
      REAL(8),DIMENSION(nf_non) :: cf_g_non,cf_l_non      
!    
      nv=0
      nf_number=nf_number_id(nv)
      istart=istart_nf(1,nf_number)      
      DO zz=1,num_valveloc                
         tt=mapping_valve(zz)             
         IF(valve_closed(tt).eq.0) CYCLE   
         DO i=1,num_valveface(zz)
            i1=n_face_valve(zz,i) 
            ir=i1-istart
            cf_g_non(ir)=0.d0
            cf_l_non(ir)=0.d0
         ENDDO
      ENDDO    
!         
      RETURN
      END SUBROUTINE valve_model_imp_end_diffusion1
!
!-------------------------------------------------------------------
!
      SUBROUTINE valve_model_imp_end_diffusion2(diag_x_nf,off_diag_x_non_i)
!      
      USE Vol_DATA     , ONLY: cell
      USE Znum_cell    , ONLY: istart_nf,nf_number_id,istart_nfs
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Zvec_geo     , ONLY: fac1_non,fac_non,sap_nf
      USE Zvec_param   , ONLY: nf_non,nf_inl
      USE Zvalve       , ONLY: num_valveloc,valve_closed,num_valveface,n_face_valve,mapping_valve
!      
      IMPLICIT NONE      
!
      INTEGER zz,tt,i1,ir,i,nf_number,nv,istart,istart0,i0,ii,kk
      REAL(8) cf_x
      REAL(8),DIMENSION(nf_non+nf_inl) :: diag_x_nf 
      REAL(8),DIMENSION(nf_non) :: off_diag_x_non_i
!    
      nv=0
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      DO zz=1,num_valveloc                
         tt=mapping_valve(zz)             
         IF(valve_closed(tt).eq.0) CYCLE   
         DO i=1,num_valveface(zz)
            i1=n_face_valve(zz,i)
            i0=istart0-istart+i1    !i0=istart0+i=istart0+(i1-istart)  
            ir=i1-istart
            ii=left_nf(i1)
            kk=right_non(ir)
!            
            cf_x =( fac1_non(ir)*cell%alphag_o(ii)*cell%rhog(ii)*cell%mdiff(ii) &
                   +fac_non(ir) *cell%alphag_o(kk)*cell%rhog(kk)*cell%mdiff(kk) )*sap_nf(i1)
            diag_x_nf(i0)=cf_x
            off_diag_x_non_i(ir)=off_diag_x_non_i(ir)+diag_x_nf(i0)
            diag_x_nf(i0)=0.d0
!            
         ENDDO
      ENDDO
!         
      RETURN
      END SUBROUTINE valve_model_imp_end_diffusion2
!
!-------------------------------------------------------------------
!
      SUBROUTINE valve_model_imp_end_diffusion3(src_x_nf,src_g_nf)
!      
      USE Vol_DATA     , ONLY: cell
      USE Znum_cell    , ONLY: istart_nf,nf_number_id,istart_nfs
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Zvec_geo     , ONLY: fac1_non,fac_non,sap_nf
      USE Zvec_param   , ONLY: nf_non,nf_inl,nf_fsw,nf_ctw,nf_chw
      USE Zvalve       , ONLY: num_valveloc,valve_closed,num_valveface,n_face_valve,mapping_valve
!      
      IMPLICIT NONE      
!
      INTEGER zz,tt,i1,ir,i,nf_number,nv,istart,istart0,i0,ii,kk
      REAL(8) cf_xe
      REAL(8),DIMENSION(nf_non+nf_inl) :: src_x_nf 
      REAL(8),DIMENSION(nf_non+nf_inl+nf_fsw+nf_ctw+nf_chw) :: src_g_nf
!    
      nv=0
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      DO zz=1,num_valveloc                
         tt=mapping_valve(zz)             
         IF(valve_closed(tt).eq.0) CYCLE   
         DO i=1,num_valveface(zz)
            i1=n_face_valve(zz,i)
            i0=istart0-istart+i1    !i0=istart0+i=istart0+(i1-istart)  
            ir=i1-istart
            ii=left_nf(i1)
            kk=right_non(ir)
!            
            cf_xe= fac1_non(ir)*cell%alphag_o(ii)*cell%rhog(ii)*cell%mdiff(ii)*(cell%ha(ii)-cell%hg(ii)) &
                  +fac_non(ir) *cell%alphag_o(kk)*cell%rhog(kk)*cell%mdiff(kk)*(cell%ha(kk)-cell%hg(kk))
            cf_xe=cf_xe*sap_nf(i1)            

            src_x_nf(i0)=0.d0 
            src_g_nf(i0)=src_g_nf(i0)-cf_xe*(cell%quala_o(kk)-cell%quala_o(ii))            
!            
         ENDDO
      ENDDO
!         
      RETURN
      END SUBROUTINE valve_model_imp_end_diffusion3
!
!-------------------------------------------------------------------
!
      SUBROUTINE valve_model_pressure_matrix1(svoli_l_non,svoli_g_non,svoli_d_non,vl_nf,vg_nf,vd_nf)
!      
      USE Zconst2      , ONLY: dt
      USE Zparam       , ONLY: ndim
      USE Zpress       , ONLY: p,dpdx
      USE Znum_cell    , ONLY: istart_nf
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Zvec_geo     , ONLY: fac1_non,fac_non,dji_x_nf,djia_nf,xn_nf
      USE Zvec_param   , ONLY: nf_flux,nf_non
      USE Zvalve       , ONLY: num_valveloc,valve_closed,num_valveface,n_face_valve,mapping_valve
!      
      IMPLICIT NONE      
!
      INTEGER zz,tt,i1,ir,i,nf_number,istart,ii,kk
      REAL(8) dp,dpi,tmp
      REAL(8),DIMENSION(nf_non) :: svoli_l_non,svoli_g_non,svoli_d_non      
      REAL(8),DIMENSION(nf_flux,ndim) :: vl_nf,vg_nf,vd_nf
!    
      nf_number=0
      istart=istart_nf(1,nf_number)
      DO zz=1,num_valveloc                
         tt=mapping_valve(zz)             
         IF(valve_closed(tt).eq.0) CYCLE   
         DO i=1,num_valveface(zz)
            i1=n_face_valve(zz,i) 
            ir=i1-istart
            ii=left_nf(i1)
            kk=right_non(ir)
            dp =p(kk)-p(ii)
            dpi= (fac1_non(ir)*dpdx(ii,1)+fac_non(ir)*dpdx(kk,1))*dji_x_nf(i1,1) &
                +(fac1_non(ir)*dpdx(ii,2)+fac_non(ir)*dpdx(kk,2))*dji_x_nf(i1,2)
            IF(ndim.eq.3) dpi=dpi+(fac1_non(ir)*dpdx(ii,3)+fac_non(ir)*dpdx(kk,3))*dji_x_nf(i1,3)
            tmp=dt/djia_nf(i1)            

            vl_nf(i1,1)=vl_nf(i1,1)+tmp*svoli_l_non(ir)*(dp-dpi)*xn_nf(i1,1)
            vg_nf(i1,1)=vg_nf(i1,1)+tmp*svoli_g_non(ir)*(dp-dpi)*xn_nf(i1,1)
            vd_nf(i1,1)=vd_nf(i1,1)+tmp*svoli_d_non(ir)*(dp-dpi)*xn_nf(i1,1)
            vl_nf(i1,2)=vl_nf(i1,2)+tmp*svoli_l_non(ir)*(dp-dpi)*xn_nf(i1,2)
            vg_nf(i1,2)=vg_nf(i1,2)+tmp*svoli_g_non(ir)*(dp-dpi)*xn_nf(i1,2)
            vd_nf(i1,2)=vd_nf(i1,2)+tmp*svoli_d_non(ir)*(dp-dpi)*xn_nf(i1,2)
            IF(ndim.eq.3) THEN
               vl_nf(i1,3)=vl_nf(i1,3)+tmp*svoli_l_non(ir)*(dp-dpi)*xn_nf(i1,3)
               vg_nf(i1,3)=vg_nf(i1,3)+tmp*svoli_g_non(ir)*(dp-dpi)*xn_nf(i1,3)
               vd_nf(i1,3)=vd_nf(i1,3)+tmp*svoli_d_non(ir)*(dp-dpi)*xn_nf(i1,3)            
            ENDIF   
!            
         ENDDO
      ENDDO           
!         
      RETURN
      END SUBROUTINE valve_model_pressure_matrix1
!
!-------------------------------------------------------------------
!    
      SUBROUTINE valve_model_pressure_matrix2(poiss_diag_nf,poiss_non_i,poiss_non_k)
!      
      USE Zbc_index    , ONLY: npb
      USE Znum_cell    , ONLY: istart_nf,nf_number_id,istart_nfs
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Zvec_param   , ONLY: nf_fluxk2,nf_non,nf_nonk
      USE Zvalve       , ONLY: num_valveloc,valve_closed,num_valveface,n_face_valve,mapping_valve,nonk_valve
      USE Zzone        , ONLY: ncell_fluid
!      
      IMPLICIT NONE      
!
      INTEGER zz,tt,i1,ir,i,nf_number,istart,ii,kk,nv,istart0,i0,ik
      REAL(8),DIMENSION(nf_fluxk2) :: poiss_diag_nf
      REAL(8),DIMENSION(nf_non) :: poiss_non_i
      REAL(8),DIMENSION(nf_nonk) :: poiss_non_k      
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
            ir=i1-istart
            ii=left_nf(i1)
!            
            IF(npb(ii).eq.0) THEN
               poiss_diag_nf(i0)=0.d0
               poiss_non_i(ir)=0.d0
            ENDIF
!            
!           Asymmetric face 
            kk=right_non(i1-istart) !i1=istart+i            
            IF(kk.le.ncell_fluid) THEN
               IF(npb(kk).eq.0) THEN
                  ik=nonk_valve(zz,i)      ! k=right_nb_k(ik) --> left_nf(k), right_non(k)  
                  poiss_diag_nf(ik)=0.d0
                  poiss_non_k(ik)=0.d0
               ENDIF
            ENDIF
         ENDDO
      ENDDO            
!         
      RETURN
      END SUBROUTINE valve_model_pressure_matrix2
!
!-------------------------------------------------------------------
! 
      SUBROUTINE valve_model_pressure_solve(iflag)
!      
      USE Zconst2      , ONLY: dt
      USE Znum_cell    , ONLY: istart_nf
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Zvalve       , ONLY: num_valveloc,valve_closed,num_valveface,n_face_valve,mapping_valve
      USE Zvec_major   , ONLY: flux_l_nf,flux_g_nf,flux_d_nf
      USE Zvec_geo     , ONLY: xn_nf,sap_nf,dji_x_nf,djia_nf,f0,f1,fac1_non,fac_non,xloc_m_non_i,xloc_m_non_k
      USE Zparam       , ONLY: ndim
      USE Zpress       , ONLY: pp,dpdx
      USE Zpress_coeff , ONLY: coefp_l,coefp_g,coefp_d
!      
      IMPLICIT NONE      
!
      INTEGER iflag
      INTEGER zz,tt,i1,ir,i,nf_number,istart,ii,kk
      REAL(8) :: dp,cf,dpx,dpi,dpj
      REAL(8) :: svolli,svolgi,svoldi
      REAL(8) :: grdPx1,grdPx2,grdPx3
!   
      nf_number=0
      istart=istart_nf(1,nf_number)      
      IF(iflag.eq.0) THEN !non_orth.eq.0
         DO zz=1,num_valveloc                
            tt=mapping_valve(zz)             
            IF(valve_closed(tt).eq.0) CYCLE   
            DO i=1,num_valveface(zz)
               i1=n_face_valve(zz,i) 
               ir=i1-istart
               ii=left_nf(i1)
               kk=right_non(ir)
!               
               cf=dt*sap_nf(i1)    
               dp=pp(kk)-pp(ii)
!               
               svolli=f1(ir)*coefp_l(ii)+f0(ir)*coefp_l(kk)
               svolgi=f1(ir)*coefp_g(ii)+f0(ir)*coefp_g(kk)
               svoldi=f1(ir)*coefp_d(ii)+f0(ir)*coefp_d(kk)
               flux_l_nf(i1)=flux_l_nf(i1)+cf*svolli*dp
               flux_g_nf(i1)=flux_g_nf(i1)+cf*svolgi*dp
               flux_d_nf(i1)=flux_d_nf(i1)+cf*svoldi*dp               
            ENDDO
         ENDDO   
!         
      ELSEIF(iflag.eq.1) THEN !non_orth.eq.1
         DO zz=1,num_valveloc                
            tt=mapping_valve(zz)             
            IF(valve_closed(tt).eq.0) CYCLE   
            DO i=1,num_valveface(zz)
               i1=n_face_valve(zz,i) 
               ir=i1-istart
               ii=left_nf(i1)
               kk=right_non(ir)
!               
               cf=dt*sap_nf(i1)    
               dp=pp(kk)-pp(ii)
               grdPx1=0.5d0*(dpdx(ii,1)+dpdx(kk,1))
               grdPx2=0.5d0*(dpdx(ii,2)+dpdx(kk,2))
               IF(ndim.eq.3) grdPx3=0.5d0*(dpdx(ii,3)+dpdx(kk,3))                   
               dpx=grdPx1*xn_nf(i1,1)+grdPx2*xn_nf(i1,2)
               IF(ndim.eq.3) dpx=dpx+grdPx3*xn_nf(i1,3)
               dp=dp+dpx*djia_nf(i1)
               dpx=grdPx1*dji_x_nf(i1,1)+grdPx2*dji_x_nf(i1,2)
               IF(ndim.eq.3) dpx=dpx+grdPx3*dji_x_nf(i1,3)
               dp=dp-dpx
!               
               svolli=fac1_non(ir)*coefp_l(ii)+fac_non(ir)*coefp_l(kk)
               svolgi=fac1_non(ir)*coefp_g(ii)+fac_non(ir)*coefp_g(kk)
               svoldi=fac1_non(ir)*coefp_d(ii)+fac_non(ir)*coefp_d(kk)
               flux_l_nf(i1)=flux_l_nf(i1)+cf*svolli*dp
               flux_g_nf(i1)=flux_g_nf(i1)+cf*svolgi*dp
               flux_d_nf(i1)=flux_d_nf(i1)+cf*svoldi*dp              
            ENDDO
         ENDDO       
!         
      ELSEIF(iflag.eq.2) THEN !non_orth.eq.2
         DO zz=1,num_valveloc                
            tt=mapping_valve(zz)             
            IF(valve_closed(tt).eq.0) CYCLE   
            DO i=1,num_valveface(zz)
               i1=n_face_valve(zz,i) 
               ir=i1-istart
               ii=left_nf(i1)
               kk=right_non(ir)
!               
               cf=dt*sap_nf(i1)    
               dp=pp(kk)-pp(ii)
               dpi=dpdx(ii,1)*xloc_m_non_i(ir,1)+dpdx(ii,2)*xloc_m_non_i(ir,2)
               IF(ndim.eq.3) dpi=dpi+dpdx(ii,3)*xloc_m_non_i(ir,3)
               dpj=dpdx(kk,1)*xloc_m_non_k(ir,1)+dpdx(kk,2)*xloc_m_non_k(ir,2)
               IF(ndim.eq.3) dpj=dpj+dpdx(kk,3)*xloc_m_non_k(ir,3)
               dp=dp+(dpj-dpi)
!               
               svolli=fac1_non(ir)*coefp_l(ii)+fac_non(ir)*coefp_l(kk)
               svolgi=fac1_non(ir)*coefp_g(ii)+fac_non(ir)*coefp_g(kk)
               svoldi=fac1_non(ir)*coefp_d(ii)+fac_non(ir)*coefp_d(kk)
               flux_l_nf(i1)=flux_l_nf(i1)+cf*svolli*dp
               flux_g_nf(i1)=flux_g_nf(i1)+cf*svolgi*dp
               flux_d_nf(i1)=flux_d_nf(i1)+cf*svoldi*dp              
            ENDDO
         ENDDO                  
          
      ENDIF
!         
      RETURN
      END SUBROUTINE valve_model_pressure_solve
!
!-------------------------------------------------------------------
!
      SUBROUTINE valve_model_scalar_energy_convection
!      
      USE Znum_cell    , ONLY: istart_nf
      USE Zvalve       , ONLY: num_valveloc,valve_closed,num_valveface,n_face_valve,mapping_valve
      USE Zvec_major   , ONLY: ecnvc_l_nf,ecnvc_g_nf,ecnvc_d_nf
!      
      IMPLICIT NONE      
!
      INTEGER zz,tt,i1,i,nf_number,istart
!    
      nf_number=0
      istart=istart_nf(1,nf_number)      
      DO zz=1,num_valveloc                
         tt=mapping_valve(zz)             
         IF(valve_closed(tt).eq.0) CYCLE   
         DO i=1,num_valveface(zz)
            i1=n_face_valve(zz,i) 
            ecnvc_l_nf(i1)=0.d0  
            ecnvc_g_nf(i1)=0.d0
            ecnvc_d_nf(i1)=0.d0    
         ENDDO
      ENDDO  
!         
      RETURN
      END SUBROUTINE valve_model_scalar_energy_convection
!
!-------------------------------------------------------------------
!
      SUBROUTINE valve_model_scalar_energy_diffusion(fluxl_diff_nf,fluxg_diff_nf)
!      
      USE Znum_cell    , ONLY: istart_nf,nf_number_id,istart_nfs
      USE Zvalve       , ONLY: num_valveloc,valve_closed,num_valveface,n_face_valve,mapping_valve
      USE Zvec_param   , ONLY: nf_non,nf_inl,nf_fsw,nf_ctw,nf_chw      
!      
      IMPLICIT NONE      
!
      INTEGER zz,tt,i1,i,nf_number,istart,nv,istart0,i0
      REAL(8),DIMENSION(nf_non+nf_inl+nf_fsw+nf_ctw+nf_chw) :: fluxl_diff_nf,fluxg_diff_nf      
!    
      nv=0
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)          
      DO zz=1,num_valveloc                
         tt=mapping_valve(zz)             
         IF(valve_closed(tt).eq.0) CYCLE   
         DO i=1,num_valveface(zz)
            i1=n_face_valve(zz,i) 
            i0=istart0-istart+i1  !i0=istart0+i=istart0+(i1-start) 
            fluxl_diff_nf(i0)=0.d0  
            fluxg_diff_nf(i0)=0.d0 
         ENDDO
      ENDDO            
!         
      RETURN
      END SUBROUTINE valve_model_scalar_energy_diffusion
!
!-------------------------------------------------------------------
!    
      SUBROUTINE valve_model_scalar_mass_convection
!      
      USE Znum_cell    , ONLY: istart_nf
      USE Zvalve       , ONLY: num_valveloc,valve_closed,num_valveface,n_face_valve,mapping_valve
      USE Zvec_major   , ONLY: liq_conv_nf,vap_conv_nf,drp_conv_nf
!      
      IMPLICIT NONE      
!
      INTEGER zz,tt,i1,i,nf_number,istart
!    
      nf_number=0
      istart=istart_nf(1,nf_number)      
      DO zz=1,num_valveloc                
         tt=mapping_valve(zz)             
         IF(valve_closed(tt).eq.0) CYCLE   
         DO i=1,num_valveface(zz)
            i1=n_face_valve(zz,i) 
            liq_conv_nf(i1)=0.d0  
            vap_conv_nf(i1)=0.d0
            drp_conv_nf(i1)=0.d0
         ENDDO
      ENDDO               
!         
      RETURN
      END SUBROUTINE valve_model_scalar_mass_convection
!
!-------------------------------------------------------------------
! 
      SUBROUTINE valve_model_scalar_mass_diffusion(fluxg_diff_nf,efluxg_diff_nf)
!      
      USE Znum_cell    , ONLY: istart_nf,nf_number_id,istart_nfs
      USE Zvalve       , ONLY: num_valveloc,valve_closed,num_valveface,n_face_valve,mapping_valve
      USE Zvec_param   , ONLY: nf_non,nf_inl
!      
      IMPLICIT NONE      
!
      INTEGER zz,tt,i1,i,nf_number,istart,nv,istart0,i0
      REAL(8) :: fluxg_diff_nf(nf_non+nf_inl)
      REAL(8) :: efluxg_diff_nf(nf_non+nf_inl)      
!    
      nv=0
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)         
      DO zz=1,num_valveloc                
         tt=mapping_valve(zz)             
         IF(valve_closed(tt).eq.0) CYCLE   
         DO i=1,num_valveface(zz)
            i1=n_face_valve(zz,i) 
            i0=istart0-istart+i1  !i0=istart0+i=istart0+(i1-start) 
            fluxg_diff_nf(i0)=0.d0  
            efluxg_diff_nf(i0)=0.d0  
         ENDDO
      ENDDO            
!         
      RETURN
      END SUBROUTINE valve_model_scalar_mass_diffusion
!
!-------------------------------------------------------------------
!    
      SUBROUTINE valve_model_scalar_work_convection
!      
      USE Znum_cell    , ONLY: istart_nf
      USE Zvalve       , ONLY: num_valveloc,valve_closed,num_valveface,n_face_valve,mapping_valve
      USE Zvec_major   , ONLY: al_conv_nf,ad_conv_nf,void_conv_nf
!      
      IMPLICIT NONE      
!
      INTEGER zz,tt,i1,i,nf_number,istart
!    
      nf_number=0
      istart=istart_nf(1,nf_number)      
      DO zz=1,num_valveloc                
         tt=mapping_valve(zz)             
         IF(valve_closed(tt).eq.0) CYCLE   
         DO i=1,num_valveface(zz)
            i1=n_face_valve(zz,i) 
            void_conv_nf(i1)=0.d0  
            al_conv_nf(i1)=0.d0
            ad_conv_nf(i1)=0.d0  
         ENDDO
      ENDDO               
!         
      RETURN
      END SUBROUTINE valve_model_scalar_work_convection
!
!-------------------------------------------------------------------
!    
      SUBROUTINE valve_model_scalar_xn_convection
!      
      USE Znum_cell    , ONLY: istart_nf
      USE Zvalve       , ONLY: num_valveloc,valve_closed,num_valveface,n_face_valve,mapping_valve
      USE Zvec_major   , ONLY: quala_conv_nf
!      
      IMPLICIT NONE      
!
      INTEGER zz,tt,i1,i,nf_number,istart
!    
      nf_number=0
      istart=istart_nf(1,nf_number)      
      DO zz=1,num_valveloc                
         tt=mapping_valve(zz)             
         IF(valve_closed(tt).eq.0) CYCLE   
         DO i=1,num_valveface(zz)
            i1=n_face_valve(zz,i) 
            quala_conv_nf(i1)=0.d0  
         ENDDO
      ENDDO               
!         
      RETURN
      END SUBROUTINE valve_model_scalar_xn_convection    
!
!-------------------------------------------------------------------
! 
      SUBROUTINE valve_model_grad_press_ls(dx_non,f_non)
!      
      USE Znum_cell    , ONLY: istart_nf
      USE Zparam       , ONLY: ndim 
      USE Zvalve       , ONLY: num_valveloc,valve_closed,num_valveface,n_face_valve,mapping_valve
      USE Zvec_param   , ONLY: nf_non
      USE Zvec_index   , ONLY: left_nf,right_non
!      
      IMPLICIT NONE      
!
      INTEGER zz,tt,i1,i,nf_number,istart,ii,kk,ir
      REAL(8) ds_non
      REAL(8) :: dx_non(nf_non,ndim)
      REAL(8) :: f_non(nf_non,ndim)   
!    
      nf_number=0
      istart=istart_nf(1,nf_number)        
      DO zz=1,num_valveloc                
         tt=mapping_valve(zz)             
         IF(valve_closed(tt).eq.0) CYCLE   
         DO i=1,num_valveface(zz)
            i1=n_face_valve(zz,i) 
            ii=left_nf(i1)
            ir=i1-istart            
            kk=right_non(ir)
            ds_non=0.d0 !s(kk)-s(ii)
            dx_non(ir,1)=0.d0 !xloc(kk,1)-xloc(ii,1)
            dx_non(ir,2)=0.d0 !xloc(kk,2)-xloc(ii,2)
            If(ndim.eq.3) dx_non(ir,3)=0.d0 !xloc(kk,3)-xloc(ii,3)
            f_non(ir,1)=dx_non(ir,1)*ds_non
            f_non(ir,2)=dx_non(ir,2)*ds_non
            If(ndim.eq.3) f_non(ir,3)=dx_non(ir,3)*ds_non  
         ENDDO
      ENDDO            
!         
      RETURN
      END SUBROUTINE valve_model_grad_press_ls  
!
!-------------------------------------------------------------------
! 
      SUBROUTINE valve_model_grad_press_ls_c_2d(b_non1,b_non2)
!      
      USE Znum_cell    , ONLY: istart_nf
      USE Zparam       , ONLY: ndim 
      USE Zvalve       , ONLY: num_valveloc,valve_closed,num_valveface,n_face_valve,mapping_valve
      USE Zvec_param   , ONLY: nf_non
      USE Zvec_index   , ONLY: left_nf,right_non
!      
      IMPLICIT NONE      
!
      INTEGER zz,tt,i1,i,nf_number,istart,ii,kk,ir
      REAL(8) dx1,dx2
      REAL(8) :: b_non1(nf_non,ndim)   
      REAL(8) :: b_non2(nf_non,ndim:ndim)      
!    
      nf_number=0
      istart=istart_nf(1,nf_number)   
      DO zz=1,num_valveloc                
         tt=mapping_valve(zz)             
         IF(valve_closed(tt).eq.0) CYCLE   
         DO i=1,num_valveface(zz)
            i1=n_face_valve(zz,i) 
            ii=left_nf(i1)
            ir=i1-istart            
            kk=right_non(ir)
            dx1=0.d0 !xloc(kk,1)-xloc(ii,1)
            dx2=0.d0 !xloc(kk,2)-xloc(ii,2)
            b_non1(ir,1)=dx1*dx1
            b_non1(ir,2)=dx1*dx2
            b_non2(ir,2)=dx2*dx2
         ENDDO
      ENDDO  
!         
      RETURN
      END SUBROUTINE valve_model_grad_press_ls_c_2d
!
!-------------------------------------------------------------------
! 
      SUBROUTINE valve_model_grad_press_ls_c_3d(b_non1,b_non2,b_non3)
!      
      USE Znum_cell    , ONLY: istart_nf
      USE Zparam       , ONLY: ndim 
      USE Zvalve       , ONLY: num_valveloc,valve_closed,num_valveface,n_face_valve,mapping_valve
      USE Zvec_param   , ONLY: nf_non
      USE Zvec_index   , ONLY: left_nf,right_non
!      
      IMPLICIT NONE      
!
      INTEGER zz,tt,i1,i,nf_number,istart,ii,kk,ir
      REAL(8) dx1,dx2,dx3
      REAL(8) :: b_non1(nf_non,ndim)
      REAL(8) :: b_non2(nf_non,2:ndim)
      REAL(8) :: b_non3(nf_non,ndim:ndim)      
!    
      nf_number=0
      istart=istart_nf(1,nf_number)   
      DO zz=1,num_valveloc                
         tt=mapping_valve(zz)             
         IF(valve_closed(tt).eq.0) CYCLE   
         DO i=1,num_valveface(zz)
            i1=n_face_valve(zz,i) 
            ii=left_nf(i1)
            ir=i1-istart
            kk=right_non(ir)
            dx1=0.d0 !xloc(kk,1)-xloc(ii,1)
            dx2=0.d0 !xloc(kk,2)-xloc(ii,2)
            dx3=0.d0 !xloc(kk,3)-xloc(ii,3)
            b_non1(ir,1)=dx1*dx1
            b_non1(ir,2)=dx1*dx2
            b_non1(ir,3)=dx1*dx3
            b_non2(ir,2)=dx2*dx2
            b_non2(ir,3)=dx2*dx3
            b_non3(ir,3)=dx3*dx3
         ENDDO
      ENDDO             
!
      RETURN
      END SUBROUTINE valve_model_grad_press_ls_c_3d