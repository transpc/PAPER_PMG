!
      SUBROUTINE momentum_diffusion(diff_liq,diff_gas,diff_drp)
!
!     This routine calculates diffusive fluxes through the cell face.
!     Diffusive fluxes are discretized using central differences.
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell            
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: np
      USE Zparam       , ONLY: ndim
      USE Zvec_param   , ONLY: nf_non,nf_inl,nf_adw,nf_fsw,nf_ctw,nf_chw,nf_nonk
      USE Znum_cell    , ONLY: istart_nf,istart_nbcon_nf,right_nb_k, &
                               nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zvec_index   , ONLY: left_nf,right_non,nbcon_nf
      USE Zb_condition , ONLY: alphab_gas,alphab_liq,alphab_drp,vb_gas,vb_liq,vb_drp,vin_gas,vin_liq,vin_drp,v_wall
      USE Zconst1      , ONLY: mdiffscheme,lowreynolds
      USE Zface        , ONLY: twall_model,laminar,Free_slip
      USE Zgradoption  , ONLY: non_orth_diff
      USE Zturb        , ONLY: wvis_liq,wvis_gas,wallnr
      USE Zvector      , ONLY: vg_o,vl_o,vd_o
      USE Zbc_index    , ONLY: vin_norm
      USE Zuserdefined , ONLY: vel_bc_profile_inl
      USE Zmodel       , ONLY:i_droplet
      USE Zturb        , ONLY:s_macroturb_source       
      USE Zvec_geo     , ONLY: xn_nf,sap_nf,sa_nf,             &
                               f1,f0,fac1_non,fac_non,fac_fsw, &
                               sad_non,dnj_non
      USE Zvec_geo     , ONLY: xn_nf
      USE Zmcp        
      USE Zrv_model    , ONLY: rv_mcp,rv_choke,rv_valve
!
      IMPLICIT NONE
!
!.....Input
      REAL(8),DIMENSION(ncell_fluid,ndim) :: diff_liq,diff_gas,diff_drp
!.....Local variables
      INTEGER :: i,k,ix
      INTEGER :: ii,kk
      INTEGER :: nv,nf_number,len,istart0,istart,istart2,i0,i1,i2
      REAL(8) :: vl2,vg2,vd2
      REAL(8) :: sv1
      REAL(8) :: dvli11,dvli12,dvli13,dvli21,dvli22,dvli23,dvli31,dvli32,dvli33
      REAL(8) :: dvgi11,dvgi12,dvgi13,dvgi21,dvgi22,dvgi23,dvgi31,dvgi32,dvgi33
      REAL(8) :: dvlj1,dvlj2,dvlj3
      REAL(8) :: dvgj1,dvgj2,dvgj3
      REAL(8) :: avisli,avisgi,avisdi
      REAL(8) :: vl,vg,vd
      REAL(8) :: f_profile
      REAL(8) :: viscl,viscg
!.....Local arrays
      REAL(8),DIMENSION(ncell_fp,ndim,ndim) :: dvldx,dvgdx
!.....Local vector arrays
      REAL(8),DIMENSION(nf_nonk+nf_non) :: avisli_non,avisgi_non,avisdi_non      
      REAL(8),DIMENSION(nf_nonk+nf_non+nf_inl+nf_adw+nf_fsw+nf_ctw+nf_chw,ndim) :: fluxl_diff_nf,fluxg_diff_nf,fluxd_diff_nf        
!
!.....Only for vv_prob="KBJ". It solves the fluid particle-based two-fluid momentum eq.
!
      IF(mdiffscheme.gt.0) THEN
         CALL momentum_diffusion_user(diff_liq,diff_gas,diff_drp)
         RETURN
      ENDIF
!
!.....Calculate velocity gradient at cell center for non-orthogonal grid
!
      IF(non_orth_diff.eq.1) THEN
         CALL grad_vel(2,vl_o,dvldx,vb_liq,vin_liq)
         CALL grad_vel(1,vg_o,dvgdx,vb_gas,vin_gas)
      ELSEIF(non_orth_diff.eq.2)THEN
         CALL grad_frink_vel(vl_o,dvldx,vb_liq,vin_liq)
         CALL grad_frink_vel(vg_o,dvgdx,vb_gas,vin_gas)
      ENDIF
      IF(np.gt.1) THEN 
         IF(non_orth_diff.eq.1 .or. non_orth_diff.eq.2)THEN
            IF(np.gt.1) CALL communicate_3d(dvldx, &
                                            dvgdx)
         ENDIF
      ENDIF
!
!.....Build summation info for non,inl,adw,fsw,ctw,chw
!
      nf_number_nb=5
      nf_number_id(-1)=-1
      nf_number_id(0)=0
      nf_number_id(1)=2
      nf_number_id(2)=4
      nf_number_id(3)=5
      nf_number_id(4)=6
      nf_number_id(5)=7
      istart_nfs(-1)=0
      istart_nfs(0)=istart_nfs(-1)+nf_nonk
      istart_nfs(1)=istart_nfs(0)+nf_non
      istart_nfs(2)=istart_nfs(1)+nf_inl
      istart_nfs(3)=istart_nfs(2)+nf_adw
      istart_nfs(4)=istart_nfs(3)+nf_fsw
      istart_nfs(5)=istart_nfs(4)+nf_ctw
      lens         =istart_nfs(5)+nf_chw    
!
!.....Interface diffusion flux
!
      nv=0
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)
         avisli_non(i)=f1(i)*cell%alphal(ii)*cell%eviscosl(ii)+f0(i)*cell%alphal(kk)*cell%eviscosl(kk)
         avisgi_non(i)=f1(i)*cell%alphag(ii)*cell%eviscosg(ii)+f0(i)*cell%alphag(kk)*cell%eviscosg(kk)
         avisdi_non(i)=f1(i)*cell%alphad(ii)*cell%eviscosd(ii)+f0(i)*cell%alphad(kk)*cell%eviscosd(kk)
      ENDDO
      DO ix=1,ndim
         DO i=1,len  
            i0=istart0+i
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            vl=vl_o(kk,ix)-vl_o(ii,ix)
            vg=vg_o(kk,ix)-vg_o(ii,ix)
            vd=vd_o(kk,ix)-vd_o(ii,ix)
            avisli=avisli_non(i)
            avisgi=avisgi_non(i)
            avisdi=avisdi_non(i)
            fluxl_diff_nf(i0,ix)=avisli*vl*sad_non(i)
            fluxg_diff_nf(i0,ix)=avisgi*vg*sad_non(i)
            fluxd_diff_nf(i0,ix)=avisdi*vd*sad_non(i)
         ENDDO
      ENDDO
!
      nv=-1
      nf_number=nf_number_id(nv)
      len=istart_nf(2,nf_number)
      DO i=1,len  
         k=right_nb_k(i)
         ii=right_non(k)
         kk=left_nf(k)
         avisli_non(i)=f1(k)*cell%alphal(kk)*cell%eviscosl(kk)+f0(k)*cell%alphal(ii)*cell%eviscosl(ii)
         avisgi_non(i)=f1(k)*cell%alphag(kk)*cell%eviscosg(kk)+f0(k)*cell%alphag(ii)*cell%eviscosg(ii)
         avisdi_non(i)=f1(k)*cell%alphad(kk)*cell%eviscosd(kk)+f0(k)*cell%alphad(ii)*cell%eviscosd(ii)         
      ENDDO      
      DO ix=1,ndim
         DO i=1,len       
            k=right_nb_k(i)
            ii=right_non(k)
            kk=left_nf(k)
            vl=vl_o(kk,ix)-vl_o(ii,ix)
            vg=vg_o(kk,ix)-vg_o(ii,ix)
            vd=vd_o(kk,ix)-vd_o(ii,ix)
            avisli=avisli_non(i)
            avisgi=avisgi_non(i)
            avisdi=avisdi_non(i)
            fluxl_diff_nf(i,ix)=avisli*vl*sad_non(k)
            fluxg_diff_nf(i,ix)=avisgi*vg*sad_non(k)
            fluxd_diff_nf(i,ix)=avisdi*vd*sad_non(k)            
         ENDDO  
      ENDDO      
!                 
      IF(non_orth_diff.gt.0)THEN
!
         nv=0
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(nv)
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)  
         DO i=1,len  
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            avisli_non(i)=f1(i)*cell%alphal(ii)*cell%eviscosl(ii)+f0(i)*cell%alphal(kk)*cell%eviscosl(kk)
            avisgi_non(i)=f1(i)*cell%alphag(ii)*cell%eviscosg(ii)+f0(i)*cell%alphag(kk)*cell%eviscosg(kk)
            avisdi_non(i)=f1(i)*cell%alphad(ii)*cell%eviscosd(ii)+f0(i)*cell%alphad(kk)*cell%eviscosd(kk)
         ENDDO         
         IF(ndim.eq.2)THEN
            DO i=1,len  
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
               avisli=avisli_non(i)
               avisgi=avisgi_non(i)
               avisdi=avisdi_non(i)
               dvli11=fac1_non(i)*dvldx(ii,1,1)+fac_non(i)*dvldx(kk,1,1)
               dvgi11=fac1_non(i)*dvgdx(ii,1,1)+fac_non(i)*dvgdx(kk,1,1)
               dvli12=fac1_non(i)*dvldx(ii,2,1)+fac_non(i)*dvldx(kk,2,1)
               dvgi12=fac1_non(i)*dvgdx(ii,2,1)+fac_non(i)*dvgdx(kk,2,1)
               dvli21=fac1_non(i)*dvldx(ii,1,2)+fac_non(i)*dvldx(kk,1,2)
               dvgi21=fac1_non(i)*dvgdx(ii,1,2)+fac_non(i)*dvgdx(kk,1,2)
               dvli22=fac1_non(i)*dvldx(ii,2,2)+fac_non(i)*dvldx(kk,2,2)
               dvgi22=fac1_non(i)*dvgdx(ii,2,2)+fac_non(i)*dvgdx(kk,2,2)
               dvlj1=dvli11*dnj_non(i,1)+dvli12*dnj_non(i,2)
               dvgj1=dvgi11*dnj_non(i,1)+dvgi12*dnj_non(i,2)
               dvlj2=dvli21*dnj_non(i,1)+dvli22*dnj_non(i,2)
               dvgj2=dvgi21*dnj_non(i,1)+dvgi22*dnj_non(i,2)
               fluxl_diff_nf(i0,1)=fluxl_diff_nf(i0,1)+avisli*dvlj1*sa_nf(i1)
               fluxg_diff_nf(i0,1)=fluxg_diff_nf(i0,1)+avisgi*dvgj1*sa_nf(i1)
               fluxl_diff_nf(i0,2)=fluxl_diff_nf(i0,2)+avisli*dvlj2*sa_nf(i1)
               fluxg_diff_nf(i0,2)=fluxg_diff_nf(i0,2)+avisgi*dvgj2*sa_nf(i1)
            ENDDO
         ELSE
            DO i=1,len  
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)         
               avisli=avisli_non(i)
               avisgi=avisgi_non(i)
               avisdi=avisdi_non(i)
               dvli11=fac1_non(i)*dvldx(ii,1,1)+fac_non(i)*dvldx(kk,1,1)
               dvgi11=fac1_non(i)*dvgdx(ii,1,1)+fac_non(i)*dvgdx(kk,1,1)
               dvli12=fac1_non(i)*dvldx(ii,2,1)+fac_non(i)*dvldx(kk,2,1)
               dvgi12=fac1_non(i)*dvgdx(ii,2,1)+fac_non(i)*dvgdx(kk,2,1)
               dvli13=fac1_non(i)*dvldx(ii,3,1)+fac_non(i)*dvldx(kk,3,1)
               dvgi13=fac1_non(i)*dvgdx(ii,3,1)+fac_non(i)*dvgdx(kk,3,1)
               dvli21=fac1_non(i)*dvldx(ii,1,2)+fac_non(i)*dvldx(kk,1,2)
               dvgi21=fac1_non(i)*dvgdx(ii,1,2)+fac_non(i)*dvgdx(kk,1,2)
               dvli22=fac1_non(i)*dvldx(ii,2,2)+fac_non(i)*dvldx(kk,2,2)
               dvgi22=fac1_non(i)*dvgdx(ii,2,2)+fac_non(i)*dvgdx(kk,2,2)
               dvli23=fac1_non(i)*dvldx(ii,3,2)+fac_non(i)*dvldx(kk,3,2)
               dvgi23=fac1_non(i)*dvgdx(ii,3,2)+fac_non(i)*dvgdx(kk,3,2)
               dvli31=fac1_non(i)*dvldx(ii,1,3)+fac_non(i)*dvldx(kk,1,3)
               dvgi31=fac1_non(i)*dvgdx(ii,1,3)+fac_non(i)*dvgdx(kk,1,3)
               dvli32=fac1_non(i)*dvldx(ii,2,3)+fac_non(i)*dvldx(kk,2,3)
               dvgi32=fac1_non(i)*dvgdx(ii,2,3)+fac_non(i)*dvgdx(kk,2,3)
               dvli33=fac1_non(i)*dvldx(ii,3,3)+fac_non(i)*dvldx(kk,3,3)
               dvgi33=fac1_non(i)*dvgdx(ii,3,3)+fac_non(i)*dvgdx(kk,3,3)
               dvlj1=dvli11*dnj_non(i,1)+dvli12*dnj_non(i,2)+dvli13*dnj_non(i,3)
               dvgj1=dvgi11*dnj_non(i,1)+dvgi12*dnj_non(i,2)+dvgi13*dnj_non(i,3)
               dvlj2=dvli21*dnj_non(i,1)+dvli22*dnj_non(i,2)+dvli23*dnj_non(i,3)
               dvgj2=dvgi21*dnj_non(i,1)+dvgi22*dnj_non(i,2)+dvgi23*dnj_non(i,3)
               dvlj3=dvli31*dnj_non(i,1)+dvli32*dnj_non(i,2)+dvli33*dnj_non(i,3)
               dvgj3=dvgi31*dnj_non(i,1)+dvgi32*dnj_non(i,2)+dvgi33*dnj_non(i,3)
               fluxl_diff_nf(i0,1)=fluxl_diff_nf(i0,1)+avisli*dvlj1*sa_nf(i1)
               fluxg_diff_nf(i0,1)=fluxg_diff_nf(i0,1)+avisgi*dvgj1*sa_nf(i1)
               fluxl_diff_nf(i0,2)=fluxl_diff_nf(i0,2)+avisli*dvlj2*sa_nf(i1)
               fluxg_diff_nf(i0,2)=fluxg_diff_nf(i0,2)+avisgi*dvgj2*sa_nf(i1)
               fluxl_diff_nf(i0,3)=fluxl_diff_nf(i0,3)+avisli*dvlj3*sa_nf(i1)
               fluxg_diff_nf(i0,3)=fluxg_diff_nf(i0,3)+avisgi*dvgj3*sa_nf(i1)
            ENDDO
         ENDIF
!
         nv=-1
         nf_number=nf_number_id(nv)
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)    
         DO i=1,len  
            k=right_nb_k(i)
            ii=right_non(k)
            kk=left_nf(k) 
            avisli_non(i)=f1(k)*cell%alphal(kk)*cell%eviscosl(kk)+f0(k)*cell%alphal(ii)*cell%eviscosl(ii)
            avisgi_non(i)=f1(k)*cell%alphag(kk)*cell%eviscosg(kk)+f0(k)*cell%alphag(ii)*cell%eviscosg(ii)
            avisdi_non(i)=f1(k)*cell%alphad(kk)*cell%eviscosd(kk)+f0(k)*cell%alphad(ii)*cell%eviscosd(ii)         
         ENDDO           
         IF(ndim.eq.2)THEN
            DO i=1,len  
               k=right_nb_k(i)
               kk=right_non(k)
               ii=left_nf(k) 
               avisli=avisli_non(i)
               avisgi=avisgi_non(i)
               avisdi=avisdi_non(i)
               dvli11=fac1_non(k)*dvldx(ii,1,1)+fac_non(k)*dvldx(kk,1,1)
               dvgi11=fac1_non(k)*dvgdx(ii,1,1)+fac_non(k)*dvgdx(kk,1,1)
               dvli12=fac1_non(k)*dvldx(ii,2,1)+fac_non(k)*dvldx(kk,2,1)
               dvgi12=fac1_non(k)*dvgdx(ii,2,1)+fac_non(k)*dvgdx(kk,2,1)
               dvli21=fac1_non(k)*dvldx(ii,1,2)+fac_non(k)*dvldx(kk,1,2)
               dvgi21=fac1_non(k)*dvgdx(ii,1,2)+fac_non(k)*dvgdx(kk,1,2)
               dvli22=fac1_non(k)*dvldx(ii,2,2)+fac_non(k)*dvldx(kk,2,2)
               dvgi22=fac1_non(k)*dvgdx(ii,2,2)+fac_non(k)*dvgdx(kk,2,2)
               dvlj1=dvli11*dnj_non(k,1)+dvli12*dnj_non(k,2)
               dvgj1=dvgi11*dnj_non(k,1)+dvgi12*dnj_non(k,2)
               dvlj2=dvli21*dnj_non(k,1)+dvli22*dnj_non(k,2)
               dvgj2=dvgi21*dnj_non(k,1)+dvgi22*dnj_non(k,2)
               fluxl_diff_nf(i,1)=fluxl_diff_nf(i,1)-avisli*dvlj1*sa_nf(k)
               fluxg_diff_nf(i,1)=fluxg_diff_nf(i,1)-avisgi*dvgj1*sa_nf(k)
               fluxl_diff_nf(i,2)=fluxl_diff_nf(i,2)-avisli*dvlj2*sa_nf(k)
               fluxg_diff_nf(i,2)=fluxg_diff_nf(i,2)-avisgi*dvgj2*sa_nf(k)
            ENDDO
         ELSE
            DO i=1,len  
               k=right_nb_k(i) 
               kk=right_non(k)
               ii=left_nf(k) 
               avisli=avisli_non(i)
               avisgi=avisgi_non(i)
               avisdi=avisdi_non(i)
               dvli11=fac1_non(k)*dvldx(ii,1,1)+fac_non(k)*dvldx(kk,1,1)
               dvgi11=fac1_non(k)*dvgdx(ii,1,1)+fac_non(k)*dvgdx(kk,1,1)
               dvli12=fac1_non(k)*dvldx(ii,2,1)+fac_non(k)*dvldx(kk,2,1)
               dvgi12=fac1_non(k)*dvgdx(ii,2,1)+fac_non(k)*dvgdx(kk,2,1)
               dvli13=fac1_non(k)*dvldx(ii,3,1)+fac_non(k)*dvldx(kk,3,1)
               dvgi13=fac1_non(k)*dvgdx(ii,3,1)+fac_non(k)*dvgdx(kk,3,1)
               dvli21=fac1_non(k)*dvldx(ii,1,2)+fac_non(k)*dvldx(kk,1,2)
               dvgi21=fac1_non(k)*dvgdx(ii,1,2)+fac_non(k)*dvgdx(kk,1,2)
               dvli22=fac1_non(k)*dvldx(ii,2,2)+fac_non(k)*dvldx(kk,2,2)
               dvgi22=fac1_non(k)*dvgdx(ii,2,2)+fac_non(k)*dvgdx(kk,2,2)
               dvli23=fac1_non(k)*dvldx(ii,3,2)+fac_non(k)*dvldx(kk,3,2)
               dvgi23=fac1_non(k)*dvgdx(ii,3,2)+fac_non(k)*dvgdx(kk,3,2)
               dvli31=fac1_non(k)*dvldx(ii,1,3)+fac_non(k)*dvldx(kk,1,3)
               dvgi31=fac1_non(k)*dvgdx(ii,1,3)+fac_non(k)*dvgdx(kk,1,3)
               dvli32=fac1_non(k)*dvldx(ii,2,3)+fac_non(k)*dvldx(kk,2,3)
               dvgi32=fac1_non(k)*dvgdx(ii,2,3)+fac_non(k)*dvgdx(kk,2,3)
               dvli33=fac1_non(k)*dvldx(ii,3,3)+fac_non(k)*dvldx(kk,3,3)
               dvgi33=fac1_non(k)*dvgdx(ii,3,3)+fac_non(k)*dvgdx(kk,3,3)
               dvlj1=dvli11*dnj_non(k,1)+dvli12*dnj_non(k,2)+dvli13*dnj_non(k,3)
               dvgj1=dvgi11*dnj_non(k,1)+dvgi12*dnj_non(k,2)+dvgi13*dnj_non(k,3)
               dvlj2=dvli21*dnj_non(k,1)+dvli22*dnj_non(k,2)+dvli23*dnj_non(k,3)
               dvgj2=dvgi21*dnj_non(k,1)+dvgi22*dnj_non(k,2)+dvgi23*dnj_non(k,3)
               dvlj3=dvli31*dnj_non(k,1)+dvli32*dnj_non(k,2)+dvli33*dnj_non(k,3)
               dvgj3=dvgi31*dnj_non(k,1)+dvgi32*dnj_non(k,2)+dvgi33*dnj_non(k,3)
               fluxl_diff_nf(i,1)=fluxl_diff_nf(i,1)-avisli*dvlj1*sa_nf(k)
               fluxg_diff_nf(i,1)=fluxg_diff_nf(i,1)-avisgi*dvgj1*sa_nf(k)
               fluxl_diff_nf(i,2)=fluxl_diff_nf(i,2)-avisli*dvlj2*sa_nf(k)
               fluxg_diff_nf(i,2)=fluxg_diff_nf(i,2)-avisgi*dvgj2*sa_nf(k)
               fluxl_diff_nf(i,3)=fluxl_diff_nf(i,3)-avisli*dvlj3*sa_nf(k)
               fluxg_diff_nf(i,3)=fluxg_diff_nf(i,3)-avisgi*dvgj3*sa_nf(k)
            ENDDO
         ENDIF               
         
      ENDIF
!
!.....fluxBC: choke model, mcp model, valve model
!      
      IF(rv_valve.eq.1.or.rv_choke.eq.1.or.rv_mcp.eq.1) CALL fluxBC_diffusion_ice(fluxl_diff_nf,fluxg_diff_nf,fluxd_diff_nf)
!
      nv=1
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      len   =istart_nf(2,nf_number)
      DO ix=1,ndim
         DO i=1,len  
            i0=istart0+i
            i1=istart+i
            i2=istart2+i
            ii=left_nf(i1)
            k=nbcon_nf(i2)
            f_profile=vel_bc_profile_inl(i)
            IF(vin_norm(k).eq.0)THEN
              vl2=vb_liq(k,ix)*f_profile
              vg2=vb_gas(k,ix)*f_profile
              vd2=vb_drp(k,ix)*f_profile    
            ELSE
              vl2=vin_liq(k)*xn_nf(i1,ix)*f_profile
              vg2=vin_gas(k)*xn_nf(i1,ix)*f_profile
              vd2=vin_drp(k)*xn_nf(i1,ix)*f_profile
            ENDIF
            sv1=sap_nf(i1)
            avisli=alphab_liq(k)*cell%eviscosl(ii)
            avisgi=alphab_gas(k)*cell%eviscosg(ii)
            avisdi=alphab_drp(k)*cell%eviscosd(ii)
            fluxl_diff_nf(i0,ix)=avisli*(vl2-vl_o(ii,ix))*sv1
            fluxg_diff_nf(i0,ix)=avisgi*(vg2-vg_o(ii,ix))*sv1
            fluxd_diff_nf(i0,ix)=avisdi*(vd2-vd_o(ii,ix))*sv1
         ENDDO   
      ENDDO   
!
      IF(Twall_Model.eq.Laminar)THEN
         DO nv=2,5
            nf_number=nf_number_id(nv)
            istart0=istart_nfs(nv)
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            IF(nf_number.eq.5) THEN
               DO ix=1,ndim
                  DO i=1,len  
                     i0=istart0+i
                     i1=istart+i
                     ii=left_nf(i1)
                     sv1=sap_nf(i1)/fac_fsw(i)
                     fluxl_diff_nf(i0,ix)=cell%alphal(ii)*cell%eviscosl(ii)*(v_wall(ix)-vl_o(ii,ix))*sv1
                     fluxg_diff_nf(i0,ix)=cell%alphag(ii)*cell%eviscosg(ii)*(v_wall(ix)-vg_o(ii,ix))*sv1
                     fluxd_diff_nf(i0,ix)=cell%alphad(ii)*cell%eviscosd(ii)*(v_wall(ix)-vd_o(ii,ix))*sv1
                  ENDDO
               ENDDO
            ELSE
               DO ix=1,ndim
                  DO i=1,len  
                     i0=istart0+i
                     i1=istart+i
                     ii=left_nf(i1)
                     sv1=sap_nf(i1)
                     fluxl_diff_nf(i0,ix)=cell%alphal(ii)*cell%eviscosl(ii)*(v_wall(ix)-vl_o(ii,ix))*sv1
                     fluxg_diff_nf(i0,ix)=cell%alphag(ii)*cell%eviscosg(ii)*(v_wall(ix)-vg_o(ii,ix))*sv1
                     fluxd_diff_nf(i0,ix)=cell%alphad(ii)*cell%eviscosd(ii)*(v_wall(ix)-vd_o(ii,ix))*sv1
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ELSE
!                     
!........LSJ: turbulence on wall is considered by wall friciton model using the hydraulic diameter.
!
         IF(s_macroturb_source.eq.'nakayama'.or.s_macroturb_source.eq.'chandesris') THEN
            DO nv=2,5
               nf_number=nf_number_id(nv)
               istart0=istart_nfs(nv)
               istart=istart_nf(1,nf_number)
               len   =istart_nf(2,nf_number)
               DO ix=1,ndim
                  DO i=1,len  
                     i0=istart0+i
                     i1=istart+i
                     ii=left_nf(i1)
                     sv1=sap_nf(i1)
                     viscl=cell%lviscosl(ii)
                     viscg=cell%lviscosg(ii)
                     fluxl_diff_nf(i0,ix)=cell%alphal(ii)*viscl*(v_wall(ix)-vl_o(ii,ix))*sv1
                     fluxg_diff_nf(i0,ix)=cell%alphag(ii)*viscg*(v_wall(ix)-vg_o(ii,ix))*sv1
                     fluxd_diff_nf(i0,ix)=cell%alphad(ii)*viscl*(v_wall(ix)-vd_o(ii,ix))*sv1
                  ENDDO   
               ENDDO   
            ENDDO   
         ELSE
            IF(lowreynolds.ge.1) THEN
               DO nv=2,5
                  nf_number=nf_number_id(nv)
                  istart0=istart_nfs(nv)
                  istart=istart_nf(1,nf_number)
                  len   =istart_nf(2,nf_number)
                  DO ix=1,ndim
                     DO i=1,len  
                        i0=istart0+i
                        i1=istart+i
                        ii=left_nf(i1)
                        sv1=sa_nf(i1)*wallnr(ii)
                        viscl=cell%eviscosl(ii)
                        viscg=cell%eviscosg(ii)
                        fluxl_diff_nf(i0,ix)=cell%alphal(ii)*viscl*(v_wall(ix)-vl_o(ii,ix))*sv1
                        fluxg_diff_nf(i0,ix)=cell%alphag(ii)*viscg*(v_wall(ix)-vg_o(ii,ix))*sv1
                        fluxd_diff_nf(i0,ix)=cell%alphad(ii)*viscl*(v_wall(ix)-vd_o(ii,ix))*sv1
                     ENDDO   
                  ENDDO   
               ENDDO   
            ELSEIF(Twall_Model.eq.Free_slip)THEN
               DO nv=2,5
                  nf_number=nf_number_id(nv)
                  istart0=istart_nfs(nv)
                  len   =istart_nf(2,nf_number)
                  DO ix=1,ndim
                     DO i=1,len  
                        i0=istart0+i
                        fluxl_diff_nf(i0,ix)=0.0d0
                        fluxg_diff_nf(i0,ix)=0.0d0
                        fluxd_diff_nf(i0,ix)=0.0d0
                     ENDDO   
                  ENDDO   
               ENDDO            
            ELSE
               DO nv=2,5
                  nf_number=nf_number_id(nv)
                  istart0=istart_nfs(nv)
                  istart=istart_nf(1,nf_number)
                  len   =istart_nf(2,nf_number)
                  DO ix=1,ndim
                     DO i=1,len  
                        i0=istart0+i
                        i1=istart+i
                        ii=left_nf(i1)
                        sv1=sa_nf(i1)*wallnr(ii)
                        viscl=wvis_liq(ii)
                        viscg=wvis_gas(ii)
                        fluxl_diff_nf(i0,ix)=cell%alphal(ii)*viscl*(v_wall(ix)-vl_o(ii,ix))*sv1
                        fluxg_diff_nf(i0,ix)=cell%alphag(ii)*viscg*(v_wall(ix)-vg_o(ii,ix))*sv1
                        fluxd_diff_nf(i0,ix)=cell%alphad(ii)*viscl*(v_wall(ix)-vd_o(ii,ix))*sv1
                     ENDDO   
                  ENDDO   
               ENDDO   
            ENDIF
         ENDIF
      ENDIF
!                
      CALL sum_nf_ndim(0,0,ncell_fluid,       &
                       fluxl_diff_nf,diff_liq, &
                       fluxg_diff_nf,diff_gas, &
                       fluxd_diff_nf,diff_drp)    
!
      IF(i_droplet.ge.1) diff_drp(:,:)=0.0d0
!
      END SUBROUTINE momentum_diffusion