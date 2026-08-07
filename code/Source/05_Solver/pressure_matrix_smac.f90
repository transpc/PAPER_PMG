!
      SUBROUTINE pressure_matrix_smac(poiss_diag,poiss_non_i,poiss_non_k,src)
!                                                                                                                                   
!     This routine sets the matrix coefficients for pressure calculation with smac option
!     flux_vol is the volume flux through the face.
!     poiss_diag - diagonal component
!     poiss - off diagonal component
!     src - source vector
!
      USE Zinterface
      USE VOL_DATA      , ONLY: cell
      USE Zzone         , ONLY: ncell_fluid
      USE Zparam        , ONLY: ndim
      USE Zvec_param    , ONLY: nf_flux1,nf_fluxk2,nf_flux,nf_nonk,nf_non,nf_mcc,nf_inl
      USE Znum_cell     , ONLY: istart_nf,istart_nbcon_nf,                 &
                                nf_number_nb,lens,nf_number_id,istart_nfs, &
                                right_nb_k
      USE Zvec_index    , ONLY: left_nf,right_non,nbcon_nf
      USE Zbc_index     , ONLY: npb
      USE Zconst2       , ONLY: dt
      USE Zb_condition  , ONLY: vb_gas,vb_liq,vb_drp,vin_gas,vin_liq,vin_drp
      USE Zconst1       , ONLY: fric_face,nd_face
      USE Zdel_scalar   , ONLY: del_eg,del_el,del_x, &
                                smac3_pres_eng
      USE Zgradoption   , ONLY: irc_damp
      USE Zmars         , ONLY: n_marsbc
      USE Zporous       , ONLY: vfporous,                  &
                                tm_mas_l,tm_mas_g,         &
                                vd_mas_l,vd_mas_g,         &
                                mixing_vane_l,             &
                                l_subchannel,l_mixing_vane
      USE Zpress        , ONLY: p,dpdx
      USE Zpress_coeff  , ONLY: coefp_g,coefp_l,coefp_d,coefm_g,coefm_l
      USE Zqvol         , ONLY: gamma,gamma_wall,h_ig,h_il
      USE Ztimecon      , ONLY: smac,iso_thermal,repeat_smac 
      USE Zuserdefined  , ONLY: vel_bc_profile_inl
      USE Zvector       , ONLY: vl_n,vg_n,vd_n,vl_f_non,vg_f_non
      USE c3com_cupid   , ONLY: i3invtbl,c3dpv
      USE Zbc_index     , ONLY: vin_norm
      USE Zvec_major    , ONLY: flux_l_nf,flux_g_nf,flux_d_nf, &
                                mflux_l_nf,mflux_g_nf,mflux_d_nf
      USE Zvec_scalar   , ONLY: arli_nf,argi_nf,ardi_nf
      USE Zvec_geo      , ONLY: xn_nf,sv_nf,svp_nf,      &
                                sap_nf,dji_x_nf,         &
                                f0,f1,fac1_non,fac_non,  &
                                perm_non,perm_out
      !OPR1000 rod-scale (EVVD)
      USE Zcoord3       , ONLY: volpr,        &
                                volp,porosity
      USE Zrv_model     , ONLY: rv_choke,rv_mcp,rv_valve
!
      IMPLICIT NONE 
!      
      INCLUDE '../10_LinkToMARS/c3com.h' 
! 
!.....Output
      REAL(8),DIMENSION(nf_non) :: poiss_non_i
      REAL(8),DIMENSION(nf_nonk) :: poiss_non_k
      REAL(8),DIMENSION(ncell_fluid) :: poiss_diag,src
!.....Local variables
      INTEGER :: i,k
      INTEGER :: ii,kk,mm
      INTEGER :: idx
      INTEGER :: nv,nf_number,len,istart0,istart,istart2,i0,i1,i2
      LOGICAL,SAVE :: flag_iso_thermal,face_body_face
      REAL(8) :: al,ag_rhog,al_rhol 
      REAL(8) :: arli,argi,ardi 
      REAL(8) :: dp_fluxl,dp_fluxg,dp_fluxd 
      REAL(8) :: del_dtvde,del_dtvdx,del_dtlde
      REAL(8) :: del_tsg
      REAL(8) :: del_tsl,del_dtvdp,del_dtldp
      REAL(8) :: c_if,c_ig,del_sv,PsP 
      REAL(8) :: svoli_l,svoli_g,svoli_d
      REAL(8) :: tmp,cf,dp
      REAL(8) :: vli1,vli2,vli3
      REAL(8) :: vgi1,vgi2,vgi3
      REAL(8) :: vdi1,vdi2,vdi3
      REAL(8) :: f_profile      
      REAL(8) :: a,b,c
      REAL(8) :: ddhfg
      REAL(8) :: poiss_diag_mcc_s
!.....Local arrays
      REAL(8),DIMENSION(ncell_fluid) :: source_liq,source_gas,source_drp
      !OPR1000 rod-scale (EVVD)
      REAL(8),DIMENSION(ncell_fluid) :: evvd,evvdg,evvdl
!.....Local vector arrays
      REAL(8),DIMENSION(nf_non) :: svoli_l_non,svoli_g_non,svoli_d_non
      REAL(8),DIMENSION(nf_fluxk2) :: poiss_diag_nf
      REAL(8),DIMENSION(nf_flux1) :: source_liq_nf,source_gas_nf,source_drp_nf
      REAL(8),DIMENSION(nf_mcc) :: source_liq_mcc,source_gas_mcc,source_drp_mcc
      REAL(8),DIMENSION(nf_flux,ndim) :: vl_nf,vg_nf,vd_nf
!.....Local vector for Mass Flux on non faces
      REAL(8),DIMENSION(nf_flux,ndim) :: ml_nf,mg_nf,md_nf
!
      flag_iso_thermal=.false.
      IF(smac.eq.3.and.iso_thermal.eq.1) flag_iso_thermal=.true.
      face_body_face=.false.
      IF(fric_face+nd_face.gt.0) face_body_face=.true.
!      
      CALL vectorize_scalar_upwind  
!
!.....Cells non
!
      nf_number=0
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      IF(vfporous.eq.0) THEN
         IF(ndim.eq.2)THEN
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
               a=f1(i)
               b=f0(i)
               vl_nf(i1,1)=a*vl_n(ii,1)+b*vl_n(kk,1)
               vg_nf(i1,1)=a*vg_n(ii,1)+b*vg_n(kk,1)
               vd_nf(i1,1)=a*vd_n(ii,1)+b*vd_n(kk,1)
               vl_nf(i1,2)=a*vl_n(ii,2)+b*vl_n(kk,2)
               vg_nf(i1,2)=a*vg_n(ii,2)+b*vg_n(kk,2)
               vd_nf(i1,2)=a*vd_n(ii,2)+b*vd_n(kk,2)
            ENDDO
         ELSE
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
               a=f1(i)
               b=f0(i)
               vl_nf(i1,1)=a*vl_n(ii,1)+b*vl_n(kk,1)
               vg_nf(i1,1)=a*vg_n(ii,1)+b*vg_n(kk,1)
               vd_nf(i1,1)=a*vd_n(ii,1)+b*vd_n(kk,1)
               vl_nf(i1,2)=a*vl_n(ii,2)+b*vl_n(kk,2)
               vg_nf(i1,2)=a*vg_n(ii,2)+b*vg_n(kk,2)
               vd_nf(i1,2)=a*vd_n(ii,2)+b*vd_n(kk,2)
               vl_nf(i1,3)=a*vl_n(ii,3)+b*vl_n(kk,3)
               vg_nf(i1,3)=a*vg_n(ii,3)+b*vg_n(kk,3)
               vd_nf(i1,3)=a*vd_n(ii,3)+b*vd_n(kk,3)
            ENDDO
         ENDIF
      ELSE
         IF(ndim.eq.2)THEN
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
               a=f1(i)*porosity(ii)
               b=f0(i)*porosity(kk)
               c=1.0d0/perm_non(i)
               vl_nf(i1,1)=(a*vl_n(ii,1)+b*vl_n(kk,1))*c
               vg_nf(i1,1)=(a*vg_n(ii,1)+b*vg_n(kk,1))*c
               vd_nf(i1,1)=(a*vd_n(ii,1)+b*vd_n(kk,1))*c
               vl_nf(i1,2)=(a*vl_n(ii,2)+b*vl_n(kk,2))*c
               vg_nf(i1,2)=(a*vg_n(ii,2)+b*vg_n(kk,2))*c
               vd_nf(i1,2)=(a*vd_n(ii,2)+b*vd_n(kk,2))*c
            ENDDO
         ELSE
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
               a=f1(i)*porosity(ii)
               b=f0(i)*porosity(kk)
               c=1.0d0/perm_non(i)
               vl_nf(i1,1)=(a*vl_n(ii,1)+b*vl_n(kk,1))*c
               vg_nf(i1,1)=(a*vg_n(ii,1)+b*vg_n(kk,1))*c
               vd_nf(i1,1)=(a*vd_n(ii,1)+b*vd_n(kk,1))*c
               vl_nf(i1,2)=(a*vl_n(ii,2)+b*vl_n(kk,2))*c
               vg_nf(i1,2)=(a*vg_n(ii,2)+b*vg_n(kk,2))*c
               vd_nf(i1,2)=(a*vd_n(ii,2)+b*vd_n(kk,2))*c
               vl_nf(i1,3)=(a*vl_n(ii,3)+b*vl_n(kk,3))*c
               vg_nf(i1,3)=(a*vg_n(ii,3)+b*vg_n(kk,3))*c
               vd_nf(i1,3)=(a*vd_n(ii,3)+b*vd_n(kk,3))*c
            ENDDO
         ENDIF
      ENDIF
!
      IF(rv_choke.eq.1) CALL choke_massflowrate1(ml_nf,mg_nf,md_nf)      
!      
      IF(face_body_face)THEN
         IF(ndim.eq.2)THEN
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
               a=fac1_non(i)*coefm_l(ii)+fac_non(i)*coefm_l(kk)
               b=fac1_non(i)*coefm_g(ii)+fac_non(i)*coefm_g(kk)
               vl_nf(i1,1)=vl_nf(i1,1)+vl_f_non(i,1)*a
               vg_nf(i1,1)=vg_nf(i1,1)+vg_f_non(i,1)*b
               vl_nf(i1,2)=vl_nf(i1,2)+vl_f_non(i,2)*a
               vg_nf(i1,2)=vg_nf(i1,2)+vg_f_non(i,2)*b
            ENDDO
         ELSE
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
               a=fac1_non(i)*coefm_l(ii)+fac_non(i)*coefm_l(kk)
               b=fac1_non(i)*coefm_g(ii)+fac_non(i)*coefm_g(kk)
               vl_nf(i1,1)=vl_nf(i1,1)+vl_f_non(i,1)*a
               vg_nf(i1,1)=vg_nf(i1,1)+vg_f_non(i,1)*b
               vl_nf(i1,2)=vl_nf(i1,2)+vl_f_non(i,2)*a
               vg_nf(i1,2)=vg_nf(i1,2)+vg_f_non(i,2)*b
               vl_nf(i1,3)=vl_nf(i1,3)+vl_f_non(i,3)*a
               vg_nf(i1,3)=vg_nf(i1,3)+vg_f_non(i,3)*b
            ENDDO
         ENDIF
      ENDIF
!
!.....Cells mcc
!
      nf_number=1
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      IF(ndim.eq.2)THEN
         DO i=1,len
            i1=istart+i
            idx=i3invtbl(i)
            vl_nf(i1,1)=c3vl(1,idx)*xn_nf(i1,1)
            vg_nf(i1,1)=c3vg(1,idx)*xn_nf(i1,1)
            vd_nf(i1,1)=c3vl(1,idx)*xn_nf(i1,1)
            vl_nf(i1,2)=c3vl(1,idx)*xn_nf(i1,2)
            vg_nf(i1,2)=c3vg(1,idx)*xn_nf(i1,2)
            vd_nf(i1,2)=c3vl(1,idx)*xn_nf(i1,2)
         ENDDO
      ELSE
         DO i=1,len
            i1=istart+i
            idx=i3invtbl(i)
            vl_nf(i1,1)=c3vl(1,idx)*xn_nf(i1,1)
            vg_nf(i1,1)=c3vg(1,idx)*xn_nf(i1,1)
            vd_nf(i1,1)=c3vl(1,idx)*xn_nf(i1,1)
            vl_nf(i1,2)=c3vl(1,idx)*xn_nf(i1,2)
            vg_nf(i1,2)=c3vg(1,idx)*xn_nf(i1,2)
            vd_nf(i1,2)=c3vl(1,idx)*xn_nf(i1,2)
            vl_nf(i1,3)=c3vl(1,idx)*xn_nf(i1,3)
            vg_nf(i1,3)=c3vg(1,idx)*xn_nf(i1,3)
            vd_nf(i1,3)=c3vl(1,idx)*xn_nf(i1,3)
         ENDDO
      ENDIF
!
!.....Cells inl
!
      nf_number=2
      istart=istart_nf(1,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      len   =istart_nf(2,nf_number)
      IF(ndim.eq.2)THEN
         DO i=1,len
            i1=istart+i
            i2=istart2+i
            k=nbcon_nf(i2)
            f_profile=vel_bc_profile_inl(i)
            IF(vin_norm(k).eq.0)THEN
               vl_nf(i1,1)=vb_liq(k,1)*f_profile
               vg_nf(i1,1)=vb_gas(k,1)*f_profile
               vd_nf(i1,1)=vb_drp(k,1)*f_profile
               vl_nf(i1,2)=vb_liq(k,2)*f_profile
               vg_nf(i1,2)=vb_gas(k,2)*f_profile
               vd_nf(i1,2)=vb_drp(k,2)*f_profile
            ELSE
               vl_nf(i1,1)=vin_liq(k)*xn_nf(i1,1)*f_profile
               vg_nf(i1,1)=vin_gas(k)*xn_nf(i1,1)*f_profile
               vd_nf(i1,1)=vin_drp(k)*xn_nf(i1,1)*f_profile
               vl_nf(i1,2)=vin_liq(k)*xn_nf(i1,2)*f_profile
               vg_nf(i1,2)=vin_gas(k)*xn_nf(i1,2)*f_profile
               vd_nf(i1,2)=vin_drp(k)*xn_nf(i1,2)*f_profile
            ENDIF
         ENDDO
      ELSE
         DO i=1,len
            i1=istart+i
            i2=istart2+i
            k=nbcon_nf(i2)
            f_profile=vel_bc_profile_inl(i)
            IF(vin_norm(k).eq.0)THEN
               vl_nf(i1,1)=vb_liq(k,1)*f_profile
               vg_nf(i1,1)=vb_gas(k,1)*f_profile
               vd_nf(i1,1)=vb_drp(k,1)*f_profile
               vl_nf(i1,2)=vb_liq(k,2)*f_profile
               vg_nf(i1,2)=vb_gas(k,2)*f_profile
               vd_nf(i1,2)=vb_drp(k,2)*f_profile
               vl_nf(i1,3)=vb_liq(k,3)*f_profile
               vg_nf(i1,3)=vb_gas(k,3)*f_profile
               vd_nf(i1,3)=vb_drp(k,3)*f_profile
            ELSE
               vl_nf(i1,1)=vin_liq(k)*xn_nf(i1,1)*f_profile
               vg_nf(i1,1)=vin_gas(k)*xn_nf(i1,1)*f_profile
               vd_nf(i1,1)=vin_drp(k)*xn_nf(i1,1)*f_profile
               vl_nf(i1,2)=vin_liq(k)*xn_nf(i1,2)*f_profile
               vg_nf(i1,2)=vin_gas(k)*xn_nf(i1,2)*f_profile
               vd_nf(i1,2)=vin_drp(k)*xn_nf(i1,2)*f_profile
               vl_nf(i1,3)=vin_liq(k)*xn_nf(i1,3)*f_profile
               vg_nf(i1,3)=vin_gas(k)*xn_nf(i1,3)*f_profile
               vd_nf(i1,3)=vin_drp(k)*xn_nf(i1,3)*f_profile
            ENDIF
         ENDDO
      ENDIF
!
!.....Cells out
!
      nf_number=3
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      IF(vfporous.eq.0) THEN
         IF(ndim.eq.2)THEN
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               vl_nf(i1,1)=vl_n(ii,1)
               vg_nf(i1,1)=vg_n(ii,1)
               vd_nf(i1,1)=vd_n(ii,1)
               vl_nf(i1,2)=vl_n(ii,2)
               vg_nf(i1,2)=vg_n(ii,2)
               vd_nf(i1,2)=vd_n(ii,2)
            ENDDO
         ELSE
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               vl_nf(i1,1)=vl_n(ii,1)
               vg_nf(i1,1)=vg_n(ii,1)
               vd_nf(i1,1)=vd_n(ii,1)
               vl_nf(i1,2)=vl_n(ii,2)
               vg_nf(i1,2)=vg_n(ii,2)
               vd_nf(i1,2)=vd_n(ii,2)
               vl_nf(i1,3)=vl_n(ii,3)
               vg_nf(i1,3)=vg_n(ii,3)
               vd_nf(i1,3)=vd_n(ii,3)
            ENDDO
         ENDIF
      ELSE
         IF(ndim.eq.2)THEN
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               a=porosity(ii)/perm_out(i)
               vl_nf(i1,1)=vl_n(ii,1)*a
               vg_nf(i1,1)=vg_n(ii,1)*a
               vd_nf(i1,1)=vd_n(ii,1)*a
               vl_nf(i1,2)=vl_n(ii,2)*a
               vg_nf(i1,2)=vg_n(ii,2)*a
               vd_nf(i1,2)=vd_n(ii,2)*a
            ENDDO
         ELSE
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               a=porosity(ii)/perm_out(i)
               vl_nf(i1,1)=vl_n(ii,1)*a
               vg_nf(i1,1)=vg_n(ii,1)*a
               vd_nf(i1,1)=vd_n(ii,1)*a
               vl_nf(i1,2)=vl_n(ii,2)*a
               vg_nf(i1,2)=vg_n(ii,2)*a
               vd_nf(i1,2)=vd_n(ii,2)*a
               vl_nf(i1,3)=vl_n(ii,3)*a
               vg_nf(i1,3)=vg_n(ii,3)*a
               vd_nf(i1,3)=vd_n(ii,3)*a
            ENDDO
         ENDIF
      ENDIF
!
      IF(ndim.eq.2)THEN
         DO nf_number=0,3
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               i1=istart+i
               vli1=vl_nf(i1,1)
               vli2=vl_nf(i1,2)
               vgi1=vg_nf(i1,1)
               vgi2=vg_nf(i1,2)
               vdi1=vd_nf(i1,1)
               vdi2=vd_nf(i1,2)
               flux_l_nf(i1)=vli1*svp_nf(i1,1)+vli2*svp_nf(i1,2)
               flux_g_nf(i1)=vgi1*svp_nf(i1,1)+vgi2*svp_nf(i1,2)
               flux_d_nf(i1)=vdi1*svp_nf(i1,1)+vdi2*svp_nf(i1,2)
            ENDDO
         ENDDO
      ELSE
         DO nf_number=0,3
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               i1=istart+i
               vli1=vl_nf(i1,1)
               vli2=vl_nf(i1,2)
               vli3=vl_nf(i1,3)
               vgi1=vg_nf(i1,1)
               vgi2=vg_nf(i1,2)
               vgi3=vg_nf(i1,3)
               vdi1=vd_nf(i1,1)
               vdi2=vd_nf(i1,2)
               vdi3=vd_nf(i1,3)
               flux_l_nf(i1)=vli1*svp_nf(i1,1)+vli2*svp_nf(i1,2)+vli3*svp_nf(i1,3)
               flux_g_nf(i1)=vgi1*svp_nf(i1,1)+vgi2*svp_nf(i1,2)+vgi3*svp_nf(i1,3)
               flux_d_nf(i1)=vdi1*svp_nf(i1,1)+vdi2*svp_nf(i1,2)+vdi3*svp_nf(i1,3)
            ENDDO
         ENDDO
      ENDIF
!      
      IF(rv_choke.eq.1) CALL choke_massflowrate2(ml_nf,mg_nf,md_nf,mflux_l_nf,mflux_g_nf,mflux_d_nf)
!
!.....fluxBC: choke model, mcp model, valve model
!                
      IF(rv_choke.eq.1.or.rv_mcp.eq.1.or.rv_valve.eq.1) CALL fluxBC_main  
!
!.....Modify the volume flux based on the Rhie-Chow scheme: smac series are limited to irc_damp=1 only.
!
      nf_number=0
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      IF(irc_damp.eq.1)THEN
          IF(ndim.eq.2)THEN
             DO i=1,len
                i1=istart+i
                ii=left_nf(i1)
                kk=right_non(i)              
                svoli_l=fac1_non(i)*coefp_l(ii)+fac_non(i)*coefp_l(kk)
                svoli_g=fac1_non(i)*coefp_g(ii)+fac_non(i)*coefp_g(kk)
                svoli_d=fac1_non(i)*coefp_d(ii)+fac_non(i)*coefp_d(kk)
                tmp=dt*sap_nf(i1)
                dp=0.0d0 
                dp=dp+(fac1_non(i)*dpdx(ii,1)+fac_non(i)*dpdx(kk,1))*dji_x_nf(i1,1)
                dp=dp+(fac1_non(i)*dpdx(ii,2)+fac_non(i)*dpdx(kk,2))*dji_x_nf(i1,2)
                dp=(p(kk)-p(ii))-dp
                dp=tmp*dp
                dp_fluxl=svoli_l*dp
                dp_fluxg=svoli_g*dp
                dp_fluxd=svoli_d*dp
                flux_l_nf(i1)=flux_l_nf(i1)-dp_fluxl
                flux_g_nf(i1)=flux_g_nf(i1)-dp_fluxg
                flux_d_nf(i1)=flux_d_nf(i1)-dp_fluxd
             ENDDO
          ELSEIF(ndim.eq.3)THEN
             DO i=1,len
                i1=istart+i
                ii=left_nf(i1)
                kk=right_non(i)              
                svoli_l=fac1_non(i)*coefp_l(ii)+fac_non(i)*coefp_l(kk)
                svoli_g=fac1_non(i)*coefp_g(ii)+fac_non(i)*coefp_g(kk)
                svoli_d=fac1_non(i)*coefp_d(ii)+fac_non(i)*coefp_d(kk)
                tmp=dt*sap_nf(i1)
                dp=0.0d0 
                dp=dp+(fac1_non(i)*dpdx(ii,1)+fac_non(i)*dpdx(kk,1))*dji_x_nf(i1,1)
                dp=dp+(fac1_non(i)*dpdx(ii,2)+fac_non(i)*dpdx(kk,2))*dji_x_nf(i1,2)
                dp=dp+(fac1_non(i)*dpdx(ii,3)+fac_non(i)*dpdx(kk,3))*dji_x_nf(i1,3)
                dp=(p(kk)-p(ii))-dp
                dp=tmp*dp
                dp_fluxl=svoli_l*dp
                dp_fluxg=svoli_g*dp
                dp_fluxd=svoli_d*dp
                flux_l_nf(i1)=flux_l_nf(i1)-dp_fluxl
                flux_g_nf(i1)=flux_g_nf(i1)-dp_fluxg
                flux_d_nf(i1)=flux_d_nf(i1)-dp_fluxd
             ENDDO
          ENDIF
      ELSEIF(irc_damp.eq.2)THEN
          IF(ndim.eq.2)THEN
             DO i=1,len
                i1=istart+i  
                ii=left_nf(i1)
                kk=right_non(i)              
                svoli_l=fac1_non(i)*coefp_l(ii)+fac_non(i)*coefp_l(kk)
                svoli_g=fac1_non(i)*coefp_g(ii)+fac_non(i)*coefp_g(kk)
                svoli_d=fac1_non(i)*coefp_d(ii)+fac_non(i)*coefp_d(kk)
                tmp=dt*sap_nf(i1)
                dp=0.0d0 
                dp=dp+(fac1_non(i)*dpdx(ii,1)+fac_non(i)*dpdx(kk,1))*sv_nf(i1,1)
                dp=dp+(fac1_non(i)*dpdx(ii,2)+fac_non(i)*dpdx(kk,2))*sv_nf(i1,2)
                dp=tmp*(p(kk)-p(ii))-dp*dt
                dp_fluxl=svoli_l*dp
                dp_fluxg=svoli_g*dp
                dp_fluxd=svoli_d*dp
                flux_l_nf(i1)=flux_l_nf(i1)-dp_fluxl
                flux_g_nf(i1)=flux_g_nf(i1)-dp_fluxg
                flux_d_nf(i1)=flux_d_nf(i1)-dp_fluxd
             ENDDO
          ELSEIF(ndim.eq.3)THEN
             DO i=1,len
                i1=istart+i  
                ii=left_nf(i1)
                kk=right_non(i)              
                svoli_l=fac1_non(i)*coefp_l(ii)+fac_non(i)*coefp_l(kk)
                svoli_g=fac1_non(i)*coefp_g(ii)+fac_non(i)*coefp_g(kk)
                svoli_d=fac1_non(i)*coefp_d(ii)+fac_non(i)*coefp_d(kk)
                tmp=dt*sap_nf(i1)
                dp=0.0d0 
                dp=dp+(fac1_non(i)*dpdx(ii,1)+fac_non(i)*dpdx(kk,1))*sv_nf(i1,1)
                dp=dp+(fac1_non(i)*dpdx(ii,2)+fac_non(i)*dpdx(kk,2))*sv_nf(i1,2)
                dp=dp+(fac1_non(i)*dpdx(ii,3)+fac_non(i)*dpdx(kk,3))*sv_nf(i1,3)
                dp=tmp*(p(kk)-p(ii))-dp*dt
                dp_fluxl=svoli_l*dp
                dp_fluxg=svoli_g*dp
                dp_fluxd=svoli_d*dp
                flux_l_nf(i1)=flux_l_nf(i1)-dp_fluxl
                flux_g_nf(i1)=flux_g_nf(i1)-dp_fluxg
                flux_d_nf(i1)=flux_d_nf(i1)-dp_fluxd
             ENDDO
          ENDIF
      ENDIF   
!
!.....fluxBC: choke model, mcp model, valve model
!           
      IF(rv_choke.eq.1.or.rv_mcp.eq.1.or.rv_valve.eq.1) CALL fluxBC_flux_update(flux_l_nf,flux_g_nf,flux_d_nf)
!
!.....Build summation info for non,mcc
!
      nf_number_nb=1
      nf_number_id(-1)=-1
      nf_number_id(0) =0
      nf_number_id(1) =1
      istart_nfs(-1)=0
      istart_nfs(0)=istart_nfs(-1)+nf_nonk
      istart_nfs(1)=istart_nfs(0) +nf_non
      lens         =istart_nfs(1) +nf_mcc
!
!.....Make pressure matrix coefficients
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
         svoli_l_non(i)=f1(i)*coefp_l(ii)+f0(i)*coefp_l(kk)
         svoli_g_non(i)=f1(i)*coefp_g(ii)+f0(i)*coefp_g(kk)
         svoli_d_non(i)=f1(i)*coefp_d(ii)+f0(i)*coefp_d(kk)
      ENDDO
!DIR$ SIMD
      DO i=1,len
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)
         arli=arli_nf(i1)
         argi=argi_nf(i1)
         ardi=ardi_nf(i1)
         tmp=dt*sap_nf(i1)
         arli=arli*svoli_l_non(i)
         argi=argi*svoli_g_non(i)
         ardi=ardi*svoli_d_non(i)
         cf=tmp*(arli/cell%rhol(ii)+argi/cell%rhog(ii)+ardi/cell%rhod(ii))
         poiss_diag_nf(i0)=-cf
         poiss_non_i(i)=cf
      ENDDO
!...........To be removed all runvv passed
      DO i=1,len
         kk=right_non(i)
         IF(npb(kk).gt.0) poiss_non_i(i)=0.d0
      ENDDO
!      
      nv=-1
      nf_number=nf_number_id(nv)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         k=right_nb_k(i)
         ii=right_non(k)
         arli=arli_nf(k)
         argi=argi_nf(k)
         ardi=ardi_nf(k)
         tmp=dt*sap_nf(k)
         arli=arli*svoli_l_non(k)
         argi=argi*svoli_g_non(k)
         ardi=ardi*svoli_d_non(k)
         cf=tmp*(arli/cell%rhol(ii)+argi/cell%rhog(ii)+ardi/cell%rhod(ii))
         poiss_diag_nf(i)=-cf
         poiss_non_k(i)=cf
      ENDDO
!...........To be removed all runvv passed
      DO i=1,len
         k=right_nb_k(i)
         kk=left_nf(k) 
         IF(npb(kk).gt.0) poiss_non_k(i)=0.d0
      ENDDO
!
!......valve model
!         
      IF(rv_valve.eq.1) CALL valve_model_pressure_matrix2(poiss_diag_nf,poiss_non_i,poiss_non_k) 
!         
!.....MARS interface
!
      nv=1
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
!DIR$ SIMD
      DO i=1,len
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)
         idx=i3invtbl(i)
         argi=c3dpv(idx,3)
         arli=c3dpv(idx,5)
         ardi=c3dpv(idx,3)
         svoli_l=c3betaf(1,idx)
         svoli_g=c3betag(1,idx)
         svoli_d=c3betaf(1,idx)
         tmp=dt*sap_nf(i1)
         cf=tmp*(arli*svoli_l/cell%rhol(ii)+argi*svoli_g/cell%rhog(ii))
         poiss_diag_mcc_s=0.0d0
         DO mm=1,n_marsbc
            poiss_diag_mcc_s=poiss_diag_mcc_s-cf*(1.0d0-c3yeta(1,idx,mm)) 
         ENDDO
         poiss_diag_nf(i0)=poiss_diag_mcc_s
      ENDDO
!
      CALL sum_nf(0,0,                      &
                  poiss_diag_nf,poiss_diag)
!
!.....Build summation info for non,mcc,inl
!
      nf_number_nb=2
      nf_number_id(0)=0
      nf_number_id(1)=1
      nf_number_id(2)=2
      istart_nfs(0)=0
      istart_nfs(1)=istart_nfs(0)+nf_non
      istart_nfs(2)=istart_nfs(1)+nf_mcc
      lens         =istart_nfs(2)+nf_inl
!
!.....Make pressure source vector non,mcc,inl
!
      DO nf_number=0,2
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len
            i1=istart+i
            ii=left_nf(i1)
            source_liq_nf(i1)=arli_nf(i1)*flux_l_nf(i1)
            source_gas_nf(i1)=argi_nf(i1)*flux_g_nf(i1)
            source_drp_nf(i1)=ardi_nf(i1)*flux_d_nf(i1)
         ENDDO
      ENDDO
!
      CALL sum_nf(0,-1,                     &
                  source_liq_nf,source_liq, &
                  source_gas_nf,source_gas, &
                  source_drp_nf,source_drp)
!
!.....Build summation info for mcc
!
      nf_number_nb=0
      nf_number_id(0)=1
      istart_nfs(0)=0
      lens         =istart_nfs(0)+nf_mcc
!
      nv=0
      nf_number=nf_number_id(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      IF(len.gt.0) THEN
!DIR$ SIMD
         DO i=1,len
            i1=istart+i
            idx=i3invtbl(i)
            argi=c3dpv(idx,3)
            arli=c3dpv(idx,5)
            ardi=c3dpv(idx,3)
            svoli_l=c3betaf(1,idx)
            svoli_g=c3betag(1,idx)
            svoli_d=c3betaf(1,idx)
            tmp=dt*sap_nf(i1)
            source_liq_mcc(i)=tmp*svoli_l*arli*c3xi(1,idx)
            source_gas_mcc(i)=tmp*svoli_g*argi*c3xi(1,idx)
            source_drp_mcc(i)=tmp*svoli_d*ardi*c3xi(1,idx)
         ENDDO
!
         CALL sum_nf(1,1,                       &
                     source_liq_mcc,source_liq, &
                     source_gas_mcc,source_gas, &
                     source_drp_mcc,source_drp)
      ENDIF
!
      DO i=1,ncell_fluid
         src(i)= source_liq(i)/cell%rhol(i) &
                +source_gas(i)/cell%rhog(i) &
                +source_drp(i)/cell%rhod(i)
      ENDDO
!
!--------------------------------------------------------
!........Complete the Poisson matrix and the source vector adding source terms from density and phase chages
!--------------------------------------------------------
!
      IF(flag_iso_thermal)THEN
         DO i=1,ncell_fluid 
            al=cell%alphal(i)+cell%alphad(i) 
            ag_rhog=cell%alphag(i)/cell%rhog(i) 
            al_rhol=al/cell%rhol(i) 
            poiss_diag(i)=poiss_diag(i)-volp(i) &
                          *((cell%drholdp(i)*al_rhol+cell%drhogdp(i)*ag_rhog)/dt)
         ENDDO 
      ELSE 
         DO i=1,ncell_fluid 
!
!...........EVVD
            IF(l_subchannel)then
               evvd(i)=0.0d0
               evvdg(i)=tm_mas_g(i)+vd_mas_g(i)
               evvdl(i)=tm_mas_l(i)+vd_mas_l(i) 
               IF(l_mixing_vane)then
                  evvdl(i)=evvdl(i)+mixing_vane_l(2,i)*volpr(i)
               ENDIF
            ENDIF      
         
            IF(gamma(i).ge.0.d0)THEN 
               ddhfg=cell%hgsat(i)-cell%hl(i)  
            ELSE 
               ddhfg=cell%hg(i)-cell%hlsat(i)  
            ENDIF 
! commented to match serial
            ddhfg=1.0d0/ddhfg
!           IF(ddhfg.ne.0)THEN
!              ddhfg=1.0d0/ddhfg
!           ELSE
!              ddhfg=1.0d8
!           ENDIF   
            PsP=cell%pps(i)/p(i) 
            c_ig=ddhfg*H_ig(i)*PsP 
            c_if=ddhfg*H_il(i)
            del_tsg=cell%ts(i)-cell%tg(i) 
            del_tsl=cell%ts(i)-cell%tl(i) 
            del_dtvdp=cell%dtsdp(i)-cell%dtgdp(i) 
            del_dtldp=cell%dtsdp(i)-cell%dtldp(i) 
            del_sv=1.d0/cell%rhol(i)-1.d0/cell%rhog(i) 
            al=cell%alphal(i)+cell%alphad(i) 
            ag_rhog=cell%alphag(i)/cell%rhog(i) 
            al_rhol=al/cell%rhol(i) 
!            
            poiss_diag(i)=poiss_diag(i)                                                                &
                                       -volp(i)*( (cell%drholdp(i)*al_rhol+cell%drhogdp(i)*ag_rhog)/dt &
                                                 -(c_ig*del_dtvdp+c_if*del_dtldp)*del_sv )
!
            IF(smac3_pres_eng.eq.0)THEN
               del_dtvde=(cell%dtsde(i)-cell%dtgde(i))*del_eg(i) 
               del_dtlde=cell%dtsde(i)*del_eg(i)-cell%dtlde(i)*del_el(i) 
               del_dtvdx=(cell%dtsdx(i)-cell%dtgdx(i))*del_x(i) 
               del_tsg=del_tsg+del_dtvde+del_dtvdx 
               del_tsl=del_tsl+del_dtlde+cell%dtsdx(i)*del_x(i)
               IF(repeat_smac)THEN
                  del_tsg=0.0d0
                  del_tsl=0.0d0
               ENDIF
               src(i)=src(i)-(                                                        &
                          volp(i)*((c_ig*del_tsg+c_if*del_tsl-gamma_wall(i))*del_sv         &
                         -(del_el(i)*cell%drholde(i)*al_rhol                                &
                         +(del_eg(i)*cell%drhogde(i)+del_x(i)*cell%drhogdx(i))*ag_rhog)/dt))
            ELSE
               IF(repeat_smac)THEN
                  del_tsg=0.0d0
                  del_tsl=0.0d0
               ENDIF
               src(i)=src(i)                                                          &
                         -volp(i)*((c_ig*del_tsg+c_if*del_tsl-gamma_wall(i))*del_sv)
            ENDIF
!            
            IF(l_subchannel)then
               src(i)=src(i)-volp(i)*(evvdg(i)/cell%rhog(i)+evvdl(i)/cell%rhol(i)) !EVVD          
            ENDIF   
         ENDDO
      ENDIF          
!
      END SUBROUTINE pressure_matrix_smac
