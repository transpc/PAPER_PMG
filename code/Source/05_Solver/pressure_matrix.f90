!
      SUBROUTINE pressure_matrix(poiss_diag,poiss_non_i,poiss_non_k,src)
!
!     This routine sets the matrix coefficients for pressure calculation.
!     flux_vol is the volume flux through the face.
!     poiss_diag - diagonal component
!     poiss - off diagonal component
!     src - source vector
!
!!!!!!!!!!!carefull not tested for all possnble irc_damp
!
      USE Zinterface
      USE VOL_DATA      , ONLY: cell
      USE Zmpi          , ONLY: jperm
      USE Zzone         , ONLY: ncell_fluid
      USE Zparam        , ONLY: ndim
      USE Zvec_param    , ONLY: nf_flux,nf_fluxk1,nf_fluxk2,nf_nonk,nf_non,nf_mcc,nf_inl
      USE Znum_cell     , ONLY: right_nb_k,istart_nf,istart_nbcon_nf,     &
                                nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zvec_index    , ONLY: left_nf,right_non,nbcon_nf
      USE Zbc_index     , ONLY: npb,vin_norm
      USE Zmars         , ONLY: n_marsbc
      USE Zare          , ONLY: ar_gas,ar_liq,ar_drp
      USE Zconst1       , ONLY: fric_face,nd_face
      USE Zconst2       , ONLY: dt,grav,gfactor
      USE Zcoord3       , ONLY: porosity
      USE Zb_condition  , ONLY: vb_liq,vb_gas,vb_drp,vin_liq,vin_gas,vin_drp
      USE Zgradoption   , ONLY: irc_damp
      USE Zporous       , ONLY: vfporous      
      USE Zpress        , ONLY: p,dpdx
      USE Zpress_coeff  , ONLY: coefp_g,coefp_l,coefp_d,coefm_g,coefm_l
      USE Zuserdefined  , ONLY: vel_bc_profile_inl
      USE Zvector       , ONLY: vl_n,vg_n,vd_n,vl_f_non,vg_f_non
      USE c3com_cupid   , ONLY: i3invtbl,i3cupid
      USE Zscalar_coeff , ONLY: sb
      USE Zscalar_coeff , ONLY: sfg6_nf,sfl6_nf,sfd6_nf,        &
                                sfg6_non_k,sfl6_non_k,sfd6_non_k
      USE Zvec_major    , ONLY: flux_l_nf,flux_g_nf,flux_d_nf,  &
                                mflux_l_nf,mflux_g_nf,mflux_d_nf
      USE Zvec_geo      , ONLY: xn_nf,sv_nf,svp_nf,dxfc_nf,     &
                                sap_nf,djia_nf,dji_x_nf,        & 
                                f0,f1,fac1_non,fac_non,         &
                                dxfc_non_k,                     &
                                perm_non,perm_out
      USE Zrv_model     , ONLY: rv_choke,rv_mcp,rv_valve
      USE Zmcp          
      USE Zvalve
!
      IMPLICIT NONE
!      
      INCLUDE '../10_LinkToMARS/c3com.h' 
!.....Output
      REAL(8),DIMENSION(nf_non) :: poiss_non_i
      REAL(8),DIMENSION(nf_nonk) :: poiss_non_k
      REAL(8),DIMENSION(ncell_fluid) :: poiss_diag,src
!.....Local variables
      INTEGER :: i,idx,mm
      INTEGER :: ii,kk,k
      INTEGER :: nv,nf_number,istart,len,istart0,istart2,i0,i1,i2
      LOGICAL :: face_body_face
      REAL(8) :: svoli_l,svoli_g,svoli_d
      REAL(8) :: tmp,cf,dp,dpi
      REAL(8) :: dp_fluxl,dp_fluxg,dp_fluxd
      REAL(8) :: vli1,vli2,vli3
      REAL(8) :: vgi1,vgi2,vgi3
      REAL(8) :: vdi1,vdi2,vdi3
      REAL(8) :: sdpe,dpe, alphag_f,alphal_f,alphad_f, rhog_f,rhol_f,rhod_f,g_gi,g_gk
      REAL(8) :: dg_g,dg_l,dg_d, gdrk,gdri,gdr, g_dk,g_di,vnormabs, g_li,g_lk
      REAL(8) :: f_profile      
      REAL(8) :: a,b,c
      REAL(8) :: a_li,a_lk,a_gi,a_gk,a_di,a_dk      
      REAL(8) :: poiss_diag_mcc_s
!.....Local vector arrays
      REAL(8),DIMENSION(nf_non) :: svoli_l_non,svoli_g_non,svoli_d_non
      REAL(8),DIMENSION(nf_flux) :: fluxt_l_nf,fluxt_g_nf,fluxt_d_nf
      REAL(8),DIMENSION(nf_fluxk2) :: poiss_diag_nf
      REAL(8),DIMENSION(nf_fluxk1) :: source_nf
      REAL(8),DIMENSION(nf_mcc) :: source_mcc
      REAL(8),DIMENSION(nf_flux,ndim) :: vl_nf,vg_nf,vd_nf
!.....Local vector for Mass Flux on non faces
      REAL(8),DIMENSION(nf_flux,ndim) :: ml_nf,mg_nf,md_nf
!
      face_body_face=.false.
      IF(fric_face+nd_face.gt.0) face_body_face=.true.
!
      ml_nf=0.d0
      mg_nf=0.d0
      md_nf=0.d0
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
!.....Discrete sources defined at cell face
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
!.....Calculate volume fluxes
!
      DO i=1,len  
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)         
         svoli_l_non(i)=f1(i)*coefp_l(ii)+f0(i)*coefp_l(kk)
         svoli_g_non(i)=f1(i)*coefp_g(ii)+f0(i)*coefp_g(kk)
         svoli_d_non(i)=f1(i)*coefp_d(ii)+f0(i)*coefp_d(kk)
      ENDDO
!
!.....Rhie-Chow's scheme
!
      IF(irc_damp.eq.1)THEN
         IF(ndim.eq.2)THEN
            DO i=1,len  
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)         
               dp =p(kk)-p(ii)
               dpi= (fac1_non(i)*dpdx(ii,1)+fac_non(i)*dpdx(kk,1))*dji_x_nf(i1,1) &
                   +(fac1_non(i)*dpdx(ii,2)+fac_non(i)*dpdx(kk,2))*dji_x_nf(i1,2)
               tmp=dt/djia_nf(i1)
!                
               vl_nf(i1,1)=vl_nf(i1,1)-tmp*svoli_l_non(i)*(dp-dpi)*xn_nf(i1,1)
               vg_nf(i1,1)=vg_nf(i1,1)-tmp*svoli_g_non(i)*(dp-dpi)*xn_nf(i1,1)
               vd_nf(i1,1)=vd_nf(i1,1)-tmp*svoli_d_non(i)*(dp-dpi)*xn_nf(i1,1)
               vl_nf(i1,2)=vl_nf(i1,2)-tmp*svoli_l_non(i)*(dp-dpi)*xn_nf(i1,2)
               vg_nf(i1,2)=vg_nf(i1,2)-tmp*svoli_g_non(i)*(dp-dpi)*xn_nf(i1,2)
               vd_nf(i1,2)=vd_nf(i1,2)-tmp*svoli_d_non(i)*(dp-dpi)*xn_nf(i1,2)
            ENDDO
         ELSE
            DO i=1,len  
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)         
               dp =p(kk)-p(ii)
               dpi= (fac1_non(i)*dpdx(ii,1)+fac_non(i)*dpdx(kk,1))*dji_x_nf(i1,1) &
                   +(fac1_non(i)*dpdx(ii,2)+fac_non(i)*dpdx(kk,2))*dji_x_nf(i1,2) &
                   +(fac1_non(i)*dpdx(ii,3)+fac_non(i)*dpdx(kk,3))*dji_x_nf(i1,3)
               tmp=dt/djia_nf(i1)
!                         
               vl_nf(i1,1)=vl_nf(i1,1)-tmp*svoli_l_non(i)*(dp-dpi)*xn_nf(i1,1)
               vg_nf(i1,1)=vg_nf(i1,1)-tmp*svoli_g_non(i)*(dp-dpi)*xn_nf(i1,1)
               vd_nf(i1,1)=vd_nf(i1,1)-tmp*svoli_d_non(i)*(dp-dpi)*xn_nf(i1,1)
               vl_nf(i1,2)=vl_nf(i1,2)-tmp*svoli_l_non(i)*(dp-dpi)*xn_nf(i1,2)
               vg_nf(i1,2)=vg_nf(i1,2)-tmp*svoli_g_non(i)*(dp-dpi)*xn_nf(i1,2)
               vd_nf(i1,2)=vd_nf(i1,2)-tmp*svoli_d_non(i)*(dp-dpi)*xn_nf(i1,2)
               vl_nf(i1,3)=vl_nf(i1,3)-tmp*svoli_l_non(i)*(dp-dpi)*xn_nf(i1,3)
               vg_nf(i1,3)=vg_nf(i1,3)-tmp*svoli_g_non(i)*(dp-dpi)*xn_nf(i1,3)
               vd_nf(i1,3)=vd_nf(i1,3)-tmp*svoli_d_non(i)*(dp-dpi)*xn_nf(i1,3)
            ENDDO
         ENDIF
      ENDIF
!      
       IF(rv_valve.eq.1) CALL valve_model_pressure_matrix1(svoli_l_non,svoli_g_non,svoli_d_non,vl_nf,vg_nf,vd_nf)
!
!.....MARS interface
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
               fluxt_l_nf(i1)=flux_l_nf(i1)
               fluxt_g_nf(i1)=flux_g_nf(i1)
               fluxt_d_nf(i1)=flux_d_nf(i1)
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
               fluxt_l_nf(i1)=flux_l_nf(i1)
               fluxt_g_nf(i1)=flux_g_nf(i1)
               fluxt_d_nf(i1)=flux_d_nf(i1)
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
!.....Modify the volume flux based on the Rhie-Chow scheme: irc_damp=1 & 2.
!
      nf_number=0
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      IF(irc_damp.eq.2)THEN
          IF(ndim.eq.2)THEN
             DO i=1,len  
                i1=istart+i
                ii=left_nf(i1)
                kk=right_non(i)              
                  a_li=fac1_non(i)*coefp_l(ii)
                  a_lk=fac_non(i) *coefp_l(kk)
                  a_gi=fac1_non(i)*coefp_g(ii)
                  a_gk=fac_non(i) *coefp_g(kk)
                  a_di=fac1_non(i)*coefp_d(ii)
                  a_dk=fac_non(i) *coefp_d(kk)
                  dp_fluxl= (a_li*dpdx(ii,1)+a_lk*dpdx(kk,1))*sv_nf(i1,1) &
                           +(a_li*dpdx(ii,2)+a_lk*dpdx(kk,2))*sv_nf(i1,2)
                  dp_fluxg= (a_gi*dpdx(ii,1)+a_gk*dpdx(kk,1))*sv_nf(i1,1) &
                           +(a_gi*dpdx(ii,2)+a_gk*dpdx(kk,2))*sv_nf(i1,2)
                  dp_fluxd= (a_di*dpdx(ii,1)+a_dk*dpdx(kk,1))*sv_nf(i1,1) &
                           +(a_di*dpdx(ii,2)+a_dk*dpdx(kk,2))*sv_nf(i1,2)
                  tmp=dt*sap_nf(i1)
                  dp=p(kk)-p(ii)
                  dp_fluxl=tmp*svoli_l_non(i)*dp-dt*dp_fluxl
                  dp_fluxg=tmp*svoli_g_non(i)*dp-dt*dp_fluxg
                  dp_fluxd=tmp*svoli_d_non(i)*dp-dt*dp_fluxd
                fluxt_l_nf(i1)=fluxt_l_nf(i1)-dp_fluxl
                fluxt_g_nf(i1)=fluxt_g_nf(i1)-dp_fluxg
                fluxt_d_nf(i1)=fluxt_d_nf(i1)-dp_fluxd
             ENDDO
          ELSEIF(ndim.eq.3)THEN
             DO i=1,len  
                i1=istart+i
                ii=left_nf(i1)
                kk=right_non(i)              
                  a_li=fac1_non(i)*coefp_l(ii)
                  a_lk=fac_non(i) *coefp_l(kk)
                  a_gi=fac1_non(i)*coefp_g(ii)
                  a_gk=fac_non(i) *coefp_g(kk)
                  a_di=fac1_non(i)*coefp_d(ii)
                  a_dk=fac_non(i) *coefp_d(kk)
                  dp_fluxl= (a_li*dpdx(ii,1)+a_lk*dpdx(kk,1))*sv_nf(i1,1) &
                           +(a_li*dpdx(ii,2)+a_lk*dpdx(kk,2))*sv_nf(i1,2) &
                           +(a_li*dpdx(ii,3)+a_lk*dpdx(kk,3))*sv_nf(i1,3)
                  dp_fluxg= (a_gi*dpdx(ii,1)+a_gk*dpdx(kk,1))*sv_nf(i1,1) &
                           +(a_gi*dpdx(ii,2)+a_gk*dpdx(kk,2))*sv_nf(i1,2) &
                           +(a_gi*dpdx(ii,3)+a_gk*dpdx(kk,3))*sv_nf(i1,3)
                  dp_fluxd= (a_di*dpdx(ii,1)+a_dk*dpdx(kk,1))*sv_nf(i1,1) &
                           +(a_di*dpdx(ii,2)+a_dk*dpdx(kk,2))*sv_nf(i1,2) &
                           +(a_di*dpdx(ii,3)+a_dk*dpdx(kk,3))*sv_nf(i1,3)
                  tmp=dt*sap_nf(i1)
                  dp=p(kk)-p(ii)
                  dp_fluxl=tmp*svoli_l_non(i)*dp-dt*dp_fluxl
                  dp_fluxg=tmp*svoli_g_non(i)*dp-dt*dp_fluxg
                  dp_fluxd=tmp*svoli_d_non(i)*dp-dt*dp_fluxd
                fluxt_l_nf(i1)=fluxt_l_nf(i1)-dp_fluxl
                fluxt_g_nf(i1)=fluxt_g_nf(i1)-dp_fluxg
                fluxt_d_nf(i1)=fluxt_d_nf(i1)-dp_fluxd
             ENDDO
          ENDIF
      ELSEIF(irc_damp.eq.3)THEN
          IF(ndim.eq.2)THEN
             DO i=1,len  
                i1=istart+i
                ii=left_nf(i1)
                kk=right_non(i)              
                  dp = p(kk) - p(ii)
                    dpi = (dpdx(ii,1)+dpdx(kk,1))*dji_x_nf(i1,1) &
                         +(dpdx(ii,2)+dpdx(kk,2))*dji_x_nf(i1,2)
                  tmp=dt*sap_nf(i1)

                  dp = dp - 0.5d0 *dpi

                  dp_fluxl = (tmp*dp)*svoli_l_non(i)
                  dp_fluxg = (tmp*dp)*svoli_g_non(i)
                  dp_fluxd = (tmp*dp)*svoli_d_non(i)
                flux_l_nf(i1)=flux_l_nf(i1)-dp_fluxl
                flux_g_nf(i1)=flux_g_nf(i1)-dp_fluxg
                flux_d_nf(i1)=flux_d_nf(i1)-dp_fluxd
             ENDDO
          ELSEIF(ndim.eq.3)THEN
             DO i=1,len  
                i1=istart+i
                ii=left_nf(i1)
                kk=right_non(i)              
                  dp = p(kk) - p(ii)
                    dpi = (dpdx(ii,1)+dpdx(kk,1))*dji_x_nf(i1,1) &
                         +(dpdx(ii,2)+dpdx(kk,2))*dji_x_nf(i1,2) &
                         +(dpdx(ii,3)+dpdx(kk,3))*dji_x_nf(i1,3)
                  tmp=dt*sap_nf(i1)

                  dp = dp - 0.5d0 *dpi

                  dp_fluxl = (tmp*dp)*svoli_l_non(i)
                  dp_fluxg = (tmp*dp)*svoli_g_non(i)
                  dp_fluxd = (tmp*dp)*svoli_d_non(i)
                flux_l_nf(i1)=flux_l_nf(i1)-dp_fluxl
                flux_g_nf(i1)=flux_g_nf(i1)-dp_fluxg
                flux_d_nf(i1)=flux_d_nf(i1)-dp_fluxd
             ENDDO
          ENDIF
!
!.....Limited Rhie-Chow's scheme
!
      ELSEIF(irc_damp.eq.4)THEN
          IF(ndim.eq.2)THEN
             DO i=1,len  
                i1=istart+i
                ii=left_nf(i1)
                kk=right_non(i)              
                  dp=p(kk)-p(ii)

                  dpi= (dpdx(ii,1)+dpdx(kk,1))*dji_x_nf(i1,1) &
                      +(dpdx(ii,2)+dpdx(kk,2))*dji_x_nf(i1,2)

                  dp = dp - 0.5d0*dpi

                  tmp=dt*sap_nf(i1)
                  dp_fluxl = tmp*svoli_l_non(i)*dp
                  dp_fluxg = tmp*svoli_g_non(i)*dp
                  dp_fluxd = tmp*svoli_d_non(i)*dp

                  ! local limiting
                  sdpe = sign(1.0d0, dp_fluxl)
                  vnormabs = abs(flux_l_nf(i1))
                  dpe = abs(dp_fluxl)/(vnormabs + 1.0d-15)
                  dp_fluxl = min( 0.03d0, dpe )*vnormabs
                  flux_l_nf(i1) = flux_l_nf(i1) - sdpe*dp_fluxl

                  sdpe = sign(1.0d0, dp_fluxg)
                  vnormabs = abs(flux_g_nf(i1))
                  dpe = abs(dp_fluxg)/(vnormabs + 1.0d-15)
                  dp_fluxg = min( 0.03d0, dpe )*vnormabs
                  flux_g_nf(i1) = flux_g_nf(i1) - sdpe*dp_fluxg

                  sdpe = sign(1.0d0, dp_fluxd)
                  vnormabs = abs(flux_d_nf(i1))
                  dpe = abs(dp_fluxd)/(vnormabs + 1.0d-15)
                  dp_fluxd = min( 0.03d0, dpe )*vnormabs
                  flux_d_nf(i1) = flux_d_nf(i1) - sdpe*dp_fluxd
             ENDDO
          ELSEIF(ndim.eq.3)THEN
             DO i=1,len  
                i1=istart+i
                ii=left_nf(i1)
                kk=right_non(i)              
                  dp=p(kk)-p(ii)

                  dpi= (dpdx(ii,1)+dpdx(kk,1))*dji_x_nf(i1,1) &
                      +(dpdx(ii,2)+dpdx(kk,2))*dji_x_nf(i1,2) &
                      +(dpdx(ii,3)+dpdx(kk,3))*dji_x_nf(i1,3)

                  dp = dp - 0.5d0*dpi

                  tmp=dt*sap_nf(i1)
                  dp_fluxl = tmp*svoli_l_non(i)*dp
                  dp_fluxg = tmp*svoli_g_non(i)*dp
                  dp_fluxd = tmp*svoli_d_non(i)*dp

                  ! local limiting
                  sdpe = sign(1.0d0, dp_fluxl)
                  vnormabs = abs(flux_l_nf(i1))
                  dpe = abs(dp_fluxl)/(vnormabs + 1.0d-15)
                  dp_fluxl = min( 0.03d0, dpe )*vnormabs
                  flux_l_nf(i1) = flux_l_nf(i1) - sdpe*dp_fluxl

                  sdpe = sign(1.0d0, dp_fluxg)
                  vnormabs = abs(flux_g_nf(i1))
                  dpe = abs(dp_fluxg)/(vnormabs + 1.0d-15)
                  dp_fluxg = min( 0.03d0, dpe )*vnormabs
                  flux_g_nf(i1) = flux_g_nf(i1) - sdpe*dp_fluxg

                  sdpe = sign(1.0d0, dp_fluxd)
                  vnormabs = abs(flux_d_nf(i1))
                  dpe = abs(dp_fluxd)/(vnormabs + 1.0d-15)
                  dp_fluxd = min( 0.03d0, dpe )*vnormabs
                  flux_d_nf(i1) = flux_d_nf(i1) - sdpe*dp_fluxd
             ENDDO
          ENDIF
!
!.....Rhie-Chow's scheme with body force
!
      ELSEIF(irc_damp.eq.5)THEN
          IF(ndim.eq.2)THEN
             DO i=1,len  
                i1=istart+i
                ii=left_nf(i1)
                kk=right_non(i)              
                  ! pressure damping term
                  dp=p(kk)-p(ii)
                  dpi  = (dpdx(ii,1)+dpdx(kk,1))*dji_x_nf(i1,1) &
                        +(dpdx(ii,2)+dpdx(kk,2))*dji_x_nf(i1,2)
                  gdri = grav(1)*dxfc_nf(i1,1)  +grav(2)*dxfc_nf(i1,2)
                  gdrk = grav(1)*dxfc_non_k(i,1)+grav(2)*dxfc_non_k(i,2)
                  gdr  = grav(1)*dji_x_nf(i1,1)   +grav(2)*dji_x_nf(i1,2)
                  dp = dp - 0.5d0*dpi

                  ! body force damping term
                  alphag_f = 1.d0/max(1.0d-8, 0.5d0*(cell%alphag(ii) + cell%alphag(kk)))
                  alphal_f = 1.d0/max(1.0d-8, 0.5d0*(cell%alphal(ii) + cell%alphal(kk)))
                  alphad_f = 1.d0/max(1.0d-8, 0.5d0*(cell%alphad(ii) + cell%alphad(kk)))

                  rhog_f = max(1.0d-8, 0.5d0*(cell%rhog(ii) + cell%rhog(kk)))
                  rhol_f = max(1.0d-8, 0.5d0*(cell%rhol(ii) + cell%rhol(kk)))
                  rhod_f = max(1.0d-8, 0.5d0*(cell%rhod(ii) + cell%rhod(kk)))

                  gdri = gdri*(gfactor(ii)*gfactor(kk))
                  g_gi =  ar_gas(ii)*alphag_f*gdri
                  g_li =  ar_liq(ii)*alphal_f*gdri
                  g_di =  ar_drp(ii)*alphad_f*gdri

                  gdrk = gdrk*(gfactor(kk)*gfactor(ii))
                  g_gk =  ar_gas(kk)*alphag_f*gdrk
                  g_lk =  ar_liq(kk)*alphal_f*gdrk
                  g_dk =  ar_drp(kk)*alphad_f*gdrk

                  dg_g = g_gi - g_gk
                  dg_l = g_li - g_lk
                  dg_d = g_di - g_dk

                  gdr = gdr*(gfactor(ii)*gfactor(kk))

                  dg_g = dg_g - rhog_f*gdr
                  dg_l = dg_l - rhol_f*gdr
                  dg_d = dg_d - rhod_f*gdr

                  tmp=dt*sap_nf(i1)
                  dp_fluxl = tmp*svoli_l_non(i)*(dp - dg_l)
                  dp_fluxg = tmp*svoli_g_non(i)*(dp - dg_g)
                  dp_fluxd = tmp*svoli_d_non(i)*(dp - dg_d)

                flux_l_nf(i1)=flux_l_nf(i1)-dp_fluxl
                flux_g_nf(i1)=flux_g_nf(i1)-dp_fluxg
                flux_d_nf(i1)=flux_d_nf(i1)-dp_fluxd
             ENDDO
          ELSEIF(ndim.eq.3)THEN
             DO i=1,len  
                i1=istart+i
                ii=left_nf(i1)
                kk=right_non(i)              
                  ! pressure damping term
                  dp=p(kk)-p(ii)
                  dpi= (dpdx(ii,1)+dpdx(kk,1))*dji_x_nf(i1,1) &
                      +(dpdx(ii,2)+dpdx(kk,2))*dji_x_nf(i1,2) & 
                      +(dpdx(ii,3)+dpdx(kk,3))*dji_x_nf(i1,3)
                  gdri=grav(1)*dxfc_nf(i1,1)  +grav(2)*dxfc_nf(i1,2)  +grav(3)*dxfc_nf(i1,3)
                  gdrk=grav(1)*dxfc_non_k(i,1)+grav(2)*dxfc_non_k(i,2)+grav(3)*dxfc_non_k(i,3)
                  gdr =grav(1)*dji_x_nf(i1,1)   +grav(2)*dji_x_nf(i1,2)   +grav(3)*dji_x_nf(i1,3)
                  dp =dp - 0.5d0*dpi

                  ! body force damping term
                  alphag_f = 1.d0/max(1.0d-8, 0.5d0*(cell%alphag(ii) + cell%alphag(kk)))
                  alphal_f = 1.d0/max(1.0d-8, 0.5d0*(cell%alphal(ii) + cell%alphal(kk)))
                  alphad_f = 1.d0/max(1.0d-8, 0.5d0*(cell%alphad(ii) + cell%alphad(kk)))

                  rhog_f = max(1.0d-8, 0.5d0*(cell%rhog(ii) + cell%rhog(kk)))
                  rhol_f = max(1.0d-8, 0.5d0*(cell%rhol(ii) + cell%rhol(kk)))
                  rhod_f = max(1.0d-8, 0.5d0*(cell%rhod(ii) + cell%rhod(kk)))

                  gdri = gdri*(gfactor(ii)*gfactor(kk))
                  g_gi =  ar_gas(ii)*alphag_f*gdri
                  g_li =  ar_liq(ii)*alphal_f*gdri
                  g_di =  ar_drp(ii)*alphad_f*gdri

                  gdrk = gdrk*(gfactor(kk)*gfactor(ii))
                  g_gk =  ar_gas(kk)*alphag_f*gdrk
                  g_lk =  ar_liq(kk)*alphal_f*gdrk
                  g_dk =  ar_drp(kk)*alphad_f*gdrk

                  dg_g = g_gi - g_gk
                  dg_l = g_li - g_lk
                  dg_d = g_di - g_dk

                  gdr = gdr*(gfactor(ii)*gfactor(kk))

                  dg_g = dg_g - rhog_f*gdr
                  dg_l = dg_l - rhol_f*gdr
                  dg_d = dg_d - rhod_f*gdr

                  tmp=dt*sap_nf(i1)
                  dp_fluxl = tmp*svoli_l_non(i)*(dp - dg_l)
                  dp_fluxg = tmp*svoli_g_non(i)*(dp - dg_g)
                  dp_fluxd = tmp*svoli_d_non(i)*(dp - dg_d)

                flux_l_nf(i1)=flux_l_nf(i1)-dp_fluxl
                flux_g_nf(i1)=flux_g_nf(i1)-dp_fluxg
                flux_d_nf(i1)=flux_d_nf(i1)-dp_fluxd
             ENDDO

                  ! local limiting
                  !sdpe = dsign(1.0d0, dp_fluxl)
                  !vnormabs = dabs(fluxvol_l(j,i))
                  !dpe = dabs(dp_fluxl)/(vnormabs + 1.0d-15)
                  !dp_fluxl = min( 0.03d0, dpe )*vnormabs
                  !fluxvol_l(j,i) = fluxvol_l(j,i) - sdpe*dp_fluxl

                  !sdpe = dsign(1.0d0, dp_fluxg)
                  !vnormabs = dabs(fluxvol_g(j,i))
                  !dpe = dabs(dp_fluxg)/(vnormabs + 1.0d-15)
                  !dp_fluxg = min( 0.03d0, dpe )*vnormabs
                  !fluxvol_g(j,i) = fluxvol_g(j,i) - sdpe*dp_fluxg

                  !sdpe = dsign(1.0d0, dp_fluxd)
                  !vnormabs = dabs(fluxvol_d(j,i))
                  !dpe = dabs(dp_fluxd)/(vnormabs + 1.0d-15)
                  !dp_fluxd = min( 0.03d0, dpe )*vnormabs
                  !fluxvol_d(j,i) = fluxvol_d(j,i) - sdpe*dp_fluxd
          ENDIF
      ENDIF  
!
!...fluxBC: choke model, mcp model, valve model
!           
      IF(rv_valve.eq.1.or.rv_choke.eq.1.or.rv_mcp.eq.1) CALL fluxBC_flux_update_ice(flux_l_nf,flux_g_nf,flux_d_nf,fluxt_l_nf,fluxt_g_nf,fluxt_d_nf)
!
!.....Build summation info for non,mcc
!
      nf_number_nb=1
      nf_number_id(-1)=-1
      nf_number_id(0)=0
      nf_number_id(1)=1
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
!DIR$ SIMD
      DO i=1,len
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)
         tmp=dt*sap_nf(i1)
         IF(npb(ii).eq.0) THEN
            cf=tmp*(svoli_l_non(i)*sfl6_nf(i1)+svoli_g_non(i)*sfg6_nf(i1)+svoli_d_non(i)*sfd6_nf(i1)) !A-1G*betag = svoli_l_non(ir)*sfl6_nf(i1)
            poiss_diag_nf(i0)=-cf
            poiss_non_i(i)=cf
         ELSE
            poiss_diag_nf(i0)=0.d0
         ENDIF
      ENDDO
!
      nv=-1
      nf_number=nf_number_id(nv)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         k=right_nb_k(i)
         ii=right_non(k)
         tmp=dt*sap_nf(k)
         IF(npb(ii).eq.0) THEN
            cf=tmp*(svoli_l_non(k)*sfl6_non_k(i)+svoli_g_non(k)*sfg6_non_k(i)+svoli_d_non(k)*sfd6_non_k(i))
            poiss_diag_nf(i)=-cf
            poiss_non_k(i)=cf
         ELSE
            poiss_diag_nf(i)=0.d0
         ENDIF
      ENDDO
!
!.....valve model
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
         IF(npb(ii).eq.0) THEN
            idx=i3invtbl(i) 
            svoli_l=-c3betaf(1,idx)
            svoli_g=-c3betag(1,idx)
            svoli_d=-c3betaf(1,idx)
            tmp=dt*sap_nf(i1)
            cf=tmp*(svoli_l*sfl6_nf(i1)+svoli_g*sfg6_nf(i1)+svoli_d*sfd6_nf(i1))
            poiss_diag_mcc_s=0.0d0
            DO mm=1,n_marsbc
               kk=i3cupid(mm)
               IF(kk.eq.jperm(ii))THEN
                  poiss_diag_mcc_s=poiss_diag_mcc_s-cf*(1.0d0-c3yeta(1,idx,mm)) 
               ENDIF
            ENDDO         
            poiss_diag_nf(i0)=poiss_diag_mcc_s
         ELSE
            poiss_diag_nf(i0)=0.d0
         ENDIF
      ENDDO
!
      CALL sum_nf(0,0,                      &
                  poiss_diag_nf,poiss_diag)
!
!.....Build summation info for non,mcc,inl
!
      nf_number_nb=2
      nf_number_id(-1)=-1
      nf_number_id(0)=0
      nf_number_id(1)=1
      nf_number_id(2)=2
      istart_nfs(-1)=0
      istart_nfs(0)=istart_nfs(-1)+nf_nonk
      istart_nfs(1)=istart_nfs(0) +nf_non
      istart_nfs(2)=istart_nfs(1) +nf_mcc
      lens         =istart_nfs(2) +nf_inl
!
      nv=-1
      nf_number=nf_number_id(nv)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         k=right_nb_k(i)
         ii=right_non(k)
         IF(npb(ii).eq.0) THEN
            source_nf(i)=-sfl6_non_k(i)*fluxt_l_nf(k) &
                         -sfg6_non_k(i)*fluxt_g_nf(k) &
                         -sfd6_non_k(i)*fluxt_d_nf(k)
         ELSE
            source_nf(i)=0.d0
         ENDIF
      ENDDO
!
      DO nv=0,2
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(nv)
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
!DIR$ SIMD
         DO i=1,len
            i0=istart0+i
            i1=istart+i
            ii=left_nf(i1)
            ii=left_nf(i1)
            IF(npb(ii).eq.0) THEN
               source_nf(i0)= sfl6_nf(i1)*fluxt_l_nf(i1) &
                             +sfg6_nf(i1)*fluxt_g_nf(i1) &
                             +sfd6_nf(i1)*fluxt_d_nf(i1)
            ELSE
               source_nf(i0)=0.d0
            ENDIF
         ENDDO
      ENDDO
!
      CALL sum_nf(0,0,          &
                  source_nf,src)
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
            ii=left_nf(i1)
            IF(npb(ii).eq.0) THEN
               idx=i3invtbl(i)
               svoli_l=-c3betaf(1,idx)
               svoli_g=-c3betag(1,idx)
               svoli_d=-c3betaf(1,idx)
               tmp=dt*sap_nf(i1)
               source_mcc(i)= tmp*svoli_l*sfl6_nf(i1)*c3xi(1,idx) &
                             +tmp*svoli_g*sfg6_nf(i1)*c3xi(1,idx) &
                             +tmp*svoli_d*sfd6_nf(i1)*c3xi(1,idx)
            ELSE
               source_mcc(i)=0.d0
            ENDIF
         ENDDO
!
         CALL sum_nf(1,1,            &
                     source_mcc,src)
!
      ENDIF
!
      DO i=1,ncell_fluid
         IF(npb(i).gt.0) CYCLE
         poiss_diag(i)=poiss_diag(i)+1.d0
         src(i)=src(i)+sb(i,6) ![A-1G+ ... ]+A-1B  for p-matrix
      ENDDO
!
      END SUBROUTINE pressure_matrix
