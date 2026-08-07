!
      SUBROUTINE momentum_diffusion_user(diff_liq,diff_gas,diff_drp)
!
!     This routine calculates diffusive fluxes through the cell face.
!     Diffusive fluxes are discretized using central differences.
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell            
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: ndim
      USE Zvec_param   , ONLY: nf_nonk,nf_non,nf_inl,nf_adw,nf_fsw,nf_ctw,nf_chw
      USE Znum_cell    , ONLY: istart_nf,istart_nbcon_nf,                 &
                               nf_number_nb,lens,nf_number_id,istart_nfs, &
                               right_nb_k
      USE Zvec_index   , ONLY: left_nf,right_non,nbcon_nf
      USE Zb_condition , ONLY: alphab_gas,alphab_liq,alphab_drp,vb_gas,vb_liq,vb_drp,vin_gas,vin_liq,vin_drp,v_wall
      USE Zconst1      , ONLY: mdiffscheme,lowreynolds
      USE Zface        , ONLY: twall_model,laminar
      USE Zturb        , ONLY: wvis_liq,wvis_gas,wallnr
      USE Zuserdefined , ONLY: vel_bc_profile_inl
      USE Zmodel       , ONLY: i_droplet
      USE Zvector      , ONLY: vg_o,vl_o,vd_o
      USE Zbc_index    , ONLY: vin_norm
      USE Zvec_geo     , ONLY: f1,f0,fac_fsw,xn_nf,sap_nf
!
      IMPLICIT NONE
!
!.....Output
      REAL(8),DIMENSION(ncell_fluid,ndim) :: diff_liq,diff_gas,diff_drp
!.....Local variables
      INTEGER :: i,k
      INTEGER :: ii,kk,ix
      INTEGER :: nv,nf_number,len,istart0,istart,istart2,i0,i1,i2
      REAL(8) :: vl11,vg11,vd11,vl21,vg21,vd21
      REAL(8) :: sv1
      REAL(8) :: avisli,avisgi,avisdi
      REAL(8) :: avisgil,avisgit,avislil,avislit
      REAL(8) :: avisgil_i,avisgit_i,avislil_i,avislit_i
      REAL(8) :: avisgil_k,avisgit_k,avislil_k,avislit_k
      REAL(8) :: vl1,vl2,vl3,vg1,vg2,vg3,vd1,vd2,vd3
      REAL(8) :: f_profile
      REAL(8) :: viscl,viscg
      REAL(8) :: a,b
!.....Local vector arrays
      REAL(8),DIMENSION(nf_nonk+nf_non+nf_inl+nf_adw+nf_fsw+nf_ctw+nf_chw,ndim) :: fluxl_diff_nf,fluxg_diff_nf
      REAL(8),DIMENSION(nf_non+nf_inl+nf_adw+nf_fsw+nf_ctw+nf_chw,ndim) :: fluxd_diff_nf
!
!.....Build summation info for all nf
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
      istart_nfs(1)=istart_nfs(0) +nf_non
      istart_nfs(2)=istart_nfs(1) +nf_inl
      istart_nfs(3)=istart_nfs(2) +nf_adw
      istart_nfs(4)=istart_nfs(3) +nf_fsw
      istart_nfs(5)=istart_nfs(4) +nf_ctw
      lens         =istart_nfs(5) +nf_chw
!
!.....Cells non
!
      IF(mdiffscheme.eq.1) THEN
         IF(ndim.eq.2)THEN
            nv=0
            nf_number=nf_number_id(nv)
            istart0=istart_nfs(nv)
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
!
               vl1=vl_o(kk,1)-vl_o(ii,1)
               vg1=vg_o(kk,1)-vg_o(ii,1)
               vd1=vd_o(kk,1)-vd_o(ii,1)
               vl2=vl_o(kk,2)-vl_o(ii,2)
               vg2=vg_o(kk,2)-vg_o(ii,2)
               vd2=vd_o(kk,2)-vd_o(ii,2)
               sv1=sap_nf(i1)
!
               a=f1(i)*cell%lviscosl(ii)+f0(i)*cell%lviscosl(kk)
               avislil_i=cell%alphal(ii)*a
               avislit_i=f1(i)*cell%tviscosl(ii)*cell%alphal(ii)+f0(i)*cell%tviscosl(kk)*cell%alphal(kk)
               fluxl_diff_nf(i0,1)=(avislil_i*vl1+avislit_i*vl1)*sv1
               fluxl_diff_nf(i0,2)=(avislil_i*vl2+avislit_i*vl2)*sv1
!
               a=f1(i)*cell%lviscosl(ii)+f0(i)*cell%lviscosl(kk) 
               avisgil_i=cell%alphag(ii)*a
               avisgit_i=f1(i)*cell%tviscosg(ii)*cell%alphag(ii)+f0(i)*cell%tviscosg(kk)*cell%alphag(kk)
               fluxg_diff_nf(i0,1)=(avisgil_i*vl1+avisgit_i*vg1)*sv1
               fluxg_diff_nf(i0,2)=(avisgil_i*vl2+avisgit_i*vg2)*sv1
!
               avisdi=f1(i)*cell%eviscosd(ii)*cell%alphad(ii)+f0(i)*cell%eviscosd(kk)*cell%alphad(kk)
               fluxd_diff_nf(i1,1)=(avisdi*vd1)*sv1
               fluxd_diff_nf(i1,2)=(avisdi*vd2)*sv1
            ENDDO
!
            nv=-1
            nf_number=nf_number_id(nv)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               k=right_nb_k(i)
               ii=right_non(k)
               kk=left_nf(k)
!
               vl1=vl_o(kk,1)-vl_o(ii,1)
               vg1=vg_o(kk,1)-vg_o(ii,1)
               vd1=vd_o(kk,1)-vd_o(ii,1)
               vl2=vl_o(kk,2)-vl_o(ii,2)
               vg2=vg_o(kk,2)-vg_o(ii,2)
               vd2=vd_o(kk,2)-vd_o(ii,2)
               sv1=sap_nf(k)
!
               a=f1(k)*cell%lviscosl(kk)+f0(k)*cell%lviscosl(ii)
               avislil_k=cell%alphal(ii)*a
               avislit_k=f1(k)*cell%tviscosl(kk)*cell%alphal(kk)+f0(k)*cell%tviscosl(ii)*cell%alphal(ii)
               fluxl_diff_nf(i,1)=(avislil_k*vl1+avislit_k*vl1)*sv1
               fluxl_diff_nf(i,2)=(avislil_k*vl2+avislit_k*vl2)*sv1
!
               a=f1(k)*cell%lviscosl(kk)+f0(k)*cell%lviscosl(ii) 
               avisgil_k=cell%alphag(ii)*a
               avisgit_k=f1(k)*cell%tviscosg(kk)*cell%alphag(kk)+f0(k)*cell%tviscosg(ii)*cell%alphag(ii)
               fluxg_diff_nf(i,1)=(avisgil_k*vl1+avisgit_k*vg1)*sv1
               fluxg_diff_nf(i,2)=(avisgil_k*vl2+avisgit_k*vg2)*sv1
            ENDDO
         ELSE
            nv=0
            nf_number=nf_number_id(nv)
            istart0=istart_nfs(nv)
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
!
               vl1=vl_o(kk,1)-vl_o(ii,1)
               vg1=vg_o(kk,1)-vg_o(ii,1)
               vd1=vd_o(kk,1)-vd_o(ii,1)
               vl2=vl_o(kk,2)-vl_o(ii,2)
               vg2=vg_o(kk,2)-vg_o(ii,2)
               vd2=vd_o(kk,2)-vd_o(ii,2)
               vl3=vl_o(kk,3)-vl_o(ii,3)
               vg3=vg_o(kk,3)-vg_o(ii,3)
               vd3=vd_o(kk,3)-vd_o(ii,3)
               sv1=sap_nf(i1)
!
               a=f1(i)*cell%lviscosl(ii)+f0(i)*cell%lviscosl(kk)
               avislil_i=cell%alphal(ii)*a
               avislit_i=f1(i)*cell%tviscosl(ii)*cell%alphal(ii)+f0(i)*cell%tviscosl(kk)*cell%alphal(kk)
               fluxl_diff_nf(i0,1)=(avislil_i*vl1+avislit_i*vl1)*sv1
               fluxl_diff_nf(i0,2)=(avislil_i*vl2+avislit_i*vl2)*sv1
               fluxl_diff_nf(i0,3)=(avislil_i*vl3+avislit_i*vl3)*sv1
!
               a=f1(i)*cell%lviscosl(ii)+f0(i)*cell%lviscosl(kk) 
               avisgil_i=cell%alphag(ii)*a
               avisgit_i=f1(i)*cell%tviscosg(ii)*cell%alphag(ii)+f0(i)*cell%tviscosg(kk)*cell%alphag(kk)
               fluxg_diff_nf(i0,1)=(avisgil_i*vl1+avisgit_i*vg1)*sv1
               fluxg_diff_nf(i0,2)=(avisgil_i*vl2+avisgit_i*vg2)*sv1
               fluxg_diff_nf(i0,3)=(avisgil_i*vl3+avisgit_i*vg3)*sv1
!
               avisdi=f1(i)*cell%eviscosd(ii)*cell%alphad(ii)+f0(i)*cell%eviscosd(kk)*cell%alphad(kk)
               fluxd_diff_nf(i1,1)=(avisdi*vd1)*sv1
               fluxd_diff_nf(i1,2)=(avisdi*vd2)*sv1
               fluxd_diff_nf(i1,3)=(avisdi*vd3)*sv1
            ENDDO
!
            nv=-1
            nf_number=nf_number_id(nv)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               k=right_nb_k(i)
               ii=right_non(k)
               kk=left_nf(k)
!
               vl1=vl_o(kk,1)-vl_o(ii,1)
               vg1=vg_o(kk,1)-vg_o(ii,1)
               vd1=vd_o(kk,1)-vd_o(ii,1)
               vl2=vl_o(kk,2)-vl_o(ii,2)
               vg2=vg_o(kk,2)-vg_o(ii,2)
               vd2=vd_o(kk,2)-vd_o(ii,2)
               vl3=vl_o(kk,3)-vl_o(ii,3)
               vg3=vg_o(kk,3)-vg_o(ii,3)
               vd3=vd_o(kk,3)-vd_o(ii,3)
               sv1=sap_nf(k)
!
               a=f1(k)*cell%lviscosl(kk)+f0(k)*cell%lviscosl(ii)
               avislil_k=cell%alphal(ii)*a
               avislit_k=f1(k)*cell%tviscosl(kk)*cell%alphal(kk)+f0(k)*cell%tviscosl(ii)*cell%alphal(ii)
               fluxl_diff_nf(i,1)=(avislil_k*vl1+avislit_k*vl1)*sv1
               fluxl_diff_nf(i,2)=(avislil_k*vl2+avislit_k*vl2)*sv1
               fluxl_diff_nf(i,3)=(avislil_k*vl3+avislit_k*vl3)*sv1
!
               a=f1(k)*cell%lviscosl(kk)+f0(k)*cell%lviscosl(ii) 
               avisgil_k=cell%alphag(ii)*a
               avisgit_k=f1(k)*cell%tviscosg(kk)*cell%alphag(kk)+f0(k)*cell%tviscosg(ii)*cell%alphag(ii)
               fluxg_diff_nf(i,1)=(avisgil_k*vl1+avisgit_k*vg1)*sv1
               fluxg_diff_nf(i,2)=(avisgil_k*vl2+avisgit_k*vg2)*sv1
               fluxg_diff_nf(i,3)=(avisgil_k*vl3+avisgit_k*vg3)*sv1
            ENDDO
         ENDIF
      ELSEIF(mdiffscheme.eq.2) THEN
         IF(ndim.eq.2)THEN
            nv=0
            nf_number=nf_number_id(nv)
            istart0=istart_nfs(nv)
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
!
               vl1=vl_o(kk,1)-vl_o(ii,1)
               vg1=vg_o(kk,1)-vg_o(ii,1)
               vd1=vd_o(kk,1)-vd_o(ii,1)
               vl2=vl_o(kk,2)-vl_o(ii,2)
               vg2=vg_o(kk,2)-vg_o(ii,2)
               vd2=vd_o(kk,2)-vd_o(ii,2)
               sv1=sap_nf(i1)
!
               a=f1(i)*cell%lviscosl(ii)+f0(i)*cell%lviscosl(kk)
               b=f1(i)*cell%tviscosl(ii)+f0(i)*cell%tviscosl(kk)
               avislil_i=cell%alphal(ii)*a
               avislit_i=cell%alphal(ii)*b
               fluxl_diff_nf(i0,1)=(avislil_i*vl1+avislit_i*vl1)*sv1
               fluxl_diff_nf(i0,2)=(avislil_i*vl2+avislit_i*vl2)*sv1
!
               avisgil_i=cell%alphag(ii)*a
               avisgit_i=cell%alphag(ii)*b
               fluxg_diff_nf(i0,1)=(avisgil_i*vl1+avisgit_i*vl1)*sv1
               fluxg_diff_nf(i0,2)=(avisgil_i*vl2+avisgit_i*vl2)*sv1
!
               avisdi=f1(i)*cell%eviscosd(ii)*cell%alphad(ii)+f0(i)*cell%eviscosd(kk)*cell%alphad(kk)
               fluxd_diff_nf(i1,1)=(avisdi*vd1)*sv1
               fluxd_diff_nf(i1,2)=(avisdi*vd2)*sv1
            ENDDO
!
            nv=-1
            nf_number=nf_number_id(nv)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               k=right_nb_k(i)
               ii=right_non(k)
               kk=left_nf(k)
!
               vl1=vl_o(kk,1)-vl_o(ii,1)
               vg1=vg_o(kk,1)-vg_o(ii,1)
               vd1=vd_o(kk,1)-vd_o(ii,1)
               vl2=vl_o(kk,2)-vl_o(ii,2)
               vg2=vg_o(kk,2)-vg_o(ii,2)
               vd2=vd_o(kk,2)-vd_o(ii,2)
               sv1=sap_nf(k)
!
               a=f1(k)*cell%lviscosl(kk)+f0(k)*cell%lviscosl(ii)
               b=f1(k)*cell%tviscosl(kk)+f0(k)*cell%tviscosl(ii)
               avislil_k=cell%alphal(ii)*a
               avislit_k=cell%alphal(ii)*b
               fluxl_diff_nf(i,1)=(avislil_k*vl1+avislit_k*vl1)*sv1
               fluxl_diff_nf(i,2)=(avislil_k*vl2+avislit_k*vl2)*sv1
!
               avisgil_k=cell%alphag(ii)*a
               avisgit_k=cell%alphag(ii)*b
               fluxg_diff_nf(i,1)=(avisgil_k*vl1+avisgit_k*vl1)*sv1
               fluxg_diff_nf(i,2)=(avisgil_k*vl2+avisgit_k*vl2)*sv1
            ENDDO
         ELSE
            nv=0
            nf_number=nf_number_id(nv)
            istart0=istart_nfs(nv)
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
!
               vl1=vl_o(kk,1)-vl_o(ii,1)
               vg1=vg_o(kk,1)-vg_o(ii,1)
               vd1=vd_o(kk,1)-vd_o(ii,1)
               vl2=vl_o(kk,2)-vl_o(ii,2)
               vg2=vg_o(kk,2)-vg_o(ii,2)
               vd2=vd_o(kk,2)-vd_o(ii,2)
               vl3=vl_o(kk,3)-vl_o(ii,3)
               vg3=vg_o(kk,3)-vg_o(ii,3)
               vd3=vd_o(kk,3)-vd_o(ii,3)
               sv1=sap_nf(i1)
!
               a=f1(i)*cell%lviscosl(ii)+f0(i)*cell%lviscosl(kk)
               b=f1(i)*cell%tviscosl(ii)+f0(i)*cell%tviscosl(kk)
               avislil_i=cell%alphal(ii)*a
               avislit_i=cell%alphal(ii)*b
               fluxl_diff_nf(i0,1)=(avislil_i*vl1+avislit_i*vl1)*sv1
               fluxl_diff_nf(i0,2)=(avislil_i*vl2+avislit_i*vl2)*sv1
               fluxl_diff_nf(i0,3)=(avislil_i*vl3+avislit_i*vl3)*sv1
!
               avisgil_i=cell%alphag(ii)*a
               avisgit_i=cell%alphag(ii)*b
               fluxg_diff_nf(i0,1)=(avisgil_i*vl1+avisgit_i*vl1)*sv1
               fluxg_diff_nf(i0,2)=(avisgil_i*vl2+avisgit_i*vl2)*sv1
               fluxg_diff_nf(i0,3)=(avisgil_i*vl3+avisgit_i*vl3)*sv1
!
               avisdi=f1(i)*cell%eviscosd(ii)*cell%alphad(ii)+f0(i)*cell%eviscosd(kk)*cell%alphad(kk)
               fluxd_diff_nf(i1,1)=(avisdi*vd1)*sv1
               fluxd_diff_nf(i1,2)=(avisdi*vd2)*sv1
               fluxd_diff_nf(i1,3)=(avisdi*vd3)*sv1
            ENDDO
!
            nv=-1
            nf_number=nf_number_id(nv)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               k=right_nb_k(i)
               ii=right_non(k)
               kk=left_nf(k)
!
               vl1=vl_o(kk,1)-vl_o(ii,1)
               vg1=vg_o(kk,1)-vg_o(ii,1)
               vd1=vd_o(kk,1)-vd_o(ii,1)
               vl2=vl_o(kk,2)-vl_o(ii,2)
               vg2=vg_o(kk,2)-vg_o(ii,2)
               vd2=vd_o(kk,2)-vd_o(ii,2)
               vl3=vl_o(kk,3)-vl_o(ii,3)
               vg3=vg_o(kk,3)-vg_o(ii,3)
               vd3=vd_o(kk,3)-vd_o(ii,3)
               sv1=sap_nf(k)
!
               a=f1(k)*cell%lviscosl(kk)+f0(k)*cell%lviscosl(ii)
               b=f1(k)*cell%tviscosl(kk)+f0(k)*cell%tviscosl(ii)
               avislil_k=cell%alphal(ii)*a
               avislit_k=cell%alphal(ii)*b
               fluxl_diff_nf(i,1)=(avislil_k*vl1+avislit_k*vl1)*sv1
               fluxl_diff_nf(i,2)=(avislil_k*vl2+avislit_k*vl2)*sv1
               fluxl_diff_nf(i,3)=(avislil_k*vl3+avislit_k*vl3)*sv1
!
               avisgil_k=cell%alphag(ii)*a
               avisgit_k=cell%alphag(ii)*b
               fluxg_diff_nf(i,1)=(avisgil_k*vl1+avisgit_k*vl1)*sv1
               fluxg_diff_nf(i,2)=(avisgil_k*vl2+avisgit_k*vl2)*sv1
               fluxg_diff_nf(i,3)=(avisgil_k*vl3+avisgit_k*vl3)*sv1
            ENDDO
         ENDIF
!
      ENDIF
!
!.....Cells inl
!
      nv=1
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      len   =istart_nf(2,nf_number)
      IF(mdiffscheme.eq.1) THEN
         DO ix=1,ndim
            DO i=1,len  
               i0=istart0+i
               i1=istart+i
               i2=istart2+i
               ii=left_nf(i1)
               k=nbcon_nf(i2)
!
               vl11=vl_o(ii,ix)
               vg11=vg_o(ii,ix)
               vd11=vd_o(ii,ix)
               avisli=alphab_liq(k)*cell%eviscosl(ii)
               avisgi=alphab_gas(k)*cell%eviscosg(ii)
               avisdi=alphab_drp(k)*cell%eviscosd(ii)
               f_profile=vel_bc_profile_inl(i)
               IF(vin_norm(k).eq.0)THEN
                  vl21=vb_liq(k,ix)*f_profile
                  vg21=vb_gas(k,ix)*f_profile
                  vd21=vb_drp(k,ix)*f_profile
               ELSE
                  vl21=vin_liq(k)*xn_nf(i1,ix)*f_profile
                  vg21=vin_gas(k)*xn_nf(i1,ix)*f_profile
                  vd21=vin_drp(k)*xn_nf(i1,ix)*f_profile
               ENDIF
               vl1=vl21-vl11
               vg1=vg21-vg11
               vd1=vd21-vd11
               sv1=sap_nf(i1)
               avislil=alphab_liq(k)*cell%lviscosl(ii)
               avislit=alphab_liq(k)*cell%tviscosl(ii)
               avisgil=alphab_gas(k)*cell%lviscosl(ii)
               avisgit=alphab_gas(k)*cell%tviscosg(ii)
               fluxl_diff_nf(i1,ix)=(avislil*vl1+avislit*vl1)*sv1
               fluxg_diff_nf(i1,ix)=(avisgil*vl1+avisgit*vg1)*sv1
               fluxd_diff_nf(i1,ix)=(avisdi*vd1)*sv1
            ENDDO
         ENDDO
      ELSEIF(mdiffscheme.eq.2) THEN
         DO ix=1,ndim
            DO i=1,len  
               i0=istart0+i
               i1=istart+i
               i2=istart2+i
               ii=left_nf(i1)
               k=nbcon_nf(i2)
!
               vl11=vl_o(ii,ix)
               vg11=vg_o(ii,ix)
               vd11=vd_o(ii,ix)
               avisli=alphab_liq(k)*cell%eviscosl(ii)
               avisgi=alphab_gas(k)*cell%eviscosg(ii)
               avisdi=alphab_drp(k)*cell%eviscosd(ii)
               f_profile=vel_bc_profile_inl(i)
               IF(vin_norm(k).eq.0)THEN
                  vl21=vb_liq(k,ix)*f_profile
                  vg21=vb_gas(k,ix)*f_profile
                  vd21=vb_drp(k,ix)*f_profile
               ELSE
                  vl21=vin_liq(k)*xn_nf(i1,ix)*f_profile
                  vg21=vin_gas(k)*xn_nf(i1,ix)*f_profile
                  vd21=vin_drp(k)*xn_nf(i1,ix)*f_profile
               ENDIF
               vl1=vl21-vl11
               vg1=vg21-vg11
               vd1=vd21-vd11
               sv1=sap_nf(i1)
               avislil=alphab_liq(k)*cell%lviscosl(ii)
               avislit=alphab_liq(k)*cell%tviscosl(ii)
               avisgil=alphab_gas(k)*cell%lviscosl(ii)
               avisgit=alphab_gas(k)*cell%tviscosg(ii)
               fluxl_diff_nf(i1,ix)=(avislil*vl1+avislit*vl1)*sv1
               fluxg_diff_nf(i1,ix)=(avisgil*vl1+avisgit*vl1)*sv1
               fluxd_diff_nf(i1,ix)=(avisdi*vd1)*sv1
            ENDDO
         ENDDO
      ENDIF
!
!.....The rest
!
      DO nv=2,5
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(nv)
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         IF(Twall_Model.eq.Laminar)THEN
            IF(nf_number.eq.5) THEN
               DO ix=1,ndim
                  DO i=1,len  
                     i0=istart0+i
                     i1=istart+i
                     ii=left_nf(i1)
                     vl1=v_wall(ix)-vl_o(ii,ix)
                     vg1=v_wall(ix)-vg_o(ii,ix)
                     vd1=v_wall(ix)-vd_o(ii,ix)
                     sv1=sap_nf(i1)/fac_fsw(i)
                     fluxl_diff_nf(i0,ix)=(cell%alphal(ii)*cell%lviscosl(ii)*vl1+cell%alphal(ii)*cell%tviscosl(ii)*vl1)*sv1
                     fluxg_diff_nf(i0,ix)=(cell%alphag(ii)*cell%lviscosl(ii)*vl1+cell%alphag(ii)*cell%tviscosg(ii)*vg1)*sv1
                     fluxd_diff_nf(i1,ix)=(cell%alphad(ii)*cell%eviscosd(ii)*vd1)*sv1
                  ENDDO
               ENDDO
            ELSE
               DO ix=1,ndim
                  DO i=1,len  
                     i0=istart0+i
                     i1=istart+i
                     ii=left_nf(i1)
                     vl1=v_wall(ix)-vl_o(ii,ix)
                     vg1=v_wall(ix)-vg_o(ii,ix)
                     vd1=v_wall(ix)-vd_o(ii,ix)
                     sv1=sap_nf(i1)
                     fluxl_diff_nf(i0,ix)=(cell%alphal(ii)*cell%lviscosl(ii)*vl1+cell%alphal(ii)*cell%tviscosl(ii)*vl1)*sv1
                     fluxg_diff_nf(i0,ix)=(cell%alphag(ii)*cell%lviscosl(ii)*vl1+cell%alphag(ii)*cell%tviscosg(ii)*vg1)*sv1
                     fluxd_diff_nf(i1,ix)=(cell%alphad(ii)*cell%eviscosd(ii)*vd1)*sv1
                  ENDDO
               ENDDO
            ENDIF
         ELSE
            IF(mdiffscheme.eq.1) THEN
               DO ix=1,ndim
                  DO i=1,len  
                     i0=istart0+i
                     i1=istart+i
                     ii=left_nf(i1)
                     vl1=v_wall(ix)-vl_o(ii,ix)
                     vg1=v_wall(ix)-vg_o(ii,ix)
                     vd1=v_wall(ix)-vd_o(ii,ix)
                     sv1=sap_nf(i1)*wallnr(ii)
                     IF(lowreynolds.ge.1) THEN
                        viscl=cell%eviscosl(ii)
                        viscg=cell%eviscosg(ii)
                     ELSE
                        viscl=wvis_liq(ii)
                        viscg=wvis_gas(ii)
                     ENDIF
!                     
                     fluxl_diff_nf(i0,ix)=(cell%alphal(ii)*viscl*vl1)*sv1
                     fluxg_diff_nf(i0,ix)=(cell%alphag(ii)*cell%lviscosl(ii)*vl1+cell%alphag(ii)*(viscg-cell%lviscosg(ii))*vg1)*sv1
                     fluxd_diff_nf(i1,ix)=(cell%alphad(ii)*viscl*vd1)*sv1
                  ENDDO
               ENDDO
            ELSEIF(mdiffscheme.eq.2) THEN
               DO ix=1,ndim
                  DO i=1,len  
                     i0=istart0+i
                     i1=istart+i
                     ii=left_nf(i1)
                     vl1=v_wall(ix)-vl_o(ii,ix)
                     vd1=v_wall(ix)-vd_o(ii,ix)
                     sv1=sap_nf(i1)*wallnr(ii)
                     IF(lowreynolds.ge.1) THEN
                        viscl=cell%eviscosl(ii)
                        viscg=cell%eviscosg(ii)
                     ELSE
                        viscl=wvis_liq(ii)
                        viscg=wvis_gas(ii)
                     ENDIF
!                     
                     fluxl_diff_nf(i0,ix)=(cell%alphal(ii)*viscl*vl1)*sv1
                     fluxg_diff_nf(i0,ix)=(cell%alphag(ii)*viscl*vl1)*sv1
                     fluxd_diff_nf(i1,ix)=(cell%alphad(ii)*viscl*vd1)*sv1
                  ENDDO
               ENDDO
             ENDIF
         ENDIF
      ENDDO
!
      CALL sum_nf_ndim(0,0,ncell_fluid,     &
                       fluxl_diff_nf,diff_liq, &
                       fluxg_diff_nf,diff_gas)
!
!.....Build summation info for non,inl,adw,fsw,ctw,chw
!
      nf_number_nb=5
      nf_number_id(0)=0
      nf_number_id(1)=2
      nf_number_id(2)=4
      nf_number_id(3)=5
      nf_number_id(4)=6
      nf_number_id(5)=7
      istart_nfs(0)=0
      istart_nfs(1)=istart_nfs(0)+nf_non
      istart_nfs(2)=istart_nfs(1)+nf_inl
      istart_nfs(3)=istart_nfs(2)+nf_adw
      istart_nfs(4)=istart_nfs(3)+nf_fsw
      istart_nfs(5)=istart_nfs(4)+nf_ctw
      lens         =istart_nfs(5)+nf_chw
!
      CALL sum_nf_ndim(0,-1,ncell_fluid,       &
                       fluxd_diff_nf,diff_drp)
!
      IF(i_droplet.ge.1) diff_drp(:,:)=0.0d0
!
      END SUBROUTINE momentum_diffusion_user
