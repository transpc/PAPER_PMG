!
      SUBROUTINE imp_mom_diffusion(diag_g,diag_l,src_g,src_l,iter,    &
                                   off_diag_g_non_i,off_diag_g_non_k, &
                                   off_diag_l_non_i,off_diag_l_non_k)
!
!     Implicit momentum diffusion      
!      Not tested for iter>1
      USE Zinterface
      USE VOL_DATA      , ONLY: cell
      USE Zmpi          , ONLY: ncell_fp
      USE Zzone         , ONLY: ncell_fluid,ncell_fluid_pad
      USE Zcore         , ONLY: np
      USE Zparam        , ONLY: ndim
      USE Zvec_param    , ONLY: nf_nonk,nf_non,nf_inl,nf_adw,nf_fsw,nf_ctw,nf_chw
      USE Znum_cell     , ONLY: istart_nf,istart_nbcon_nf,                 &
                                nf_number_nb,lens,nf_number_id,istart_nfs, &
                                right_nb_k
      USE Zvec_index    , ONLY: left_nf,right_non,nbcon_nf
      USE Zb_condition  , ONLY: alphab_gas,alphab_liq,vb_gas,vb_liq,vin_gas,vin_liq,v_wall
      USE Zconst1       , ONLY: lowreynolds
      USE Zface         , ONLY: twall_model,laminar,Free_slip
      USE Zgradoption   , ONLY: non_orth_diff
      USE Zturb         , ONLY: wvis_liq,wvis_gas,wallnr
      USE Zvector       , ONLY: vl_n,vg_n,vg_o,vl_o
      USE Zuserdefined  , ONLY: vel_bc_profile_inl
      USE Zturb         , ONLY: s_macroturb_source      
      USE Zbc_index     , ONLY: vin_norm
      USE Zvec_geo      , ONLY: xn_nf,sap_nf,sa_nf,       &
                                fac1_non,fac_non,fac_fsw, &
                                djir_non,dnj_non
      USE Zrv_model     , ONLY: rv_mcp,rv_choke,rv_valve
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: iter
!.....Output
      REAL(8),DIMENSION(ncell_fluid) :: diag_g,diag_l
      REAL(8),DIMENSION(ncell_fluid_pad,ndim) :: src_g,src_l
      REAL(8),DIMENSION(nf_non) :: off_diag_g_non_i,off_diag_l_non_i
      REAL(8),DIMENSION(nf_nonk) :: off_diag_g_non_k,off_diag_l_non_k
!.....Local variables
      INTEGER :: i,k
      INTEGER :: ii,kk
      INTEGER :: nv,nf_number,istart0,istart,len,istart2,i0,i1,i2
      REAL(8) :: dvli1,dvli2,dvli3
      REAL(8) :: dvgi1,dvgi2,dvgi3
      REAL(8) :: cf_g,cf_l,sv1
      REAL(8) :: vl1,vl2,vl3,vg1,vg2,vg3
      REAL(8) :: dvlj,dvgj
      REAL(8) :: viscl,viscg      
      REAL(8) :: f_profile
!.....Local arrays
      REAL(8),DIMENSION(ncell_fp,ndim,ndim) :: dvldx,dvgdx
!.....Local vector arrays
      REAL(8),DIMENSION(nf_non+nf_inl+nf_adw+nf_fsw+nf_ctw+nf_chw) :: diag_g_nf,diag_l_nf,cf_g_nf,cf_l_nf
      REAL(8),DIMENSION(nf_non+nf_inl+nf_adw+nf_fsw+nf_ctw+nf_chw,ndim) :: src_g_nf,src_l_nf
!
!.....Calculate velocity gradient at cell center for non-orthogonal grid
!
      IF(iter.eq.1) THEN
         IF(non_orth_diff.eq.1)THEN
            CALL grad_vel(2,vl_o,dvldx,vb_liq,vin_liq)
            CALL grad_vel(1,vg_o,dvgdx,vb_gas,vin_gas)
         ELSEIF(non_orth_diff.eq.2)THEN
            CALL grad_frink_vel(vl_o,dvldx,vb_liq,vin_liq)
            CALL grad_frink_vel(vg_o,dvgdx,vb_gas,vin_gas)
         ENDIF
         IF(np.gt.1) THEN
            IF(non_orth_diff.eq.1 .or. non_orth_diff.eq.2)THEN
               CALL communicate_3d(dvldx, &
                                   dvgdx)
            ENDIF
         ENDIF
      ENDIF
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
!.....Cells non
!
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
         cf_g_nf(i0)=(fac1_non(i)*cell%alphag(ii)*cell%eviscosg(ii)+fac_non(i)*cell%alphag(kk)*cell%eviscosg(kk))*sa_nf(i1)
         cf_l_nf(i0)=(fac1_non(i)*cell%alphal(ii)*cell%eviscosl(ii)+fac_non(i)*cell%alphal(kk)*cell%eviscosl(kk))*sa_nf(i1)
      ENDDO
!
!.....Cells inl
!
      nv=1
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i0=istart0+i
         i1=istart+i
         i2=istart2+i
         ii=left_nf(i1)
         k=nbcon_nf(i2)
         cf_g_nf(i0)=alphab_gas(k)*cell%eviscosg(ii)*sap_nf(i1)
         cf_l_nf(i0)=alphab_liq(k)*cell%eviscosl(ii)*sap_nf(i1)
      ENDDO
!
      IF(Twall_Model.eq.Laminar)THEN
         DO nv=2,5
            nf_number=nf_number_id(nv)
            istart0=istart_nfs(nv)
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len  
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               sv1=sap_nf(i1)
               IF(nf_number.eq.5) sv1=sap_nf(i1)/fac_fsw(i)
               cf_g_nf(i0)=cell%alphag(ii)*cell%eviscosg(ii)*sv1
               cf_l_nf(i0)=cell%alphal(ii)*cell%eviscosl(ii)*sv1
            ENDDO
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
               DO i=1,len  
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  sv1=sap_nf(i1)
                  viscl=cell%lviscosl(ii)
                  viscg=cell%lviscosg(ii)
                  cf_l_nf(i0)=cell%alphal(ii)*viscl*sv1
                  cf_g_nf(i0)=cell%alphag(ii)*viscg*sv1
               ENDDO
            ENDDO
         ELSE
            IF(lowreynolds.ge.1) THEN
               DO nv=2,5
                  nf_number=nf_number_id(nv)
                  istart0=istart_nfs(nv)
                  istart=istart_nf(1,nf_number)
                  len   =istart_nf(2,nf_number)
                  DO i=1,len  
                     i0=istart0+i
                     i1=istart+i
                     ii=left_nf(i1)
                     sv1=sa_nf(i1)*wallnr(ii)
                     viscl=cell%eviscosl(ii)
                     viscg=cell%eviscosg(ii)
                     cf_l_nf(i0)=cell%alphal(ii)*viscl*sv1
                     cf_g_nf(i0)=cell%alphag(ii)*viscg*sv1
                  ENDDO
               ENDDO
            ELSEIF(Twall_Model.eq.Free_slip)THEN
               DO nv=2,5
                  nf_number=nf_number_id(nv)
                  istart0=istart_nfs(nv)
                  len   =istart_nf(2,nf_number)
                  DO i=1,len  
                     i0=istart0+i
                     cf_l_nf(i0)=0.0d0
                     cf_g_nf(i0)=0.0d0
                  ENDDO
               ENDDO         
            ELSE
               DO nv=2,5
                  nf_number=nf_number_id(nv)
                  istart0=istart_nfs(nv)
                  istart=istart_nf(1,nf_number)
                  len   =istart_nf(2,nf_number)
                  DO i=1,len  
                     i0=istart0+i
                     i1=istart+i
                     ii=left_nf(i1)
                     sv1=sa_nf(i1)*wallnr(ii)
                     viscl=wvis_liq(ii)
                     viscg=wvis_gas(ii)
                     cf_l_nf(i0)=cell%alphal(ii)*viscl*sv1
                     cf_g_nf(i0)=cell%alphag(ii)*viscg*sv1
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
      ENDIF
!
      IF(iter.eq.1)THEN
!
!........Cells non
!
         nv=0
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(nv)
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
!
!........Non-orthogonal grid contribution
!
         IF(non_orth_diff.gt.0)THEN
            IF(ndim.eq.2) THEN
               DO i=1,len  
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  kk=right_non(i)
!
                  cf_g=cf_g_nf(i0)
                  cf_l=cf_l_nf(i0)
                  dvli1=fac1_non(i)*dvldx(ii,1,1)+fac_non(i)*dvldx(kk,1,1)
                  dvgi1=fac1_non(i)*dvgdx(ii,1,1)+fac_non(i)*dvgdx(kk,1,1)
                  dvli2=fac1_non(i)*dvldx(ii,2,1)+fac_non(i)*dvldx(kk,2,1)
                  dvgi2=fac1_non(i)*dvgdx(ii,2,1)+fac_non(i)*dvgdx(kk,2,1)
                  dvlj=dvli1*dnj_non(i,1)+dvli2*dnj_non(i,2)
                  dvgj=dvgi1*dnj_non(i,1)+dvgi2*dnj_non(i,2)
                  src_l_nf(i0,1)=cf_l*dvlj
                  src_g_nf(i0,1)=cf_g*dvgj
!
                  dvli1=fac1_non(i)*dvldx(ii,1,2)+fac_non(i)*dvldx(kk,1,2)
                  dvgi1=fac1_non(i)*dvgdx(ii,1,2)+fac_non(i)*dvgdx(kk,1,2)
                  dvli2=fac1_non(i)*dvldx(ii,2,2)+fac_non(i)*dvldx(kk,2,2)
                  dvgi2=fac1_non(i)*dvgdx(ii,2,2)+fac_non(i)*dvgdx(kk,2,2)
                  dvlj=dvli1*dnj_non(i,1)+dvli2*dnj_non(i,2)
                  dvgj=dvgi1*dnj_non(i,1)+dvgi2*dnj_non(i,2)
                  src_l_nf(i0,2)=cf_l*dvlj
                  src_g_nf(i0,2)=cf_g*dvgj
!
                  cf_g=cf_g*djir_non(i)
                  cf_l=cf_l*djir_non(i)
                  diag_g_nf(i0)=cf_g
                  diag_l_nf(i0)=cf_l
                  off_diag_g_non_i(i)=-cf_g
                  off_diag_l_non_i(i)=-cf_l
               ENDDO
               DO i=1,nf_nonk
                  k=right_nb_k(i)
                  off_diag_g_non_k(i)=-diag_g_nf(k)
                  off_diag_l_non_k(i)=-diag_l_nf(k)
               ENDDO
            ELSE
               DO i=1,len  
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  kk=right_non(i)
!
                  cf_g=cf_g_nf(i0)
                  cf_l=cf_l_nf(i0)
                  dvli1=fac1_non(i)*dvldx(ii,1,1)+fac_non(i)*dvldx(kk,1,1)
                  dvgi1=fac1_non(i)*dvgdx(ii,1,1)+fac_non(i)*dvgdx(kk,1,1)
                  dvli2=fac1_non(i)*dvldx(ii,2,1)+fac_non(i)*dvldx(kk,2,1)
                  dvgi2=fac1_non(i)*dvgdx(ii,2,1)+fac_non(i)*dvgdx(kk,2,1)
                  dvli3=fac1_non(i)*dvldx(ii,3,1)+fac_non(i)*dvldx(kk,3,1)
                  dvgi3=fac1_non(i)*dvgdx(ii,3,1)+fac_non(i)*dvgdx(kk,3,1)
                  dvlj=dvli1*dnj_non(i,1)+dvli2*dnj_non(i,2)+dvli3*dnj_non(i,3)
                  dvgj=dvgi1*dnj_non(i,1)+dvgi2*dnj_non(i,2)+dvgi3*dnj_non(i,3)
                  src_l_nf(i0,1)=cf_l*dvlj
                  src_g_nf(i0,1)=cf_g*dvgj
!
                  dvli1=fac1_non(i)*dvldx(ii,1,2)+fac_non(i)*dvldx(kk,1,2)
                  dvgi1=fac1_non(i)*dvgdx(ii,1,2)+fac_non(i)*dvgdx(kk,1,2)
                  dvli2=fac1_non(i)*dvldx(ii,2,2)+fac_non(i)*dvldx(kk,2,2)
                  dvgi2=fac1_non(i)*dvgdx(ii,2,2)+fac_non(i)*dvgdx(kk,2,2)
                  dvli3=fac1_non(i)*dvldx(ii,3,2)+fac_non(i)*dvldx(kk,3,2)
                  dvgi3=fac1_non(i)*dvgdx(ii,3,2)+fac_non(i)*dvgdx(kk,3,2)
                  dvlj=dvli1*dnj_non(i,1)+dvli2*dnj_non(i,2)+dvli3*dnj_non(i,3)
                  dvgj=dvgi1*dnj_non(i,1)+dvgi2*dnj_non(i,2)+dvgi3*dnj_non(i,3)
                  src_l_nf(i0,2)=cf_l*dvlj
                  src_g_nf(i0,2)=cf_g*dvgj
!
                  dvli1=fac1_non(i)*dvldx(ii,1,3)+fac_non(i)*dvldx(kk,1,3)
                  dvgi1=fac1_non(i)*dvgdx(ii,1,3)+fac_non(i)*dvgdx(kk,1,3)
                  dvli2=fac1_non(i)*dvldx(ii,2,3)+fac_non(i)*dvldx(kk,2,3)
                  dvgi2=fac1_non(i)*dvgdx(ii,2,3)+fac_non(i)*dvgdx(kk,2,3)
                  dvli3=fac1_non(i)*dvldx(ii,3,3)+fac_non(i)*dvldx(kk,3,3)
                  dvgi3=fac1_non(i)*dvgdx(ii,3,3)+fac_non(i)*dvgdx(kk,3,3)
                  dvlj=dvli1*dnj_non(i,1)+dvli2*dnj_non(i,2)+dvli3*dnj_non(i,3)
                  dvgj=dvgi1*dnj_non(i,1)+dvgi2*dnj_non(i,2)+dvgi3*dnj_non(i,3)
                  src_l_nf(i0,3)=cf_l*dvlj
                  src_g_nf(i0,3)=cf_g*dvgj
!
                  cf_g=cf_g*djir_non(i)
                  cf_l=cf_l*djir_non(i)
                  diag_g_nf(i0)=cf_g
                  diag_l_nf(i0)=cf_l
                  off_diag_g_non_i(i)=-cf_g
                  off_diag_l_non_i(i)=-cf_l
               ENDDO
               DO i=1,nf_nonk
                  k=right_nb_k(i)
                  off_diag_g_non_k(i)=-diag_g_nf(k)
                  off_diag_l_non_k(i)=-diag_l_nf(k)
               ENDDO
            ENDIF 
         ELSE
            DO i=1,len  
               i0=istart0+i
               i1=istart+i
!
               cf_g=cf_g_nf(i0)
               cf_l=cf_l_nf(i0)
!
               cf_g=cf_g*djir_non(i)
               cf_l=cf_l*djir_non(i)
               diag_g_nf(i0)=cf_g
               diag_l_nf(i0)=cf_l
               off_diag_g_non_i(i)=-cf_g
               off_diag_l_non_i(i)=-cf_l
            ENDDO
            DO i=1,nf_nonk
               k=right_nb_k(i)
               off_diag_g_non_k(i)=-diag_g_nf(k)
               off_diag_l_non_k(i)=-diag_l_nf(k)
            ENDDO
            IF(ndim.eq.2) THEN
               DO i=1,len
                  i0=istart0+i
                  i1=istart+i
                  src_g_nf(i0,1)=0.d0
                  src_l_nf(i0,1)=0.d0
                  src_g_nf(i0,2)=0.d0
                  src_l_nf(i0,2)=0.d0
               ENDDO
            ELSE
               DO i=1,len
                  i0=istart0+i
                  i1=istart+i
                  src_g_nf(i0,1)=0.d0
                  src_l_nf(i0,1)=0.d0
                  src_g_nf(i0,2)=0.d0
                  src_l_nf(i0,2)=0.d0
                  src_g_nf(i0,3)=0.d0
                  src_l_nf(i0,3)=0.d0
               ENDDO
            ENDIF
         ENDIF 
!    
!........fluxBC: choke model, mcp model, valve model
!         
         IF(rv_valve.eq.1.or.rv_choke.eq.1.or.rv_mcp.eq.1) CALL fluxBC_diffusion_smac3(diag_l_nf,diag_g_nf,off_diag_l_non_i,off_diag_g_non_i,off_diag_l_non_k,off_diag_g_non_k,diag_l,diag_g,src_l,src_g,iter)
!
!........Cells inl
!
         nv=1
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(nv)
         istart=istart_nf(1,nf_number)
         istart2=istart_nbcon_nf(nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len  
            i0=istart0+i
            i1=istart+i
            cf_g=cf_g_nf(i0)
            cf_l=cf_l_nf(i0)
            diag_g_nf(i0)=cf_g
            diag_l_nf(i0)=cf_l
         ENDDO
         IF(ndim.eq.2) THEN
            DO i=1,len
               i0=istart0+i
               i1=istart+i
               i2=istart2+i
               k=nbcon_nf(i2)
               f_profile=vel_bc_profile_inl(i)
               IF(vin_norm(k).eq.0)THEN
                  vg1=vb_gas(k,1)*f_profile
                  vl1=vb_liq(k,1)*f_profile
                  vg2=vb_gas(k,2)*f_profile
                  vl2=vb_liq(k,2)*f_profile
               ELSE
                  vg1=vin_gas(k)*xn_nf(i1,1)*f_profile
                  vl1=vin_liq(k)*xn_nf(i1,1)*f_profile
                  vg2=vin_gas(k)*xn_nf(i1,2)*f_profile
                  vl2=vin_liq(k)*xn_nf(i1,2)*f_profile
               ENDIF
               cf_g=cf_g_nf(i0)
               cf_l=cf_l_nf(i0)
               src_g_nf(i0,1)=cf_g*vg1
               src_l_nf(i0,1)=cf_l*vl1
               src_g_nf(i0,2)=cf_g*vg2
               src_l_nf(i0,2)=cf_l*vl2
            ENDDO
         ELSE
            DO i=1,len
               i0=istart0+i
               i1=istart+i
               i2=istart2+i
               k=nbcon_nf(i2)
               f_profile=vel_bc_profile_inl(i)
               IF(vin_norm(k).eq.0)THEN
                  vg1=vb_gas(k,1)*f_profile
                  vl1=vb_liq(k,1)*f_profile
                  vg2=vb_gas(k,2)*f_profile
                  vl2=vb_liq(k,2)*f_profile
                  vg3=vb_gas(k,3)*f_profile
                  vl3=vb_liq(k,3)*f_profile
               ELSE
                  vg1=vin_gas(k)*xn_nf(i1,1)*f_profile
                  vl1=vin_liq(k)*xn_nf(i1,1)*f_profile
                  vg2=vin_gas(k)*xn_nf(i1,2)*f_profile
                  vl2=vin_liq(k)*xn_nf(i1,2)*f_profile
                  vg3=vin_gas(k)*xn_nf(i1,3)*f_profile
                  vl3=vin_liq(k)*xn_nf(i1,3)*f_profile
               ENDIF
               cf_g=cf_g_nf(i0)
               cf_l=cf_l_nf(i0)
               src_g_nf(i0,1)=cf_g*vg1
               src_l_nf(i0,1)=cf_l*vl1
               src_g_nf(i0,2)=cf_g*vg2
               src_l_nf(i0,2)=cf_l*vl2
               src_g_nf(i0,3)=cf_g*vg3
               src_l_nf(i0,3)=cf_l*vl3
            ENDDO
         ENDIF
!
!........Walls: adiabatic, constant temperature, constant heat flux, fluid-solid interface
!
         DO nv=2,5
            nf_number=nf_number_id(nv)
            istart0=istart_nfs(nv)
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len  
               i0=istart0+i
               i1=istart+i
               cf_g=cf_g_nf(i0)
               cf_l=cf_l_nf(i0)
               diag_g_nf(i0)=cf_g
               diag_l_nf(i0)=cf_l
            ENDDO
            IF(ndim.eq.2) THEN
               DO i=1,len
                  i0=istart0+i
                  i1=istart+i
                  cf_g=cf_g_nf(i0)
                  cf_l=cf_l_nf(i0)
                  src_g_nf(i0,1)=cf_g*v_wall(1)
                  src_l_nf(i0,1)=cf_l*v_wall(1)
                  src_g_nf(i0,2)=cf_g*v_wall(2)
                  src_l_nf(i0,2)=cf_l*v_wall(2)
               ENDDO
            ELSE
               DO i=1,len
                  i0=istart0+i
                  i1=istart+i
                  cf_g=cf_g_nf(i0)
                  cf_l=cf_l_nf(i0)
                  src_g_nf(i0,1)=cf_g*v_wall(1)
                  src_l_nf(i0,1)=cf_l*v_wall(1)
                  src_g_nf(i0,2)=cf_g*v_wall(2)
                  src_l_nf(i0,2)=cf_l*v_wall(2)
                  src_g_nf(i0,3)=cf_g*v_wall(3)
                  src_l_nf(i0,3)=cf_l*v_wall(3)
               ENDDO
            ENDIF
         ENDDO
!
         CALL sum_nf(1,1,              &
                     diag_g_nf,diag_g, &
                     diag_l_nf,diag_l)
         CALL sum_nf_ndim(1,1,ncell_fluid_pad, &
                          src_g_nf,src_g,      &
                          src_l_nf,src_l)
!
      ELSE  ! iter
!
!........Cells non
!
         nv=0
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(nv)
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len  
            i0=istart0+i
            i1=istart+i
!
            cf_g=cf_g_nf(i0)
            cf_l=cf_l_nf(i0)
!
            cf_g=cf_g*djir_non(i)
            cf_l=cf_l*djir_non(i)
            diag_g_nf(i0)=cf_g
            diag_l_nf(i0)=cf_l
            off_diag_g_non_i(i)=-cf_g
            off_diag_l_non_i(i)=-cf_l
         ENDDO
         DO i=1,nf_nonk
            k=right_nb_k(i)
            off_diag_g_non_k(i)=-diag_g_nf(k)
            off_diag_l_non_k(i)=-diag_l_nf(k)
         ENDDO
         IF(ndim.eq.2) THEN
            DO i=1,len
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
               src_g_nf(i0,1)=diag_g_nf(i0)*(vg_n(ii,1)-vg_n(kk,1))
               src_l_nf(i0,1)=diag_l_nf(i0)*(vl_n(ii,1)-vl_n(kk,1))
               src_g_nf(i0,2)=diag_g_nf(i0)*(vg_n(ii,2)-vg_n(kk,2))
               src_l_nf(i0,2)=diag_l_nf(i0)*(vl_n(ii,2)-vl_n(kk,2))
            ENDDO
         ELSE
            DO i=1,len
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
               src_g_nf(i0,1)=diag_g_nf(i0)*(vg_n(ii,1)-vg_n(kk,1))
               src_l_nf(i0,1)=diag_l_nf(i0)*(vl_n(ii,1)-vl_n(kk,1))
               src_g_nf(i0,2)=diag_g_nf(i0)*(vg_n(ii,2)-vg_n(kk,2))
               src_l_nf(i0,2)=diag_l_nf(i0)*(vl_n(ii,2)-vl_n(kk,2))
               src_g_nf(i0,3)=diag_g_nf(i0)*(vg_n(ii,3)-vg_n(kk,3))
               src_l_nf(i0,3)=diag_l_nf(i0)*(vl_n(ii,3)-vl_n(kk,3))
            ENDDO
         ENDIF
!
!........fluxBC: choke model, mcp model, valve model
!         
!         CALL fluxBC_diffusion_smac3(diag_l_nf,diag_g_nf,off_diag_l_non_i,off_diag_g_non_i,off_diag_l_non_k,off_diag_g_non_k,diag_l,diag_g,src_l,src_g,iter) !not yet (iter_mom is only 1.)
!
!........Cells inl
!
         nv=1
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(nv)
         istart=istart_nf(1,nf_number)
         istart2=istart_nbcon_nf(nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len  
            i0=istart0+i
            i1=istart+i
            i2=istart2+i
            ii=left_nf(i1)
            k=nbcon_nf(i2)
            cf_g=alphab_gas(k)*cell%eviscosg(ii)*sap_nf(i1)
            cf_l=alphab_liq(k)*cell%eviscosl(ii)*sap_nf(i1)
            diag_g_nf(i0)=cf_g
            diag_l_nf(i0)=cf_l
         ENDDO
         IF(ndim.eq.2) THEN
            DO i=1,len
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               src_g_nf(i0,1)=diag_g_nf(i0)*vg_n(ii,1)
               src_l_nf(i0,1)=diag_l_nf(i0)*vl_n(ii,1)
               src_g_nf(i0,2)=diag_g_nf(i0)*vg_n(ii,2)
               src_l_nf(i0,2)=diag_l_nf(i0)*vl_n(ii,2)
            ENDDO
         ELSE
            DO i=1,len
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               src_g_nf(i0,1)=diag_g_nf(i0)*vg_n(ii,1)
               src_l_nf(i0,1)=diag_l_nf(i0)*vl_n(ii,1)
               src_g_nf(i0,2)=diag_g_nf(i0)*vg_n(ii,2)
               src_l_nf(i0,2)=diag_l_nf(i0)*vl_n(ii,2)
               src_g_nf(i0,3)=diag_g_nf(i0)*vg_n(ii,3)
               src_l_nf(i0,3)=diag_l_nf(i0)*vl_n(ii,3)
            ENDDO
         ENDIF
!
!........Walls: adiabatic, constant temperature, constant heat flux, fluid-solid interface
!
         DO nv=2,5
            nf_number=nf_number_id(nv)
            istart0=istart_nfs(nv)
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len  
               i0=istart0+i
               i1=istart+i
               cf_g=cf_g_nf(i0)
               cf_l=cf_l_nf(i0)
               diag_g_nf(i0)=cf_g
               diag_l_nf(i0)=cf_l
            ENDDO
         ENDDO
         IF(ndim.eq.2) THEN
            DO i=1,len
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               src_g_nf(i0,1)=diag_g_nf(i0)*vg_n(ii,1)
               src_l_nf(i0,1)=diag_l_nf(i0)*vl_n(ii,1)
               src_g_nf(i0,2)=diag_g_nf(i0)*vg_n(ii,2)
               src_l_nf(i0,2)=diag_l_nf(i0)*vl_n(ii,2)
            ENDDO
         ELSE
            DO i=1,len
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               src_g_nf(i0,1)=diag_g_nf(i0)*vg_n(ii,1)
               src_l_nf(i0,1)=diag_l_nf(i0)*vl_n(ii,1)
               src_g_nf(i0,2)=diag_g_nf(i0)*vg_n(ii,2)
               src_l_nf(i0,2)=diag_l_nf(i0)*vl_n(ii,2)
               src_g_nf(i0,3)=diag_g_nf(i0)*vg_n(ii,3)
               src_l_nf(i0,3)=diag_l_nf(i0)*vl_n(ii,3)
            ENDDO
         ENDIF
!
         CALL sum_nf(1,1,              &
                     diag_g_nf,diag_g, &
                     diag_l_nf,diag_l)
         CALL sum_nf_ndim(1,1,ncell_fluid_pad, &
                          src_g_nf,src_g,      &
                          src_l_nf,src_l)
!
      ENDIF ! iter
!
      END SUBROUTINE imp_mom_diffusion
