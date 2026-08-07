!
      SUBROUTINE imp_mom_convection(diag_g,src_g,diag_l,src_l,iter,    &
                                    off_diag_g_non_i,off_diag_l_non_i, &
                                    off_diag_g_non_k,off_diag_l_non_k)
!
!     Implicit momentum convection
! Not tested for iter>1
      USE Zinterface
      USE Zzone         , ONLY: ncell_fluid,ncell_fluid_pad
      USE Zparam        , ONLY: ndim
      USE Zvec_param    , ONLY: nf_nonk,nf_non,nf_mcc,nf_inl,nf_flux1,nf_fluxk1
      USE Zvec_index    , ONLY: left_nf,right_non,nbcon_nf
      USE Znum_cell     , ONLY: istart_nf,istart_nbcon_nf,                  &
                                nf_number_nb,lens,nf_number_id,istart_nfs , &
                                right_nb_k
      USE Z2nd_order    , ONLY: mom_conv_2nd
      USE Zare          , ONLY: ar_gas,ar_liq
      USE Zb_condition  , ONLY: alphab_gas,alphab_liq,vb_gas,vb_liq,vin_gas,vin_liq,rhob_liq,rhob_gas
      USE Zbc_index     , ONLY: vin_norm
      USE Zuserdefined  , ONLY: vel_bc_profile_inl
      USE Zvector       , ONLY: vl_n,vg_n
      USE Zvec_major    , ONLY: flux_l_nf,flux_g_nf
      USE Zvec_scalar   , ONLY: arli_nf,argi_nf
      USE c3com_cupid   , ONLY: i3invtbl
      USE Zvec_geo      , ONLY: xn_nf
      USE Zrv_choke
      USE Zmcp   
      USE Zrv_model     , ONLY: rv_mcp,rv_choke,rv_valve
!
      IMPLICIT NONE
!
      INCLUDE '../10_LinkToMARS/c3com.h' 
!.....Input
      INTEGER :: iter
!.....Output
      REAL(8),DIMENSION(ncell_fluid) :: diag_g,diag_l
      REAL(8),DIMENSION(ncell_fluid_pad,ndim) :: src_g,src_l
      REAL(8),DIMENSION(nf_non) :: off_diag_g_non_i,off_diag_l_non_i
      REAL(8),DIMENSION(nf_nonk) :: off_diag_g_non_k,off_diag_l_non_k
!.....Local variables
      INTEGER :: i,k
      INTEGER :: ii,kk,idx,ix
      INTEGER :: nv,nf_number,istart0,istart,len,istart2,i0,i1,i2
      REAL(8) :: f_profile
      REAL(8) :: xn
      REAL(8) :: vgb,vlb
!.....Local vector arrays
      REAL(8),DIMENSION(nf_non) :: diag_g1_non,diag_l1_non
      REAL(8),DIMENSION(nf_fluxk1) :: diag_g_nf,diag_l_nf
      REAL(8),DIMENSION(nf_flux1,ndim) :: src_g_nf,src_l_nf
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
         diag_l_nf(i0)  =-MIN(flux_l_nf(i1),0.d0)*ar_liq(kk)
         diag_g_nf(i0)  =-MIN(flux_g_nf(i1),0.d0)*ar_gas(kk)
         off_diag_l_non_i(i)=off_diag_l_non_i(i)-diag_l_nf(i0)
         off_diag_g_non_i(i)=off_diag_g_non_i(i)-diag_g_nf(i0)
      ENDDO
!
      IF(mom_conv_2nd.gt.0) THEN
         DO i=1,len
            i0=istart0+i
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            IF(flux_l_nf(i1).lt.0)THEN
               diag_l1_non(i)  =-flux_l_nf(i1)*ar_liq(kk)
            ELSE
               diag_l1_non(i)  =-flux_l_nf(i1)*ar_liq(ii)
            ENDIF
            IF(flux_g_nf(i1).lt.0)THEN
               diag_g1_non(i)  =-flux_g_nf(i1)*ar_gas(kk)
            ELSE
               diag_g1_non(i)  =-flux_g_nf(i1)*ar_gas(ii)
            ENDIF
         ENDDO
      ENDIF
!
      IF(iter.gt.1) THEN
         DO ix=1,ndim
            DO i=1,len
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
               src_l_nf(i1,ix)=diag_l_nf(i0)*(vl_n(ii,ix)-vl_n(kk,ix))
               src_g_nf(i1,ix)=diag_g_nf(i0)*(vg_n(ii,ix)-vg_n(kk,ix))
            ENDDO
         ENDDO
      ENDIF
!
      nv=-1
      nf_number=nf_number_id(nv)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         k=right_nb_k(i)
         kk=left_nf(k)
         diag_l_nf(i)=MAX(flux_l_nf(k),0.d0)*ar_liq(kk)
         diag_g_nf(i)=MAX(flux_g_nf(k),0.d0)*ar_gas(kk)
         off_diag_l_non_k(i)=off_diag_l_non_k(i)-diag_l_nf(i)
         off_diag_g_non_k(i)=off_diag_g_non_k(i)-diag_g_nf(i)
      ENDDO
!      
!.....fluxBC: choke model, mcp model, valve model
!      
      IF(rv_valve.eq.1.or.rv_choke.eq.1.or.rv_mcp.eq.1) CALL fluxBC_convection_smac3(diag_g_nf,diag_l_nf,off_diag_l_non_i,off_diag_g_non_i,off_diag_l_non_k,off_diag_g_non_k,src_l,src_g)      
!
!.....Cells mcc
! 
      nv=1
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)
         idx=i3invtbl(i)
         diag_l_nf(i0)=-MIN(flux_l_nf(i1),0.d0)*argi_nf(i1)
         diag_g_nf(i0)=-MIN(flux_g_nf(i1),0.d0)*arli_nf(i1)
      ENDDO
      DO ix=1,ndim
         DO i=1,len
            i0=istart0+i
            i1=istart+i
            idx=i3invtbl(i)
            src_l_nf(i1,ix)=diag_l_nf(i0)*c3vl(1,idx)*xn_nf(i1,ix)
            src_g_nf(i1,ix)=diag_g_nf(i0)*c3vg(1,idx)*xn_nf(i1,ix)
         ENDDO
      ENDDO
!
!.....Cells inl
!
      nv=2
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
         diag_l_nf(i0)=-MIN(flux_l_nf(i1),0.d0)*(alphab_liq(k)*rhob_liq(k))
         diag_g_nf(i0)=-MIN(flux_g_nf(i1),0.d0)*(alphab_gas(k)*rhob_gas(k))
      ENDDO
      IF(iter.eq.1) then
         DO ix=1,ndim
            DO i=1,len
               i0=istart0+i
               i1=istart+i
               i2=istart2+i
               ii=left_nf(i1)
               k=nbcon_nf(i2)
               f_profile=vel_bc_profile_inl(i)
               IF(vin_norm(k).eq.0)THEN
                  vlb=vb_liq(k,ix)*f_profile
                  vgb=vb_gas(k,ix)*f_profile
               ELSE
                  xn=xn_nf(i1,ix)*f_profile
                  vlb=vin_liq(k)*xn
                  vgb=vin_gas(k)*xn
               ENDIF
               src_l_nf(i1,ix)=diag_l_nf(i0)*vlb
               src_g_nf(i1,ix)=diag_g_nf(i0)*vgb
            ENDDO
         ENDDO
      ELSE
         DO ix=1,ndim
            DO i=1,len
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               src_l_nf(i1,ix)=diag_l_nf(i0)*vl_n(ii,ix)
               src_g_nf(i1,ix)=diag_g_nf(i0)*vg_n(ii,ix)
            ENDDO
         ENDDO
      ENDIF
!
      CALL sum_nf(1,0,              &
                  diag_g_nf,diag_g, &
                  diag_l_nf,diag_l)
!
      IF(iter.eq.1) then
!
!........Build summation info for mcc,inl
!
         nf_number_nb=1
         nf_number_id(0)=1
         nf_number_id(1)=2
         istart_nfs(0)=nf_non
         istart_nfs(1)=istart_nfs(0)+nf_mcc
         lens         =istart_nfs(1)+nf_inl
!
         CALL sum_nf_ndim(1,-1,ncell_fluid_pad, &
                          src_g_nf,src_g,       &
                          src_l_nf,src_l)
!
         IF(mom_conv_2nd.gt.0) CALL mom_2nd_conv_imp(diag_g1_non,diag_l1_non, &
                                                     src_g,src_l)
      ELSE ! iter
!
!........Build summation info for non,mcc,inl
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
         CALL sum_nf_ndim(1,-1,ncell_fluid_pad, &
                          src_g_nf,src_g,       &
                          src_l_nf,src_l)
      ENDIF ! iter
!
      END SUBROUTINE imp_mom_convection
