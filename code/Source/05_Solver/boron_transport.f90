!
      SUBROUTINE boron_transport
!
!     This routine solves the boron transport equation
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell
      USE Zzone        , ONLY: ncell_fluid
      USE Zconst2      , ONLY: dt
      USE Zvec_param   , ONLY: nf_non,nf_mcc,nf_inl,nf_out,nf_flux
      USE Znum_cell    , ONLY: istart_nf,                                &
                               nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zvec_index   , ONLY: left_nf
      USE Zimplicit    , ONLY: imp_boron_trans
      USE Ztimecon     , ONLY: alpha_min
      USE Zare         , ONLY: ar_liq
      USE Zcoord3      , ONLY: vol
      USE Zvec_major   , ONLY: flux_l_nf,lbor_conv_nf
!
      IMPLICIT NONE
!
!.....Local variable
      INTEGER :: nv,nf_number,len,istart0,istart,i0,i1
      INTEGER :: i,ii
      REAL(8) :: boron_flux
      REAL(8) :: a,b
!.....Local arrays
      REAL(8),DIMENSION(ncell_fluid) :: boron_fluxv
!.....Local vector arrays
      REAL(8),DIMENSION(nf_flux) :: boron_flux_nf
!
!.....Implicit boron transport
!
      IF(imp_boron_trans.eq.1)THEN
         ar_liq(:)=cell%alphal(:)*cell%rhol(:) 
         CALL imp_boron_transport    
         RETURN
      ENDIF
!
!.....Calculate boron convection
!
      CALL boron_convection
!
!.....Build summation info for non,mcc,inl,out
!
      nf_number_nb=3
      nf_number_id(0)=0
      nf_number_id(1)=1
      nf_number_id(2)=2
      nf_number_id(3)=3
      istart_nfs(0)=0
      istart_nfs(1)=istart_nfs(0)+nf_non
      istart_nfs(2)=istart_nfs(1)+nf_mcc
      istart_nfs(3)=istart_nfs(2)+nf_inl
      lens         =istart_nfs(3)+nf_out
!
!.....non,mcc,inl,out
!
      DO nv=0,nf_number_nb
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(nv)
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len  
            i0=istart0+i
            i1=istart+i
            ii=left_nf(i1)
            boron_flux_nf(i0)=-lbor_conv_nf(i1)*flux_l_nf(i1)
         ENDDO
      ENDDO
      CALL sum_nf(0,-1,                      &
                  boron_flux_nf,boron_fluxv)
!   
!.....Update boron concentration
!
      DO i=1,ncell_fluid
         IF(cell%alphag(i).gt.1.0-alpha_min)THEN
            cell%cboron(i)=cell%cboron_o(i)
         ELSE
            boron_flux=boron_fluxv(i)
            a=1.d0-cell%alphag_o(i)
            b=1.d0-cell%alphag(i)
            cell%cboron(i)= (a*cell%rhol_o(i)*cell%cboron_o(i) + dt/vol(i)*boron_flux) &
                           /(b*cell%rhol(i))
         ENDIF
      ENDDO
!
      END SUBROUTINE boron_transport
!
!---------------------------------------------------------------------------------------------------
!
      SUBROUTINE imp_boron_transport
!
!     This routine calculates implicit boron transport
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell
      USE Zmpi         , ONLY: au,ju_a,ia_a
      USE Zzone        , ONLY: ncell_fluid
      USE Zvec_param   , ONLY: nf_nonk,nf_non,nf_mcc,nf_inl,nf_fluxk1
      USE Zconst2      , ONLY: dt
      USE Zvec_index   , ONLY: left_nf,right_non,nbcon_nf
      USE Znum_cell    , ONLY: istart_nf,istart_nbcon_nf, &
                               nf_number_nb,lens,nf_number_id,istart_nfs, &
                               right_nb_k
      USE Zcoord3      , ONLY: volr
      USE Zboron       , ONLY: cboronb_liq
      USE Zare         , ONLY: ar_liq
      USE Zvec_major   , ONLY: flux_l_nf
      USE Zb_condition , ONLY: alphab_liq,rhob_liq
      USE Zimplicit    , ONLY: eps_imp_boron,max_iter_boron
      USE Z2nd_order   , ONLY: boron_conv_2nd
      USE c3com_cupid  , ONLY: i3invtbl,mcdirect,c3dpv
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,k,idx
      INTEGER :: ii,jj,kk
      INTEGER :: nv,nf_number,len,istart0,istart,istart2,i0,i1,i2
      REAL(8) :: vt
      REAL(8) :: ar2,arv1,arv2,arb2,arbv2
!.....Local arrays
      REAL(8),DIMENSION(ncell_fluid) :: diag_cb,src_cb
!.....Local vector array
      REAL(8),DIMENSION(nf_non) :: off_diag_non_i
      REAL(8),DIMENSION(nf_nonk) :: off_diag_non_k
      REAL(8),DIMENSION(nf_fluxk1) :: diag_nf
      REAL(8),DIMENSION(nf_mcc+nf_inl) :: src_nf
!
!
!.....Setup matrix
!
      vt=1.d0/dt
      DO i=1,ncell_fluid
         diag_cb(i)=ar_liq(i)*vt
         src_cb(i)=ar_liq(i)*vt*cell%cboron_o(i)
      ENDDO
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
!.....Computing Cell
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
         arv2=MIN(flux_l_nf(i1),0.d0)*ar_liq(kk)
         diag_nf(i0)=-arv2*volr(ii)
         off_diag_non_i(i)= arv2*volr(ii)
      ENDDO
!
      nv=-1
      nf_number=nf_number_id(nv)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         k=right_nb_k(i)
         ii=right_non(k)
         kk=left_nf(k)
         arv1=MAX(flux_l_nf(k),0.d0)*ar_liq(kk)
         diag_nf(i)=arv1*volr(ii)
         off_diag_non_k(i)=-arv1*volr(ii)
      ENDDO
!
!.....MARS interface
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
!         
         idx=i3invtbl(i)
         IF(mcdirect(idx).lt.0)THEN !from mars to cupid
            ar2=c3dpv(i,5)
            arb2=c3dpv(i,8)
         ELSE
            ar2=ar_liq(ii)
         ENDIF
!
         IF(flux_l_nf(i1).lt.0.d0) THEN !from mars
            arv2=flux_l_nf(i1)*ar2
            arbv2=flux_l_nf(i1)*arb2
            diag_nf(i0)=-arv2*volr(ii)
         ELSE
            diag_nf(i0)=0.d0
         ENDIF         
      ENDDO       
!
!.....Inlet
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
         ar2=alphab_liq(k)*rhob_liq(k)
         IF(flux_l_nf(i1).lt.0.d0) THEN
            arv2=flux_l_nf(i1)*ar2
            diag_nf(i0)=-arv2*volr(ii)
         ELSE
            diag_nf(i0)=0.d0
         ENDIF
      ENDDO
!
      CALL sum_nf(1,0,             &
                  diag_nf,diag_cb)
!
!.....Build summation info for mcc,inl
!
      nf_number_nb=1
      nf_number_id(0)=1
      nf_number_id(1)=2
      istart_nfs(0)=0
      istart_nfs(1)=istart_nfs(0)+nf_mcc
      lens         =istart_nfs(1)+nf_inl
!
!.....MARS interface
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
!         
         idx=i3invtbl(i)
         IF(mcdirect(idx).lt.0)THEN !from mars to cupid
            ar2=c3dpv(i,5)
            arb2=c3dpv(i,8)
         ELSE
            ar2=ar_liq(ii)
         ENDIF
!
         IF(flux_l_nf(i1).lt.0.d0) THEN !from mars
            arv2=flux_l_nf(i1)*ar2
            arbv2=flux_l_nf(i1)*arb2
            src_nf(i0)=-arbv2*volr(ii)
         ELSE
            src_nf(i0)=0.d0
         ENDIF         
      ENDDO       
!
!........Inlet
!                          
      nv=1
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      DO i=1,len  
         i0=istart0+i
         i1=istart+i
         i2=istart2+i
         ii=left_nf(i1)
         k=nbcon_nf(i2)
         ar2=alphab_liq(k)*rhob_liq(k)
         IF(flux_l_nf(i1).lt.0.d0) THEN
            arv2=flux_l_nf(i1)*ar2
            src_nf(i0)=-arv2*cboronb_liq(k)*volr(ii)
         ELSE
            src_nf(i0)=0.d0
         ENDIF
      ENDDO
!
      CALL sum_nf(1,-1,          &
                  src_nf,src_cb)
!
      DO i=1,ncell_fluid
         IF(cell%alphal(i).lt.1.d-15) THEN
            diag_cb(i)=1.d0
            src_cb(i)=cell%cboron_o(i)
         ENDIF         
      ENDDO
!     
      IF(boron_conv_2nd.eq.1) THEN
         CALL boron_2nd_conv_imp(src_cb)
      ENDIF
!		
!.....Build directly solverCSR  array here
!		
      CALL csr_build_a(diag_cb,off_diag_non_i,off_diag_non_k)
      DO i=1,ncell_fluid
         IF(cell%alphal(i).lt.1.d-15) THEN
            DO jj=ia_a(i),ju_a(i)-1
               au(jj)=0.d0
            ENDDO
            DO jj=ju_a(i)+1,ia_a(i+1)-1
               au(jj)=0.d0
            ENDDO
         ENDIF
      ENDDO
!
      CALL csr_cg_solvers_scalar(diag_cb,src_cb,cell%cboron,eps_imp_boron,max_iter_boron)
!
      END SUBROUTINE imp_boron_transport
