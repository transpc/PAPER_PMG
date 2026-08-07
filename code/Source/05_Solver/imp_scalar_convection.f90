!
      SUBROUTINE imp_scalar_convection(diag_g,diag_l,diag_x,src_g,src_l,src_x,             &
                                       off_diag_g_non_i,off_diag_l_non_i,off_diag_x_non_i, &
                                       off_diag_g_non_k,off_diag_l_non_k,off_diag_x_non_k, &
                                       iter)      
!
!     Implicit scalar convection
!
      USE Zinterface
      USE VOL_DATA      , ONLY: cell
      USE Zzone         , ONLY: ncell_fluid
      USE Zcore         , ONLY: np
      USE Zvec_param    , ONLY: nf_nonk,nf_non,nf_mcc,nf_inl,nf_fluxk1
      USE Znum_cell     , ONLY: istart_nf,istart_nbcon_nf,        &
                                nf_number_nb,lens,                &
                                right_nb_k,                       &
                                istart_nfs,nf_number_id,istart_nf
      USE Zvec_index    , ONLY: left_nf,right_non,nbcon_nf
      USE Z2nd_order    , ONLY: eng_conv_2nd, qula_conv_2nd
      USE Zare          , ONLY: ar_gas,ar_liq
      USE Zb_condition  , ONLY: alphab_gas,alphab_liq,rhob_liq,rhob_gas
      USE Zvec_major    , ONLY: flux_l_nf,flux_g_nf
      USE Zb_condition  , ONLY: eb_gas,eb_liq,qualab
      USE Zvec_scalar   , ONLY: arli_nf,argi_nf
      USE c3com_cupid   , ONLY: i3invtbl,c3dpv
!
      IMPLICIT NONE
!      
      INCLUDE '../10_LinkToMARS/c3com.h'
!
!.....Input
      INTEGER iter
!.....Output
      REAL(8),DIMENSION(nf_non) :: off_diag_g_non_i,off_diag_l_non_i,off_diag_x_non_i
      REAL(8),DIMENSION(nf_nonk) :: off_diag_g_non_k,off_diag_l_non_k,off_diag_x_non_k
      REAL(8),DIMENSION(ncell_fluid) :: diag_l,diag_g,diag_x, &
                                        src_l,src_g,src_x
!.....Local variables
      INTEGER :: i,k
      INTEGER :: ii,kk,idx
      INTEGER :: nv,nf_number,istart0,istart,len,istart2,i0,i1,i2
      REAL(8) :: arvl,arvg
      REAL(8) :: arvl0,arvl1
      REAL(8) :: arvg0,arvg1
!.....Local vector arrays 
      REAL(8),DIMENSION(nf_non) :: diag_l1_non,diag_g1_non
      REAL(8),DIMENSION(nf_fluxk1) :: diag_l_nf,diag_g_nf,diag_x_nf, &
                                      src_l_nf,src_g_nf,src_x_nf
!
      IF(np.gt.1) CALL communicate_1d(cell%drhogde, &
                                      cell%drholde, &
                                      cell%drhogdx)
      IF(np.gt.1) CALL communicate_1d(cell%eg_o, &
                                      cell%el_o, &
                                      cell%quala_o) 
!                                                                                               
!.....Scalar convection is not considered for wall, symmetry, and MARS code boundaries
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
         arvl0=-MIN(flux_l_nf(i1),0.d0)*ar_liq(kk)
         arvl1= MAX(flux_l_nf(i1),0.d0)*ar_liq(ii)
         diag_l_nf(i0)=arvl0
         IF(iter.eq.1) THEN
            src_l_nf(i0)=arvl0*(cell%el_o(kk)-cell%el_o(ii))
         ELSE
            src_l_nf(i0)=0.d0
         ENDIF
         off_diag_l_non_i(i)=-arvl0
!
         arvg0=-MIN(flux_g_nf(i1),0.d0)*ar_gas(kk)
         arvg1= MAX(flux_g_nf(i1),0.d0)*ar_gas(ii)
         diag_g_nf(i0)=arvg0
         diag_x_nf(i0)=arvg0
         IF(iter.eq.1) THEN
            src_g_nf(i0)=arvg0*(cell%eg_o(kk)-cell%eg_o(ii))
            src_x_nf(i0)=arvg0*(cell%quala_o(kk)-cell%quala_o(ii))
         ELSE
            src_g_nf(i0)=0.d0
            src_x_nf(i0)=0.d0
         ENDIF
         off_diag_g_non_i(i)=-arvg0
         off_diag_x_non_i(i)=-arvg0
      ENDDO
      IF(iter.eq.1)THEN
         IF(eng_conv_2nd.gt.0) THEN
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
               IF(flux_l_nf(i1).lt.0.d0)THEN
                  arvl=-flux_l_nf(i1)*ar_liq(kk)
                  diag_l1_non(i)=arvl
               ELSE
                  arvl=flux_l_nf(i1)*ar_liq(ii)
                  diag_l1_non(i)=-arvl !see sum_nf_nz_k_6v
               ENDIF
            ENDDO
         ENDIF
         IF(eng_conv_2nd.gt.0 .or. qula_conv_2nd.gt.0) THEN
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
               IF(flux_g_nf(i1).lt.0)THEN
                  arvg=-flux_g_nf(i1)*ar_gas(kk)
                  diag_g1_non(i)=arvg
               ELSE
                  arvg=flux_g_nf(i1)*ar_gas(ii)
                  diag_g1_non(i)=-arvg !see sum_nf_nz_k_6v
               ENDIF
            ENDDO
         ENDIF
      ENDIF
!
      nv=-1
      nf_number=nf_number_id(nv)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         k=right_nb_k(i)
         ii=right_non(k)
         kk=left_nf(k)
!
         arvl1= MAX(flux_l_nf(k),0.d0)*ar_liq(kk)
         arvg1= MAX(flux_g_nf(k),0.d0)*ar_gas(kk)
         diag_l_nf(i)=arvl1
         IF(iter.eq.1) THEN
            src_l_nf(i)=arvl1*(cell%el_o(kk)-cell%el_o(ii)) 
         ELSE
            src_l_nf(i)=0.d0
         ENDIF
!
         diag_g_nf(i)=arvg1
         diag_x_nf(i)=arvg1
         IF(iter.eq.1) THEN
            src_g_nf(i)=arvg1*(cell%eg_o(kk)-cell%eg_o(ii))       
            src_x_nf(i)=arvg1*(cell%quala_o(kk)-cell%quala_o(ii)) 
         ELSE
            src_g_nf(i)=0.d0
            src_x_nf(i)=0.d0
         ENDIF
!
         off_diag_l_non_k(i)=-arvl1
         off_diag_g_non_k(i)=-arvg1
         off_diag_x_non_k(i)=-arvg1
      ENDDO
!                                                                                               
!.....Source term contribution for MARS boundary                                         
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
         diag_l_nf(i0)=-flux_l_nf(i1)*arli_nf(i1)
         diag_g_nf(i0)=-flux_g_nf(i1)*argi_nf(i1)
         diag_x_nf(i0)=-flux_g_nf(i1)*argi_nf(i1)
         IF(iter.eq.1) THEN
            src_l_nf(i0)=-diag_l_nf(i0)*cell%el_o(ii)-flux_l_nf(i1)*c3dpv(idx,2)
            src_g_nf(i0)=-diag_g_nf(i0)*cell%eg_o(ii)-flux_g_nf(i1)*c3dpv(idx,1)
            src_x_nf(i0)=-diag_g_nf(i0)*cell%quala_o(ii)-flux_g_nf(i1)*c3dpv(idx,6)
         ELSE
            src_l_nf(i0)=0.0d0
            src_g_nf(i0)=0.0d0
            src_x_nf(i0)=0.0d0
         ENDIF
      ENDDO
!                                                                                               
!.....Source term contribution for inlet boundary                                         
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
!     
         arvl=flux_l_nf(i1)*alphab_liq(k)*rhob_liq(k)
         arvg=flux_g_nf(i1)*alphab_gas(k)*rhob_gas(k)
!
         diag_l_nf(i0)=-arvl
         diag_g_nf(i0)=-arvg
         diag_x_nf(i0)=-arvg
         IF(iter.eq.1) THEN
            src_l_nf(i0)=arvl*cell%el_o(ii)-arvl*eb_liq(k)
            src_g_nf(i0)=arvg*cell%eg_o(ii)-arvg*eb_gas(k)
            src_x_nf(i0)=arvg*cell%quala_o(ii)-arvg*qualab(k)
         ELSE
            src_l_nf(i0)=0.0d0
            src_g_nf(i0)=0.0d0
            src_x_nf(i0)=0.0d0
         ENDIF
      ENDDO
!
      CALL sum_nf(1,0,              &
                  diag_l_nf,diag_l, &
                  diag_g_nf,diag_g, &
                  diag_x_nf,diag_x, &
                  src_l_nf,src_l,   &
                  src_g_nf,src_g,   &
                  src_x_nf,src_x)
!
!.....2nd order convection
!
      IF(iter.eq.1)THEN
         IF(eng_conv_2nd.gt.0) CALL energy_2nd_conv_imp(diag_l1_non,diag_g1_non, &
                                                        src_l,src_g)      
         IF(qula_conv_2nd.gt.0) CALL quality_2nd_conv_imp(cell%quala_o,diag_g1_non,src_x) !'mom_2nd_conv'
      ENDIF   
!
      END SUBROUTINE imp_scalar_convection
