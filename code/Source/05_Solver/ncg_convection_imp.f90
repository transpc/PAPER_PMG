!
      SUBROUTINE ncg_convection_imp(diag,src,qn,nc,               &
                                    off_diag_non_i,off_diag_non_k)
!
!     Implicit scalar convection
!
      USE Zinterface
      USE VOL_DATA      , ONLY: cell
      USE Zmpi          , ONLY: ncell_fp
      USE Zzone         , ONLY: ncell_fluid
      USE Znum_cell     , ONLY: istart_nf,istart_nbcon_nf,                 &
                                nf_number_nb,lens,nf_number_id,istart_nfs, &
                                right_nb_k
      USE Zvec_param    , ONLY: nf_nonk,nf_non,nf_inl
      USE Zvec_index    , ONLY: left_nf,right_non,nbcon_nf
      USE Z2nd_order    , ONLY: ncg_conv_2nd
      USE Zare          , ONLY: ar_gas
      USE Zb_condition  , ONLY: alphab_gas,rhob_gas,qualab
      USE Zncg          , ONLY: qn_nvin,ncg_diff
      USE Zvec_major    , ONLY: flux_g_nf
      USE Zvec_geo      , ONLY: fac1_non,fac_non,sap_nf
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: nc
      REAL(8),DIMENSION(ncell_fp) :: qn
!.....Output
      REAL(8),DIMENSION(nf_non) :: off_diag_non_i
      REAL(8),DIMENSION(nf_nonk) :: off_diag_non_k
      REAL(8),DIMENSION(ncell_fluid) :: diag,src
!.....Local variables
      INTEGER :: i,k
      INTEGER :: ii,kk
      INTEGER :: nv,nf_number,len,istart0,istart,istart2,i0,i1,i2
      REAL(8) :: ardi
!.....Local vector array s
      REAL(8),DIMENSION(nf_nonk+nf_non+nf_inl) :: diag_nf
      REAL(8),DIMENSION(nf_non) :: temp_non,ardi_non
      REAL(8),DIMENSION(nf_inl) :: src_inl
!
      nf_number=0
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)
         IF(flux_g_nf(i1).lt.0.0d0)THEN
            temp_non(i)=flux_g_nf(i1)*ar_gas(kk)*cell%quala(kk)
            off_diag_non_i(i)=temp_non(i)
         ELSE
            temp_non(i)=flux_g_nf(i1)*ar_gas(ii)*cell%quala(ii)
            off_diag_non_i(i)=0.d0
         ENDIF
      ENDDO
!
      IF(ncg_conv_2nd.gt.0) CALL mult_ncg_2nd_conv_imp(qn,temp_non,src)
!
!.....Build summation info for non,inl
!
      nf_number_nb=1
      nf_number_id(-1)=-1
      nf_number_id(0)=0
      nf_number_id(1)=2
      istart_nfs(-1)=0
      istart_nfs(0)=istart_nfs(-1)+nf_nonk
      istart_nfs(1)=istart_nfs(0)+nf_non
      lens         =istart_nfs(1)+nf_inl
!
!.....Computing cells
!
      nv=0
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      IF(ncg_diff.gt.0)THEN
         DO i=1,len
            i0=istart0+i
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            ardi=fac1_non(i)*ar_gas(ii)*cell%mdiff(ii)+fac_non(i)*ar_gas(kk)*cell%mdiff(kk)
            ardi_non(i)=ardi*sap_nf(i1)
            off_diag_non_i(i)=off_diag_non_i(i)-ardi_non(i)*cell%quala(kk)
            IF(flux_g_nf(i1).lt.0.d0)THEN
               diag_nf(i0)=-temp_non(i)+ardi_non(i)*cell%quala(kk)
            ELSE
               diag_nf(i0)=ardi_non(i)*cell%quala(kk)
            ENDIF
         ENDDO
      ELSE
         DO i=1,len
            i0=istart0+i
            i1=istart+i
            IF(flux_g_nf(i1).lt.0.d0)THEN
               diag_nf(i0)=-temp_non(i)
            ELSE
               diag_nf(i0)=0.d0
            ENDIF
         ENDDO
      ENDIF
!
      nv=-1
      nf_number=nf_number_id(nv)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         k=right_nb_k(i)
         IF(flux_g_nf(k).lt.0.0d0)THEN
            off_diag_non_k(i)=0.d0
         ELSE
            off_diag_non_k(i)=-temp_non(k)
         ENDIF
      ENDDO
      IF(ncg_diff.gt.0)THEN
         DO i=1,len
            k=right_nb_k(i)
            kk=left_nf(k)
            off_diag_non_k(i)=off_diag_non_k(i)-ardi_non(k)*cell%quala(kk)
            IF(flux_g_nf(k).lt.0.d0)THEN
               diag_nf(i)=ardi_non(k)*cell%quala(kk)
            ELSE
               diag_nf(i)=temp_non(k)+ardi_non(k)*cell%quala(kk)
            ENDIF
         ENDDO
      ELSE
         DO i=1,len
            k=right_nb_k(i)
            IF(flux_g_nf(k).lt.0.d0)THEN
               diag_nf(i)=0.d0
            ELSE
               diag_nf(i)=temp_non(k)
            ENDIF
         ENDDO
      ENDIF
!
!.....Source term contribution for inlet boundary                                         
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
         diag_nf(i0)=-flux_g_nf(i1)*alphab_gas(k)*rhob_gas(k)*qualab(k)
      ENDDO
!
      CALL sum_nf(1,0,          &
                  diag_nf,diag)
!
!.....Build summation info for inl
!
      nf_number_nb=0
      nf_number_id(0)=2
      istart_nfs(0)=0
      lens         =istart_nfs(0)+nf_inl
!
      nv=0
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
         src_inl(i)=-flux_g_nf(i1)*alphab_gas(k)*rhob_gas(k)*qualab(k) &
                    *qn_nvin(k,nc)
      ENDDO
!
      CALL sum_nf(1,-1,      &
                  src_inl,src)
!
      END SUBROUTINE ncg_convection_imp
