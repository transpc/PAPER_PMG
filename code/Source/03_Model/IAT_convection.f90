!
      SUBROUTINE IAT_convection
!
!     This routine calculates IAT source by convection
!
      USE Zinterface
      USE Zcore       , ONLY: np
      USE Zvec_param  , ONLY: nf_non,nf_inl,nf_out
      USE Znum_cell   , ONLY: istart_nf,istart_nbcon_nf, &
                              nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zvec_index  , ONLY: left_nf,right_non,nbcon_nf
      USE Ziat        , ONLY: ia_b,ia_conv,ia_old
      USE Zvec_major  , ONLY: flux_g_nf
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,k
      INTEGER :: ii,kk
      INTEGER :: nv,nf_number,len,istart0,istart,istart2,i0,i1,i2
      REAL(8) :: area1,area2
!.....Local vector arrays
      REAL(8) :: flux_conv_nf(nf_non+nf_inl+nf_out)
!
      IF(np.gt.1) CALL communicate_1d(ia_old) 
!
!.....Build summation info for non,inl,out
!
      nf_number_nb=2
      nf_number_id(0)=0
      nf_number_id(1)=2
      nf_number_id(2)=3
      istart_nfs(0)=0
      istart_nfs(1)=istart_nfs(0)+nf_non
      istart_nfs(2)=istart_nfs(1)+nf_inl
      lens         =istart_nfs(2)+nf_out
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
         area1=ia_old(ii)
         area2=ia_old(kk)
         flux_conv_nf(i0)=MIN(flux_g_nf(i1),0.d0)*area2+MAX(flux_g_nf(i1),0.d0)*area1
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
         area1=ia_old(ii)
         area2=ia_b(k)
         flux_conv_nf(i0)=MIN(flux_g_nf(i1),0.d0)*area2+MAX(flux_g_nf(i1),0.d0)*area1
      ENDDO
!               
!.....Cells out
!
      nv=2
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len 
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)
         area1=ia_old(ii)
         flux_conv_nf(i0)=flux_g_nf(i1)*area1
      ENDDO
!
      CALL sum_nf(0,-1,                 &
                  flux_conv_nf,ia_conv)
!
      END SUBROUTINE IAT_convection
