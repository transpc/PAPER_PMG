      SUBROUTINE grad_conv(s,dsdx)
!
!     This routine calculates the components of the gradient based on the gauss theorem.
!
      USE Zinterface
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: np
      USE Zparam       , ONLY: ndim
      USE Zvec_param   , ONLY: nf_non,nf_mcc,nf_inl,nf_out,nf_adw,nf_fsw,nf_ctw,nf_chw,nf_sym,nf_tot
      USE Znum_cell    , ONLY: istart_nf, &
                               nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zcoord3      , ONLY: volr
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Zvec_geo     , ONLY: sv_nf,fac1_non,fac_non
!
      IMPLICIT NONE
!
!.....Input
      REAL(8),DIMENSION(ncell_fp) :: s
!.....Output
      REAL(8),DIMENSION(ncell_fp,ndim) :: dsdx
!.....Local variables
      INTEGER :: i,ix
      INTEGER :: ii,kk
      INTEGER :: nv,nf_number,istart,len,i1 
!.....Local vector arrays
      REAL(8),DIMENSION(nf_tot,ndim) :: fie_nf
!
!.....Build summation info for non,inl
!
      nf_number_nb=8
      nf_number_id(0)=0
      nf_number_id(1)=1
      nf_number_id(2)=2
      nf_number_id(3)=3
      nf_number_id(4)=4
      nf_number_id(5)=5
      nf_number_id(6)=6
      nf_number_id(7)=7
      nf_number_id(8)=8
      istart_nfs(0)=0
      istart_nfs(1)=istart_nfs(0)+nf_non
      istart_nfs(2)=istart_nfs(1)+nf_mcc
      istart_nfs(3)=istart_nfs(2)+nf_inl
      istart_nfs(4)=istart_nfs(3)+nf_out
      istart_nfs(5)=istart_nfs(4)+nf_adw
      istart_nfs(6)=istart_nfs(5)+nf_fsw
      istart_nfs(7)=istart_nfs(6)+nf_ctw
      istart_nfs(8)=istart_nfs(7)+nf_chw
      lens         =istart_nfs(8)+nf_sym
!
!.....Computing cell
!
      nv=0
      nf_number=nf_number_id(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO ix=1,ndim
         DO i=1,len  
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            fie_nf(i1,ix)=(fac1_non(i)*s(ii)+fac_non(i)*s(kk))*sv_nf(i1,ix)
         ENDDO
      ENDDO
!
!.....The rest
!
      DO nv=1,8
         nf_number=nf_number_id(nv)
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO ix=1,ndim
            DO i=1,len  
               i1=istart+i
               ii=left_nf(i1)
               fie_nf(i1,ix)=s(ii)*sv_nf(i1,ix)
            ENDDO
         ENDDO
      ENDDO
!
      CALL sum_nf_ndim(0,-1,ncell_fp, &
                       fie_nf,dsdx)
!
      DO ix=1,ndim
         DO i=1,ncell_fluid
            dsdx(i,ix)=dsdx(i,ix)*volr(i)
         ENDDO
      ENDDO
!
      IF(np.gt.1) CALL communicate_2d(dsdx)
!
      END SUBROUTINE grad_conv
