!
      SUBROUTINE boron_2nd_conv_imp(src_cb)
!
!     This routine calculates convective mass fluxes through the cell face
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: np
      USE Zparam       , ONLY: ndim
      USE Zvec_param   , ONLY: nf_non
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Znum_cell    , ONLY: istart_nf, &
                               nf_number_nb,lens,nf_number_id,istart_nfs
      USE Z2nd_order   , ONLY: boron_conv_2nd,limiter_boron
      USE Zcoord3      , ONLY: volr
      USE Zvec_major   , ONLY: flux_l_nf
      USE Zare         , ONLY: ar_liq
      USE Zvec_geo     , ONLY: dxfc_nf,dxfc_non_k
!
      IMPLICIT NONE
!
!.....Output
      REAL(8) :: src_cb(ncell_fluid)
!.....Local variables
      INTEGER :: ii,kk
      INTEGER :: i
      INTEGER :: nv,nf_number,istart,len,i1
      REAL(8) :: dx11,dx12,dx13
      REAL(8) :: dx21,dx22,dx23
      REAL(8) :: dq,ar_l
!.....Local arrays
      REAL(8) :: src_cb1(ncell_fluid)
      REAL(8) :: q(ncell_fp)
      REAL(8) :: dqdx(ncell_fp,ndim)
!.....Local vector arrays
      REAL(8) :: src_xs_non(nf_non)
      
      q(:)=cell%cboron(:)
!
!.....Gradient for 2nd order interpolation
!
      IF(boron_conv_2nd.eq.2)THEN
!
!........Frink method
!
!
      ELSE
!
!........Green-Gauss (Default)
!
         CALL grad_conv(q,dqdx)
!
      ENDIF
!
!.....Limiter of the gradient
!
      IF(limiter_boron.gt.0) THEN
         CALL grad_limiter(q,dqdx,limiter_boron)
!
!.....MPI communication
!
         IF(np.gt.1) CALL communicate_2d(dqdx)
      ENDIF
!
!.....Build summation info for non,mcc,inl,out
!
      nf_number_nb=0
      nf_number_id(0)=0
      istart_nfs(0)=0
      lens         =istart_nfs(0)+nf_non
!
!.....Cells non
!
      nv=0
      nf_number=nf_number_id(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      IF(ndim.eq.2)THEN
         DO i=1,len
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            dx11=dxfc_nf(i1,1)
            dx12=dxfc_nf(i1,2)
            dx21=dxfc_non_k(i,1)
            dx22=dxfc_non_k(i,2)
            IF(flux_l_nf(i).lt.0)THEN
               ar_l=ar_liq(kk)
               dq=dqdx(kk,1)*dx21+dqdx(kk,2)*dx22
               src_xs_non(i)=-ar_l*flux_l_nf(i)*dq
            ELSE
               ar_l=ar_liq(ii)
               dq=dqdx(ii,1)*dx11+dqdx(ii,2)*dx12
               src_xs_non(i)=-ar_l*flux_l_nf(i)*dq
            ENDIF
         ENDDO
      ELSE
         DO i=1,len
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            dx11=dxfc_nf(i1,1)
            dx12=dxfc_nf(i1,2)
            dx13=dxfc_nf(i1,3)
            dx21=dxfc_non_k(i,1)
            dx22=dxfc_non_k(i,2)
            dx23=dxfc_non_k(i,3)
            IF(flux_l_nf(i).lt.0)THEN
               ar_l=ar_liq(kk)
               dq=dqdx(kk,1)*dx21+dqdx(kk,2)*dx22+dqdx(kk,3)*dx23
               src_xs_non(i)=-ar_l*flux_l_nf(i)*dq
            ELSE
               ar_l=ar_liq(ii)
               dq=dqdx(ii,1)*dx11+dqdx(ii,2)*dx12+dqdx(ii,3)*dx13
               src_xs_non(i)=-ar_l*flux_l_nf(i)*dq
            ENDIF
         ENDDO
      ENDIF
!
      DO i=1,ncell_fluid
!        src_cb1(i)=0.d0
      ENDDO
      CALL sum_nf(0,-1,               &
                  src_xs_non,src_cb1)
!
      DO i=1,ncell_fluid
         src_cb(i)=src_cb(i)+src_cb1(i)*volr(i)
      ENDDO
!       
      END SUBROUTINE boron_2nd_conv_imp
