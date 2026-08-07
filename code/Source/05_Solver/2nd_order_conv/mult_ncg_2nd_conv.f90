!
      SUBROUTINE mult_ncg_2nd_conv(qn,temp_non)
!
!     This routine calculates convective mass fluxes through the cell face
!
      USE Zinterface
      USE Zmpi         , ONLY: ncell_fp
      USE Zcore        , ONLY: np
      USE Zparam       , ONLY: ndim
      USE Z2nd_order   , ONLY: ncg_conv_2nd,limiter_qn
      USE Zvec_param   , ONLY: nf_non
      USE Znum_cell    , ONLY: istart_nf
      USE Zvec_major   , ONLY: flux_g_nf
      USE Zvec_geo     , ONLY: dxfc_nf,         &
                               dxfc_non_k
      USE Zvec_index   , ONLY: left_nf,right_non
!
      IMPLICIT NONE
!
!.....Input
      REAL(8) :: qn(ncell_fp)
!.....Output
      REAL(8) :: temp_non(nf_non)
!.....Local variables
      INTEGER :: i
      INTEGER :: ii,kk
      INTEGER :: nf_number,istart,len,i1
      REAL(8) :: dq,dx1,dx2,dx3
!.....Local arrays 
      REAL(8) :: dqdx(ncell_fp,ndim)
!
!.....Gradient for 2nd order interpolation
!
      IF(ncg_conv_2nd.eq.2)THEN
!
!........Frink method
!
!
      ELSE
!
!........Green-Gauss (Default)
!
         CALL grad_conv(qn,dqdx)
!
      ENDIF
!
!.....Limiter of the gradient
!
      IF(limiter_qn.gt.0) THEN
         CALL grad_limiter(qn,dqdx,limiter_qn)
         IF(np.gt.1) CALL communicate_2d(dqdx)
      ENDIF
!
!.....Computing cells
!
      nf_number=0
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      IF(ndim.eq.2)THEN
         DO i=1,len
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
!
            IF(flux_g_nf(i1).lt.0.0d0)THEN
               dx1=dxfc_non_k(i,1)
               dx2=dxfc_non_k(i,2)
               dq=dqdx(kk,1)*dx1+dqdx(kk,2)*dx2
               temp_non(i)=temp_non(i)+dq*flux_g_nf(i1)
            ELSEIF(flux_g_nf(i1).gt.0.0d0)THEN
               dx1=dxfc_nf(i1,1)
               dx2=dxfc_nf(i1,2)
               dq=dqdx(ii,1)*dx1+dqdx(ii,2)*dx2
               temp_non(i)=temp_non(i)+dq*flux_g_nf(i1)
            ENDIF
         ENDDO
      ELSE
         DO i=1,len
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
!
            IF(flux_g_nf(i1).lt.0.0d0)THEN
               dx1=dxfc_non_k(i,1)
               dx2=dxfc_non_k(i,2)
               dx3=dxfc_non_k(i,3)
               dq=dqdx(kk,1)*dx1+dqdx(kk,2)*dx2+dqdx(kk,3)*dx3
               temp_non(i)=temp_non(i)+dq*flux_g_nf(i1)
            ELSEIF(flux_g_nf(i1).gt.0.0d0)THEN
               dx1=dxfc_nf(i1,1)
               dx2=dxfc_nf(i1,2)
               dx3=dxfc_nf(i1,3)
               dq=dqdx(ii,1)*dx1+dqdx(ii,2)*dx2+dqdx(ii,3)*dx3
               temp_non(i)=temp_non(i)+dq*flux_g_nf(i1)
            ENDIF
         ENDDO
      ENDIF
!       
      END SUBROUTINE mult_ncg_2nd_conv
!
      SUBROUTINE mult_ncg_2nd_conv1(qn,temp_non)
!
!     This routine calculates convective mass fluxes through the cell face
!
      USE Zinterface
      USE Zmpi         , ONLY: ncell_fp
      USE Zcore        , ONLY: np
      USE Zparam       , ONLY: ndim
      USE Z2nd_order   , ONLY: ncg_conv_2nd,limiter_qn
      USE Zvec_param   , ONLY: nf_non
      USE Znum_cell    , ONLY: istart_nf
      USE Zvec_major   , ONLY: flux_g_nf
      USE Zvec_geo     , ONLY: dxfc_nf,        &
                               dxfc_non_k
      USE Zvec_index   , ONLY: left_nf,right_non
!
      IMPLICIT NONE
!
!     input
      REAL(8) :: qn(ncell_fp)
!     output
      REAL(8) :: temp_non(nf_non)
!     local variables
      INTEGER :: i
      INTEGER :: ii,kk
      INTEGER :: nf_number,istart,len,i1
      REAL(8) :: dq,dx1,dx2,dx3
!     local arrays 
      REAL(8) :: dqdx(ncell_fp,ndim)
!
!.....Gradient for 2nd order interpolation
!
      IF(ncg_conv_2nd.eq.2)THEN
!
!........Frink method
!
!
      ELSE
!
!........Green-Gauss (Default)
!
         CALL grad_conv(qn,dqdx)
!
      ENDIF
!
!.....Limiter of the gradient
!
      IF(limiter_qn.gt.0) THEN
         CALL grad_limiter(qn,dqdx,limiter_qn)
         IF(np.gt.1) CALL communicate_2d(dqdx)
      ENDIF
!
!.....Computing cells
!
      nf_number=0
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      IF(ndim.eq.2)THEN
         DO i=1,len
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
!
            IF(flux_g_nf(i1).lt.0.0d0)THEN
               dx1=dxfc_non_k(i,1)
               dx2=dxfc_non_k(i,2)
               dq=dqdx(kk,1)*dx1+dqdx(kk,2)*dx2
               temp_non(i)=temp_non(i)+dq
            ELSEIF(flux_g_nf(i1).gt.0.0d0)THEN
               dx1=dxfc_nf(i1,1)
               dx2=dxfc_nf(i1,2)
               dq=dqdx(ii,1)*dx1+dqdx(ii,2)*dx2
               temp_non(i)=temp_non(i)+dq
            ENDIF
         ENDDO
      ELSE
         DO i=1,len
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
!
            IF(flux_g_nf(i1).lt.0.0d0)THEN
               dx1=dxfc_non_k(i,1)
               dx2=dxfc_non_k(i,2)
               dx3=dxfc_non_k(i,3)
               dq=dqdx(kk,1)*dx1+dqdx(kk,2)*dx2+dqdx(kk,3)*dx3
               temp_non(i)=temp_non(i)+dq
            ELSEIF(flux_g_nf(i1).gt.0.0d0)THEN
               dx1=dxfc_nf(i1,1)
               dx2=dxfc_nf(i1,2)
               dx3=dxfc_nf(i1,3)
               dq=dqdx(ii,1)*dx1+dqdx(ii,2)*dx2+dqdx(ii,3)*dx3
               temp_non(i)=temp_non(i)+dq
            ENDIF
         ENDDO
      ENDIF
!       
      END SUBROUTINE mult_ncg_2nd_conv1
