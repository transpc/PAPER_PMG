!
      SUBROUTINE mult_ncg_2nd_conv_imp(qn,temp_non,src)
!
!     This routine calculates convective mass fluxes through the cell face
!
      USE Zinterface
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: np
      USE Zparam       , ONLY: ndim
      USE Zvec_param   , ONLY: nf_non
      USE Znum_cell    , ONLY: istart_nf, &
                               nf_number_nb,lens,nf_number_id,istart_nfs
      USE Z2nd_order   , ONLY: ncg_conv_2nd,limiter_qn
      USE Zvec_geo     , ONLY: dxfc_nf,    &
                               dxfc_non_k
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Zvec_major   , ONLY: flux_g_nf
!
      IMPLICIT NONE
!
!.....Input
      REAL(8),DIMENSION(ncell_fp) :: qn
      REAL(8),DIMENSION(nf_non) :: temp_non
!.....Output
      REAL(8),DIMENSION(ncell_fluid) :: src
!.....Local variable
      INTEGER :: i
      INTEGER :: ii,kk
      INTEGER :: nv,nf_number,istart,len,i1
!
      REAL(8) :: dq
      REAL(8) :: dx1,dx2,dx3
!.....Local arrays
      REAL(8),DIMENSION(ncell_fp,ndim) :: dqdx
!.....Local vector arrays
      REAL(8),DIMENSION(nf_non) :: dq_non
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
!.....Build summation info for non
!
      nf_number_nb=0
      nf_number_id(0)=0
      istart_nfs(0)=0
      lens=istart_nfs(1)+nf_non
!
!.....Computing cells
!
      nv=0
      nf_number=nf_number_id(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      IF(ndim.eq.2) THEN
         DO i=1,len
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            IF(flux_g_nf(i1).lt.0.d0)THEN
               dx1=dxfc_non_k(i,1)
               dx2=dxfc_non_k(i,2)
               dq=dqdx(kk,1)*dx1+dqdx(kk,2)*dx2
            ELSE
               dx1=dxfc_nf(i1,1)
               dx2=dxfc_nf(i1,2)
               dq=dqdx(ii,1)*dx1+dqdx(ii,2)*dx2
            ENDIF
            dq_non(i)=-temp_non(i)*dq 
         ENDDO
      ELSE
         DO i=1,len
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            IF(flux_g_nf(i1).lt.0.d0)THEN
               dx1=dxfc_non_k(i,1)
               dx2=dxfc_non_k(i,2)
               dx3=dxfc_non_k(i,3)
               dq=dqdx(kk,1)*dx1+dqdx(kk,2)*dx2+dqdx(kk,3)*dx3
            ELSE
               dx1=dxfc_nf(i1,1)
               dx2=dxfc_nf(i1,2)
               dx3=dxfc_nf(i1,3)
               dq=dqdx(ii,1)*dx1+dqdx(ii,2)*dx2+dqdx(ii,3)*dx3
            ENDIF
            dq_non(i)=-temp_non(i)*dq 
         ENDDO
      ENDIF
!
      CALL sum_nf(1,-1,       &
                  dq_non,src)
!       
      END SUBROUTINE mult_ncg_2nd_conv_imp
