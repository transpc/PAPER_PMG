!
      SUBROUTINE boron_2nd_conv
!
!     This routine calculates convective mass fluxes through the cell face
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell
      USE Zmpi         , ONLY: ncell_fp
      USE Zcore        , ONLY: np
      USE Zparam       , ONLY: ndim
      USE Z2nd_order   , ONLY: limiter_boron,boron_conv_2nd
      USE Zare         , ONLY: ar_liq
      USE Znum_cell    , ONLY: istart_nf
      USE Zvec_geo     , ONLY: dxfc_nf,         &
                               dxfc_non_k
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Zvec_major   , ONLY: flux_l_nf,lbor_conv_nf
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i
      INTEGER :: ii,kk
      INTEGER :: nf_number,istart,len,i1
      REAL(8) :: ar1_liq,ar2_liq
      REAL(8) :: dq1,dq2
      REAL(8) :: dx1,dx2,dx3
!.....Local arrays
      REAL(8) :: dqdx(ncell_fp,ndim)
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
         IF(np.gt.1) CALL communicate_1d(cell%cboron)
         CALL grad_conv(cell%cboron,dqdx)
!
      ENDIF
!
!.....Limiter of the gradient
!
      IF(limiter_boron.gt.0) THEN
         CALL grad_limiter(cell%cboron,dqdx,limiter_boron)
         IF(np.gt.1) CALL communicate_2d(dqdx)
      ENDIF
!
!.....Computing cells
!
!........Apply 2nd order upwind
!         
      nf_number=0
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      IF(ndim.eq.2)THEN
         DO i=1,len  
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            IF    (flux_l_nf(i1).gt.0.d0)THEN
               dx1=dxfc_nf(i1,1)
               dx2=dxfc_nf(i1,2)
               dq1=dqdx(ii,1)*dx1+dqdx(ii,2)*dx2
               ar1_liq=ar_liq(ii)
               lbor_conv_nf(i1)=lbor_conv_nf(i1)+ar1_liq*dq1
            ELSEIF(flux_l_nf(i1).lt.0.d0)THEN
               dx1=dxfc_non_k(i,1)
               dx2=dxfc_non_k(i,2)
               dq2=dqdx(kk,1)*dx1+dqdx(kk,2)*dx2
               ar2_liq=ar_liq(kk)
               lbor_conv_nf(i1)=lbor_conv_nf(i1)+ar2_liq*dq2
            ENDIF
         ENDDO
      ELSE
         DO i=1,len  
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            IF    (flux_l_nf(i1).gt.0.d0)THEN
               dx1=dxfc_nf(i1,1)
               dx2=dxfc_nf(i1,2)
               dx3=dxfc_nf(i1,3)
               dq1=dqdx(ii,1)*dx1+dqdx(ii,2)*dx2+dqdx(ii,3)*dx3
               ar1_liq=ar_liq(ii)
               lbor_conv_nf(i1)=lbor_conv_nf(i1)+ar1_liq*dq1
            ELSEIF(flux_l_nf(i1).lt.0.d0)THEN
               dx1=dxfc_non_k(i,1)
               dx2=dxfc_non_k(i,2)
               dx3=dxfc_non_k(i,3)
               dq2=dqdx(kk,1)*dx1+dqdx(kk,2)*dx2+dqdx(kk,3)*dx3
               ar2_liq=ar_liq(kk)
               lbor_conv_nf(i1)=lbor_conv_nf(i1)+ar2_liq*dq2
            ENDIF
         ENDDO
      ENDIF
!       
      END SUBROUTINE boron_2nd_conv
