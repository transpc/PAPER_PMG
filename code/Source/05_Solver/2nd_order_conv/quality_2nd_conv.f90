!
      SUBROUTINE quality_2nd_conv
!
!     This routine calculates convective mass fluxes through the cell face
!
      USE Zinterface
      USE VOL_DATA      , ONLY: cell
      USE Zmpi         , ONLY: ncell_fp
      USE Zcore        , ONLY: np
      USE Zparam       , ONLY: ndim
      USE Z2nd_order   , ONLY: qula_conv_2nd,limiter_xn
      USE Zare         , ONLY: ar_gas
      USE Znum_cell    , ONLY: istart_nf
      USE Ztimecon     , ONLY: alpha_min
      USE Zvec_geo     , ONLY: dxfc_nf,         &
                               dxfc_non_k
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Zvec_major   , ONLY: flux_g_nf,quala_conv_nf
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i
      INTEGER :: ii,kk
      INTEGER :: nf_number,istart,len,i1
      REAL(8) :: quala_2nd_conv,ar1_qul,ar2_qul,dq1,dq2
      REAL(8) :: dx1,dx2,dx3
!.....Local arrays
      REAL(8) :: dqdx(ncell_fp,ndim)
!
!.....Gradient for 2nd order interpolation
!
      IF(qula_conv_2nd.eq.2)THEN
!
!........Frink method
!
!
      ELSE
!
!........Green-Gauss (Default)
!
         CALL grad_conv(cell%quala,dqdx)
!
      ENDIF
!
!.....Limiter of the gradient
!
      IF(limiter_xn.gt.0) THEN
         CALL grad_limiter(cell%quala,dqdx,limiter_xn)
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
!
!........Define owner cell values
!........Assume zero cell values when the phase fraction is less than alpha_min
!
            IF(cell%alphag_o(ii).le.alpha_min) THEN
               ar1_qul=0.d0
            ELSE
               dx1=dxfc_nf(i1,1)
               dx2=dxfc_nf(i1,2)
               dq1=dqdx(ii,1)*dx1+dqdx(ii,2)*dx2
               ar1_qul=ar_gas(ii)*dq1
            ENDIF 
!
!........Define neighbor cell values
!........Assume zero cell values when the phase fraction is less than alpha_min
!
            IF(cell%alphag_o(kk).le.alpha_min) THEN
               ar2_qul=0.d0 
            ELSE
               dx1=dxfc_non_k(i,1)
               dx2=dxfc_non_k(i,2)
               dq2=dqdx(kk,1)*dx1+dqdx(kk,2)*dx2
               ar2_qul=ar_gas(kk)*dq2
            ENDIF 
!
!........Apply 2nd order upwind
!
            IF    (flux_g_nf(i1).gt.0.d0)THEN
               quala_2nd_conv=ar1_qul
               quala_conv_nf(i1)=quala_conv_nf(i1)+quala_2nd_conv
            ELSEIF(flux_g_nf(i1).lt.0.d0)THEN
               quala_2nd_conv=ar2_qul
               quala_conv_nf(i1)=quala_conv_nf(i1)+quala_2nd_conv
            ENDIF
         ENDDO
      ELSE
         DO i=1,len  
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
!
!
!........Define owner cell values
!........Assume zero cell values when the phase fraction is less than alpha_min
!
            IF(cell%alphag_o(ii).le.alpha_min) THEN
               ar1_qul=0.d0
            ELSE
               dx1=dxfc_nf(i1,1)
               dx2=dxfc_nf(i1,2)
               dx3=dxfc_nf(i1,3)
               dq1=dqdx(ii,1)*dx1+dqdx(ii,2)*dx2+dqdx(ii,3)*dx3
               ar1_qul=ar_gas(ii)*dq1
            ENDIF 
!
!........Define neighbor cell values
!........Assume zero cell values when the phase fraction is less than alpha_min
!
            IF(cell%alphag_o(kk).le.alpha_min) THEN
               ar2_qul=0.d0 
            ELSE
               dx1=dxfc_non_k(i,1)
               dx2=dxfc_non_k(i,2)
               dx3=dxfc_non_k(i,3)
               dq2=dqdx(kk,1)*dx1+dqdx(kk,2)*dx2+dqdx(kk,3)*dx3
               ar2_qul=ar_gas(kk)*dq2
            ENDIF 
!
!........Apply 2nd order upwind
!
            IF    (flux_g_nf(i1).gt.0.d0)THEN
               quala_2nd_conv=ar1_qul
            ELSEIF(flux_g_nf(i1).lt.0.d0)THEN
               quala_2nd_conv=ar2_qul
            ELSE
               quala_2nd_conv=(ar1_qul+ar2_qul)*0.5d0
            ENDIF
            quala_conv_nf(i1)=quala_conv_nf(i1)+quala_2nd_conv
         ENDDO
      ENDIF
!       
      END SUBROUTINE quality_2nd_conv
