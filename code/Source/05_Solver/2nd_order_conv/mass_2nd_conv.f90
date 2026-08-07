!
      SUBROUTINE mass_2nd_conv
!
!     This routine calculates convective mass fluxes through the cell face
!
      USE Zinterface
      USE VOL_DATA
      USE Zmpi         , ONLY: ncell_fp
      USE Zcore        , ONLY: np
      USE Zparam       , ONLY: ndim
      USE Z2nd_order   , ONLY: mass_conv_2nd,limiter_mass
      USE Znum_cell    , ONLY: istart_nf
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Zvec_major   , ONLY: flux_l_nf,flux_g_nf,liq_conv_nf,vap_conv_nf
      USE Zvec_geo     , ONLY: dxfc_nf,    &
                               dxfc_non_k
!
      IMPLICIT NONE
!.....Local variables
      INTEGER :: i
      INTEGER :: ii,kk
      INTEGER :: nf_number,istart,len,i1
!
      REAL(8) :: liq_2nd_conv,vap_2nd_conv,da1,da2
      REAL(8) :: dx11,dx12,dx13
      REAL(8) :: dx21,dx22,dx23
!.....Local arrays
      REAL(8) :: smin(ncell_fp),smax(ncell_fp)
      REAL(8) :: dadx(ncell_fp,ndim)
!
!.....Gradient for 2nd order interpolation
!
      IF(mass_conv_2nd.eq.2)THEN
!
!........Frink method
!
!
      ELSE
!
!........Green-Gauss (Default)
!
         CALL grad_conv(cell%alphag_o,dadx)
!
      ENDIF
!
!.....Limiter of the gradient
!
      IF(limiter_mass.eq.1.or.limiter_mass.eq.2) THEN
         CALL grad_limiter(cell%alphag_o,dadx,limiter_mass)
         IF(np.gt.1) CALL communicate_2d(dadx)
      ENDIF
!
!.....MIN-MAX values
!
      IF(limiter_mass.eq.3) THEN
         DO i=1,ncell_fp
            smin(i)=cell%alphag_o(i)
            smax(i)=cell%alphag_o(i)
         ENDDO
!
!........Obtain min and max value among neighboring cells
!
         nf_number=0
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len  
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            smax(ii)=DMAX1(smax(ii),cell%alphag_o(kk))
            smin(ii)=DMIN1(smin(ii),cell%alphag_o(kk))
            smax(kk)=DMAX1(smax(kk),cell%alphag_o(ii))
            smin(kk)=DMIN1(smin(kk),cell%alphag_o(ii))
         ENDDO
         IF(np.gt.1) CALL communicate_1d(smin, &
                                         smax)
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
            dx11=dxfc_nf(i1,1)
            dx12=dxfc_nf(i1,2)
            dx21=dxfc_non_k(i,1)
            dx22=dxfc_non_k(i,2)
!
            da1=dadx(ii,1)*dx11+dadx(ii,2)*dx12
            da2=dadx(kk,1)*dx21+dadx(kk,2)*dx22
!
!........MIN-MAX filter
!
            IF(limiter_mass.eq.3)THEN
               da1=DMAX1(da1,smin(ii)-cell%alphag_o(ii))
               da1=DMIN1(da1,smax(ii)-cell%alphag_o(ii))
               da2=DMAX1(da2,smin(kk)-cell%alphag_o(kk))
               da2=DMIN1(da2,smax(kk)-cell%alphag_o(kk))
            ENDIF
!
!........Apply 2nd order upwind
!
            IF    (flux_g_nf(i1).gt.0.d0)THEN
               vap_2nd_conv=cell%rhog_o(ii)*da1
               vap_conv_nf(i1)=vap_conv_nf(i1)+vap_2nd_conv
            ELSEIF(flux_g_nf(i1).lt.0.d0)THEN
               vap_2nd_conv=cell%rhog_o(kk)*da2
               vap_conv_nf(i1)=vap_conv_nf(i1)+vap_2nd_conv
            ENDIF
!
            IF    (flux_l_nf(i1).gt.0.d0)THEN
               liq_2nd_conv=cell%rhol_o(ii)*da1
               liq_conv_nf(i1)=liq_conv_nf(i1)-liq_2nd_conv
            ELSEIF(flux_l_nf(i1).lt.0.d0)THEN
               liq_2nd_conv=cell%rhol_o(kk)*da2
               liq_conv_nf(i1)=liq_conv_nf(i1)-liq_2nd_conv
            ENDIF
         ENDDO
      ELSE
         DO i=1,len  
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
!
            dx11=dxfc_nf(i1,1)
            dx12=dxfc_nf(i1,2)
            dx13=dxfc_nf(i1,3)
            dx21=dxfc_non_k(i,1)
            dx22=dxfc_non_k(i,2)
            dx23=dxfc_non_k(i,3)
!
            da1=dadx(ii,1)*dx11+dadx(ii,2)*dx12+dadx(ii,3)*dx13
            da2=dadx(kk,1)*dx21+dadx(kk,2)*dx22+dadx(kk,3)*dx23
!
!........MIN-MAX filter
!
            IF(limiter_mass.eq.3)THEN
               da1=DMAX1(da1,smin(ii)-cell%alphag_o(ii))
               da1=DMIN1(da1,smax(ii)-cell%alphag_o(ii))
               da2=DMAX1(da2,smin(kk)-cell%alphag_o(kk))
               da2=DMIN1(da2,smax(kk)-cell%alphag_o(kk))
            ENDIF
!
!........Apply 2nd order upwind
!
            IF    (flux_g_nf(i1).gt.0.d0)THEN
               vap_2nd_conv=cell%rhog_o(ii)*da1
               vap_conv_nf(i1)=vap_conv_nf(i1)+vap_2nd_conv
            ELSEIF(flux_g_nf(i1).lt.0.d0)THEN
               vap_2nd_conv=cell%rhog_o(kk)*da2
               vap_conv_nf(i1)=vap_conv_nf(i1)+vap_2nd_conv
            ENDIF
!
            IF    (flux_l_nf(i1).gt.0.d0)THEN
               liq_2nd_conv=cell%rhol_o(ii)*da1
               liq_conv_nf(i1)=liq_conv_nf(i1)-liq_2nd_conv
            ELSEIF(flux_l_nf(i1).lt.0.d0)THEN
               liq_2nd_conv=cell%rhol_o(kk)*da2
               liq_conv_nf(i1)=liq_conv_nf(i1)-liq_2nd_conv
            ENDIF
         ENDDO
      ENDIF
!       
      END SUBROUTINE mass_2nd_conv
