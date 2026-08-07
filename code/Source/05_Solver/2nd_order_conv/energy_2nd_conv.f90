!
      SUBROUTINE energy_2nd_conv
!
!     This routine calculates convective energy fluxes through the cell face
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell
      USE Zmpi         , ONLY: ncell_fp
      USE Zcore        , ONLY: np
      USE Zparam       , ONLY: ndim
      USE Z2nd_order   , ONLY: eng_conv_2nd,limiter_eng
      USE Zare         , ONLY: ar_liq,ar_drp,ar_gas
      USE Znum_cell    , ONLY: istart_nf
      USE Zvec_geo     , ONLY: dxfc_nf,         &
                               dxfc_non_k
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Zvec_major   , ONLY: flux_l_nf,flux_g_nf,flux_d_nf,ecnvc_l_nf,ecnvc_g_nf,ecnvc_d_nf
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i
      INTEGER :: ii,kk
      INTEGER :: nf_number,istart,len,i1
      REAL(8) :: vap_2nd_conv,da1g,da2g
      REAL(8) :: liq_2nd_conv,da1l,da2l
      REAL(8) :: drp_2nd_conv,da1d,da2d
      REAL(8) :: dx11,dx12,dx13
      REAL(8) :: dx21,dx22,dx23
!.....Local arrays
!     REAL(8) :: smin(ncell_fp),smax(ncell_fp)
      REAL(8) :: dagdx(ncell_fp,ndim),daldx(ncell_fp,ndim),daddx(ncell_fp,ndim)
!
!.....Gradient for 2nd order interpolation
!
      IF(eng_conv_2nd.eq.2)THEN
!
!........Frink method
!
!
      ELSE
!
!........Green-Gauss (Default)
!
         CALL grad_conv(cell%eg_o,dagdx)
         CALL grad_conv(cell%el_o,daldx)
         cell%ed_o=cell%el_o
         CALL grad_conv(cell%ed_o,daddx)
!
      ENDIF
!
!.....Limiter of the gradient
!
      IF(limiter_eng.eq.1.or.limiter_eng.eq.2) THEN
         CALL grad_limiter(cell%eg_o,dagdx,limiter_eng)
         CALL grad_limiter(cell%el_o,daldx,limiter_eng)
         CALL grad_limiter(cell%ed_o,daddx,limiter_eng)
         IF(np.gt.1) CALL communicate_2d(dagdx, &
                                         daldx, &
                                         daddx)
      ENDIF
!
!.....MIN-MAX values
!
!      IF(limiter_eng.eq.3) THEN
!         DO i=1,ncell_fp
!            smin(i)=cell%eg_o(i)
!            smax(i)=cell%eg_o(i)
!         ENDDO
!!
!!........Obtain min and max value among neighboring cells
!!
!         nf_number=0
!         istart=istart_nf(1,nf_number)
!         len   =istart_nf(2,nf_number)
!         DO i=1,len  
!            ii=left_nf(i)
!            kk=right_non(i)
!            smax(ii)=DMAX1(smax(ii),cell%eg_o(ii))
!            smin(ii)=DMIN1(smin(ii),cell%eg_o(ii))
!            smax(kk)=DMAX1(smax(kk),cell%eg_o(kk))
!            smin(kk)=DMIN1(smin(kk),cell%eg_o(kk))
!         ENDDO
!         IF(np.gt.1)THEN
!            CALL communicate_1d(smin)
!            CALL communicate_1d(smax)
!         ENDIF
!      ENDIF
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
            da1g=dagdx(ii,1)*dx11+dagdx(ii,2)*dx12
            da1l=daldx(ii,1)*dx11+daldx(ii,2)*dx12
            da1d=daddx(ii,1)*dx11+daddx(ii,2)*dx12
            da2g=dagdx(kk,1)*dx21+dagdx(kk,2)*dx22
            da2l=daldx(kk,1)*dx21+daldx(kk,2)*dx22
            da2d=daddx(kk,1)*dx21+daddx(kk,2)*dx22
!
!........MIN-MAX filter
!
!         IF(limiter_eng.eq.3)THEN
!            da1=DMAX1(da1,smin(ii)-cell%eg_o(ii))
!            da1=DMIN1(da1,smax(ii)-cell%eg_o(ii))
!            da2=DMAX1(da2,smin(kk)-cell%eg_o(kk))
!            da2=DMIN1(da2,smax(kk)-cell%eg_o(kk))
!         ENDIF
!
!........Apply 2nd order upwind
!
            IF    (flux_g_nf(i1).gt.0.d0)THEN
               vap_2nd_conv=ar_gas(ii)*da1g
               ecnvc_g_nf(i1)=ecnvc_g_nf(i1)+vap_2nd_conv
            ELSEIF(flux_g_nf(i1).lt.0.d0)THEN
               vap_2nd_conv=ar_gas(kk)*da2g
               ecnvc_g_nf(i1)=ecnvc_g_nf(i1)+vap_2nd_conv
            ENDIF
!
            IF    (flux_l_nf(i1).gt.0.d0)THEN
               liq_2nd_conv=ar_liq(ii)*da1l
               ecnvc_l_nf(i1)=ecnvc_l_nf(i1)+liq_2nd_conv
            ELSEIF(flux_l_nf(i1).lt.0.d0)THEN
               liq_2nd_conv=ar_liq(kk)*da2l
               ecnvc_l_nf(i1)=ecnvc_l_nf(i1)+liq_2nd_conv
            ENDIF
!
            IF    (flux_d_nf(i1).gt.0.d0)THEN
               drp_2nd_conv=ar_drp(ii)*da1d
               ecnvc_d_nf(i1)=ecnvc_d_nf(i1)+drp_2nd_conv
            ELSEIF(flux_d_nf(i1).lt.0.d0)THEN
               drp_2nd_conv=ar_drp(kk)*da2d
               ecnvc_d_nf(i1)=ecnvc_d_nf(i1)+drp_2nd_conv
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
            da1g=dagdx(ii,1)*dx11+dagdx(ii,2)*dx12+dagdx(ii,3)*dx13
            da1l=daldx(ii,1)*dx11+daldx(ii,2)*dx12+daldx(ii,3)*dx13
            da1d=daddx(ii,1)*dx11+daddx(ii,2)*dx12+daddx(ii,3)*dx13
            da2g=dagdx(kk,1)*dx21+dagdx(kk,2)*dx22+dagdx(kk,3)*dx23
            da2l=daldx(kk,1)*dx21+daldx(kk,2)*dx22+daldx(kk,3)*dx23
            da2d=daddx(kk,1)*dx21+daddx(kk,2)*dx22+daddx(kk,3)*dx23
!
!........MIN-MAX filter
!
!         IF(limiter_eng.eq.3)THEN
!            da1=DMAX1(da1,smin(ii)-cell%eg_o(ii))
!            da1=DMIN1(da1,smax(ii)-cell%eg_o(ii))
!            da2=DMAX1(da2,smin(kk)-cell%eg_o(kk))
!            da2=DMIN1(da2,smax(kk)-cell%eg_o(kk))
!         ENDIF
!
!........Apply 2nd order upwind
!
            IF    (flux_g_nf(i1).gt.0.d0)THEN
               vap_2nd_conv=ar_gas(ii)*da1g
               ecnvc_g_nf(i1)=ecnvc_g_nf(i1)+vap_2nd_conv
            ELSEIF(flux_g_nf(i1).lt.0.d0)THEN
               vap_2nd_conv=ar_gas(kk)*da2g
               ecnvc_g_nf(i1)=ecnvc_g_nf(i1)+vap_2nd_conv
            ENDIF
!
            IF    (flux_l_nf(i1).gt.0.d0)THEN
               liq_2nd_conv=ar_liq(ii)*da1l
               ecnvc_l_nf(i1)=ecnvc_l_nf(i1)+liq_2nd_conv
            ELSEIF(flux_l_nf(i1).lt.0.d0)THEN
               liq_2nd_conv=ar_liq(kk)*da2l
               ecnvc_l_nf(i1)=ecnvc_l_nf(i1)+liq_2nd_conv
            ENDIF
!
            IF    (flux_d_nf(i1).gt.0.d0)THEN
               drp_2nd_conv=ar_drp(ii)*da1d
               ecnvc_d_nf(i1)=ecnvc_d_nf(i1)+drp_2nd_conv
            ELSEIF(flux_d_nf(i1).lt.0.d0)THEN
               drp_2nd_conv=ar_drp(kk)*da2d
               ecnvc_d_nf(i1)=ecnvc_d_nf(i1)+drp_2nd_conv
            ENDIF
         ENDDO
      ENDIF
!       
      END SUBROUTINE energy_2nd_conv
