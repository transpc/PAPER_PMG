!
      SUBROUTINE energy_2nd_conv_imp(diag_l_non,diag_g_non, &
                                     src_l,src_g)
!
!     This routine calculates convective energy fluxes through the cell face
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: np
      USE Zparam       , ONLY: ndim
      USE Zvec_param   , ONLY: nf_non
      USE Znum_cell    , ONLY: istart_nf, &
                               nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Z2nd_order   , ONLY: eng_conv_2nd,limiter_eng
!     USE Zcoord1      , ONLY: xloc
!     USE Zdecoupled   , ONLY: al_min_c,ag_min_c
      USE Zvec_major   , ONLY: flux_l_nf,flux_g_nf
      USE Zvec_geo     , ONLY: dxfc_nf,dxfc_non_k
!
      IMPLICIT NONE
!
!.....Input
      REAL(8),DIMENSION(nf_non) :: diag_l_non,diag_g_non
!.....Output
      REAL(8),DIMENSION(ncell_fluid) :: src_l,src_g
!.....Local variables
      INTEGER :: i
      INTEGER :: ii,kk
      INTEGER :: nv,nf_number,istart,len,i1
      REAL(8) :: dag,dal
      REAL(8) :: dx1,dx2,dx3
!.....Local arrays
!     REAL(8) :: smin(ncell_fp),smax(ncell_fp)
      REAL(8),DIMENSION(ncell_fp,ndim) :: dagdx,daldx
!.....Local vector arrays
      REAL(8),DIMENSION(nf_non) :: srcl_non,srcg_non
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
!
      ENDIF
!
!.....Limiter of the gradient
!
      IF(limiter_eng.eq.1.or.limiter_eng.eq.2) THEN
         CALL grad_limiter(cell%eg_o,dagdx,limiter_eng)
         CALL grad_limiter(cell%el_o,daldx,limiter_eng)
         IF(np.gt.1) CALL communicate_2d(dagdx, &
                                         daldx)
      ENDIF
!!
!!.....MIN-MAX values
!!
!      IF(limiter_mass.eq.3) THEN
!         DO i=1,ncell_fp
!            smin(i)=cell%alphag_o(i)
!            smax(i)=cell%alphag_o(i)
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
!            smax(ii)=DMAX1(smax(ii),cell%alphag_o(ii))
!            smin(ii)=DMIN1(smin(ii),cell%alphag_o(ii))
!            smax(kk)=DMAX1(smax(kk),cell%alphag_o(kk))
!            smin(kk)=DMIN1(smin(kk),cell%alphag_o(kk))
!         ENDDO
!         IF(np.gt.1)THEN
!            CALL communicate_1d(smin)
!            CALL communicate_1d(smax)
!         ENDIF
!      ENDIF
!
!.....Computing cells
!
!      nf_number=0
!      istart=istart_nf(1,nf_number)
!      len   =istart_nf(2,nf_number)
!      DO i=1,len  
!         i1=istart+i
!         ii=left_nf(i1)
!         kk=right_non(i)
!!
!         dx1(1)=xfcx_nf(i1)-xloc(ii,1)
!         dx2(1)=xfcx_nf(i1)-xloc(kk,1)
!         dx1(2)=xfcy_nf(i1)-xloc(ii,2)
!         dx2(2)=xfcy_nf(i1)-xloc(kk,2)
!         IF(ndim.eq.3)THEN
!            dx1(3)=xfcz_nf(i1)-xloc(ii,3)
!            dx2(3)=xfcz_nf(i1)-xloc(kk,3)
!         ENDIF
!!
!         da1g=0.0d0
!         da2g=0.0d0
!         da1l=0.0d0
!         da2l=0.0d0
!         DO ix=1,ndim
!            da1g=da1g+dagdx(ix,ii)*dx1(ix)
!            da2g=da2g+dagdx(ix,kk)*dx2(ix)
!            da1l=da1l+daldx(ix,ii)*dx1(ix)
!            da2l=da2l+daldx(ix,kk)*dx2(ix)
!         ENDDO
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
!         IF    (flux_g_nf(i1).gt.0.d0)THEN
!            vap_2nd_conv=diag_g_non(i1)*da1g
!         ELSEIF(flux_g_nf(i1).lt.0.d0)THEN
!            vap_2nd_conv=diag_g_non(i1)*da2g
!         ELSE
!            vap_2nd_conv=0.0d0
!         ENDIF
!!
!         src_g(ii)=src_g(ii)+vap_2nd_conv
!!
!         IF    (flux_l_nf(i1).gt.0.d0)THEN
!            liq_2nd_conv=diag_l_non(i1)*da1l
!         ELSEIF(flux_l_nf(i1).lt.0.d0)THEN
!            liq_2nd_conv=diag_l_non(i1)*da2l
!         ELSE
!            liq_2nd_conv=0.0d0
!         ENDIF
!!
!         src_l(ii)=src_l(ii)+liq_2nd_conv
!!
!      ENDDO
!
!.....Build summation info for non
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
      IF(ndim.eq.2) THEN
         DO i=1,len
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            IF(flux_l_nf(i1).lt.0.d0)THEN
               dx1=dxfc_non_k(i,1)
               dx2=dxfc_non_k(i,2)
               dal=daldx(kk,1)*dx1+daldx(kk,2)*dx2
               srcl_non(i)=diag_l_non(i)*dal
            ELSE
               dx1=dxfc_nf(i1,1)
               dx2=dxfc_nf(i1,2)
               dal=daldx(ii,1)*dx1+daldx(ii,2)*dx2
               srcl_non(i)=diag_l_non(i)*dal
            ENDIF
            IF(flux_g_nf(i1).lt.0)THEN
               dx1=dxfc_non_k(i,1)
               dx2=dxfc_non_k(i,2)
               dag=dagdx(kk,1)*dx1+dagdx(kk,2)*dx2
               srcg_non(i)=diag_g_non(i)*dag
            ELSE
               dx1=dxfc_nf(i1,1)
               dx2=dxfc_nf(i1,2)
               dag=dagdx(ii,1)*dx1+dagdx(ii,2)*dx2
               srcg_non(i)=diag_g_non(i)*dag
            ENDIF
         ENDDO
      ELSE
         DO i=1,len
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            IF(flux_l_nf(i1).lt.0.d0)THEN
               dx1=dxfc_non_k(i,1)
               dx2=dxfc_non_k(i,2)
               dx3=dxfc_non_k(i,3)
               dal=daldx(kk,1)*dx1+daldx(kk,2)*dx2+daldx(kk,3)*dx3
               srcl_non(i)=diag_l_non(i)*dal
            ELSE
               dx1=dxfc_nf(i1,1)
               dx2=dxfc_nf(i1,2)
               dx3=dxfc_nf(i1,3)
               dal=daldx(ii,1)*dx1+daldx(ii,2)*dx2+daldx(ii,3)*dx3
               srcl_non(i)=diag_l_non(i)*dal
            ENDIF
            IF(flux_g_nf(i1).lt.0)THEN
               dx1=dxfc_non_k(i,1)
               dx2=dxfc_non_k(i,2)
               dx3=dxfc_non_k(i,3)
               dag=dagdx(kk,1)*dx1+dagdx(kk,2)*dx2+dagdx(kk,3)*dx3
               srcg_non(i)=diag_g_non(i)*dag
            ELSE
               dx1=dxfc_nf(i1,1)
               dx2=dxfc_nf(i1,2)
               dx3=dxfc_nf(i1,3)
               dag=dagdx(ii,1)*dx1+dagdx(ii,2)*dx2+dagdx(ii,3)*dx3
               srcg_non(i)=diag_g_non(i)*dag
            ENDIF
         ENDDO
      ENDIF
!       
      CALL sum_nf(1,-1,           &
                  srcg_non,src_g, &
                  srcl_non,src_l)
!        
      END SUBROUTINE energy_2nd_conv_imp
