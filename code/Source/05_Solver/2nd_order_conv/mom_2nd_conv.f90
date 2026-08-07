!
      SUBROUTINE mom_2nd_conv(fluxl_c_nf,fluxg_c_nf,fluxd_c_nf)
!
!     This routine calculates convective mass fluxes through the cell face
!
      USE Zinterface
      USE VOL_DATA      , ONLY: cell
      USE Zmpi         , ONLY: ncell_fp
      USE Zparam       , ONLY: ndim
      USE Zcore        , ONLY: np
      USE Zvec_param   , ONLY: nf_flux
      USE Znum_cell    , ONLY: istart_nf
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Z2nd_order   , ONLY: mom_conv_2nd,limiter_mom
      USE Zare         , ONLY: ar_gas,ar_liq,ar_drp
      USE Zb_condition , ONLY: vb_liq,vin_liq,vb_gas,vin_gas,vb_drp,vin_drp
      USE Ztimecon     , ONLY: alpha_min
      USE Zvector      , ONLY: vl_o,vg_o,vd_o
      USE Zvec_geo     , ONLY: dxfc_nf,      &
                               dxfc_non_k
      USE Zvec_major   , ONLY: flux_g_nf,flux_l_nf,flux_d_nf
!
      IMPLICIT NONE
!
!.....Output
      REAL(8) :: fluxl_c_nf(nf_flux,ndim),fluxg_c_nf(nf_flux,ndim),fluxd_c_nf(nf_flux,ndim)
!.....Local variables
      INTEGER :: i
      INTEGER :: ii,kk
      INTEGER :: nf_number,istart,len,i1
!
      REAL(8) :: dvl11,dvl12,dvl13,dvg11,dvg12,dvg13,dvd11,dvd12,dvd13
      REAL(8) :: dvl21,dvl22,dvl23,dvg21,dvg22,dvg23,dvd21,dvd22,dvd23
      REAL(8) :: dx11,dx12,dx13
      REAL(8) :: dx21,dx22,dx23
      REAL(8) :: a_l,a_g,a_d
      REAL(8) :: b_l,b_g,b_d
!.....Local arrays
      REAL(8) :: dvldx(ncell_fp,ndim,ndim),dvgdx(ncell_fp,ndim,ndim),dvddx(ncell_fp,ndim,ndim)
!
!.....Gradient for 2nd order interpolation
!
      IF(mom_conv_2nd.eq.2)THEN
!
!........Frink method
!
!
      ELSE
!
!........Green-Gauss (Default)
!
         CALL grad_vel(2,vl_o,dvldx,vb_liq,vin_liq)
         CALL grad_vel(1,vg_o,dvgdx,vb_gas,vin_gas)
         CALL grad_vel(2,vd_o,dvddx,vb_drp,vin_drp)
!
      ENDIF
      IF(np.gt.1) CALL communicate_3d(dvldx, &
                                      dvgdx, &
                                      dvddx)
!
!.....Limiter of the gradient
!
      IF(limiter_mom.gt.0)THEN
         CALL grad_limiter(vl_o(1,1),dvldx(1,1,1),limiter_mom)
         CALL grad_limiter(vl_o(1,2),dvldx(1,1,2),limiter_mom)
         CALL grad_limiter(vg_o(1,1),dvgdx(1,1,1),limiter_mom)
         CALL grad_limiter(vg_o(1,2),dvgdx(1,1,2),limiter_mom)
         IF(np.gt.1) CALL communicate_3d(dvldx, &
                                         dvgdx, &
                                         dvddx)
      ENDIF
!
      IF(ndim.eq.2) THEN
!
!........Computing cells
!
         nf_number=0
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
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
!...........Assume zero cell values when the phase fraction is less than alpha_min
!
            IF(cell%alphal_o(ii).le.alpha_min) THEN
               dvl11=0.d0
               dvl12=0.d0
            ELSE
               dvl11=dvldx(ii,1,1)*dx11+dvldx(ii,2,1)*dx12
               dvl12=dvldx(ii,1,2)*dx11+dvldx(ii,2,2)*dx12
            ENDIF
            IF(cell%alphag_o(ii).le.alpha_min) THEN
               dvg11=0.d0
               dvg12=0.d0
            ELSE
               dvg11=dvgdx(ii,1,1)*dx11+dvgdx(ii,2,1)*dx12
               dvg12=dvgdx(ii,1,2)*dx11+dvgdx(ii,2,2)*dx12
            ENDIF
            IF(cell%alphad_o(ii).le.alpha_min) THEN
               dvd11=0.d0
               dvd12=0.d0
            ELSE
               dvd11=dvddx(ii,1,1)*dx11+dvddx(ii,2,1)*dx12
               dvd12=dvddx(ii,1,2)*dx11+dvddx(ii,2,2)*dx12
            ENDIF
!
            IF(cell%alphal_o(kk).le.alpha_min) THEN
               dvl21=0.d0
               dvl22=0.d0
            ELSE
               dvl21=dvldx(kk,1,1)*dx21+dvldx(kk,2,1)*dx22
               dvl22=dvldx(kk,1,2)*dx21+dvldx(kk,2,2)*dx22
            ENDIF
            IF(cell%alphag_o(kk).le.alpha_min) THEN
               dvg21=0.d0
               dvg22=0.d0
            ELSE
               dvg21=dvgdx(kk,1,1)*dx21+dvgdx(kk,2,1)*dx22
               dvg22=dvgdx(kk,1,2)*dx21+dvgdx(kk,2,2)*dx22
            ENDIF
            IF(cell%alphad_o(kk).le.alpha_min) THEN
               dvd21=0.d0
               dvd22=0.d0
            ELSE
               dvd21=dvddx(kk,1,1)*dx21+dvddx(kk,2,1)*dx22
               dvd22=dvddx(kk,1,2)*dx21+dvddx(kk,2,2)*dx22
            ENDIF
!
!........Apply 2nd order upwind
!
            a_l=min(flux_l_nf(i1),0.d0)
            b_l=max(flux_l_nf(i1),0.d0)
            a_g=min(flux_g_nf(i1),0.d0)
            b_g=max(flux_g_nf(i1),0.d0)
            a_d=min(flux_d_nf(i1),0.d0)
            b_d=max(flux_d_nf(i1),0.d0)
            fluxl_c_nf(i1,1)=fluxl_c_nf(i1,1)+(a_l*ar_liq(kk)*dvl21+b_l*ar_liq(ii)*dvl11)
            fluxl_c_nf(i1,2)=fluxl_c_nf(i1,2)+(a_l*ar_liq(kk)*dvl22+b_l*ar_liq(ii)*dvl12)
            fluxg_c_nf(i1,1)=fluxg_c_nf(i1,1)+(a_g*ar_gas(kk)*dvg21+b_g*ar_gas(ii)*dvg11)
            fluxg_c_nf(i1,2)=fluxg_c_nf(i1,2)+(a_g*ar_gas(kk)*dvg22+b_g*ar_gas(ii)*dvg12)
            fluxd_c_nf(i1,1)=fluxd_c_nf(i1,1)+(a_d*ar_drp(kk)*dvd21+b_d*ar_drp(ii)*dvd11)
            fluxd_c_nf(i1,2)=fluxd_c_nf(i1,2)+(a_d*ar_drp(kk)*dvd22+b_d*ar_drp(ii)*dvd12)
         ENDDO
!
      ELSE
!
!........Computing cells
!
         nf_number=0
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
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
!...........Assume zero cell values when the phase fraction is less than alpha_min
!
            IF(cell%alphal_o(ii).le.alpha_min) THEN
               dvl11=0.d0
               dvl12=0.d0
               dvl13=0.d0
            ELSE
               dvl11=dvldx(ii,1,1)*dx11+dvldx(ii,2,1)*dx12+dvldx(ii,3,1)*dx13
               dvl12=dvldx(ii,1,2)*dx11+dvldx(ii,2,2)*dx12+dvldx(ii,3,2)*dx13
               dvl13=dvldx(ii,1,3)*dx11+dvldx(ii,2,3)*dx12+dvldx(ii,3,3)*dx13
            ENDIF
            IF(cell%alphag_o(ii).le.alpha_min) THEN
               dvg11=0.d0
               dvg12=0.d0
               dvg13=0.d0
            ELSE
               dvg11=dvgdx(ii,1,1)*dx11+dvgdx(ii,2,1)*dx12+dvgdx(ii,3,1)*dx13
               dvg12=dvgdx(ii,1,2)*dx11+dvgdx(ii,2,2)*dx12+dvgdx(ii,3,2)*dx13
               dvg13=dvgdx(ii,1,3)*dx11+dvgdx(ii,2,3)*dx12+dvgdx(ii,3,3)*dx13
            ENDIF
            IF(cell%alphad_o(ii).le.alpha_min) THEN
               dvd11=0.d0
               dvd12=0.d0
               dvd13=0.d0
            ELSE
               dvd11=dvddx(ii,1,1)*dx11+dvddx(ii,2,1)*dx12+dvddx(ii,3,1)*dx13
               dvd12=dvddx(ii,1,2)*dx11+dvddx(ii,2,2)*dx12+dvddx(ii,3,2)*dx13
               dvd13=dvddx(ii,1,3)*dx11+dvddx(ii,3,2)*dx12+dvddx(ii,3,3)*dx13
            ENDIF
!
            IF(cell%alphal_o(kk).le.alpha_min) THEN
               dvl21=0.d0
               dvl22=0.d0
               dvl23=0.d0
            ELSE
               dvl21=dvldx(kk,1,1)*dx21+dvldx(kk,2,1)*dx22+dvldx(kk,3,1)*dx23
               dvl22=dvldx(kk,1,2)*dx21+dvldx(kk,2,2)*dx22+dvldx(kk,3,2)*dx23
               dvl23=dvldx(kk,1,3)*dx21+dvldx(kk,2,3)*dx22+dvldx(kk,3,3)*dx23
            ENDIF
            IF(cell%alphag_o(kk).le.alpha_min) THEN
               dvg21=0.d0
               dvg22=0.d0
               dvg23=0.d0
            ELSE
               dvg21=dvgdx(kk,1,1)*dx21+dvgdx(kk,2,1)*dx22+dvgdx(kk,3,1)*dx23
               dvg22=dvgdx(kk,1,2)*dx21+dvgdx(kk,2,2)*dx22+dvgdx(kk,3,2)*dx23
               dvg23=dvgdx(kk,1,3)*dx21+dvgdx(kk,2,3)*dx22+dvgdx(kk,3,2)*dx23
            ENDIF
            IF(cell%alphad_o(kk).le.alpha_min) THEN
               dvd21=0.d0
               dvd22=0.d0
               dvd23=0.d0
            ELSE
               dvd21=dvddx(kk,1,1)*dx21+dvddx(kk,2,1)*dx22+dvddx(kk,3,1)*dx23
               dvd22=dvddx(kk,1,2)*dx21+dvddx(kk,2,2)*dx22+dvddx(kk,3,2)*dx23
               dvd23=dvddx(kk,1,3)*dx21+dvddx(kk,2,3)*dx22+dvddx(kk,3,3)*dx23
            ENDIF
!
!........Apply 2nd order upwind
!
            a_l=min(flux_l_nf(i1),0.d0)
            b_l=max(flux_l_nf(i1),0.d0)
            a_g=min(flux_g_nf(i1),0.d0)
            b_g=max(flux_g_nf(i1),0.d0)
            a_d=min(flux_d_nf(i1),0.d0)
            b_d=max(flux_d_nf(i1),0.d0)
            fluxl_c_nf(i1,1)=fluxl_c_nf(i1,1)+(a_l*ar_liq(kk)*dvl21+b_l*ar_liq(ii)*dvl11)
            fluxl_c_nf(i1,2)=fluxl_c_nf(i1,2)+(a_l*ar_liq(kk)*dvl22+b_l*ar_liq(ii)*dvl12)
            fluxl_c_nf(i1,3)=fluxl_c_nf(i1,3)+(a_l*ar_liq(kk)*dvl23+b_l*ar_liq(ii)*dvl13)
            fluxg_c_nf(i1,1)=fluxg_c_nf(i1,1)+(a_g*ar_gas(kk)*dvg21+b_g*ar_gas(ii)*dvg11)
            fluxg_c_nf(i1,2)=fluxg_c_nf(i1,2)+(a_g*ar_gas(kk)*dvg22+b_g*ar_gas(ii)*dvg12)
            fluxg_c_nf(i1,3)=fluxg_c_nf(i1,3)+(a_g*ar_gas(kk)*dvg23+b_g*ar_gas(ii)*dvg13)
            fluxd_c_nf(i1,1)=fluxd_c_nf(i1,1)+(a_d*ar_drp(kk)*dvd21+b_d*ar_drp(ii)*dvd11)
            fluxd_c_nf(i1,2)=fluxd_c_nf(i1,2)+(a_d*ar_drp(kk)*dvd22+b_d*ar_drp(ii)*dvd12)
            fluxd_c_nf(i1,3)=fluxd_c_nf(i1,3)+(a_d*ar_drp(kk)*dvd23+b_d*ar_drp(ii)*dvd13)
      ENDDO
!       
      ENDIF
!       
      END SUBROUTINE mom_2nd_conv
