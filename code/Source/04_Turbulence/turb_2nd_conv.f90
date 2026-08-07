!
      SUBROUTINE turb_2nd_conv(ke,dp,ar,flux_nf,conv_k_non,conv_e_non)
!
!     This routine calculates convective mass fluxes through the cell face
!
      USE Zinterface
      USE VOL_DATA
      USE Zmpi         , ONLY: ncell_fp
      USE Zcore        , ONLY: np
      USE Zparam       , ONLY: ndim
      USE Z2nd_order   , ONLY: turb_conv_2nd,limiter_turb
      USE Znum_cell    , ONLY: istart_nf
      USE Zvec_param   , ONLY: nf_non,nf_flux
      USE Zvec_geo     , ONLY: dxfc_nf,    &
                               dxfc_non_k
      USE Zvec_index   , ONLY: left_nf,right_non
!
      IMPLICIT NONE
!
!.....Input
      REAL(8),DIMENSION(ncell_fp) :: ke,dp,ar
      REAL(8),DIMENSION(nf_flux) :: flux_nf
!.....Output
      REAL(8),DIMENSION(nf_non) :: conv_k_non,conv_e_non
!.....Local variables
      INTEGER :: i
      INTEGER :: ii,kk
      INTEGER :: nf_number,istart,len,i1
!
      REAL(8) :: k_2nd_conv,e_2nd_conv,dk1,dk2,de1,de2
      REAL(8) :: dx11,dx12,dx13
      REAL(8) :: dx21,dx22,dx23
!.....Local arrays
!     REAL(8) :: kmin(ncell_fp),kmax(ncell_fp),emin(ncell_fp),emax(ncell_fp)
      REAL(8),DIMENSION(ncell_fp,ndim)  :: dkdx,dedx
!
!.....Gradient for 2nd order interpolation
!
      IF(turb_conv_2nd.eq.2)THEN
!
!........Frink method
!
      ELSE
!
!........Green-Gauss (Default)
!
         CALL grad_conv(ke,dkdx)
         CALL grad_conv(dp,dedx)
!
      ENDIF
!
!.....Limiter of the gradient
!
      IF(limiter_turb.eq.1.or.limiter_turb.eq.2) THEN
         CALL grad_limiter(ke,dkdx,limiter_turb)
         CALL grad_limiter(dp,dedx,limiter_turb)
         IF(np.gt.1) CALL communicate_2d(dkdx, &
                                         dedx)
      ENDIF
!!
!!.....MIN-MAX values
!!
!      IF(limiter_turb.eq.3) THEN
!         DO i=1,ncell_fp
!            kmin(i)=ke(i)
!            kmax(i)=ke(i)
!            emin(i)=dp(i)
!            emax(i)=dp(i)            
!         ENDDO
!!
!!........Obtain min and max value among neighboring cells
!!
!         nf_number=0
!         istart=istart_nf(1,nf_number)
!         len   =istart_nf(2,nf_number)
!         DO i=1,len  
!            i1=istart+i 
!            ii=left_nf(i1)
!            kk=right_non(i)
!            kmax(ii)=DMAX1(kmax(ii),ke(ii))
!            kmin(ii)=DMIN1(kmin(ii),ke(ii))
!            kmax(kk)=DMAX1(kmax(kk),ke(kk))
!            kmin(kk)=DMIN1(kmin(kk),ke(kk))
!            emax(ii)=DMAX1(emax(ii),dp(ii))
!            emin(ii)=DMIN1(emin(ii),dp(ii))
!            emax(kk)=DMAX1(emax(kk),dp(kk))
!            emin(kk)=DMIN1(emin(kk),dp(kk))            
!         ENDDO
!         IF(np.gt.1)THEN
!            CALL communicate_1b(kmin)
!            CALL communicate_1b(kmax)
!            CALL communicate_1b(emin)
!            CALL communicate_1b(emax)            
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
            dk1=dkdx(ii,1)*dx11+dkdx(ii,2)*dx12
            de1=dedx(ii,1)*dx11+dedx(ii,2)*dx12
            dk2=dkdx(kk,1)*dx21+dkdx(kk,2)*dx22
            de2=dedx(kk,1)*dx21+dedx(kk,2)*dx22
!
!!
!!........MIN-MAX filter
!!
!         IF(limiter_turb.eq.3)THEN
!            dk1=DMAX1(dk1,kmin(ii)-ke(ii))
!            dk1=DMIN1(dk1,kmax(ii)-ke(ii))
!            dk2=DMAX1(dk2,kmin(kk)-ke(kk))
!            dk2=DMIN1(dk2,kmax(kk)-ke(kk))
!            de1=DMAX1(de1,emin(ii)-dp(ii))
!            de1=DMIN1(de1,emax(ii)-dp(ii))
!            de2=DMAX1(de2,emin(kk)-dp(kk))
!            de2=DMIN1(de2,emax(kk)-dp(kk))            
!         ENDIF
!
!........Apply 2nd order upwind
!
            IF    (flux_nf(i1).gt.0.d0)THEN
               k_2nd_conv=ar(ii)*dk1*flux_nf(i1)
               conv_k_non(i)=conv_k_non(i)+k_2nd_conv             ! minus sign is already considered in 'Turb_ke_calc_liq.f90: Line 268'
            ELSEIF(flux_nf(i1).lt.0.d0)THEN
               k_2nd_conv=ar(kk)*dk2*flux_nf(i1)
               conv_k_non(i)=conv_k_non(i)+k_2nd_conv             ! minus sign is already considered in 'Turb_ke_calc_liq.f90: Line 268'
            ENDIF
!
            IF    (flux_nf(i1).gt.0.d0)THEN
               e_2nd_conv=ar(ii)*de1*flux_nf(i1)
               conv_e_non(i)=conv_e_non(i)+e_2nd_conv             ! minus sign is already considered in 'Turb_ke_calc_liq.f90: Line 268'
            ELSEIF(flux_nf(i1).lt.0.d0)THEN
               e_2nd_conv=ar(kk)*de2*flux_nf(i1)
               conv_e_non(i)=conv_e_non(i)+e_2nd_conv             ! minus sign is already considered in 'Turb_ke_calc_liq.f90: Line 268'
            ENDIF
!         
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
            dk1=dkdx(ii,1)*dx11+dkdx(ii,2)*dx12+dkdx(ii,3)*dx13
            de1=dedx(ii,1)*dx11+dedx(ii,2)*dx12+dedx(ii,3)*dx13
            dk2=dkdx(kk,1)*dx21+dkdx(kk,2)*dx22+dkdx(kk,3)*dx23
            de2=dedx(kk,1)*dx21+dedx(kk,2)*dx22+dedx(kk,3)*dx23
!
!!
!!........MIN-MAX filter
!!
!         IF(limiter_turb.eq.3)THEN
!            dk1=DMAX1(dk1,kmin(ii)-ke(ii))
!            dk1=DMIN1(dk1,kmax(ii)-ke(ii))
!            dk2=DMAX1(dk2,kmin(kk)-ke(kk))
!            dk2=DMIN1(dk2,kmax(kk)-ke(kk))
!            de1=DMAX1(de1,emin(ii)-dp(ii))
!            de1=DMIN1(de1,emax(ii)-dp(ii))
!            de2=DMAX1(de2,emin(kk)-dp(kk))
!            de2=DMIN1(de2,emax(kk)-dp(kk))            
!         ENDIF
!
!........Apply 2nd order upwind
!
            IF    (flux_nf(i1).gt.0.d0)THEN
               k_2nd_conv=ar(ii)*dk1*flux_nf(i1)
               conv_k_non(i)=conv_k_non(i)+k_2nd_conv             ! minus sign is already considered in 'Turb_ke_calc_liq.f90: Line 268'
            ELSEIF(flux_nf(i1).lt.0.d0)THEN
               k_2nd_conv=ar(kk)*dk2*flux_nf(i1)
               conv_k_non(i)=conv_k_non(i)+k_2nd_conv             ! minus sign is already considered in 'Turb_ke_calc_liq.f90: Line 268'
            ENDIF
!
            IF    (flux_nf(i1).gt.0.d0)THEN
               e_2nd_conv=ar(ii)*de1*flux_nf(i1)
               conv_e_non(i)=conv_e_non(i)+e_2nd_conv             ! minus sign is already considered in 'Turb_ke_calc_liq.f90: Line 268'
            ELSEIF(flux_nf(i1).lt.0.d0)THEN
               e_2nd_conv=ar(kk)*de2*flux_nf(i1)
               conv_e_non(i)=conv_e_non(i)+e_2nd_conv             ! minus sign is already considered in 'Turb_ke_calc_liq.f90: Line 268'
            ENDIF
!         
         ENDDO
      ENDIF
!       
      END SUBROUTINE turb_2nd_conv
