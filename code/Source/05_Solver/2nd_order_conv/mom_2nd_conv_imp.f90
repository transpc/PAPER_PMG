!
      SUBROUTINE mom_2nd_conv_imp(diag_g_non,diag_l_non, &
                                  src_g,src_l)
!
!     This routine calculates 2nd oder terms of the momentum convection
!
      USE Zinterface
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid_pad
      USE Zcore        , ONLY: np
      USE Zparam       , ONLY: ndim
      USE Zvec_param   , ONLY: nf_non
      USE Znum_cell    , ONLY: istart_nf,                                &
                               nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Z2nd_order   , ONLY: mom_conv_2nd,limiter_mom
      USE Zb_condition , ONLY: vb_liq,vin_liq,vb_gas,vin_gas
      USE Zvector      , ONLY: vl_o,vg_o
      USE Zvec_geo     , ONLY: dxfc_nf,dxfc_non_k
      USE Zvec_major   , ONLY: flux_l_nf,flux_g_nf
!
      IMPLICIT NONE
!
!.....Input
      REAL(8),DIMENSION(nf_non) :: diag_g_non,diag_l_non
!.....Output
      REAL(8),DIMENSION(ncell_fluid_pad,ndim) :: src_g,src_l
!.....Local variables
      INTEGER :: i,ii,kk,ix
      INTEGER :: nv,nf_number,istart,len,i1
      REAL(8) :: dx1,dx2,dx3
      REAL(8) :: dvlx,dvly,dvlz,dvgx,dvgy,dvgz
!.....Local vector arrays
      REAL(8),DIMENSION(ncell_fp,ndim,ndim) :: dvldx,dvgdx
      REAL(8),DIMENSION(nf_non,ndim) :: src_l_non(nf_non,ndim),src_g_non
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
!
      ENDIF
      IF(np.gt.1) CALL communicate_3d(dvldx, &
                                      dvgdx)
!
!.....Limiter of the gradient
!
      IF(limiter_mom.gt.0)THEN
         DO ix=1,ndim
            CALL grad_limiter(vl_o(1,ix),dvldx(1,1,ix),limiter_mom)
            CALL grad_limiter(vg_o(1,ix),dvgdx(1,1,ix),limiter_mom)
         ENDDO
         IF(np.gt.1) CALL communicate_3d(dvldx, &
                                         dvgdx)
      ENDIF
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
            IF(flux_l_nf(i1).lt.0)THEN
               dx1=dxfc_non_k(i,1)
               dx2=dxfc_non_k(i,2)
               dvlx=dvldx(kk,1,1)*dx1+dvldx(kk,2,1)*dx2
               dvly=dvldx(kk,1,2)*dx1+dvldx(kk,2,2)*dx2
               src_l_non(i,1)=diag_l_non(i)*dvlx
               src_l_non(i,2)=diag_l_non(i)*dvly
            ELSEIF(flux_l_nf(i1).gt.0)THEN
               dx1=dxfc_nf(i1,1)
               dx2=dxfc_nf(i1,2)
               dvlx=dvldx(ii,1,1)*dx1+dvldx(ii,2,1)*dx2
               dvly=dvldx(ii,1,2)*dx1+dvldx(ii,2,2)*dx2
               src_l_non(i,1)=diag_l_non(i)*dvlx
               src_l_non(i,2)=diag_l_non(i)*dvly
            ELSE
               src_l_non(i,1)=0.d0
               src_l_non(i,2)=0.d0
         ENDIF
            IF(flux_g_nf(i1).lt.0)THEN
               dx1=dxfc_non_k(i,1)
               dx2=dxfc_non_k(i,2)
               dvgx=dvgdx(kk,1,1)*dx1+dvgdx(kk,2,1)*dx2
               dvgy=dvgdx(kk,1,2)*dx1+dvgdx(kk,2,2)*dx2
               src_g_non(i,1)=diag_g_non(i)*dvgx
               src_g_non(i,2)=diag_g_non(i)*dvgy
            ELSEIF(flux_g_nf(i1).gt.0)THEN
               dx1=dxfc_nf(i1,1)
               dx2=dxfc_nf(i1,2)
               dvgx=dvgdx(ii,1,1)*dx1+dvgdx(ii,2,1)*dx2
               dvgy=dvgdx(ii,1,2)*dx1+dvgdx(ii,2,2)*dx2
               src_g_non(i,1)=diag_g_non(i)*dvgx
               src_g_non(i,2)=diag_g_non(i)*dvgy
            ELSE
               src_g_non(i,1)=0.d0
               src_g_non(i,2)=0.d0
            ENDIF
         ENDDO
      ELSE
         DO i=1,len
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            IF(flux_l_nf(i1).lt.0)THEN
               dx1=dxfc_non_k(i,1)
               dx2=dxfc_non_k(i,2)
               dx3=dxfc_non_k(i,3)
               dvlx=dvldx(kk,1,1)*dx1+dvldx(kk,2,1)*dx2+dvldx(kk,3,1)*dx3
               dvly=dvldx(kk,1,2)*dx1+dvldx(kk,2,2)*dx2+dvldx(kk,3,2)*dx3
               dvlz=dvldx(kk,1,3)*dx1+dvldx(kk,2,3)*dx2+dvldx(kk,3,3)*dx3
               src_l_non(i,1)=diag_l_non(i)*dvlx
               src_l_non(i,2)=diag_l_non(i)*dvly
               src_l_non(i,3)=diag_l_non(i)*dvlz
            ELSEIF(flux_l_nf(i1).gt.0)THEN
               dx1=dxfc_nf(i1,1)
               dx2=dxfc_nf(i1,2)
               dx3=dxfc_nf(i1,3)
               dvlx=dvldx(ii,1,1)*dx1+dvldx(ii,2,1)*dx2+dvldx(ii,3,1)*dx3
               dvly=dvldx(ii,1,2)*dx1+dvldx(ii,2,2)*dx2+dvldx(ii,3,2)*dx3
               dvlz=dvldx(ii,1,3)*dx1+dvldx(ii,2,3)*dx2+dvldx(ii,3,3)*dx3
               src_l_non(i,1)=diag_l_non(i)*dvlx
               src_l_non(i,2)=diag_l_non(i)*dvly
               src_l_non(i,3)=diag_l_non(i)*dvlz
            ELSE
               src_l_non(i,1)=0.d0
               src_l_non(i,2)=0.d0
               src_l_non(i,3)=0.d0
            ENDIF
            IF(flux_g_nf(i1).lt.0)THEN
               dx1=dxfc_non_k(i,1)
               dx2=dxfc_non_k(i,2)
               dx3=dxfc_non_k(i,3)
               dvgx=dvgdx(kk,1,1)*dx1+dvgdx(kk,2,1)*dx2+dvgdx(kk,3,1)*dx3
               dvgy=dvgdx(kk,1,2)*dx1+dvgdx(kk,2,2)*dx2+dvgdx(kk,3,2)*dx3
               dvgz=dvgdx(kk,1,3)*dx1+dvgdx(kk,2,3)*dx2+dvgdx(kk,3,3)*dx3
               src_g_non(i,1)=diag_g_non(i)*dvgx
               src_g_non(i,2)=diag_g_non(i)*dvgy
               src_g_non(i,3)=diag_g_non(i)*dvgz
            ELSEIF(flux_g_nf(i1).gt.0)THEN
               dx1=dxfc_nf(i1,1)
               dx2=dxfc_nf(i1,2)
               dx3=dxfc_nf(i1,3)
               dvgx=dvgdx(ii,1,1)*dx1+dvgdx(ii,2,1)*dx2+dvgdx(ii,3,1)*dx3
               dvgy=dvgdx(ii,1,2)*dx1+dvgdx(ii,2,2)*dx2+dvgdx(ii,3,2)*dx3
               dvgz=dvgdx(ii,1,3)*dx1+dvgdx(ii,2,3)*dx2+dvgdx(ii,3,3)*dx3
               src_g_non(i,1)=diag_g_non(i)*dvgx
               src_g_non(i,2)=diag_g_non(i)*dvgy
               src_g_non(i,3)=diag_g_non(i)*dvgz
            ELSE
               src_g_non(i,1)=0.d0
               src_g_non(i,2)=0.d0
               src_g_non(i,3)=0.d0
            ENDIF
         ENDDO
      ENDIF
!
      CALL sum_nf_ndim(1,-1,ncell_fluid_pad, &
                       src_g_non,src_g,      &
                       src_l_non,src_l)
!       
      END SUBROUTINE mom_2nd_conv_imp
