!
      SUBROUTINE turb_2nd_conv_imp(turb_o,flux_nf,src_ke,diag_non)
!
!     This routine calculates convective mass fluxes through the cell face
!
      USE Zinterface
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: np
      USE Zparam       , ONLY: ndim
      USE Zvec_param   , ONLY: nf_non,nf_flux
      USE Znum_cell    , ONLY: istart_nf,                                 &
                                nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Z2nd_order   , ONLY: turb_conv_2nd,limiter_turb
      USE Zcoord3      , ONLY: volr      
      USE Zvec_geo     , ONLY: dxfc_nf,dxfc_non_k
!
      IMPLICIT NONE
!
!.....Input
      REAL(8) :: turb_o(ncell_fp),src_ke(ncell_fluid)
      REAL(8) :: flux_nf(nf_flux)
      REAL(8) :: diag_non(nf_non)
!.....Local variables 
      INTEGER :: i
      INTEGER :: nv,nf_number,istart,len,i1
      INTEGER :: ii,kk
      REAL(8) :: dx1,dx2,dx3
      REAL(8) :: do_non
!.....Local arrays
      REAL(8) :: src_ke0(ncell_fluid)
      REAL(8) :: dodx(ncell_fp,ndim)
!.....Local vector arrays
      REAL(8) :: src_non(nf_non)
!
!.....Gradient for 2nd order interpolation
!
      IF(turb_conv_2nd.eq.2)THEN
!
!........Frink method
!
!
      ELSE
!
!........Green-Gauss (Default)
!
         CALL grad_conv(turb_o,dodx)
!
      ENDIF
!
!.....Limiter of the gradient
!
      IF(limiter_turb.eq.1.or.limiter_turb.eq.2) THEN
         CALL grad_limiter(turb_o,dodx,limiter_turb)
         IF(np.gt.1) CALL communicate_2d(dodx)
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
      IF(ndim.eq.2)THEN
         DO i=1,len
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            IF(flux_nf(i1).lt.0)THEN
               dx1=dxfc_non_k(i,1)
               dx2=dxfc_non_k(i,2)
               do_non=dodx(kk,1)*dx1+dodx(kk,2)*dx2
               src_non(i)=diag_non(i)*do_non
            ELSE
               dx1=dxfc_nf(i1,1)
               dx2=dxfc_nf(i1,2)
               do_non=dodx(ii,1)*dx1+dodx(ii,2)*dx2
               src_non(i)=diag_non(i)*do_non
            ENDIF
         ENDDO
      ELSE
         DO i=1,len
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            IF(flux_nf(i1).lt.0)THEN
               dx1=dxfc_non_k(i,1)
               dx2=dxfc_non_k(i,2)
               dx3=dxfc_non_k(i,3)
               do_non=dodx(kk,1)*dx1+dodx(kk,2)*dx2+dodx(kk,3)*dx3
               src_non(i)=diag_non(i)*do_non
            ELSE
               dx1=dxfc_nf(i1,1)
               dx2=dxfc_nf(i1,2)
               dx3=dxfc_nf(i1,3)
               do_non=dodx(ii,1)*dx1+dodx(ii,2)*dx2+dodx(ii,3)*dx3
               src_non(i)=diag_non(i)*do_non
            ENDIF
         ENDDO
      ENDIF
!
      CALL sum_nf(0,-1,            &
                  src_non,src_ke0)
!
      DO i=1,ncell_fluid
         src_ke(i)=src_ke(i)+src_ke0(i)*volr(i)
      ENDDO
!       
      END SUBROUTINE turb_2nd_conv_imp
