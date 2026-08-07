!
      SUBROUTINE mass_2nd_conv_imp(fl_1_non,fl_2_non,fg_1_non,fg_2_non,src)
!
!     This routine calculates convective mass fluxes through the cell face
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: np
      USE Zparam       , ONLY: ndim
      USE Znum_cell    , ONLY: istart_nf
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Z2nd_order   , ONLY: mass_conv_2nd,limiter_mass
      USE Zvec_geo     , ONLY: dxfc_nf,         &
                               dxfc_non_k
      USE Zvec_param   , ONLY: nf_non
!
      IMPLICIT NONE
!
!.....Input
      REAL(8) :: fl_1_non(nf_non),fl_2_non(nf_non)
      REAL(8) :: fg_1_non(nf_non),fg_2_non(nf_non)
!.....Output
      REAL(8) :: src(ncell_fluid)
!.....Local variables
      INTEGER :: i
      INTEGER :: ii,kk
      INTEGER :: nf_number,istart,len,i1
      REAL(8) :: da11,da21
      REAL(8) :: da12,da22
      REAL(8) :: dx11,dx12,dx13
      REAL(8) :: dx21,dx22,dx23
!.....Local arrays
      REAL(8) :: smin(ncell_fp),smax(ncell_fp)
      REAL(8) :: dadx(ncell_fp,ndim)
!.....Local vector arrays
      REAL(8) :: src0_non_i(nf_non),src0_non_k(nf_non)
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
            smax(ii)=MAX(smax(ii),cell%alphag_o(kk))
            smin(ii)=MIN(smin(ii),cell%alphag_o(kk))
            smax(kk)=MAX(smax(kk),cell%alphag_o(ii))
            smin(kk)=MIN(smin(kk),cell%alphag_o(ii))
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
            da11=dadx(ii,1)*dx11+dadx(ii,2)*dx12
            da21=dadx(kk,1)*dx21+dadx(kk,2)*dx22
            da12=da21
            da22=da11
!
!...........MIN-MAX filter
!
            IF(limiter_mass.eq.3)THEN
               da11=MAX(da11,smin(ii)-cell%alphag_o(ii))
               da11=MIN(da11,smax(ii)-cell%alphag_o(ii))
               da21=MAX(da21,smin(kk)-cell%alphag_o(kk))
               da21=MIN(da21,smax(kk)-cell%alphag_o(kk))
!
               da12=MAX(da12,smin(ii)-cell%alphag_o(ii))
               da12=MIN(da12,smax(ii)-cell%alphag_o(ii))
               da22=MAX(da22,smin(kk)-cell%alphag_o(kk))
               da22=MIN(da22,smax(kk)-cell%alphag_o(kk))
            ENDIF
            src0_non_i(i)= (fg_1_non(i)*da11+fg_2_non(i)*da21)  &
                          +(fl_1_non(i)*da11+fl_2_non(i)*da21)
            src0_non_k(i)= (fg_1_non(i)*da12+fg_2_non(i)*da22)  &
                          +(fl_1_non(i)*da12+fl_2_non(i)*da22)
         ENDDO
      ELSE
         DO i=1,len
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
!
            dx11=dxfc_nf(i,1)
            dx12=dxfc_nf(i,2)
            dx13=dxfc_nf(i,3)
            dx21=dxfc_non_k(i,1)
            dx22=dxfc_non_k(i,2)
            dx23=dxfc_non_k(i,3)
            da11=dadx(ii,1)*dx11+dadx(ii,2)*dx12+dadx(ii,3)*dx13
            da21=dadx(kk,1)*dx21+dadx(kk,2)*dx22+dadx(kk,3)*dx23
            da12=da21
            da22=da11
!
!...........MIN-MAX filter
!
            IF(limiter_mass.eq.3)THEN
               da11=MAX(da11,smin(ii)-cell%alphag_o(ii))
               da11=MIN(da11,smax(ii)-cell%alphag_o(ii))
               da21=MAX(da21,smin(kk)-cell%alphag_o(kk))
               da21=MIN(da21,smax(kk)-cell%alphag_o(kk))
!
               da12=MAX(da12,smin(ii)-cell%alphag_o(ii))
               da12=MIN(da12,smax(ii)-cell%alphag_o(ii))
               da22=MAX(da22,smin(kk)-cell%alphag_o(kk))
               da22=MIN(da22,smax(kk)-cell%alphag_o(kk))
            ENDIF
            src0_non_i(i)= (fg_1_non(i)*da11+fg_2_non(i)*da21)  &
                          +(fl_1_non(i)*da11+fl_2_non(i)*da21)
            src0_non_k(i)= (fg_1_non(i)*da12+fg_2_non(i)*da22)  &
                          +(fl_1_non(i)*da12+fl_2_non(i)*da22)
         ENDDO
      ENDIF
!
      CALL mass_2nd_conv_imp_sum(src0_non_i,src0_non_k, &
                                 src) 
!       
      END SUBROUTINE mass_2nd_conv_imp
!
      SUBROUTINE mass_2nd_conv_imp_sum(src0_non_i,src0_non_k, &
                                       s1) 
!
      USE VOL_DATA     , ONLY: cell
      USE Zzone        , ONLY: ncell_fluid
      USE Zvec_param   , ONLY: nf_non
      USE Znum_cell    , ONLY: istart_nb1,               &
                               ia_nb,icell_nb,right_nb_k
      USE Zdecoupled   , ONLY: al_min_c,ag_min_c
!
      IMPLICIT NONE
!  
!.....Input
      REAL(8) :: src0_non_i(nf_non),src0_non_k(nf_non)
!.....Output
      REAL(8) :: s1(ncell_fluid)
!.....Local variables
      INTEGER :: i,k,nb
      INTEGER :: ii
      INTEGER :: nf_number,len,istart1,i1
      REAL(8) :: s1_s
!
!.....Partial sum src for non_k
!
      nf_number=-1
      len   =istart_nb1(2,nf_number)
      DO nb=1,len
         ii=icell_nb(nb)
         IF(cell%alphag_o(ii).ge.ag_min_c .or. cell%alphal_o(ii).lt.al_min_c) THEN
            s1_s=0.d0
!DIR$ NOVECTOR
            DO i=ia_nb(nb),ia_nb(nb+1)-1
               k=right_nb_k(i)
               s1_s=s1_s+src0_non_k(k)
            ENDDO
            s1(ii)=s1(ii)-s1_s
         ENDIF
      ENDDO
!
!.....Partial sum src for non
!
      nf_number=0
      istart1=istart_nb1(1,nf_number)
      len    =istart_nb1(2,nf_number)
      DO nb=1,len
         i1=istart1+nb
         ii=icell_nb(i1)
         IF(cell%alphag_o(ii).ge.ag_min_c .or. cell%alphal_o(ii).lt.al_min_c) THEN
            s1_s=0.d0
!DIR$ NOVECTOR
            DO i=ia_nb(i1),ia_nb(i1+1)-1
               s1_s=s1_s+src0_non_i(i)
            ENDDO
            s1(ii)=s1(ii)+s1_s
         ENDIF
      ENDDO
!
      END SUBROUTINE mass_2nd_conv_imp_sum
