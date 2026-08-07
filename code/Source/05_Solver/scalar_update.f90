!
      SUBROUTINE scalar_update(ncvgd)
!
!     This routine updates scalar variables after the pressure calculation
!
      USE Zinterface
      USE VOL_DATA      , ONLY: cell
      USE Zzone         , ONLY: ncell_fluid
      USE Zcore         , ONLY: np
      USE Ztimecon      , ONLY: dt_opt
      USE Zvec_param    , ONLY: nf_nonk,nf_non,nf_mcc,nf_inl,nf_fluxk1
      USE Znum_cell     , ONLY: right_nb_k,istart_nf,                     &
                                nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zbc_index     , ONLY: npb,icell_type
      USE Zvec_index    , ONLY: right_non,left_nf
      USE Zvec_major    , ONLY: flux_l_nf,flux_g_nf,flux_d_nf
      USE Zdel_scalar   , ONLY: del_eg,del_el,del_x
      USE Zscalar_coeff , ONLY: l_th_equil
      USE Zuserdefined  , ONLY: udfl_update_scalar
      USE Zconst1       , ONLY: wconden,lwconden_alphal0
      USE Zscalar_coeff , ONLY: sfg_nf,sfl_nf,sfd_nf,         &
                                sfg_non_k,sfl_non_k,sfd_non_k
      USE Zscalar_coeff , ONLY: sb
      USE Zcheck_scalar , ONLY: eps_eng,eps_vol
!
      IMPLICIT NONE
!
!.....Input
      LOGICAL :: ncvgd(3)
!.....Local variables
      INTEGER :: i,k
      INTEGER :: ii
      INTEGER :: nv,nf_number,istart0,istart,len,i0,i1
!      REAL(8) :: eps_eng=1.d-1
!      REAL(8) :: eps_vol=1.d-6
!.....Local arrays
      INTEGER :: ncvgdi(3),l_tmp(3)
      REAL(8),DIMENSION(ncell_fluid) :: eng_g,eng_l,alpha_g,alpha_l,alpha_d,qu
!.....Local vector arrays
      REAL(8),DIMENSION(nf_fluxk1) :: eng_g_nf,eng_l_nf,alpha_g_nf,alpha_d_nf,qu_nf
!
!.....Apply user subroutine when udfl_update_scalar=.true.
!
      IF(udfl_update_scalar)THEN
         CALL udfn_update_scalar(ncvgdi)
         ncvgd(1)=.false. 
         ncvgd(2)=.false. 
         ncvgd(3)=.false. 
         IF(ncvgdi(1).gt.0) ncvgd(1)=.true. 
         IF(ncvgdi(2).gt.0) ncvgd(2)=.true. 
         IF(ncvgdi(3).gt.0) ncvgd(3)=.true. 
         RETURN
      ENDIF
!
!.....Initialize re-calculation flags
!
      ncvgdi(:)=0
!
!.....Build summation info for non,mcc,inl
!
      nf_number_nb=2
      nf_number_id(-1)=-1
      nf_number_id(0)=0
      nf_number_id(1)=1
      nf_number_id(2)=2
      istart_nfs(-1)=0
      istart_nfs(0)=istart_nfs(-1)+nf_nonk
      istart_nfs(1)=istart_nfs(0) +nf_non
      istart_nfs(2)=istart_nfs(1) +nf_mcc
      lens         =istart_nfs(2) +nf_inl
!
!...........Convection terms
!
      nv=-1
      nf_number=nf_number_id(nv)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         k=right_nb_k(i)
         ii=right_non(k) 
         IF(npb(ii).eq.0) THEN
            eng_g_nf(i)  =-(sfg_non_k(i,1)*flux_g_nf(k)+sfl_non_k(i,1)*flux_l_nf(k)+sfd_non_k(i,1)*flux_d_nf(k))
            eng_l_nf(i)  =-(sfg_non_k(i,2)*flux_g_nf(k)+sfl_non_k(i,2)*flux_l_nf(k)+sfd_non_k(i,2)*flux_d_nf(k))
            alpha_g_nf(i)=-(sfg_non_k(i,3)*flux_g_nf(k)+sfl_non_k(i,3)*flux_l_nf(k)+sfd_non_k(i,3)*flux_d_nf(k))
            alpha_d_nf(i)=-(sfg_non_k(i,4)*flux_g_nf(k)+sfl_non_k(i,4)*flux_l_nf(k)+sfd_non_k(i,4)*flux_d_nf(k))
            qu_nf(i)     =-(sfg_non_k(i,5)*flux_g_nf(k)+sfl_non_k(i,5)*flux_l_nf(k)+sfd_non_k(i,5)*flux_d_nf(k))
         ELSE
            eng_g_nf(i)  =0.d0
            eng_l_nf(i)  =0.d0
            alpha_g_nf(i)=0.d0
            alpha_d_nf(i)=0.d0
            qu_nf(i)     =0.d0
         ENDIF
      ENDDO
!
      DO nv=0,2
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(nv)
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len
            i0=istart0+i
            i1=istart+i
            ii=left_nf(i1)
            IF(npb(ii).eq.0) THEN
               eng_g_nf(i0)  =sfg_nf(i1,1)*flux_g_nf(i1)+sfl_nf(i1,1)*flux_l_nf(i1)+sfd_nf(i1,1)*flux_d_nf(i1)
               eng_l_nf(i0)  =sfg_nf(i1,2)*flux_g_nf(i1)+sfl_nf(i1,2)*flux_l_nf(i1)+sfd_nf(i1,2)*flux_d_nf(i1)
               alpha_g_nf(i0)=sfg_nf(i1,3)*flux_g_nf(i1)+sfl_nf(i1,3)*flux_l_nf(i1)+sfd_nf(i1,3)*flux_d_nf(i1)
               alpha_d_nf(i0)=sfg_nf(i1,4)*flux_g_nf(i1)+sfl_nf(i1,4)*flux_l_nf(i1)+sfd_nf(i1,4)*flux_d_nf(i1)
               qu_nf(i0)     =sfg_nf(i1,5)*flux_g_nf(i1)+sfl_nf(i1,5)*flux_l_nf(i1)+sfd_nf(i1,5)*flux_d_nf(i1)
            ELSE
               eng_g_nf(i0)  =0.d0
               eng_l_nf(i0)  =0.d0
               alpha_g_nf(i0)=0.d0
               alpha_d_nf(i0)=0.d0
               qu_nf(i0)     =0.d0
            ENDIF
         ENDDO
      ENDDO
!
      CALL sum_nf(0,0,                &
                  eng_g_nf  ,eng_g,   &
                  eng_l_nf  ,eng_l,   &
                  alpha_g_nf,alpha_g, &
                  alpha_d_nf,alpha_d, &
                  qu_nf     ,qu)
!
      DO i=1,ncell_fluid
         IF(npb(i).ne.0) cycle
         eng_g(i)  =eng_g(i)  +sb(i,1)
         eng_l(i)  =eng_l(i)  +sb(i,2)
         alpha_g(i)=alpha_g(i)+cell%alphag_o(i)+sb(i,3)
         alpha_d(i)=alpha_d(i)+cell%alphad_o(i)+sb(i,4)
         qu(i)     =qu(i)     +sb(i,5)
      ENDDO
!
      DO i=1,ncell_fluid
         IF(npb(i).ne.0) cycle
!
!........Limit energy change less than eps*old_value
!
         IF (ABS(eng_l(i)).gt.eps_eng*cell%el_o(i)) ncvgdi(1)=1
         IF (ABS(eng_g(i)).gt.eps_eng*cell%eg_o(i)) ncvgdi(2)=1
!
!........Limit void fraction error less than eps_vol
!
         IF (alpha_g(i)      .lt.-eps_vol) ncvgdi(3)=1
         IF (alpha_g(i)-1.d0 .gt. eps_vol) ncvgdi(3)=1
!
!...........Update gas and liquid internal energy
!
         eng_g(i)=cell%eg_o(i)+eng_g(i)
         eng_l(i)=cell%el_o(i)+eng_l(i)
!
!........Limit gas and liquid fraction between 0.0 and 1.0
!
         alpha_g(i)=MIN(MAX(0.d0,alpha_g(i)),1.d0)
         alpha_d(i)=MIN(MAX(0.d0,alpha_d(i)),1.d0)
         alpha_l(i)=1.d0-alpha_g(i)-alpha_d(i)
!
      ENDDO
      ncvgd(1)=.false. 
      ncvgd(2)=.false. 
      ncvgd(3)=.false. 
      IF(ncvgdi(1).gt.0) ncvgd(1)=.true. 
      IF(ncvgdi(2).gt.0) ncvgd(2)=.true. 
      IF(ncvgdi(3).gt.0) ncvgd(3)=.true. 
!
!.....Convergence test
!
      l_tmp(1)=ncvgdi(1)+ncvgdi(2)+ncvgdi(3)
      l_tmp(2)=ncvgdi(1)+ncvgdi(2)
      l_tmp(3)=ncvgdi(3)
      IF(np.gt.1) CALL allreducei_i(l_tmp,3)
!
!.....Update scalar variables if the calculation has converged
!.....Otherwise retrune for re-calculation
!
      IF (((l_tmp(1).gt.0) .and. (dt_opt.eq.2.or.dt_opt.eq.3)) .or.  &
          ((l_tmp(2).gt.0) .and. (dt_opt.eq.4.or.dt_opt.eq.5)) .or.  &
          ((l_tmp(3).gt.0) .and. (dt_opt.eq.6.or.dt_opt.eq.7)))THEN    
         RETURN
      ELSE
         DO i=1,ncell_fluid
            IF(npb(i).ne.0) cycle
               cell%alphag(i)=alpha_g(i)
               cell%alphal(i)=alpha_l(i)
               cell%alphad(i)=alpha_d(i)
               cell%el(i)=eng_l(i)
               cell%eg(i)=eng_g(i)
               cell%quala(i)=cell%quala_o(i)+qu(i)
               IF(l_th_equil.gt.0)THEN
                  IF(cell%alphal(i).le.1.d-5) cell%el(i)=cell%elsat(i)
                  IF(cell%alphag(i).le.1.d-5) cell%eg(i)=cell%egsat(i)
               ELSE
                  IF(cell%alphal(i).le.1.d-8) cell%el(i)=cell%elsat(i)
                  IF(cell%alphag(i).le.1.d-8) cell%eg(i)=cell%egsat(i)
               ENDIF
               cell%quala(i)=MIN(MAX(0.d0,cell%quala(i)),1.d0)
          ENDDO
      ENDIF
!
      IF(wconden.ne.0.and.lwconden_alphal0)THEN
         DO i=1,ncell_fluid
            IF(icell_type(i).eq.1) THEN
               cell%alphag(i)=1.d0
               cell%alphal(i)=0.d0
               cell%alphad(i)=0.d0
            ENDIF
         ENDDO
      ENDIF
!
      CALL ncg_transport
!
!.....Update material properties
!
      CALL property_calc(1)
!
!.....Vapor generation rate

      CALL set_vapor_generation
!
!.....Save scalar changes
!
      DO i=1,ncell_fluid
         del_eg(i)=cell%eg(i)-cell%eg_o(i)
         del_el(i)=cell%el(i)-cell%el_o(i)
         del_x(i)=cell%quala(i)-cell%quala_o(i)
      ENDDO
!
      END SUBROUTINE scalar_update
