!
      SUBROUTINE udfn_subchannel_scalar_conv_mixingvane
!
!     This routine calculates energy convective fluxes through the cell face
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell            
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: ndim
      USE Zvec_param   , ONLY: nf_nonk,nf_non
      USE Zvec_index   , ONLY: left_nf,right_non,jneigh_nf
      USE Znum_cell    , ONLY: istart_nf,i_neigh, &
                               nf_number_nb,lens, &
                               right_nb_k,istart_nfs,nf_number_id,istart_nf
      USE Zbc_index    , ONLY: npb
      USE Zconst1      , ONLY: vv_prob
      USE Zcoord3      , ONLY: volp
      USE Zporous      , ONLY: mixing_vane_l
      USE Ztimecon     , ONLY: alpha_min
      USE Zvector      , ONLY: vl_o
      USE Zvec_geo     , ONLY: xn_nf
      USE Zporous      , ONLY: mv_fac,chn_type
      !
      IMPLICIT NONE
!      
!      INCLUDE 'c3com.h'
!
!.....Local variables
      INTEGER :: i,j,k,ix,j0
      INTEGER :: nv,nf_number,istart0,istart,len,i0,i1,ii,kk
      REAL(8) :: mvM_conv_up0,mvE_conv_up0
      REAL(8) :: MVx1_l_m,MVy1_l_m
      REAL(8) :: MVx2_l_m,MVy2_l_m
      REAL(8) :: MVx1_l_e,MVy1_l_e
      REAL(8) :: MVx2_l_e,MVy2_l_e
      REAL(8) :: F_latconv
!.....Local arrays
      REAL(8),DIMENSION(ncell_fluid) :: MVx1_l_mv,MVy1_l_mv,        &
                                        MVx1_l_ev,MVy1_l_ev,        &
                                        mvM_conv_up_s,mvE_conv_up_s
!.....Local vector arrays
      REAL(8),DIMENSION(nf_non) :: mv_fac_non
      REAL(8),DIMENSION(nf_nonk+nf_non) :: mvM_conv_up,mvE_conv_up
!
      mixing_vane_l=0.d0

      IF(vv_prob.eq.'single_assem'         .or. &
         vv_prob.eq.'single_assem_new'     .or. &
         vv_prob.eq.'PSBT_5x5'                   )THEN
         F_latconv=0.27d0 !original
      ELSEIF(vv_prob.eq.'APR1400_fullcore'           .or. &
             vv_prob.eq.'APR1400_fullcore_modmesh01' .or. &
             vv_prob.eq.'OPR1000_fullcore'           .or. & 
             vv_prob.eq.'OPR1000_quarter_core'       .or. &
             vv_prob.eq.'OPR1000_single_assem'             )THEN
         F_latconv=0.27d0 !original
      ENDIF

      DO i=1,ncell_fluid
         IF(npb(i).ne.0) CYCLE
         mvM_conv_up0=0.d0
         mvE_conv_up0=0.d0
! 
!........Define cell values
!
         MVx1_l_mv(i)=cell%vfwl_x(i)*volp(i)/(F_latconv*vl_o(i,ndim))              ! mass
         MVy1_l_mv(i)=cell%vfwl_y(i)*volp(i)/(F_latconv*vl_o(i,ndim))
         MVx1_l_ev(i)=cell%vfwl_x(i)*volp(i)/(F_latconv*vl_o(i,ndim))*cell%hl(i)   ! 180205 KSB
         MVy1_l_ev(i)=cell%vfwl_y(i)*volp(i)/(F_latconv*vl_o(i,ndim))*cell%hl(i)
!
!........Assume zero cell values when the phase fraction is less than alpha_min
!
         IF(cell%alphag_o(i) .le. alpha_min) THEN
!            are1_gas = 0.d0
         ENDIF
         IF(cell%alphal_o(i) .le. alpha_min) THEN
!            are1_liq = 0.d0
!            are1_drp = 0.d0  
            MVx1_l_mv(i)=0.d0
            MVy1_l_mv(i)=0.d0
            MVx1_l_ev(i)=0.d0
            MVy1_l_ev(i)=0.d0
         ENDIF
      ENDDO       
!           
!......Vectorize mv_fac ==>must be CSR format
!
!
!......Build summation info for non,inl
!
       nf_number_nb=0
       nf_number_id(-1)=-1
       nf_number_id(0)=0
       istart_nfs(0)=istart_nfs(-1)+nf_nonk
       lens         =istart_nfs(0) +nf_non
!
      nv=0
      nf_number=nf_number_id(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         i1=istart+i
         ii=left_nf(i1)
         j0=i_neigh(ii)-1
         j=jneigh_nf(i1)+j0
         mv_fac_non(i1)=mv_fac(j)
      ENDDO
!
      nv=0
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)
         MVx2_l_m=MVx1_l_mv(kk)
         MVy2_l_m=MVy1_l_mv(kk)
         MVx2_l_e=MVx1_l_ev(kk)
         MVy2_l_e=MVy1_l_ev(kk)
!
!........Assume zero cell values when the phase fraction is less than alpha_min
!
         IF (cell%alphal_o(kk) .le. alpha_min) THEN
             MVx2_l_m=0.d0
             MVy2_l_m=0.d0
             MVx2_l_e=0.d0
             MVy2_l_e=0.d0
         ENDIF
         DO ix=1,2
            IF(mv_fac_non(i)*xn_nf(i1,ix).gt.0.d0)THEN
               IF(ix.eq.1)THEN
                  mvM_conv_up(i0)=MVx1_l_m
                  mvE_conv_up(i0)=MVx1_l_e
               ELSEIF(ix.eq.2)THEN
                  mvM_conv_up(i0)=MVy1_l_m
                  mvE_conv_up(i0)=MVy1_l_e
               ENDIF
            ELSEIF(mv_fac_non(i)*xn_nf(i1,ix).lt.0.d0)THEN
               IF(ix.eq.1)THEN
                  mvM_conv_up(i0)=MVx2_l_m
                  mvE_conv_up(i0)=MVx2_l_e
               ELSEIF(ix.eq.2)THEN
                  mvM_conv_up(i0)=MVy2_l_m
                  mvE_conv_up(i0)=MVy2_l_e
               ENDIF
            ENDIF
         ENDDO
      ENDDO
!
      nv=-1
      nf_number=nf_number_id(nv)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         k=right_nb_k(i)
         ii=right_non(k)
         MVx2_l_m=MVx1_l_mv(ii)
         MVy2_l_m=MVy1_l_mv(ii)
         MVx2_l_e=MVx1_l_ev(ii)
         MVy2_l_e=MVy1_l_ev(ii)
!
!........Assume zero cell values when the phase fraction is less than alpha_min
!
         IF (cell%alphal_o(ii) .le. alpha_min) THEN
             MVx2_l_m=0.d0
             MVy2_l_m=0.d0
             MVx2_l_e=0.d0
             MVy2_l_e=0.d0
         ENDIF
!
            MVx2_l_m=MVx1_l_mv(ii)
            MVy2_l_m=MVy1_l_mv(ii)
            MVx2_l_e=MVx1_l_ev(ii)
            MVy2_l_e=MVy1_l_ev(ii)
            IF(cell%alphal_o(ii) .le. alpha_min) THEN
               MVx2_l_m=0.d0
               MVy2_l_m=0.d0
               MVx2_l_e=0.d0
               MVy2_l_e=0.d0
            ENDIF
            DO ix=1,2
               IF(-mv_fac_non(k)*xn_nf(k,ix).gt.0.d0)THEN
                  IF(ix.eq.1)THEN
                     mvM_conv_up(i)=MVx1_l_m
                     mvE_conv_up(i)=MVx1_l_e
                  ELSEIF(ix.eq.2)THEN
                     mvM_conv_up(i)=MVy1_l_m
                     mvE_conv_up(i)=MVy1_l_e
                  ENDIF
               ELSEIF(-mv_fac_non(k)*xn_nf(k,ix).lt.0.d0)THEN
                  IF(ix.eq.1)THEN
                     mvM_conv_up(i)=MVx2_l_m
                     mvE_conv_up(i)=MVx2_l_e
                  ELSEIF(ix.eq.2)THEN
                     mvM_conv_up(i)=MVy2_l_m
                     mvE_conv_up(i)=MVy2_l_e
                  ENDIF
               ENDIF
            ENDDO
      ENDDO
!
      CALL sum_nf(0,0,                       &
                  mvM_conv_up,mvM_conv_up_s, &
                  mvE_conv_up,mvE_conv_up_s)

      DO i=1,ncell_fluid
         IF(npb(i).ne.0) CYCLE
         mixing_vane_l(2,i)=mvM_conv_up_s(i)
         mixing_vane_l(3,i)=mvE_conv_up_s(i)
      ENDDO
      IF(vv_prob.eq.'APR1400_fullcore'        .or. &
         vv_prob.eq.'OPR1000_fullcore'        .or. &
         vv_prob.eq.'OPR1000_quarter_core'    .or. &
         vv_prob.eq.'OPR1000_single_assem'          )THEN
         DO i=1,ncell_fluid
            IF(npb(i).ne.0) CYCLE
               IF(chn_type(i).gt.6) THEN
                  mixing_vane_l(2,i)=0.d0
                  mixing_vane_l(3,i)=0.d0
               ENDIF
         ENDDO
      ENDIF
!
      END SUBROUTINE udfn_subchannel_scalar_conv_mixingvane
