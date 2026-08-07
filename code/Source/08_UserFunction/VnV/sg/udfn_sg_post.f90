!
      SUBROUTINE udfn_sg_post
!
      USE Zinterface
      USE Zmpi          , ONLY: ncell_fp
      USE Zzone         , ONLY: ncell_fluid
      USE Zcore         , ONLY: np,myrank
      USE Zvec_param    , ONLY: nf_nonk,nf_non,nf_inl 
      USE Znum_cell     , ONLY: istart_nf,istart_nbcon_nf,                   &
                                nf_number_nb,lens,nf_number_id,istart_nfs,   &
                                right_nb_k,istart_nfs,nf_number_id,istart_nf
      USE Zvec_index    , ONLY: left_nf,right_non,nbcon_nf
      USE Ztimecon      , ONLY: time
      USE Zare          , ONLY: ar_gas,ar_liq
      USE Zcoord1       , ONLY: xloc
      USE Zb_condition  , ONLY: rhob_liq
      USE Zsg           , ONLY: idc,z_econ,dz_fw
      USE Zvec_major    , ONLY: flux_l_nf,flux_g_nf
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,k,ii,kk
      INTEGER :: nv,nf_number,istart,len,istart0,istart2,i0,i1,i2
      LOGICAL,SAVE :: initial=.true.
      REAL(8) :: flow_cd_cold,flow_cd_hot,fw_econ,fw_dc,fw_tot
!.....Local arrays
      REAL(8) :: tmp(4)
      REAL(8),DIMENSION(ncell_fluid) :: flow_cd_coldav,flow_cd_hotav, &
                                        fw_econav,fw_dcav
!.....Local vector arrays
      REAL(8),DIMENSION(nf_nonk+nf_non) :: flow_cd_colda_non,flow_cd_hota_non
      REAL(8),DIMENSION(nf_inl) :: fw_econa_inl,fw_dca_inl
!.....Local allocatable arrays
      INTEGER,DIMENSION(:),ALLOCATABLE,SAVE :: idc_h,idc_c
!
      IF(initial) THEN
         initial=.false.
         IF(myrank.eq.0) OPEN(999,file='VD13_sg_ref.dat')
!
         ALLOCATE(idc_h(ncell_fp),idc_c(ncell_fp))
         DO i=1,ncell_fluid
            idc_h(i)=0
            idc_c(i)=0
            IF(idc(i).eq.2)THEN
               IF(xloc(i,3).lt.dz_fw.and.xloc(i,1).lt.0.d0) idc_h(i)=1
               IF(xloc(i,3).gt.(z_econ-dz_fw).and.xloc(i,3).lt.z_econ.and.xloc(i,1).gt.0.d0) idc_c(i)=1
            ENDIF
         ENDDO
! 
         IF(np.gt.1) CALL communicate_1d_int(idc,   &
                                             idc_h, &
                                             idc_c)
!
      ENDIF
!
!.....Build summation info for non
!
      nf_number_nb=0
      nf_number_id(-1)=-1
      nf_number_id(0)=0
      istart_nfs(-1)=0
      istart_nfs(0)=istart_nfs(-1)+nf_nonk
      lens         =istart_nfs(0) +nf_non
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
         IF(idc_c(ii).gt.0 .and. idc(kk).eq.1) THEN
            flow_cd_colda_non(i0)=ar_liq(ii)*flux_l_nf(i1)+ar_gas(ii)*flux_g_nf(i1)
         ENDIF
         IF(idc_h(ii).gt.0 .and. idc(kk).eq.1) THEN
            flow_cd_hota_non(i0)=ar_liq(ii)*flux_l_nf(i1)+ar_gas(ii)*flux_g_nf(i1)
         ENDIF
      ENDDO
!
      nv=-1
      nf_number=nf_number_id(nv)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         k=right_nb_k(i)
         ii=right_non(k)
         kk=left_nf(k)
         IF(idc_c(ii).gt.0 .and. idc(kk).eq.1) THEN
            flow_cd_colda_non(i)=ar_liq(ii)*flux_l_nf(k)+ar_gas(ii)*flux_g_nf(k)
         ENDIF
         IF(idc_h(ii).gt.0 .and. idc(kk).eq.1) THEN
            flow_cd_hota_non(i)=ar_liq(ii)*flux_l_nf(k)+ar_gas(ii)*flux_g_nf(k)
         ENDIF
      ENDDO
      CALL sum_nf0_idc(flow_cd_colda_non,flow_cd_coldav,idc_c, &
                       flow_cd_hota_non,flow_cd_hotav,idc_h)
!
      flow_cd_cold=0.d0
      flow_cd_hot=0.d0
      DO i=1,ncell_fluid
         IF(idc_c(i).gt.0) THEN
            flow_cd_cold=flow_cd_cold+flow_cd_coldav(i)
         ENDIF
         IF(idc_h(i).gt.0) THEN
            flow_cd_hot=flow_cd_hot+flow_cd_hotav(i)
         ENDIF
      ENDDO
!
!.....Build summation info inl
!
      nf_number_nb=0
      nf_number_id(0)=2
      istart_nfs(0)=0
      lens         =istart_nfs(0)+nf_inl
!
      nv=0
      nf_number=nf_number_id(nv)
      istart=istart_nf(1,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i1=istart+i
         i2=istart2+i
         ii=left_nf(i1)
         k=nbcon_nf(i2)
         IF(k.eq.1)THEN
            fw_econa_inl(i)=rhob_liq(k)*flux_l_nf(i1)
         ELSEIF(k.eq.2)THEN
            fw_dca_inl(i)=rhob_liq(k)*flux_l_nf(i1)
         ENDIF
      ENDDO
!
      CALL sum_nf2_m(fw_econa_inl,fw_econav, &
                     fw_dca_inl,fw_dcav)
!
      fw_econ=0.d0
      fw_dc  =0.d0
      DO i=1,ncell_fluid
         fw_econ=fw_econ+fw_econav(i)
         fw_dc  =fw_dc  +fw_dcav(i)
      ENDDO
!
      IF(np.gt.1)THEN
         tmp(1)=flow_cd_cold
         tmp(2)=flow_cd_hot
         tmp(3)=fw_econ
         tmp(4)=fw_dc
         CALL allreducei_r(tmp,4)
         flow_cd_cold=tmp(1)
         flow_cd_hot=tmp(2)
         fw_econ=tmp(3)
         fw_dc=tmp(4)
      ENDIF
      fw_tot=fw_econ+fw_dc
      IF(myrank.eq.0) WRITE(999,10) time,flow_cd_cold/fw_tot,flow_cd_hot/fw_tot,(flow_cd_cold+flow_cd_hot)/fw_tot
   10 FORMAT(4(e22.15,1x))
!
      END SUBROUTINE udfn_sg_post
