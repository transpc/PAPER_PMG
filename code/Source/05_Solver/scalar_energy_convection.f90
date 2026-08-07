!
      SUBROUTINE scalar_energy_convection
!
!     This routine calculates energy convective fluxes through the cell face
!
      USE VOL_DATA                 
      USE Z2nd_order   , ONLY: eng_conv_2nd
      USE Zare         , ONLY: are_liq,are_gas,are_drp
      USE Zb_condition , ONLY: alphab_liq,alphab_gas,alphab_drp,eb_liq,eb_gas,rhob_liq,rhob_gas
      USE Ztimecon     , ONLY: alpha_min
      USE Zmodel       , ONLY: i_droplet
      USE c3com_cupid  , ONLY: i3invtbl,mcdirect,c3dpv,mcgdirect
!
      USE Znum_cell    , ONLY: istart_nf,istart_nbcon_nf
      USE Zvec_index   , ONLY: left_nf,right_non,nbcon_nf
      USE Zvec_major   , ONLY: flux_l_nf,flux_g_nf,flux_d_nf,     &
                               ecnvc_l_nf,ecnvc_g_nf,ecnvc_d_nf
      USE Zrv_model    , ONLY: rv_valve
!
      IMPLICIT NONE
!      
      INCLUDE '../10_LinkToMARS/c3com.h'
!
!     local variables
      INTEGER :: i,k
      INTEGER :: ii,kk
      INTEGER :: idx
      INTEGER :: nf_number,istart,len,istart2,i1,i2
!
      REAL(8) :: ali_tmp,agi_tmp,alk_tmp,agk_tmp
      REAL(8) :: arel_i,areg_i,ared_i
      REAL(8) :: arel_k,areg_k,ared_k
!
!.....Computing Cells
!
!====> non
      nf_number=0
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)
         ali_tmp=cell%alphal_o(ii)
         agi_tmp=cell%alphag_o(ii)
         alk_tmp=cell%alphal_o(kk)
         agk_tmp=cell%alphag_o(kk)
         arel_i=are_liq(ii)
         areg_i=are_gas(ii)
         ared_i=are_drp(ii)
         arel_k=are_liq(kk)
         areg_k=are_gas(kk)
         ared_k=are_drp(kk)
!                        
         IF(agi_tmp.le.alpha_min) areg_i=0.0d0
         IF(ali_tmp.le.alpha_min) THEN
            arel_i=0.0d0
            ared_i=0.0d0
         ENDIF
         IF(agk_tmp.le.alpha_min) areg_k=0.0d0
         IF(alk_tmp.le.alpha_min)THEN
            arel_k=0.0d0
            ared_k=0.0d0
         ENDIF
!
         IF    (flux_l_nf(i1).gt.0.0d0)THEN
            ecnvc_l_nf(i1)=arel_i
         ELSEIF(flux_l_nf(i1).lt.0.0d0)THEN
            ecnvc_l_nf(i1)=arel_k
         ELSE
            ecnvc_l_nf(i1)=(arel_i+arel_k)*0.5d0
         ENDIF
         IF    (flux_g_nf(i1).gt.0.0d0)THEN
            ecnvc_g_nf(i1)=areg_i
         ELSEIF(flux_g_nf(i1).lt.0.0d0)THEN
            ecnvc_g_nf(i1)=areg_k
         ELSE
            ecnvc_g_nf(i1)=(areg_i+areg_k)*0.5d0
         ENDIF
         IF    (flux_d_nf(i1).gt.0.0d0)THEN
            ecnvc_d_nf(i1)=ared_i
         ELSEIF(flux_d_nf(i1).lt.0.0d0)THEN
            ecnvc_d_nf(i1)=ared_k
         ELSE
            ecnvc_d_nf(i1)=(ared_i+ared_k)*0.5d0
         ENDIF
!
!..............User defined energy convection of droplets
!   
         IF(i_droplet.ge.1) ecnvc_d_nf(i1)=0.d0
      ENDDO
!
!.....valve model
!        
      IF(rv_valve.eq.1) CALL valve_model_scalar_energy_convection
!
!.....2nd order convection
!
      IF(eng_conv_2nd.gt.0) CALL energy_2nd_conv
!
!.....Inlet
!
      nf_number=2
      istart=istart_nf(1,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i1=istart+i
         i2=istart2+i
         ii=left_nf(i1)
         k=nbcon_nf(i2)
         ali_tmp=cell%alphal_o(ii)
         agi_tmp=cell%alphag_o(ii)
!                  
         arel_i=are_liq(ii)
         areg_i=are_gas(ii)
         ared_i=are_drp(ii)
         arel_k=alphab_liq(k)*rhob_liq(k)*eb_liq(k)
         areg_k=alphab_gas(k)*rhob_gas(k)*eb_gas(k)
         ared_k=alphab_drp(k)*rhob_liq(k)*eb_liq(k)
!                        
         IF(agi_tmp.le.alpha_min) areg_i=0.0d0
         IF(ali_tmp.le.alpha_min) THEN
            arel_i=0.0d0
            ared_i=0.0d0
         ENDIF
!
         IF    (flux_l_nf(i1).gt.0.0d0)THEN
            ecnvc_l_nf(i1)=arel_i
         ELSEIF(flux_l_nf(i1).lt.0.0d0)THEN
            ecnvc_l_nf(i1)=arel_k
         ELSE
            ecnvc_l_nf(i1)=(arel_i+arel_k)*0.5d0
         ENDIF
         IF    (flux_g_nf(i1).gt.0.0d0)THEN
            ecnvc_g_nf(i1)=areg_i
         ELSEIF(flux_g_nf(i1).lt.0.0d0)THEN
            ecnvc_g_nf(i1)=areg_k
         ELSE
            ecnvc_g_nf(i1)=(areg_i+areg_k)*0.5d0
         ENDIF
         IF    (flux_d_nf(i1).gt.0.0d0)THEN
            ecnvc_d_nf(i1)=ared_i
         ELSEIF(flux_d_nf(i1).lt.0.0d0)THEN
            ecnvc_d_nf(i1)=ared_k
         ELSE
            ecnvc_d_nf(i1)=(ared_i+ared_k)*0.5d0
         ENDIF
!
!........User defined energy convection of droplets
!  
         IF(i_droplet.ge.1) ecnvc_d_nf(i1)=0.d0
      ENDDO
!
!.....Outlet
!
      nf_number=3
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i1=istart+i
         ii=left_nf(i1)
         ali_tmp=cell%alphal_o(ii)
         agi_tmp=cell%alphag_o(ii)
!
         arel_i=are_liq(ii)
         areg_i=are_gas(ii)
         ared_i=are_drp(ii)
!
         IF(agi_tmp.le.alpha_min) areg_i=0.0d0
         IF(ali_tmp.le.alpha_min) THEN
            arel_i=0.0d0
            ared_i=0.0d0
         ENDIF
         ecnvc_l_nf(i1)=arel_i
         ecnvc_g_nf(i1)=areg_i
         ecnvc_d_nf(i1)=ared_i
!
!........User defined energy convection of droplets
!   
         IF(i_droplet.ge.1) ecnvc_d_nf(i1)=0.d0
      ENDDO
!
!.....MARS interface
!
      nf_number=1
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i1=istart+i
         ii=left_nf(i1)
         ali_tmp=cell%alphal_o(ii)
         agi_tmp=cell%alphag_o(ii)
!
         arel_i=are_liq(ii)
         areg_i=are_gas(ii)
         ared_i=are_drp(ii)
!
         IF(agi_tmp.le.alpha_min) areg_i=0.0d0
         IF(ali_tmp.le.alpha_min) THEN
            arel_i=0.0d0
            ared_i=0.0d0
         ENDIF
!         
         idx=i3invtbl(i)
         IF(mcdirect(idx).lt.0)THEN
            arel_k=c3dpv(idx,2)
            ared_k=c3dpv(idx,7)
         ELSE
            arel_k=are_liq(ii)
            ared_k=are_drp(ii)
         ENDIF
         IF(mcgdirect(idx).lt.0)THEN
            areg_k=c3dpv(idx,1)
         ELSE
            areg_k=are_gas(ii)
         ENDIF
!
         IF    (flux_l_nf(i1).gt.0.0d0)THEN
            ecnvc_l_nf(i1)=arel_i
         ELSEIF(flux_l_nf(i1).lt.0.0d0)THEN
            ecnvc_l_nf(i1)=arel_k
         ELSE
            ecnvc_l_nf(i1)=(arel_i+arel_k)*0.5d0
         ENDIF
         IF    (flux_g_nf(i1).gt.0.0d0)THEN
            ecnvc_g_nf(i1)=areg_i
         ELSEIF(flux_g_nf(i1).lt.0.0d0)THEN
            ecnvc_g_nf(i1)=areg_k
         ELSE
            ecnvc_g_nf(i1)=(areg_i+areg_k)*0.5d0
         ENDIF
         IF    (flux_d_nf(i1).gt.0.0d0)THEN
            ecnvc_d_nf(i1)=ared_i
         ELSEIF(flux_d_nf(i1).lt.0.0d0)THEN
            ecnvc_d_nf(i1)=ared_k
         ELSE
            ecnvc_d_nf(i1)=(ared_i+ared_k)*0.5d0
         ENDIF
!
!........User defined energy convection of droplets
!   
         IF(i_droplet.ge.1) ecnvc_d_nf(i1)=0.d0
      ENDDO
!      
      RETURN
      END SUBROUTINE scalar_energy_convection
