!
      SUBROUTINE scalar_mass_convection
!
!     This routine calculates convective mass fluxes through the cell face
!
      USE VOL_DATA
      USE Z2nd_order   , ONLY: mass_conv_2nd
      USE Zare         , ONLY: ar_gas,ar_liq,ar_drp
      USE Zb_condition , ONLY: alphab_gas,alphab_liq,alphab_drp,rhob_gas,rhob_liq
      USE c3com_cupid  , ONLY: i3invtbl,mcdirect,c3dpv,mcgdirect
      USE Znum_cell    , ONLY: istart_nf,istart_nbcon_nf
      USE Ztimecon     , ONLY: alpha_min
      USE Zmodel       , ONLY: i_droplet
      USE Zvec_index   , ONLY: left_nf,right_non,nbcon_nf
      USE Zvec_major   , ONLY: flux_l_nf,flux_g_nf,flux_d_nf,liq_conv_nf,vap_conv_nf,drp_conv_nf
      USE Zrv_model    , ONLY: rv_valve       
!
      IMPLICIT NONE
!      
      INCLUDE '../10_LinkToMARS/c3com.h' 
!
      INTEGER :: ii,kk
      INTEGER :: i,k,idx
      INTEGER :: nf_number,istart,len,istart2,i1,i2
!
      REAL(8) :: ar1_gas,ar2_gas
      REAL(8) :: ar1_liq,ar2_liq
      REAL(8) :: ar1_drp,ar2_drp
      REAL(8) :: vap_conv_up,liq_conv_up,drp_conv_up
!
!.....Computing cells
!
      nf_number=0
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)
!
!........Define owner cell values
!
         ar1_gas=ar_gas(ii)
         ar1_liq=ar_liq(ii)
         ar1_drp=ar_drp(ii)
!
!...........Assume zero cell values when the phase fraction is less than alpha_min
!
         IF(cell%alphag_o(ii).le.alpha_min) ar1_gas=0.d0
         IF(cell%alphal_o(ii).le.alpha_min) ar1_liq=0.0d0
         IF(cell%alphad_o(ii).le.alpha_min) ar1_drp=0.0d0
!
!........Define neighbor cell values
!
         ar2_gas=ar_gas(kk)
         ar2_liq=ar_liq(kk)
         ar2_drp=ar_drp(kk)
!
!........Assume zero cell values when the phase fraction is less than alpha_min
!
         IF(cell%alphag_o(kk).le.alpha_min) ar2_gas = 0.d0
         IF(cell%alphal_o(kk).le.alpha_min) ar2_liq = 0.0d0
         IF(cell%alphad_o(kk).le.alpha_min) ar2_drp=0.0d0
!
!........Apply 1st order upwind
!
         IF    (flux_g_nf(i1).gt.0.d0)THEN
            vap_conv_up=ar1_gas
         ELSEIF(flux_g_nf(i1).lt.0.d0)THEN
            vap_conv_up=ar2_gas
         ELSE
             vap_conv_up=(ar1_gas+ar2_gas)*0.5d0
         ENDIF
         IF    (flux_l_nf(i1).gt.0.d0)THEN
            liq_conv_up=ar1_liq
         ELSEIF(flux_l_nf(i1).lt.0.d0)THEN
            liq_conv_up=ar2_liq
         ELSE
            liq_conv_up=(ar1_liq+ar2_liq)*0.5d0
         ENDIF
         IF    (flux_d_nf(i1).gt.0.d0)THEN
            drp_conv_up=ar1_drp
         ELSEIF(flux_d_nf(i1).lt.0.d0)THEN
            drp_conv_up=ar2_drp
         ELSE
            drp_conv_up=(ar1_drp+ar2_drp)*0.5d0
         ENDIF
!         
         liq_conv_nf(i1)  =liq_conv_up
         vap_conv_nf(i1)  =vap_conv_up
         drp_conv_nf(i1)  =drp_conv_up
!
!........User defined mass convection of droplets
!         
         IF(i_droplet.ge.1) drp_conv_nf(i1)=0.d0
      ENDDO
!
!....valve model
!        
      IF(rv_valve.eq.1) CALL valve_model_scalar_mass_convection
!
!.....2nd order convection
!
      IF(mass_conv_2nd.gt.0) CALL mass_2nd_conv
!       
!.....Inlet
!
      nf_number=2
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      DO i=1,len  
         i1=istart+i
         i2=istart2+i
         ii=left_nf(i1)
         k=nbcon_nf(i2)
!
         ar1_gas=ar_gas(ii)
         ar1_liq=ar_liq(ii)
         ar1_drp=ar_drp(ii)
!
         IF(cell%alphag_o(ii).le.alpha_min) ar1_gas=0.d0
         IF(cell%alphal_o(ii).le.alpha_min) ar1_liq=0.0d0
         IF(cell%alphad_o(ii).le.alpha_min) ar1_drp=0.0d0
!
         ar2_gas=alphab_gas(k)*rhob_gas(k)
         ar2_liq=alphab_liq(k)*rhob_liq(k)
         ar2_drp=alphab_drp(k)*rhob_liq(k)
!
         IF    (flux_g_nf(i1).gt.0.d0)THEN
            vap_conv_up=ar1_gas
         ELSEIF(flux_g_nf(i1).lt.0.d0)THEN
            vap_conv_up=ar2_gas
         ELSE
            vap_conv_up=(ar1_gas+ar2_gas)*0.5d0
         ENDIF
         IF    (flux_l_nf(i1).gt.0.d0)THEN
            liq_conv_up=ar1_liq
         ELSEIF(flux_l_nf(i1).lt.0.d0)THEN
            liq_conv_up=ar2_liq
         ELSE
            liq_conv_up=(ar1_liq+ar2_liq)*0.5d0
         ENDIF
         IF    (flux_d_nf(i1).gt.0.d0)THEN
            drp_conv_up=ar1_drp
         ELSEIF(flux_d_nf(i1).lt.0.d0)THEN
            drp_conv_up=ar2_drp
         ELSE
            drp_conv_up=(ar1_drp+ar2_drp)*0.5d0
         ENDIF
!         
         liq_conv_nf(i1)  =liq_conv_up
         vap_conv_nf(i1)  =vap_conv_up
         drp_conv_nf(i1)  =drp_conv_up
!
         IF(i_droplet.ge.1) drp_conv_nf(i1)=0.d0
!
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
!
         ar1_gas=ar_gas(ii)
         ar1_liq=ar_liq(ii)
         ar1_drp=ar_drp(ii)
!
         IF(cell%alphag_o(ii).le.alpha_min) ar1_gas=0.d0
         IF(cell%alphal_o(ii).le.alpha_min) ar1_liq=0.0d0
         IF(cell%alphad_o(ii).le.alpha_min) ar1_drp=0.0d0
!
         ar2_gas=ar_gas(ii)
         ar2_liq=ar_liq(ii)
         ar2_drp=ar_drp(ii)
!
         IF(cell%alphag_o(ii).le.alpha_min) ar2_gas=0.d0
         IF(cell%alphal_o(ii).le.alpha_min) ar2_liq=0.0d0
         IF(cell%alphad_o(ii).le.alpha_min) ar2_drp=0.0d0
!
         IF    (flux_g_nf(i1).gt.0.d0)THEN
            vap_conv_up=ar1_gas
         ELSEIF(flux_g_nf(i1).lt.0.d0)THEN
            vap_conv_up=ar2_gas
         ELSE
             vap_conv_up=(ar1_gas+ar2_gas)*0.5d0
         ENDIF
         IF    (flux_l_nf(i1).gt.0.d0)THEN
            liq_conv_up=ar1_liq
         ELSEIF(flux_l_nf(i1).lt.0.d0)THEN
            liq_conv_up=ar2_liq
         ELSE
            liq_conv_up=(ar1_liq+ar2_liq)*0.5d0
         ENDIF
         IF    (flux_d_nf(i1).gt.0.d0)THEN
            drp_conv_up=ar1_drp
         ELSEIF(flux_d_nf(i1).lt.0.d0)THEN
            drp_conv_up=ar2_drp
         ELSE
            drp_conv_up=(ar1_drp+ar2_drp)*0.5d0
         ENDIF
!         
         liq_conv_nf(i1)  =liq_conv_up
         vap_conv_nf(i1)  =vap_conv_up
         drp_conv_nf(i1)  =drp_conv_up
!         
         IF(i_droplet.ge.1) drp_conv_nf(i1)=0.d0
!
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
!
         ar1_gas=ar_gas(ii)
         ar1_liq=ar_liq(ii)
         ar1_drp=ar_drp(ii)
!
         IF(cell%alphag_o(ii).le.alpha_min) ar1_gas=0.d0
         IF(cell%alphal_o(ii).le.alpha_min) ar1_liq=0.0d0
         IF(cell%alphad_o(ii).le.alpha_min) ar1_drp=0.0d0
!
         idx=i3invtbl(i)
         IF(mcdirect(idx).lt.0)THEN
            ar2_liq=c3dpv(idx,5)
            ar2_drp=c3dpv(idx,4)
         ELSE
            ar2_liq=ar_liq(ii)
            ar2_drp=ar_drp(ii)
         ENDIF
         IF(mcgdirect(idx).lt.0)THEN
            ar2_gas=c3dpv(idx,3)
         ELSE
            ar2_gas=ar_gas(ii)
         ENDIF
!
         IF    (flux_g_nf(i1).gt.0.d0)THEN
            vap_conv_up=ar1_gas
         ELSEIF(flux_g_nf(i1).lt.0.d0)THEN
            vap_conv_up=ar2_gas
         ELSE
            vap_conv_up=(ar1_gas+ar2_gas)*0.5d0
         ENDIF
         IF    (flux_l_nf(i1).gt.0.d0)THEN
            liq_conv_up=ar1_liq
         ELSEIF(flux_l_nf(i1).lt.0.d0)THEN
            liq_conv_up=ar2_liq
         ELSE
            liq_conv_up=(ar1_liq+ar2_liq)*0.5d0
         ENDIF
         IF    (flux_d_nf(i1).gt.0.d0)THEN
            drp_conv_up=ar1_drp
         ELSEIF(flux_d_nf(i1).lt.0.d0)THEN
            drp_conv_up=ar2_drp
         ELSE
            drp_conv_up=(ar1_drp+ar2_drp)*0.5d0
         ENDIF
!
         liq_conv_nf(i1)=liq_conv_up
         vap_conv_nf(i1)=vap_conv_up
         drp_conv_nf(i1)=drp_conv_up
!
         IF(i_droplet.ge.1) drp_conv_nf(i1)=0.d0
!
      ENDDO
!
      RETURN
      END SUBROUTINE scalar_mass_convection
