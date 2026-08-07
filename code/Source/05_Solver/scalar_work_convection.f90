!
      SUBROUTINE scalar_work_convection
!
!     This routine calculates convective mass fluxes through the cell face
!
      USE VOL_DATA
      USE Zb_condition , ONLY: alphab_gas,alphab_liq,alphab_drp
      USE c3com_cupid  , ONLY: i3invtbl,mcdirect,mcgdirect
      USE Ztimecon     , ONLY: alpha_min
      USE Znum_cell    , ONLY: istart_nf,istart_nbcon_nf
      USE Zvec_index   , ONLY: left_nf,right_non,nbcon_nf
      USE Zvec_major   , ONLY: flux_l_nf,flux_g_nf,flux_d_nf,al_conv_nf,ad_conv_nf,void_conv_nf
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
      REAL(8) :: a1_gas,a2_gas
      REAL(8) :: a1_liq,a2_liq
      REAL(8) :: a1_drp,a2_drp
      REAL(8) :: void_conv_up,al_conv_up,ad_conv_up
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
         a1_gas=cell%alphag(ii)
         a1_liq=cell%alphal(ii)
         a1_drp=cell%alphad(ii)
!
!...........Assume zero cell values when the phase fraction is less than alpha_min
!
         IF(cell%alphag_o(ii).le.alpha_min) a1_gas=0.d0
         IF(cell%alphal_o(ii).le.alpha_min) a1_liq=0.d0
         IF(cell%alphad_o(ii).le.alpha_min) a1_drp=0.d0
!
!........Define neighbor cell values
!
         a2_gas=cell%alphag(kk)
         a2_liq=cell%alphal(kk)
         a2_drp=cell%alphad(kk)
!
!........Assume zero cell values when the phase fraction is less than alpha_min
!
         IF(cell%alphag_o(kk).le.alpha_min) a2_gas = 0.d0
         IF(cell%alphal_o(kk).le.alpha_min) a2_liq = 0.d0
         IF(cell%alphad_o(kk).le.alpha_min) a2_drp=0.d0
!
!........Apply 1st order upwind
!
         IF    (flux_g_nf(i1).gt.0.d0)THEN
            void_conv_up=a1_gas
         ELSEIF(flux_g_nf(i1).lt.0.d0)THEN
            void_conv_up=a2_gas
         ELSE
            void_conv_up=(a1_gas+a2_gas)*0.5d0
         ENDIF
         IF    (flux_l_nf(i1).gt.0.d0)THEN
            al_conv_up=a1_liq
         ELSEIF(flux_l_nf(i1).lt.0.d0)THEN
            al_conv_up=a2_liq
         ELSE
            al_conv_up=(a1_liq+a2_liq)*0.5d0
         ENDIF
         IF    (flux_d_nf(i1).gt.0.d0)THEN
            ad_conv_up=a1_drp
         ELSEIF(flux_d_nf(i1).lt.0.d0)THEN
            ad_conv_up=a2_drp
         ELSE
            ad_conv_up=(a1_drp+a2_drp)*0.5d0
         ENDIF
!         
         void_conv_nf(i1) =void_conv_up
         al_conv_nf(i1)   =al_conv_up
         ad_conv_nf(i1)   =ad_conv_up
!
      ENDDO
!
!.....valve model
!       
       IF(rv_valve.eq.1) CALL valve_model_scalar_work_convection
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
         a1_gas=cell%alphag(ii)
         a1_liq=cell%alphal(ii)
         a1_drp=cell%alphad(ii)
!
         IF(cell%alphag_o(ii).le.alpha_min) a1_gas=0.d0
         IF(cell%alphal_o(ii).le.alpha_min) a1_liq=0.d0
         IF(cell%alphad_o(ii).le.alpha_min) a1_drp=0.d0
!
         a2_gas=alphab_gas(k)
         a2_liq=alphab_liq(k)
         a2_drp=alphab_drp(k) 
!
         IF    (flux_g_nf(i1).gt.0.d0)THEN
            void_conv_up=a1_gas
         ELSEIF(flux_g_nf(i1).lt.0.d0)THEN
            void_conv_up=a2_gas
         ELSE
             void_conv_up=(a1_gas+a2_gas)*0.5d0
         ENDIF
         IF    (flux_l_nf(i1).gt.0.d0)THEN
            al_conv_up=a1_liq
         ELSEIF(flux_l_nf(i1).lt.0.d0)THEN
            al_conv_up=a2_liq
         ELSE
            al_conv_up=(a1_liq+a2_liq)*0.5d0
         ENDIF
         IF    (flux_d_nf(i1).gt.0.d0)THEN
            ad_conv_up=a1_drp
         ELSEIF(flux_d_nf(i1).lt.0.d0)THEN
            ad_conv_up=a2_drp
         ELSE
            ad_conv_up=(a1_drp+a2_drp)*0.5d0
         ENDIF
!         
         void_conv_nf(i1) =void_conv_up
         al_conv_nf(i1)   =al_conv_up
         ad_conv_nf(i1)   =ad_conv_up
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
         a1_gas=cell%alphag(ii)
         a1_liq=cell%alphal(ii)
         a1_drp=cell%alphad(ii)
!
         IF(cell%alphag_o(ii).le.alpha_min) a1_gas=0.d0
         IF(cell%alphal_o(ii).le.alpha_min) a1_liq=0.d0
         IF(cell%alphad_o(ii).le.alpha_min) a1_drp=0.d0
!
         a2_gas=cell%alphag(ii)
         a2_liq=cell%alphal(ii)
         a2_drp=cell%alphad(ii)
!
         IF(cell%alphag_o(ii).le.alpha_min) a2_gas=0.d0
         IF(cell%alphal_o(ii).le.alpha_min) a2_liq=0.d0
         IF(cell%alphad_o(ii).le.alpha_min) a2_drp=0.d0
!
         IF    (flux_g_nf(i1).gt.0.d0)THEN
            void_conv_up=a1_gas
         ELSEIF(flux_g_nf(i1).lt.0.d0)THEN
            void_conv_up=a2_gas
         ELSE
            void_conv_up=(a1_gas+a2_gas)*0.5d0
         ENDIF
         IF    (flux_l_nf(i1).gt.0.d0)THEN
            al_conv_up=a1_liq
         ELSEIF(flux_l_nf(i1).lt.0.d0)THEN
            al_conv_up=a2_liq
         ELSE
            al_conv_up=(a1_liq+a2_liq)*0.5d0
         ENDIF
         IF    (flux_d_nf(i1).gt.0.d0)THEN
            ad_conv_up=a1_drp
         ELSEIF(flux_d_nf(i1).lt.0.d0)THEN
            ad_conv_up=a2_drp
         ELSE
            ad_conv_up=(a1_drp+a2_drp)*0.5d0
         ENDIF
!         
         void_conv_nf(i1)=void_conv_up
         al_conv_nf(i1)=al_conv_up
         ad_conv_nf(i1)=ad_conv_up
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
         a1_gas=cell%alphag(ii)
         a1_liq=cell%alphal(ii)
         a1_drp=cell%alphad(ii)
!
         IF(cell%alphag_o(ii).le.alpha_min) a1_gas=0.d0
         IF(cell%alphal_o(ii).le.alpha_min) a1_liq=0.d0
         IF(cell%alphad_o(ii).le.alpha_min) a1_drp=0.d0
!
         idx=i3invtbl(i)
         IF(mcdirect(idx).lt.0)THEN
            a2_liq=c3rtp(1,idx,1)
            a2_drp=0.0d0
         ELSE
            a2_liq=cell%alphal(ii)
            a2_drp=cell%alphad(ii)
         ENDIF
         IF(mcgdirect(idx).lt.0)THEN
            a2_gas=c3rtp(1,idx,2)
         ELSE
            a2_gas=cell%alphag(ii)
         ENDIF
!
         IF    (flux_g_nf(i1).gt.0.d0)THEN
            void_conv_up=a1_gas
         ELSEIF(flux_g_nf(i1).lt.0.d0)THEN
            void_conv_up=a2_gas
         ELSE
            void_conv_up=(a1_gas+a2_gas)*0.5d0
         ENDIF
         IF    (flux_l_nf(i1).gt.0.d0)THEN
            al_conv_up=a1_liq
         ELSEIF(flux_l_nf(i1).lt.0.d0)THEN
            al_conv_up=a2_liq
         ELSE
            al_conv_up=(a1_liq+a2_liq)*0.5d0
         ENDIF
         IF    (flux_d_nf(i1).gt.0.d0)THEN
            ad_conv_up=a1_drp
         ELSEIF(flux_d_nf(i1).lt.0.d0)THEN
            ad_conv_up=a2_drp
         ELSE
            ad_conv_up=(a1_drp+a2_drp)*0.5d0
         ENDIF
!
         void_conv_nf(i1)=void_conv_up
         al_conv_nf(i1)=al_conv_up
         ad_conv_nf(i1)=ad_conv_up
!
      ENDDO
!
      RETURN
      END SUBROUTINE scalar_work_convection
