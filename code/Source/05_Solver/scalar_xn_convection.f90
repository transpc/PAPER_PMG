!
      SUBROUTINE scalar_xn_convection
!
!     This routine calculates convective mass fluxes through the cell face
!
      USE VOL_DATA
      USE Z2nd_order   , ONLY: qula_conv_2nd
      USE Zare         , ONLY: ar_gas
      USE Zb_condition , ONLY: alphab_gas,rhob_gas,qualab
      USE c3com_cupid  , ONLY: i3invtbl,c3dpv,mcgdirect
      USE Znum_cell    , ONLY: istart_nf,istart_nbcon_nf
      USE Ztimecon     , ONLY: alpha_min
      USE Zvec_index   , ONLY: left_nf,right_non,nbcon_nf
      USE Zvec_major   , ONLY: flux_g_nf,quala_conv_nf
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
      REAL(8) :: quala_conv_up
      REAL(8) :: ar1_qul,ar2_qul
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
         ar1_qul=ar_gas(ii)*cell%quala(ii)
!
!...........Assume zero cell values when the phase fraction is less than alpha_min
!
         IF(cell%alphag_o(ii).le.alpha_min) ar1_qul=0.d0
!
!........Define neighbor cell values
!
         ar2_qul=ar_gas(kk)*cell%quala(kk)
!
!........Assume zero cell values when the phase fraction is less than alpha_min
!
         IF(cell%alphag_o(kk).le.alpha_min) ar2_qul = 0.d0
!
!........Apply 1st order upwind
!
         IF    (flux_g_nf(i1).gt.0.d0)THEN
            quala_conv_up=ar1_qul
         ELSEIF(flux_g_nf(i1).lt.0.d0)THEN
            quala_conv_up=ar2_qul
         ELSE
            quala_conv_up=(ar1_qul+ar2_qul)*0.5d0
         ENDIF
!         
         quala_conv_nf(i1)=quala_conv_up
!
      ENDDO
!
!......valve model
!       
      IF(rv_valve.eq.1) CALL valve_model_scalar_xn_convection
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
         ar1_qul=ar_gas(ii)*cell%quala(ii)
!
         IF(cell%alphag_o(ii).le.alpha_min) ar1_qul=0.d0
!
         ar2_qul=alphab_gas(k)*rhob_gas(k)*qualab(k)
!
         IF    (flux_g_nf(i1).gt.0.d0)THEN
            quala_conv_up=ar1_qul
         ELSEIF(flux_g_nf(i1).lt.0.d0)THEN
            quala_conv_up=ar2_qul
         ELSE
             quala_conv_up=(ar1_qul+ar2_qul)*0.5d0
         ENDIF
!         
         quala_conv_nf(i1)=quala_conv_up
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
         ar1_qul=ar_gas(ii)*cell%quala(ii)
!
         IF(cell%alphag_o(ii).le.alpha_min) ar1_qul=0.d0
!
         ar2_qul=ar_gas(ii)*cell%quala(ii)
!
         IF(cell%alphag_o(ii).le.alpha_min) ar2_qul=0.d0
!
         IF    (flux_g_nf(i1).gt.0.d0)THEN
            quala_conv_up=ar1_qul
         ELSEIF(flux_g_nf(i1).lt.0.d0)THEN
            quala_conv_up=ar2_qul
         ELSE
            quala_conv_up=(ar1_qul+ar2_qul)*0.5d0
         ENDIF
!         
         quala_conv_nf(i1)=quala_conv_up
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
         ar1_qul=ar_gas(ii)*cell%quala(ii)
!
         IF(cell%alphag_o(ii).le.alpha_min) ar1_qul=0.d0
!
         idx=i3invtbl(i)
         IF(mcgdirect(idx).lt.0)THEN
            ar2_qul=c3dpv(idx,6)
         ELSE
            ar2_qul=ar_gas(ii)*cell%quala(ii)
         ENDIF
!
         IF    (flux_g_nf(i1).gt.0.d0)THEN
            quala_conv_up=ar1_qul
         ELSEIF(flux_g_nf(i1).lt.0.d0)THEN
            quala_conv_up=ar2_qul
         ELSE
             quala_conv_up=(ar1_qul+ar2_qul)*0.5d0
         ENDIF
!
         quala_conv_nf(i1)=quala_conv_up
!
      ENDDO
!
      IF(qula_conv_2nd.gt.0) CALL quality_2nd_conv
!
      RETURN
      END SUBROUTINE scalar_xn_convection
