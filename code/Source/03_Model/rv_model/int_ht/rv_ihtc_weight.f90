!
      SUBROUTINE rv_ihtc_weight(i)
!
      USE VOL_DATA        , ONLY: cell
      USE Zrv_weight      , ONLY: liquid,gas,                                    &
                                  bubbly,slug,churn,annular,mpr,                 &
                                  invann,invchn,invslg,mist,mpo,                 &
                                  wf_liquid,wf_gas,                              &
                                  wf_bubbly,wf_slug,wf_churn,wf_annular,wf_MPR,  &
                                  wf_invann,wf_invchn,wf_invslg,wf_mist,wf_MPO,  &
                                  wf_sum,                                        &
                                  pfnrgj
!
      IMPLICIT NONE
!
      REAL(8),PARAMETER :: alpha_am=0.9999d0        
!.....Local variables
      INTEGER :: i
      REAL(8) :: wf_sumr
      REAL(8) :: alphag
      REAL(8) :: regime_weight !function
!
!.....Initialize mode
!
      liquid =0
      gas    =0      
      bubbly =0
      slug   =0
      churn  =0
      annular=0
      mpr    =0
      invann =0
      invchn =0
      invslg =0
      mist   =0
      mpo    =0
!
!.....Initialize weights
!
      wf_liquid =0.d0
      wf_gas    =0.d0
      wf_bubbly =0.d0
      wf_slug   =0.d0
      wf_churn  =0.d0
      wf_annular=0.d0
      wf_MPR    =0.d0
      wf_invann =0.d0
      wf_invchn =0.d0
      wf_invslg =0.d0
      wf_mist   =0.d0
      wf_MPO    =0.d0
!
!.....Initialize gas fraction
!
      alphag=MIN(1.d0,MAX(cell%alphag(i),0.d0))
!      pfnrgj=z_trans(i)
!   
!.....Weighting factor
!
      wf_liquid =regime_weight(alphag,-1.d0           ,0.d0)                 !single-phase liquid
      wf_gas    =regime_weight(alphag,1.d0            ,2.d0)                 !single-phase gas
      wf_bubbly =regime_weight(alphag,0.d0            ,cell%alpha_bs(i))     !pre-CHF
      wf_slug   =regime_weight(alphag,cell%alpha_bs(i),cell%alpha_de(i))     !pre-CHF
      wf_churn  =regime_weight(alphag,cell%alpha_de(i),cell%alpha_sa(i))     !pre-CHF
      wf_annular=regime_weight(alphag,cell%alpha_sa(i),alpha_am)             !pre-CHF
      wf_MPR    =regime_weight(alphag,alpha_am        ,1.d0)                 !pre-CHF
      wf_invann =regime_weight(alphag,0.d0            ,cell%alpha_bs(i))     !post-dry
      wf_invchn =regime_weight(alphag,cell%alpha_bs(i),cell%alpha_cd(i))     !post-dry
      wf_invslg =regime_weight(alphag,cell%alpha_cd(i),cell%alpha_sa(i))     !post-dry
      wf_mist   =regime_weight(alphag,cell%alpha_sa(i),alpha_am)             !post-dry
      wf_MPO    =regime_weight(alphag,alpha_am        ,1.d0)                 !post-dry
!      
!.....Weighting factor control for pre-CHF or pos-dryout
!
      pfnrgj=cell%wf_dry(i)
      IF(pfnrgj.eq.0.d0)THEN      !Pre-CHF
         wf_invann=0.d0
         wf_invchn=0.d0
         wf_invslg=0.d0
         wf_mist  =0.d0
         wf_MPO   =0.d0 
      ELSEIF(pfnrgj.eq.1.d0)THEN !Post-DRY
         wf_bubbly =0.d0
         wf_slug   =0.d0
         wf_churn  =0.d0
         wf_annular=0.d0
         wf_MPR    =0.d0    
      ENDIF  
!      
!.....Weighting factor control for inter-regime smoothing (PIK)
!
      IF(wf_churn.gt.0.d0.and.wf_churn.lt.0.5d0)THEN
         IF(wf_slug.gt.0.5d0.and.wf_slug.lt.1.d0)THEN
            wf_churn=0.d0
            wf_slug =1.d0
         ELSEIF(wf_annular.gt.0.5d0.and.wf_annular.lt.1.d0)THEN
            wf_churn  =0.d0
            wf_annular=1.d0
         ENDIF
      ENDIF     
      IF(wf_churn.ge.0.5d0.and.wf_churn.lt.1.d0)THEN
         IF(wf_slug.gt.0.d0.and.wf_slug.le.0.5d0)THEN
            wf_churn=(wf_churn-0.5d0)*2.d0
            wf_slug =wf_slug*2.d0
         ELSEIF(wf_annular.ge.0.d0.and.wf_annular.le.0.5d0)THEN
            wf_churn  =(wf_churn-0.5d0)*2.d0
            wf_annular=wf_annular*2.d0
         ENDIF
      ENDIF                  
      IF(wf_invchn.gt.0.d0.and.wf_invchn.lt.0.5d0)THEN
         IF(wf_invslg.gt.0.5d0.and.wf_invslg.lt.1.d0)THEN
            wf_invchn=0.d0
            wf_invslg=1.d0
         ELSEIF(wf_invann.gt.0.5d0.and.wf_invann.lt.1.d0)THEN
            wf_invchn=0.d0
            wf_invann=1.d0
         ENDIF
      ENDIF 
      IF(wf_invchn.ge.0.5d0.and.wf_invchn.lt.1.d0)THEN
         IF(wf_invslg.gt.0.d0.and.wf_invslg.le.0.5d0)THEN
            wf_invchn=(wf_invchn-0.5d0)*2.d0
            wf_invslg=wf_invslg*2.d0
         ELSEIF(wf_invann.gt.0.d0.and.wf_invann.le.0.5d0)THEN
            wf_invchn=(wf_invchn-0.5d0)*2.d0
            wf_invann=wf_invann*2.d0
         ENDIF
      ENDIF
!
!.....Weighting factor sum
!
      wf_sum= (wf_liquid+wf_bubbly+wf_slug  +wf_churn +wf_annular    +wf_MPR+wf_gas)*(1.d0-pfnrgj) &
             +(wf_liquid+wf_invann+wf_invchn+wf_invslg+wf_mist+wf_MPO       +wf_gas)*pfnrgj
!
!.....Weighting factor update
!
      wf_sumr=1.d0/wf_sum
      wf_liquid =wf_liquid *wf_sumr
      wf_gas    =wf_gas    *wf_sumr
      wf_bubbly =wf_bubbly *wf_sumr
      wf_slug   =wf_slug   *wf_sumr
      wf_churn  =wf_churn  *wf_sumr
      wf_annular=wf_annular*wf_sumr
      wf_MPR    =wf_MPR    *wf_sumr
      wf_invann =wf_invann *wf_sumr
      wf_invchn =wf_invchn *wf_sumr
      wf_invslg =wf_invslg *wf_sumr
      wf_mist   =wf_mist   *wf_sumr
      wf_MPO    =wf_MPO    *wf_sumr
! 
!.....Flow regime mode determination
!
      IF(wf_liquid .gt.0.d0.and.wf_liquid  .le.1.d0) liquid =1
      IF(wf_gas    .gt.0.d0.and.wf_gas     .le.1.d0) gas    =1
      IF(wf_bubbly .gt.0.d0.and.wf_bubbly  .le.1.d0) bubbly =1
      IF(wf_slug   .gt.0.d0.and.wf_slug    .le.1.d0) slug   =1
      IF(wf_churn  .gt.0.d0.and.wf_churn   .le.1.d0) churn  =1
      IF(wf_annular.gt.0.d0.and.wf_annular .le.1.d0) annular=1
      IF(wf_MPR    .gt.0.d0.and.wf_MPR     .le.1.d0) MPR    =1
      IF(wf_invann .gt.0.d0.and.wf_invann  .le.1.d0) invann =1
      IF(wf_invchn .gt.0.d0.and.wf_invchn  .le.1.d0) invchn =1
      IF(wf_invslg .gt.0.d0.and.wf_invslg  .le.1.d0) invslg =1
      IF(wf_mist   .gt.0.d0.and.wf_mist    .le.1.d0) mist   =1
      IF(wf_MPO    .gt.0.d0.and.wf_mist    .le.1.d0) MPO    =1
!      
      END SUBROUTINE rv_ihtc_weight
!
!------------------------------------------------------------------------
!     FUNCTION: regime_weight
!------------------------------------------------------------------------
      FUNCTION regime_weight(alphag,alp_left,alp_right) RESULT(wff)
!      
      IMPLICIT NONE
!      
      REAL(8) :: alphag,alp_left,alp_right
      REAL(8) :: wl,wr,wff
      REAL(8),PARAMETER :: vtran_width=100.d0   !< transition width
!
      wl=0.5d0+0.5d0*(alphag-alp_left)*vtran_width
      wl=MAX(0.d0,MIN(wl,1.d0))
      wr=0.5d0+0.5d0*(alphag-alp_right)*vtran_width
      wr=MAX(0.d0,MIN(wr,1.d0))
      wff=wl-wr
!      
      END FUNCTION regime_weight      

