!
      SUBROUTINE int_htc
!
!     This routine caclulates heat transfer coefficients at the interface
!
      USE Vol_DATA, only: cell
      USE Zconst1        , ONLY: mhtc
      USE Zdel_scalar    , ONLY: limit_iht_opt,ag_min_hig,al_min_hil,qu_min_hgf
      USE Zmpi           , ONLY: ncell_fp
      USE Zqvol          , ONLY: h_ig,h_il,h_gf
      USE Ztimecon       , ONLY: alpha_min
      USE Zuserdefined   , ONLY: udfl_calc_HTC_int_i
      USE Zzone          , ONLY: ncell_fluid
!
      IMPLICIT NONE
!
      INTEGER i
!   
      REAL(8),ALLOCATABLE::ag(:),al(:),agi(:),ali(:)
!
      ALLOCATE(ag(ncell_fp),al(ncell_fp),agi(ncell_fp),ali(ncell_fp))
!
!.....When 'rv heat transfer model' is off
!
!
!.....Predict void fraction when smac=2 or 4 (ice-ice 2step or smac-ice 2step)
!
      DO i=1,ncell_fluid
         agi(i)=cell%alphag(i)
         ali(i)=cell%alphal(i)
      ENDDO
!     
!.....Use of user-defined model 
!
      IF(udfl_calc_HTC_int_i)THEN
         CALL udfn_calc_HTC_int_i(agi,ali)      
      ELSE   
!     
!.....Use of simple model (assumption of linear equation)
!
         IF(mHTC.eq.0)THEN
            CALL int_htc_simple_model(agi,ali)
!     
!........Use of Ranz&Marshall
!         
         ELSEIF(mHTC.eq.1)THEN
            CALL int_htc_simple_topology(agi)
!     
!........Use of Ranz&Marshall&STAR
!         
         ELSE  
            CALL int_htc_full_topology(agi)
         ENDIF
!         
      ENDIF   
!
      DO i=1,ncell_fluid
         IF(cell%alphal(i).lt.alpha_min) h_gf(i)=0.0d0
      ENDDO                      
!
      IF(limit_iht_opt.gt.0)THEN
         DO i=1,ncell_fluid
!
            IF(cell%alphag(i).lt.ag_min_hig)THEN
                h_ig(i)=0.0d0
                h_gf(i)=0.0d0
            ELSE
               IF(cell%quala(i).lt.qu_min_hgf) h_gf(i)=0.0d0
            ENDIF
            IF(cell%alphal(i).lt.al_min_hil) THEN
               h_il(i)=0.0d0 
               h_gf(i)=0.0d0
            ENDIF
            IF(cell%quala(i).gt.1.0d0-qu_min_hgf) THEN
               h_ig(i)=0.0d0 
            ENDIF 
!
         ENDDO
      ENDIF
!
      DEALLOCATE(ag,al,agi,ali)
!
      RETURN
      END SUBROUTINE int_htc
