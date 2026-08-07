!
      SUBROUTINE initialize_time_variant_inlet
!
!     Time-Variant inlet velocity boundary
!
      USE Zbc_index    , ONLY: nvin
      USE Zboron       , ONLY: cboronb_liq,cboronb
      USE Zconst1      , ONLY: vv_prob
      USE Ztimecon     , ONLY: time
!
      IMPLICIT NONE
!
      INTEGER i
!      
!.....boron_trans
!
      IF (vv_prob.eq.'boron_trans') THEN
         IF(time.ge.5.d0)THEN
            DO i=1,nvin
               cboronb_liq(i)=cboronb(i)
            ENDDO
          ELSE
             DO i=1,nvin
                cboronb_liq(i)=0.d0
             ENDDO
          ENDIF       
      ENDIF    
!      
      RETURN
      END SUBROUTINE initialize_time_variant_inlet
