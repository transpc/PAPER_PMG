      SUBROUTINE udfn_erg_diff
!
!     User-defined energy diffusion terms
!       
      USE VOL_DATA                 
      USE Zconst1      , ONLY: vv_prob
      USE Zenergy_diff , ONLY: ediff_liq,ediff_gas
      USE Zzone        , ONLY: ncell_fluid
!
      IMPLICIT NONE
!      
      INTEGER i
!
      IF(vv_prob.eq.'PAFS-POOL') THEN
         DO i=1,ncell_fluid
            IF(cell%regime(i).eq.11) ediff_gas(i)=0.0d0
            IF(cell%regime(i).eq.13) ediff_liq(i)=0.0d0
         ENDDO
      ENDIF
!
      IF(vv_prob.eq.'fluidic_device') THEN
         DO i=1,ncell_fluid
            ediff_gas(i)=0.0d0
            ediff_liq(i)=0.0d0
         ENDDO
      ENDIF
!
      IF(vv_prob.eq.'OPR1000_fullvessel_1x1' .or. &
         vv_prob.eq.'opr1000_rv'             .or. &
         vv_prob.eq.'apr1400_rv'                   )THEN 
         DO i=1,ncell_fluid
            ediff_gas(i)=0.0d0
            ediff_liq(i)=0.0d0
         ENDDO
      ENDIF
!
      IF(vv_prob.eq.'APR1400_fullcore'                     .or. &
         vv_prob.eq.'APR1400_fullcore_modmesh01'           .or. &
         vv_prob.eq.'APR1400_fullcore_modmesh02_rv'        .or. &
         vv_prob.eq.'OPR1000_fullcore_modmesh02_rv'        .or. &
         vv_prob.eq.'OPR1000_fullcore_modmesh02_rv_vessel' .or. &
         vv_prob.eq.'OPR1000_single_assem'                       )THEN
         DO i=1,ncell_fluid
            ediff_gas(i)=0.0d0
            ediff_liq(i)=0.0d0
         ENDDO
      ENDIF
!
      IF(vv_prob.eq.'opr1000_rv_lbloca'                   )THEN 
         DO i=1,ncell_fluid
            ediff_gas(i)=0.0d0
            ediff_liq(i)=0.0d0
         ENDDO
      ENDIF 
      
      IF(vv_prob.eq.'ST2-CT-01'.or.vv_prob.eq.'ST2-CT-02'.or. &
         vv_prob.eq.'ST2-CT-03')THEN
         DO i=1,ncell_fluid
            ediff_gas(i)=0.0d0
            ediff_liq(i)=0.0d0
         ENDDO
      ENDIF
!
      RETURN
      END SUBROUTINE udfn_erg_diff
