!
      SUBROUTINE set_vapor_generation
!
!     This routine calculate vapor generation rate
!
      USE VOL_DATA              
      USE Zpress    , ONLY: p
      USE Zqvol     , ONLY: gamma,h_ig,h_il
      USE Ztimecon  , ONLY: alpha_min
      USE Zzone     , ONLY: ncell_fluid
      USE Zrv_model    , ONLY: free_model,rv_ht_i
!
      IMPLICIT NONE
!
      INTEGER i
!
      REAL(8) hi_gas,hi_liq,PsP,gg
!
!.....Get recent interfacial heat transfer coefficients
!
      IF(free_model)CALL int_htc
      CALL int_swap(2)
      IF(rv_ht_i.gt.0)CALL rv_int_ht
      CALL int_swap(22)
!
      DO i=1,ncell_fluid
         PsP=cell%pps(i)/p(i)
         gg=-(H_ig(i)*PsP*(cell%ts(i)-cell%tg(i))+H_il(i)*(cell%ts(i)-cell%tl(i)))
         IF(gg.ge.0.d0)THEN
            hi_gas=cell%hgsat(i)
            hi_liq=cell%hl(i)
         ELSE
            hi_gas=cell%hg(i)
            hi_liq=cell%hlsat(i)
         ENDIF
         gamma(i)=gg/(hi_gas-hi_liq)
         IF(cell%alphag(i).le.alpha_min.and.gamma(i).lt.0.d0) gamma(i)=0.d0
         IF(cell%alphag(i).ge.1.d0-1.0d0*alpha_min.and.gamma(i).gt.0.d0) gamma(i)=0.d0
      ENDDO
!
      RETURN
      END SUBROUTINE set_vapor_generation
