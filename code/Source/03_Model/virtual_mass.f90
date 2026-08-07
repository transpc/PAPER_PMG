!
      SUBROUTINE virtual_mass(alphag_mi,alphal_mi,alphad_mi,rhog_mi,rhol_mi,rhod_mi, &
                               ar_gas_i,ar_liq_i,ar_drp_i,ag_i,al_i,ad_i,avm_gl_i,avm_gd_i,ncell)
!     
!     This calculates mixture, virtual mass and slip coefficient
!
      IMPLICIT NONE
!      
      REAL(8) alphag_mi(*),alphal_mi(*),alphad_mi(*)
      REAL(8) rhog_mi(*),rhol_mi(*),rhod_mi(*)
      REAL(8) ar_liq_i(*),ar_gas_i(*),ar_drp_i(*)
      REAL(8) ag_i(*),al_i(*),ad_i(*)
      REAL(8) avm_gl_i(*),avm_gd_i(*)
      INTEGER ncell
!      
      INTEGER i
      REAL(8) alpha
      REAL(8) cvm_gl_i,cvm_gd_i,c_vm
      REAL(8) rhovm_gl_i,rhovm_gd_i
!
      DO i=1,ncell
!
      alpha=(alphag_mi(i)+alphal_mi(i))
      IF(alpha.gt.0.d0)alpha=alphag_mi(i)/alpha
      cvm_gl_i=c_vm(alpha)
!
      alpha=alphag_mi(i)+alphad_mi(i)
      IF(alpha.gt.0.d0)alpha=alphag_mi(i)/alpha
      cvm_gd_i=c_vm(alpha)
!
      ar_gas_i(i)=alphag_mi(i)*rhog_mi(i)
      ar_liq_i(i)=alphal_mi(i)*rhol_mi(i)
      ar_drp_i(i)=alphad_mi(i)*rhod_mi(i)
!
      rhovm_gl_i=cvm_gl_i*(ar_liq_i(i)+ar_gas_i(i))        
      rhovm_gd_i=cvm_gd_i*(ar_drp_i(i)+ar_gas_i(i))        
      ag_i(i)=alphag_mi(i)*(rhog_mi(i)+alphal_mi(i)*rhovm_gl_i+alphad_mi(i)*rhovm_gd_i)
      al_i(i)=alphal_mi(i)*(rhol_mi(i)+alphag_mi(i)*rhovm_gl_i)
      ad_i(i)=alphad_mi(i)*(rhod_mi(i)+alphag_mi(i)*rhovm_gd_i)
!
      avm_gl_i(i)=alphag_mi(i)*alphal_mi(i)*rhovm_gl_i
      avm_gd_i(i)=alphag_mi(i)*alphad_mi(i)*rhovm_gd_i
!
      ENDDO
      CALL vectorize_scalar_upwind
!
      RETURN
      END SUBROUTINE virtual_mass
!      
!
      FUNCTION c_vm(alpha)
!
!     This function calculates virtual mass coefficient
!
      IMPLICIT NONE
!
      REAL(8) alpha,c_vm
!
      IF(alpha.lt.0.5d0)THEN
         c_vm=(0.5d0+alpha)/(1.d0-alpha)
      ELSE
         c_vm=(1.5d0-alpha)/alpha
      ENDIF
!
      RETURN 
      END FUNCTION c_vm      
   
