!
      SUBROUTINE scalar_reset
!
!     This routine resets scalar values to that of previous time step
!
      USE VOL_DATA 
      USE SOLID_DATA, ONLY: solid                   
      USE Zzone     , ONLY: ncell_fluid,ncell_cond
      USE Zconst1   , ONLY: iat,iturb
      USE Ziat      , ONLY: ia,ia_old
      USE Zpress    , ONLY: p,pp
      USE Zturb     , ONLY: turb_ke,turb_ke_o,turb_dp,turb_dp_o
!
      USE Zbc_index , ONLY: npb
!
      USE Zvec_param    , ONLY: nf_flux
      USE Zvec_major    , ONLY: flux_l_nf,flux_g_nf,flux_d_nf,       &
                                flux_l_nf_o,flux_g_nf_o,flux_d_nf_o
!
      IMPLICIT NONE 
!
      INTEGER i,i1
!
      DO i=1,ncell_fluid
         IF(npb(i).eq.0) p(i)=p(i)-pp(i)
         cell%el(i)=cell%el_o(i)
         cell%eg(i)=cell%eg_o(i)
         cell%alphal(i)=cell%alphal_o(i)
         cell%alphad(i)=cell%alphad_o(i)
         cell%alphag(i)=cell%alphag_o(i)
         cell%quala(i)=cell%quala_o(i)
      ENDDO
!      
      DO i1=1,nf_flux
         flux_l_nf(i1)=flux_l_nf_o(i1)
         flux_g_nf(i1)=flux_g_nf_o(i1)
         flux_d_nf(i1)=flux_d_nf_o(i1)
      ENDDO
!
      CALL property_calc(1)
!
      CALL set_vapor_generation
!
      IF(iturb.ge.0)THEN
         DO i=1,ncell_fluid
            turb_ke(i)=turb_ke_o(i)
            turb_dp(i)=turb_dp_o(i)
         ENDDO         
      ENDIF
      IF(iat.gt.0)THEN       
         DO i=1,ncell_fluid
            ia(i)=ia_old(i)
         ENDDO
      ENDIF
!
      IF(ncell_cond.gt.0) solid%tsol(:)=solid%tsol_o(:)
!
      RETURN
      END SUBROUTINE scalar_reset
