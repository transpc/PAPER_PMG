!
      SUBROUTINE calc_scalar(ncvgd,pcnv)
!
!     This routine discretizes and solves the linearized  energy and mass 
!     equations explicitly or implicitly(smac=3 & implicit option=on).
!     imp_mom_conv(vg,vl,vd),imp_mom_diff(vg,vl,vd),imp_alpha(ag,al,ad),
!     imp_scalar_conv(eg,el,Xn),imp_scalar_diff(eg,el)
!
      USE Zinterface
      USE VOL_DATA      , ONLY: cell
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: np,myrank
      USE Zbicg        , ONLY: pbcgind,pbcgind_max            
      USE Zare         , ONLY: ar_gas,ar_liq,ar_drp,are_gas,are_liq,are_drp
      USE Zconst1      , ONLY: mboron
      USE Zconst2      , ONLY: dt,dt_old,dtr
      USE Zimplicit    , ONLY: imp_alpha
      USE Zio_unit     , ONLY: unit_log
      USE Zncg         , ONLY: ncg_diff
      USE Ztimecon     , ONLY: smac,iso_thermal,repeat_smac,time
      USE Zrv_model    , ONLY: free_model,rv_ht_i
!
      IMPLICIT NONE
!     
!.....Local variables
      INTEGER :: i
      LOGICAL :: ncvgd(3),pcnv,exp1,exp2
!
      DO i=1,ncell_fluid
         are_liq(i)=ar_liq(i)*cell%el(i)
         are_gas(i)=ar_gas(i)*cell%eg(i)
         are_drp(i)=ar_drp(i)*cell%el(i)
      ENDDO
!
!.....Communicate are_liq,are_gas for energy diffusion and convection
!
      IF(np.gt.1) CALL communicate_1d(are_liq,        &
                                      are_gas,        &
                                      cell%condl,     &
                                      cell%condg,     &
                                      cell%rhog,      &
                                      cell%rhol,      &
                                      cell%eg,cell%el)
!
      IF(smac.eq.3) CALL pressure_solve(0,0,pcnv)
!
!.....Return if pressure correction has not been converged
!
      IF(.not.pcnv) RETURN
!
      exp1=.true.
      exp2=.true.
!
      IF(smac.eq.3) THEN
         !!!IF(imp_scalar_diff.gt.0) exp1=.false.
         !!!IF(imp_scalar_conv.gt.0) exp2=.false.
      ENDIF
!
!.....Explicit energy diffusion
!
      IF(exp1) CALL scalar_energy_diffusion
!
!.....Energy convection
!
      IF(exp2) CALL scalar_energy_convection
!
!.....Non-condensable gas convection
! 
      IF(exp2) CALL scalar_xn_convection
!
!.....Pressure work convection
! 
      CALL scalar_work_convection
!
!.....Mass convection
!
      CALL scalar_mass_convection
!
!.....Implement scalar_mass_diffusion 2015.07.29 JHLee (SNU)
!
      IF(ncg_diff.gt.0) then
         IF(np.gt.1) CALL communicate_1d(cell%mdiff,     &
                                         cell%ha,cell%hg)
         CALL scalar_mass_diffusion
      ENDIF
!      
!.....Update IHTC for scalar calculation
!
      IF(smac.eq.3) THEN
         IF(imp_alpha.eq.0)THEN
            CALL continuity_for_smac3
         ELSE
            CALL continuity_for_smac3_imp
         ENDIF
         IF(free_model)CALL int_htc
         CALL int_swap(2)
         IF(rv_ht_i.gt.0)CALL rv_int_ht
         CALL int_swap(22)
      ENDIF
!
!.....Set scalar matrix
!
      IF(smac.ne.3) CALL scalar_matrix
!
!.....Calculate cell pressure
!
      IF(smac.eq.1) CALL pressure_solve(0,0,pcnv)
      IF(smac.eq.0) CALL pressure_solve(1,1,pcnv)
!
      IF(pbcgind.gt.0)THEN
         IF(myrank.eq.0)WRITE(*,*)'Iteration fails in fluid!!!'
      ENDIF
! 
      IF(pbcgind.gt.0.and.pbcgind.le.pbcgind_max)RETURN
!
!.....Update material properties
!
      IF(smac.eq.3) THEN
         IF(iso_thermal.eq.0)THEN
            CALL scalar_update_smac3
            IF(repeat_smac)THEN
               CALL scalar_reset
               time=time-dt
               dt=dt_old*0.5d0
               dtr=1.0d0/dt
               time=time+dt
               dt_old=dt
               IF(myrank.eq.0) WRITE(*,*) '#### dt is reduced by half ####'
               IF(myrank.eq.0) WRITE(unit_log,*)'#### dt is reduced by half ####'
               IF(dt.lt.1.0d-15) STOP '### dt is less than 1.0d-15 ###'
               RETURN
            ENDIF
         ELSE
            CALL property_calc(1)
         ENDIF
      ELSE
         CALL scalar_update(ncvgd)
      ENDIF
!
!.....Calculate boron transport
!   
      IF(mboron.eq.1) CALL boron_transport
!
      END SUBROUTINE calc_scalar
