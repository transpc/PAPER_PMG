!
      SUBROUTINE turb_ke_gas
!
!     This routine calculates turbulence liquid viscosity based 
!     on k-e model.
!
      USE Zinterface
      USE VOL_DATA , ONLY: cell
      USE Zzone    , ONLY: ncell_fluid
      USE Zcore    , ONLY: np
      USE Zparam   , ONLY: cmu
      USE Zconst1  , ONLY: restart,iturb
      USE Zface    , ONLY: Kepsilon_real
      USE Zndforce , ONLY: d_bfc                 
      USE Zturb    , ONLY: turb_keg,turb_dpg,turb_keg_o,turb_dpg_o,  &
                            pro_keg,tauwg,utaug,yplusg,veltg,strn_keg,ustar_keg,w_real_keg,cmug_real 
      USE Zvector  , ONLY: vg_o
     
!
      IMPLICIT NONE
!      
!.....Local variables
      INTEGER :: i
      LOGICAL, SAVE:: initial=.true.
!
!...Initialize turbulence kinetic energy and dissipation
!
      IF(initial.and.restart.eq.0) THEN
         CALL init_ke(turb_keg,turb_dpg,vg_o)
         initial=.false.
      ENDIF
!
!.....Shift solutions
!
      DO i=1,ncell_fluid
         turb_keg_o(i)=turb_keg(i)
         turb_dpg_o(i)=turb_dpg(i)
      ENDDO
!
!.....Calculate utau and yplus
!
      IF(iturb.ne.Kepsilon_real)THEN
         DO i=1,ncell_fluid
           utaug(i)=cmu**0.25d0*DSQRT(turb_keg(i))
           yplusg(i)=cell%rhog(i)*utaug(i)*d_bfc(i)/cell%lviscosg(i) !next
         ENDDO      
      ELSE
         DO i=1,ncell_fluid
           utaug(i)=cmug_real(i)**0.25d0*DSQRT(turb_keg(i))
           yplusg(i)=cell%rhog(i)*utaug(i)*d_bfc(i)/cell%lviscosg(i) !next
         ENDDO       
      ENDIF      
!
!.....Calculate tangential velocity and wall sheer stress
!
      CALL turb_ke_vt_tauw_gas(tauwg,veltg,vg_o,turb_keg_o,yplusg)
!
!.....Communicate velt
!
      IF(np.gt.1) CALL communicate_1d(veltg)         
!
!.....Caculate kinetic energy production
!
      IF(iturb.ne.Kepsilon_real)THEN
         CALL turb_ke_product(2,pro_keg,tauwg,strn_keg)
      ELSE
         CALL turb_ke_product_real(2,pro_keg,tauwg,strn_keg,ustar_keg,w_real_keg)
      ENDIF        
!
!.....Solve k-e equation
!      
      CALL turb_ke_calc_gas
!
      END SUBROUTINE turb_ke_gas
