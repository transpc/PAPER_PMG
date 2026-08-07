!
      SUBROUTINE turb_ke_liq
!
!     This routine calculates turbulence liquid viscosity based 
!     on k-e model.
!
      USE Zinterface
      USE VOL_DATA , ONLY: cell
      USE Zzone    , ONLY: ncell_fluid
      USE Zcore    , ONLY: np           
      USE Zparam   , ONLY: cmu      
      USE Zndforce , ONLY: d_bfc      
      USE Zturb    , ONLY: turb_ke,turb_dp,turb_ke_o,turb_dp_o,                            &
                           pro_ke,tauw,utau,yplus,velt,strn_ke,ustar_ke,w_real_ke,cmu_real  
      USE Zvector  , ONLY: vl_o
      USE Zconst1  , ONLY: restart,iturb
      USE Zface    , ONLY: Kepsilon_real
!
      IMPLICIT NONE
!      
!.....Local variables
      INTEGER :: i
      LOGICAL,SAVE :: initial=.true.
!
!.....Initialize turbulence kinetic energy and dissipation
!
      IF(initial.and.restart.eq.0) THEN
         CALL init_ke(turb_ke,turb_dp,vl_o)
         initial=.false.
      ENDIF
!
!.....Shift solutions
!
      DO i=1,ncell_fluid
         turb_ke_o(i)=turb_ke(i)
         turb_dp_o(i)=turb_dp(i)
      ENDDO
!
!.....Calculate utau and yplus
!  
      IF(iturb.ne.Kepsilon_real)THEN
         DO i=1,ncell_fluid
           utau(i)=cmu**0.25d0*DSQRT(turb_ke(i))
           yplus(i)=cell%rhol(i)*utau(i)*d_bfc(i)/cell%lviscosl(i) !next
         ENDDO      
      ELSE
         DO i=1,ncell_fluid
           utau(i)=cmu_real(i)**0.25d0*DSQRT(turb_ke(i))
           yplus(i)=cell%rhol(i)*utau(i)*d_bfc(i)/cell%lviscosl(i) !next
         ENDDO       
      ENDIF
!
!.....Calculate tangential velocity and wall sheer stress
!
      CALL turb_ke_vt_tauw_liq(tauw,velt,vl_o,turb_ke_o,yplus)
!
!.....Communicate velt
!
      IF(np.gt.1) CALL communicate_1d(velt)      
!
!.....Caculate kinetic energy production
!
      IF(iturb.ne.Kepsilon_real)THEN
         CALL turb_ke_product(1,pro_ke,tauw,strn_ke)
      ELSE
         CALL turb_ke_product_real(1,pro_ke,tauw,strn_ke,ustar_ke,w_real_ke)      
      ENDIF         
!
!.....Solve k-e equation
!
      CALL turb_ke_calc_liq
!
      END SUBROUTINE turb_ke_liq
