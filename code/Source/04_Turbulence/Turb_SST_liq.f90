!
      SUBROUTINE turb_SST_liq
!
!     This routine calculates turbulence liquid viscosity based 
!     on k-e model.
!
      USE Zinterface
      USE VOL_DATA , ONLY: cell
      USE Zzone    , ONLY: ncell_fluid,ncell_fluid_all
      USE Zcore    , ONLY: np
      USE Zparam   , ONLY: cmu,clog  
      USE Zconst1  , ONLY: restart
      USE Zndforce , ONLY: d_bfc,c_bface      
      USE Zturb    , ONLY: turb_ke,turb_dp,turb_ke_o,turb_dp_o,  &
                           pro_ke,utau,yplus,velt,strn_ke
      USE Zvector  , ONLY: vl_o
!
      IMPLICIT NONE
!      
      INTEGER i
!
      REAL(8),ALLOCATABLE::utau_all(:)
!      
      LOGICAL, SAVE:: initial
!      
      DATA initial/.true./
      ALLOCATE(utau_all(ncell_fluid_all))
      utau_all(:)=0.d0
!
!.....Initialize turbulence kinetic energy and dissipation
!
      IF(initial.and.restart.eq.0) THEN
         CALL init_SST(turb_ke,turb_dp,vl_o)
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
!.....Calculate tangential velocity and wall sheer stress,utau
!
      CALL turb_SST_vt_tauw_liq(velt,vl_o)      
!
!.....Calculate yplus (using updated utau)
!
      IF(np.gt.1)then
         CALL allgatherv_r(utau,utau_all,ncell_fluid,ncell_fluid_all,0)
      ELSE
         utau_all(:)=utau(:)
      ENDIF
!      
      DO i=1,ncell_fluid
!
!        utau(i)=cmu**0.25d0*DSQRT(turb_ke(i))                   !!! CYJ k-w 수정필요???
!        tauw(i)=cell%lviscosl(i)*DABS(vt(i))
         !yplus(i)=cell%rhol(i)*utau(i)*d_bfc(i)/cell%lviscosl(i) !next
        yplus(i)=cell%rhol(i)*utau_all(c_bface(i))*d_bfc(i)/cell%lviscosl(i) !next
      ENDDO  
!
!.....Communicate velt
!
      IF(np.gt.1) CALL communicate_1d(velt)      
!
!.....Caculate kinetic energy production
!
      CALL turb_SST_product(1,pro_ke,strn_ke)
      
!
!.....Solve k-e equation
!
      CALL turb_SST_calc_liq
!
      DEALLOCATE(utau_all)
!
      END SUBROUTINE turb_SST_liq
