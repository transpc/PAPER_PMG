!
      SUBROUTINE turb_SST_vis_gas
!
!     This routine calculates liquid and gas turbulent viscosity 
!     using the solutions of k-e equations
!
      USE VOL_DATA     , ONLY: cell
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: cmu,ke_small,cappa,clog
      USE Zbc_index    , ONLY: icell_type
      USE Zconst1      , ONLY: turb_phase
      USE Zface        , ONLY: gas_only
      USE Zturb        , ONLY: turb_keg,turb_dpg,wvis_gas,f_b2,strn_keg
!
      IMPLICIT NONE
!
!     local variables
      INTEGER :: i
!      
      REAL(8) :: tv_max,tv_mult_max
      REAL(8) :: alpha_l
!
      DATA tv_mult_max/2.0d3/
!
!.....Additional constants - CYJ k-w
!
      alpha_l=5.0d0/9.0d0
!
!.....Effective viscosity and conductivity
!
      DO i=1,ncell_fluid
!
!........Turbulence viscosity
!
         IF(turb_dpg(i).gt.ke_small) THEN
!
!...........Calculation of invariant measure of the strain rate for k-w model
!
            cell%tviscosg(i)=cell%rhog(i)*alpha_l*turb_keg(i)/DMAX1(alpha_l*turb_dpg(i),f_b2(i)*strn_keg(i))              !CYJ k-w
         ELSE
            cell%tviscosg(i)=0.0d0
         ENDIF
!
!........Control the maximum turbulent viscosity by tv_max
!
         tv_max=cell%lviscosg(i)*tv_mult_max
         cell%tviscosg(i)=DMIN1(cell%tviscosg(i),tv_max)
!         
!........Effective viscosity
!         
         cell%eviscosg(i)=cell%lviscosg(i)+cell%tviscosg(i)
         IF(icell_type(i).eq.1)wvis_gas(i)=cell%eviscosg(i)
!
!........Liquid viscosity
!
         IF(turb_phase.eq.gas_only) THEN
            cell%tviscosl(i)=cell%tviscosg(i)*cell%rhol(i)/cell%rhog(i)
            cell%eviscosl(i)=cell%eviscosg(i)*cell%rhol(i)/cell%rhog(i)
         ENDIF
!
      ENDDO
!
      RETURN
      END SUBROUTINE turb_SST_vis_gas

