!
      SUBROUTINE turb_SST_gas
!
!     This routine calculates turbulence liquid viscosity based 
!     on k-e model.
!
      USE Zinterface
      USE VOL_DATA , ONLY: cell
      USE Zzone    , ONLY: ncell_fluid,ncell_fluid_all
      USE Zcore    , ONLY: np
      USE Zparam   , ONLY: cmu
      USE Zconst1  , ONLY: restart       
      USE Zndforce , ONLY: d_bfc,c_bface                 
      USE Zturb    , ONLY: turb_keg,turb_dpg,turb_keg_o,turb_dpg_o,  &
                            pro_keg,utaug,yplusg,veltg,strn_keg !,walln
      USE Zvector  , ONLY: vg_o
!
      IMPLICIT NONE
!      
!.....Local variables
      LOGICAL,SAVE :: initial2=.true.
      INTEGER :: i
!.....Local allocatablearrays
      REAL(8),DIMENSION(:),ALLOCATABLE :: utaug_all
!      
      ALLOCATE(utaug_all(ncell_fluid_all))
      utaug_all(:)=0.d0      
!
!...Initialize turbulence kinetic energy and dissipation
!
      IF(initial2.and.restart.eq.0) THEN
         CALL init_SST(turb_keg,turb_dpg,vg_o)
         initial2=.false.
      ENDIF
!
!.....Shift solutions
!
      DO i=1,ncell_fluid
         turb_keg_o(i)=turb_keg(i)
         turb_dpg_o(i)=turb_dpg(i)
      ENDDO
!
!.....Calculate tangential velocity and wall sheer stress
!
      CALL turb_SST_vt_tauw_gas(veltg,vg_o)
!
!.....Calculate yplus (using updated utau)
!
      IF(np.gt.1)THEN
         CALL allgatherv_r(utaug,utaug_all,ncell_fluid,ncell_fluid_all,0)
      ELSE
         utaug_all(:)=utaug(:)
      ENDIF
      DO i=1,ncell_fluid
        yplusg(i)=cell%rhog(i)*utaug_all(c_bface(i))*d_bfc(i)/cell%lviscosg(i) 
      ENDDO        
!
!.....Communicate velt
!
      IF(np.gt.1) CALL communicate_1d(veltg)         
!
!.....Caculate kinetic energy production
!
      CALL turb_SST_product(2,pro_keg,strn_keg)
!
!.....Solve k-e equation
!      
      CALL turb_SST_calc_gas
!
      DEALLOCATE(utaug_all)
!
      END SUBROUTINE turb_SST_gas

