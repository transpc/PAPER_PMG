!
      SUBROUTINE rv_allocate_var_2d
!
!     This routine allocates variables containing the information on models.
!
      USE Zwall_HTC    , ONLY: mode,Mul_o,Mul1d_o,  &
                                h_index,chfr,hmode_rv,twall_rv,dt_sat_rv,chf_rv
      USE Zrv_hts_2d   , ONLY: hlr_2f,hgr_2f,hstr_2f,hspr_2f,tlr_2f,tgr_2f,tstr_2f,tspr_2f,&
                                wet_b,wet_t,ztop_q,nrod_2d,wet_bi,wet_ti,nz0_2d
      USE Zrv_ncell    , ONLY: ncell_fluid_core,ncell_fuel_rod,nz_fine
      !OPR1000 rod-scale
      USE Zporous , ONLY: l_subchannel
      USE Zwall_HTC    , ONLY: gamma_wall_rod
      USE Zwall_HTC,ONLY:qflux_l0,qflux_g0,qf0,qf1,qg0,qg1,gw0,gw1    
      USE Zmpi,ONLY:ncell_fp       
!
      USE Zrv_hts_2d  , ONLY:twall_fuel
      USE Zzone       , ONLY:ncell_fluid
!
      IMPLICIT NONE
!
      INTEGER n,m,k,i,nz
!
      n=ncell_fluid_core
      m=ncell_fuel_rod
      k=nrod_2d
      nz=nz0_2d*nz_fine
!
!.....Zwall_HTC
!
      IF(n.gt.0) THEN
         ALLOCATE(mode(n),h_index(n))
         ALLOCATE(hmode_rv(n,5),twall_rv(n,5),dt_sat_rv(n),chf_rv(n))
      ELSE
         ALLOCATE(mode(1),h_index(1))
         ALLOCATE(hmode_rv(1,5),twall_rv(1,5),dt_sat_rv(1),chf_rv(1))
      ENDIF
      mode(:)=0
      h_index(:)=0 
      hmode_rv(:,:)=0
      twall_rv(:,:)=0.0d0            
      dt_sat_rv(:)=0.0d0
      chf_rv(:)=0.0d0           
      ALLOCATE(Mul_o(n),Mul1d_o(n))
      Mul_o(:)=0.0d0   
      Mul1d_o(:)=0.0d0  
      ALLOCATE(chfr(n))
      chfr(:)=0.0d0         
!
!.....Fuel rod
!
      ALLOCATE(hlr_2f(m),hgr_2f(m),hstr_2f(m),hspr_2f(m),tlr_2f(m),tgr_2f(m),tstr_2f(m),tspr_2f(m))
      hlr_2f(:)=0.0d0
      hgr_2f(:)=0.0d0
      hstr_2f(:)=0.0d0
      hspr_2f(:)=0.0d0
      tlr_2f(:)=0.0d0
      tgr_2f(:)=0.0d0
      tstr_2f(:)=0.0d0
      tspr_2f(:)=0.0d0
      ALLOCATE(wet_b(k),wet_t(k),ztop_q(k))
      wet_b(:)=0.0d0
      wet_t(:)=0.0d0
      ztop_q(:)=0.0d0
      ALLOCATE(wet_bi(k),wet_ti(k))
      wet_bi(:)=0
      wet_ti(:)=0
!
      !OPR1000 rod-scale
      IF(l_subchannel)then
         ALLOCATE(gamma_wall_rod(m))
         gamma_wall_rod=0.0d0

         ALLOCATE(qflux_l0(m),qflux_g0(m))
         ALLOCATE(qf0(ncell_fp),qf1(ncell_fp))
         ALLOCATE(qg0(ncell_fp),qg1(ncell_fp))
         ALLOCATE(gw0(ncell_fp),gw1(ncell_fp))

         do i=1,m
            qflux_l0(i)=0.0d0
            qflux_g0(i)=0.0d0
         enddo
         do i=1,ncell_fp
            qf0(i)=0.0d0
            qf1(i)=0.0d0
            qg0(i)=0.0d0
            qg1(i)=0.0d0
            gw0(i)=0.0d0
            gw1(i)=0.0d0
         enddo
      ENDIF
!
!.....Fuel surface temperature getting from FRAPTRAN
!      
      ALLOCATE(twall_fuel(ncell_fluid))      
      twall_fuel(:)=0.d0
!      
      RETURN
      END SUBROUTINE rv_allocate_var_2d
