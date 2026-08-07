!
      SUBROUTINE rv_wall_HT_2d
!
!     This routine calculates wall heat transfer coefficient and heat flux to each phase according to 'mode'
!
!     mode 0  -- air-water mixture convection                            
!     mode 1  -- single phase liquid convection (Natural convection)
!     mode 2  -- single phase liquid convection (Forced convection)
!     mode 3  -- subcooled nucleate boiling                              
!     mode 4  -- saturated nucleate boiling                              
!     mode 5  -- subcooled transition film boiling                       
!     mode 6  -- saturated transition film boiling                       
!     mode 7  -- subcooled film boiling                                  
!     mode 8  -- saturated film boiling                                  
!     mode 9  -- single phase vapor convection                           
!     mode 10 -- single phase liquid convection (p .ge. pcritical)                          
!     mode 11 -- condensation                         
!     Add 20 to mode number if quala .gt. 1.0e-9 
!     Add 40 to mode number if Reflood=1
!
      USE VOL_DATA     , ONLY: cell
      USE Zcore        , ONLY: np
      USE STM_TBL_cupid , ONLY: st_tbl,    &
                                nt,ndxstd, &
                                pcrit
      USE Zparam       , ONLY: pi
      USE Zqvol        , ONLY: qporous_gas,qporous_liq
      USE Zvector      , ONLY: vl_n,vg_n
      USE Zwall_HTC    , ONLY: mode,f_direc,c_direc1,c_direc2,                              &
                               reflood,base_bundle,zqf,zqf_top,h_bundle,hmode_rv,           &
                               chfr,dt_sat_rv,chf_rv,k_grid,dis_grid,                       &
                               f_pp_axial,horiz_chf,angle_horiz,twall_rv,                   &
                               qflux_t,qflux_l,qflux_g,HTC_tl,HTC_tst,HTC_tg,HTC_tgp,       &
                               mflux_liqa,mflux_gasa,mflux_tota,tw,                         &
                               chf,tsat_t,dt_sat,hfg,sat_hfp,h_liq,sigma,hfg_p,qual_eq,vfg, &
                               mul_o,incnd
      USE Zrv_hts_2d   , ONLY: nr_2d,nz0_2d,twall_fuel 
      USE Zrv_hts_2d   , ONLY: hlr_2f,hgr_2f,hstr_2f,hspr_2f,tlr_2f,tgr_2f,tstr_2f,tspr_2f, &
                                t_fuel,z_fuel,wet_b,wet_t,ztop_q,dz_fuel,nrod_2d,           &
                                ht_area_fuel0,ht_geo_2d,ri_2d,dz_fuel0,wet_bi,wet_ti
      USE Zrv_ncell    , ONLY: ncell_fluid_core,ncell_fuel_rod,cupid_cell_hts2d,nz_fuel_rod, &
                                channel_cell_hts2d,n_channel_fluid,ncell_fluid_core,         &
                                cupid_cell_channel,nrod_fuel_rod,nz_fine,num_fuel_rod
      !OPR1000 rod-scale
      USE Zporous      , ONLY: l_subchannel   
      USE Zio_unit     , ONLY: unit_log                             
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,k,m,n,nz,nrod
      LOGICAL :: err
      LOGICAL,SAVE :: initial=.true.
      REAL(8) :: x_flow,h_gas,h_mix
      REAL(8) :: qe,qc,ep
      REAL(8) :: zi_2d_m
      REAL(8) pcritical
!
!.....Local arrays
      REAL(8) s(36)
      INTEGER,DIMENSION(:),ALLOCATABLE :: index_qfb,index_qft
      REAL(8),DIMENSION(:),ALLOCATABLE :: qflux_liq,qflux_gas
!      
      !OPR1000 rod-scale
      IF(l_subchannel)then
         CALL rv_wall_HT_2d_rod
         RETURN
      ENDIF
      pcritical=pcrit
      ep=1.0d-10 
!
!.....Initialization of variables for Reflood calculation
!
      ALLOCATE(qflux_liq(ncell_fluid_core),qflux_gas(ncell_fluid_core))
      qflux_liq(:)=0.0d0
      qflux_gas(:)=0.0d0
      chf=0.0d0
!
      IF(reflood.eq.1)THEN
         ALLOCATE(index_qfb(nrod_2d),index_qft(nrod_2d))
         index_qfb(:)=0
         index_qft(:)=0
      ENDIF
!
!.....Save wall tempertures calculated in Heat Structure module at the old time step. This variable is used only for Post-processing not for calculation.
!
      twall_rv(:,1)=0.0d0
      DO i=1,ncell_fuel_rod        
         m=channel_cell_hts2d(i)
         twall_rv(m,1)=twall_rv(m,1)+t_fuel(i,nr_2d)
         twall_rv(m,2)=twall_rv(m,2)+t_fuel(i,1) 
      ENDDO
      twall_rv(:,1)=twall_rv(:,1)/nz_fine
      twall_rv(:,2)=twall_rv(:,2)/nz_fine 
!
!.....ht_area_fuel0 for various dz  
!
      IF(initial) then
         ALLOCATE(ht_area_fuel0(ncell_fuel_rod))
         ALLOCATE(dz_fuel(ncell_fuel_rod))
         ht_area_fuel0=0.0d0
         dz_fuel=0.0d0
         DO i=1,ncell_fuel_rod
            k=cupid_cell_hts2d(i)
            m=channel_cell_hts2d(i)
            n=n_channel_fluid(m)
            nz=nz_fuel_rod(i)
            dz_fuel(i)=dz_fuel0(nz) !/nz_fine
            IF(ht_geo_2d.eq.1) THEN
               ht_area_fuel0(i)=4.0d0*ri_2d(nr_2d)*ri_2d(nr_2d)*dz_fuel(i)*DBLE(num_fuel_rod(n))
            ELSEIF(ht_geo_2d.eq.2) THEN
!               ht_area_fuel0(i)=2.0d0*pi*ri_2d(nr_2d)*dz_fuel(i)*DBLE(num_fuel_rod(n))
               ht_area_fuel0(i)=2.0d0*pi*ri_2d(nr_2d)*dz_fuel0(nz)*DBLE(num_fuel_rod(n))
            ENDIF
         ENDDO
         initial=.false.
      ENDIF
!
!.....Main Loop for Wall HTC Calculation
!
      DO i=1,ncell_fuel_rod
         k=channel_cell_hts2d(i)
         m=cupid_cell_hts2d(i)
         n=n_channel_fluid(k)
         nz=nz_fuel_rod(i)
         nrod=nrod_fuel_rod(i)
!
!........Load wall temp. from heat structure
!
         tw=t_fuel(i,nr_2d) 
!         
         twall_fuel(m)=tw    
! 
!........Set bagic properties
!             
         tsat_t=cell%tst(m)   
         dt_sat=tw-tsat_t
         qe=0.0d0
         qc=0.0d0
!
         hfg=cell%hgsat(m)-cell%hlsat(m)
         sat_hfp=cell%hlsat(m)    
         h_liq=cell%el(m)+cell%p(m)/cell%rhol(m)  
         sigma=cell%sigma(m)                  
!
!........Calculate thermal equilibrium quality(qual_eq), liquid enthalpy  on total pressure (sat_hfp), heat of vaporization on total pressure (hfg_p)
!
         mflux_gasa=DMAX1(DABS(vg_n(m,f_direc)*cell%rhog(m)*cell%alphag(m)),0.001d0) 
         mflux_liqa=DMAX1(DABS(vl_n(m,f_direc)*cell%rhol(m)*cell%alphal(m)),0.001d0) 
         IF(cell%quala(m).gt.1.d-9)THEN
            s(2)=cell%p(m)
            IF(s(2).gt.pcritical)THEN
               WRITE(*,*) 'Pressure exceeds saturation bound'
               WRITE(unit_log,*) 'Pressure exceeds saturation bound'
               s(2)=DMIN1(pcritical,s(2))
            ENDIF
            s(9)=0.0d0
            CALL sth2x2_cupid(st_tbl(ndxstd),st_tbl(ndxstd+nt),s,err)
            hfg_p=s(16)-s(15)
            sat_hfp=s(15)
            x_flow=(mflux_gasa+0.01*cell%alphag(m)*cell%rhog(m))/(mflux_liqa+mflux_gasa+0.01*cell%rhom(m))
            h_gas=cell%eg(m)+cell%p(m)/cell%rhog(m)
            IF(s(2).lt.pcritical)THEN
               h_mix=h_liq+x_flow*(h_gas-h_liq)
               qual_eq=(h_mix-sat_hfp)/hfg_p
            ELSE
               qual_eq=cell%quals(m)
            ENDIF         
         ELSE
            hfg_p=hfg
            qual_eq=cell%quals(m)
         ENDIF
!
!........Set mass flux and relative velocity in major direction (1-D base)
!   
         mflux_tota=mflux_gasa+mflux_liqa 
         vfg=DABS(vg_n(m,f_direc)-vl_n(m,f_direc)) 
!
!........Calculate 'Reflood parameters'
!
         IF(reflood.eq.1)THEN
            zi_2d_m=z_fuel(nz)-base_bundle
            zi_2d_m=zi_2d_m+0.5d0*dz_fuel0(nz)
            ztop_q(nrod)=h_bundle-wet_t(nrod)
!
!            zqf=zi_2d_m-wet_b(nrod)
            zqf=DMAX1(0.0d0,zi_2d_m-wet_b(nrod))  !  yjm
!            
            zqf_top=DMAX1(0.d0,ztop_q(nrod)-zi_2d_m)
            !zqf_top=ztop_q(nrod)-zi_2d_m
         ENDIF      
!
!-----------------------------------------------------------------------------------------      
!      
!........Mode Selection
!
         IF(cell%quala(m).gt.0.99999999d0)THEN
            mode(k)=0                                 ! Air-Water
            CALL single_phase_HTC(m,mode(k))
         ELSEIF(cell%p(m).gt.pcritical)THEN
            mode(k)=10                                ! Critical fluid
            CALL single_phase_HTC(m,mode(k))
         ELSEIF(tw.lt.cell%ts(m)-0.001d0)THEN
            IF(cell%alphag(m).lt.0.1d0)THEN
               mode(k)=2                              ! Liquid 1-phase
               !!!!cyj: CHF calculation to calculate CHF margin, not for HTC
            CALL single_phase_HTC(m,mode(k))
            ELSEIF(cell%quala(m).gt.0.999d0)THEN
               mode(k)=0                              ! Air-Water
            CALL single_phase_HTC(m,mode(k))
            ELSEIF(tw.gt.cell%tl(m).and.cell%alphag(m).lt.0.999d0)THEN
               mode(k)=2                              ! Liquid 1-phase
            CALL single_phase_HTC(m,mode(k))
            ELSE
               mode(k)=11                             ! Condensation    
               incnd=0
               CALL condensation_HTC(m,mode(k),ht_area_fuel0(i),qc)
            ENDIF
         ELSEIF(dt_sat.le.0.0d0)THEN
            mode(k)=0                                  ! Air-Water
            CALL single_phase_HTC(m,mode(k))
         ELSEIF(tw.lt.cell%tl(m))THEN
            mode(k)=0                                  ! Air-Water
            CALL single_phase_HTC(m,mode(k))
         ELSE
            IF(cell%alphag(m).ge.0.999d0)THEN
               mode(k)=9                                ! Gas 1-phase
               CALL single_phase_HTC(m,mode(k))
            ELSEIF(dt_sat.gt.600.0d0)THEN
               !mode(i)=7~8                            ! Film boiling
               CALL CHF_calc(m)
               IF(reflood.eq.1)THEN
                  CALL trans_film_reflood_HTC(m,mode(k))
               ELSE
                  CALL trans_film_boiling_HTC(m,mode(k))
               ENDIF 
               CALL subcooled_boiling(m,mul_o(k),ht_area_fuel0(i),qe)              ! Calculate Gamma_wall !!!
            ELSEIF(dt_sat.gt.100.0d0)THEN 
               !mode(i)=5~8                              ! Transient or Film boiling
               CALL CHF_calc(m)
               IF(reflood.eq.1)THEN
                  CALL trans_film_reflood_HTC(m,mode(k))
               ELSE
                  CALL trans_film_boiling_HTC(m,mode(k))
               ENDIF
                  CALL subcooled_boiling(m,mul_o(k),ht_area_fuel0(i),qe)              ! Calculate Gamma_wall !!!
            ELSE            
               CALL CHF_calc(m)
               CALL nucl_boiling_HTC(m,mode(k))
               IF(qflux_t.ge.chf)THEN
                  !mode(i)=5~8                            ! Transient or Film boiling
                  IF(reflood.eq.1)THEN
                     CALL trans_film_reflood_HTC(m,mode(k))
                  ELSE
                     CALL trans_film_boiling_HTC(m,mode(k))
                  ENDIF 
                  CALL subcooled_boiling(m,mul_o(k),ht_area_fuel0(i),qe)              ! Calculate Gamma_wall !!!
               ELSEIF(qflux_t.gt.0.0d0)THEN
                  !mode(i)=3,4
                  IF(cell%tl(m).lt.cell%tst(m))THEN
                     mode(k)=3                            ! Subcooled Boiling
                  ELSE
                     mode(k)=4                            ! Nucleate Boling
                  ENDIF
                  CALL subcooled_boiling(m,mul_o(k),ht_area_fuel0(i),qe)              ! Calculate Gamma_wall !!!
               ELSE
                  mode(k)=1                               ! Do nothing in MARS!!!
               ENDIF 
!
            ENDIF
         ENDIF !End of mode selection
!
!........Energy Partitioning !!!cyj
!                      
         qflux_liq(k)=qflux_liq(k)+(qflux_l-qe)*ht_area_fuel0(i)
         qflux_gas(k)=qflux_gas(k)+(qflux_g+qc)*ht_area_fuel0(i)
!         
         IF(chf.ne.0) chfr(k)=qflux_t/chf     
!       
!........Calculate 'Quenching Front,QF'
!
         IF(reflood.eq.1)THEN
            IF(mode(k).ge.5.and.mode(k).ne.10)THEN   !shrink
               IF(zi_2d_m.lt.wet_b(nrod)-ep.and.dt_sat.gt.40.d0)THEN
                  index_qfb(nrod)=index_qfb(nrod)-1
               ENDIF
               IF(zi_2d_m.gt.(ztop_q(nrod)+ep).and.dt_sat.gt.40.d0)THEN
                  index_qft(nrod)=index_qft(nrod)-1
               ENDIF
            ELSE                                     !advance
               IF(zi_2d_m.gt.wet_b(nrod)+ep)THEN        
                  IF(zi_2d_m.le.(wet_b(nrod)+dz_fuel0(nz)+ep))THEN                                      !+1.0d-6
                     IF(zi_2d_m.lt.ztop_q(nrod)-ep)THEN     
!                        IF(dt_sat.lt.40.d0)index_qfb(nrod)=index_qfb(nrod)+1                     
                        index_qfb(nrod)=index_qfb(nrod)+1   !  yjm
                     ENDIF
                  ELSEIF((zi_2d_m+dz_fuel0(nz)+ep).ge.ztop_q(nrod))THEN
                     IF(zi_2d_m.gt.wet_b(nrod)+dz_fuel0(nz)+ep)THEN  
                        IF(ztop_q(nrod).ge.zi_2d_m-ep)THEN
                           IF(zi_2d_m.le.ztop_q(nrod)-ep)THEN
                              IF(dt_sat.lt.40.d0)THEN
                                 index_qft(nrod)=index_qft(nrod)+1
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
            ztop_q(nrod)=h_bundle-wet_t(nrod)
            zqf=zi_2d_m-wet_b(nrod)
            zqf_top=ztop_q(nrod)-zi_2d_m
         ENDIF
!
!........Save calculated HTC in HS array
!
         hlr_2f(i)=HTC_tl
         hgr_2f(i)=HTC_tg
         hstr_2f(i)=HTC_tst
         hspr_2f(i)=HTC_tgp
!
         tlr_2f(i)=cell%tl(m)
         tgr_2f(i)=cell%tg(m)
         tstr_2f(i)=cell%tst(m)
         tspr_2f(i)=cell%ts(m)
!
         dt_sat_rv(k)=dt_sat
         chf_rv(k)=chf                
! 
      ENDDO
!
      DO k=1,ncell_fluid_core
         m=cupid_cell_channel(k)
         qporous_liq(m)=qporous_liq(m)+qflux_liq(k)
         qporous_gas(m)=qporous_gas(m)+qflux_gas(k)
      ENDDO
!
      DEALLOCATE(qflux_liq,qflux_gas)
!
      hmode_rv(:,1)=mode(:)
!
      IF(reflood.eq.1)THEN
         IF(np.gt.1) THEN
            CALL allreducei_i(index_qfb,nrod_2d)
            CALL allreducei_i(index_qft,nrod_2d)
         ENDIF
!
!.....When different dz, dz_fuel(i), is defined,
!     Checek array size !!!
!       - dz_fuel0(nz0_2d*nz_fine)
!       - dz_fuel (ncell_fuel_rod)
!     
         DO i=1, nrod_2d
            wet_bi(i)=wet_bi(i)+index_qfb(i)
            wet_bi(i)=MIN(MAX(0,wet_bi(i)),nz0_2d*nz_fine)
            IF(wet_bi(i).eq.0)THEN
               wet_b(i)=0.0d0
            ELSE
               wet_b(i)=z_fuel(wet_bi(i)) 
            ENDIF
            wet_ti(i)=wet_ti(i)+index_qft(i)     
            wet_ti(i)=MIN(MAX(0,wet_ti(i)),nz0_2d*nz_fine)      
            IF(wet_ti(i).eq.0)THEN
               wet_t(i)=0.0d0
               ztop_q(i)=h_bundle
            ELSE
               wet_t(i)=h_bundle-z_fuel(wet_ti(i))
               ztop_q(i)=z_fuel(wet_ti(i))
            ENDIF             
         ENDDO
!
!         DO i=1, nrod_2d
!            wet_b(i)=wet_b(i)+index_qfb(i)*dz_fuel0(i)
!            wet_b(i)=DMIN1(DMAX1(0.0d0,wet_b(i)),h_bundle)
!!
!            wet_t(i)=wet_t(i)+index_qft(i)*dz_fuel0(i)
!            wet_t(i)=DMIN1(DMAX1(0.0d0,wet_t(i)),h_bundle)
!            ztop_q(i)=h_bundle-wet_t(i)
!         ENDDO
         DEALLOCATE(index_qfb,index_qft)
      ENDIF
!
      END SUBROUTINE rv_wall_HT_2d
