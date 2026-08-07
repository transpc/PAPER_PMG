!
      SUBROUTINE rv_wall_HT_1d(rfluid)
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
!     mode 10 -- single phase liquid convection (p .ge. pcrit)                          
!     mode 11 -- condensation                         
!     add 20 to mode number if quala .gt. 1.0e-9 
!
      USE VOL_DATA       , ONLY: cell
      USE STM_TBL_cupid  , ONLY: st_tbl,    &
                                 nt,ndxstd, &
                                 pcrit
      USE Zqvol        , ONLY: qporous_gas,qporous_liq
      USE Zrv_hts_1d   , ONLY: ncell_hts_1d,nr_1d,t_hts_1d,cupid_cell_1d,ig_hts_1d,         &
                               hll_1d,hstl_1d,hspl_1d,hgl_1d,tll_1d,tstl_1d,tspl_1d,tgl_1d, &
                               hlr_1d,hstr_1d,hspr_1d,hgr_1d,tlr_1d,tstr_1d,tspr_1d,tgr_1d, &
                               ht_area_left_1d,ht_area_right_1d,bcl_1d,bcr_1d
      USE Zwall_HTC    , ONLY: mode_1d,f_direc,c_direc1,c_direc2,                          &
                               chfr,k_grid,dis_grid,                                       &
                               f_pp_axial,horiz_chf,angle_horiz,mul1d_o,                   &
                               qflux_t,qflux_l,qflux_g,HTC_tl,HTC_tst,HTC_tg,HTC_tgp,      &
                               mflux_liqa,mflux_gasa,mflux_tota,tw,                        &
                               chf,tsat_t,dt_sat,hfg,sat_hfp,h_liq,sigma,hfg_p,qual_eq,vfg
      USE Zvector      , ONLY: vl_n,vg_n
      USE Zio_unit     , ONLY: unit_log
!
      IMPLICIT NONE
!      
      INTEGER i,k,m,ig,rfluid
!
      LOGICAL err      
!
      REAL(8) s(36),qe,qc
      REAL(8) x_flow,h_gas,h_mix     
      REAL(8) ht_area_1d
!      
!.....pcrit is constant read in stread.f90 why change?
!     pcrit=22.4d6
      chf=0.0d0
!
!.....Main Loop for Wall HTC Calculation
!
!      IF(rfluid.eq.1)THEN
!         qporous_liq(:)=0.0d0 !apr1400_lbloca_porous
!         qporous_gas(:)=0.0d0 !apr1400_lbloca_porous
!      ENDIF
!
!.....Main Loop for Wall HTC Calculation
!
      DO k=1,ncell_hts_1d
!
         i=cupid_cell_1d(k)
         m=i
         ig=ig_hts_1d(k)
!
!........Load wall temp. from heat structure
!        
         IF(rfluid.eq.1)THEN
            tw=t_hts_1d(k,1)
            ht_area_1d=ht_area_left_1d(k)
            IF(bcl_1d(ig).gt.3)CYCLE
         ELSEIF(rfluid.eq.2)THEN
            tw=t_hts_1d(k,nr_1d)
            ht_area_1d=ht_area_right_1d(k)
            IF(bcr_1d(ig).gt.3)CYCLE
         ENDIF
!
!........Set bagic properties
!             
         tsat_t=cell%tst(i)   
         dt_sat=tw-tsat_t
         qe=0.0d0
         qc=0.0d0
!
         hfg=cell%hgsat(i)-cell%hlsat(i)
         sat_hfp=cell%hlsat(i)    
         h_liq=cell%el(i)+cell%p(i)/cell%rhol(i)  
         sigma=cell%sigma(i)                  
!
!........Calculate thermal equilibrium quality(qual_eq), liquid enthalpy  on total pressure (sat_hfp), heat of vaporization on total pressure (hfg_p)
!
         mflux_gasa=MAX(ABS(vg_n(m,f_direc)*cell%rhog(m)*cell%alphag(m)),0.001d0) 
         mflux_liqa=MAX(ABS(vl_n(m,f_direc)*cell%rhol(m)*cell%alphal(m)),0.001d0) 
         IF(cell%quala(i).gt.1.d-9)THEN
            s(2)=cell%p(i)
            IF(s(2).gt.pcrit)THEN
               WRITE(*,*) 'Pressure exceeds saturation bound'
               WRITE(unit_log,*) 'Pressure exceeds saturation bound'
               s(2)=MIN(pcrit,s(2))
            ENDIF
            s(9)=0.0d0
            CALL sth2x2_cupid(st_tbl(ndxstd),st_tbl(ndxstd+nt),s,err)
            hfg_p=s(16)-s(15)
            sat_hfp=s(15)
            x_flow=(mflux_gasa+0.01*cell%alphag(i)*cell%rhog(i))/(mflux_liqa+mflux_gasa+0.01*cell%rhom(i))
            h_gas=cell%eg(i)+cell%p(i)/cell%rhog(i)
            IF(s(2).lt.pcrit)THEN
               h_mix=h_liq+x_flow*(h_gas-h_liq)
               qual_eq=(h_mix-sat_hfp)/hfg_p
            ELSE
               qual_eq=cell%quals(i)
            ENDIF         
         ELSE
            hfg_p=hfg
            qual_eq=cell%quals(i)
         ENDIF
!
!...........Set mass flux and relative velocity in major direction (1-D base)
!   
         mflux_tota=mflux_gasa+mflux_liqa 
         vfg=ABS(vg_n(i,f_direc)-vl_n(i,f_direc))     
!
!--------------------------------------------------------------------------------------------      
!      
!...........Mode Selection
!
         IF(cell%quala(i).gt.0.99999999d0)THEN
            mode_1d(k)=0                                 ! Air-Water
            CALL single_phase_HTC(m,mode_1d(k))
         ELSEIF(cell%p(i).ge.pcrit)THEN
            mode_1d(k)=10                                ! Critical fluid
            CALL single_phase_HTC(m,mode_1d(k))
         ELSEIF(tw.lt.cell%ts(i))THEN
            IF(cell%alphag(i).lt.0.1d0)THEN
               mode_1d(k)=2                              ! Liquid 1-phase
               !!!!cyj: CHF calculation to calculate CHF margin, not for HTC
               CALL single_phase_HTC(m,mode_1d(k))
            ELSEIF(cell%quala(i).gt.0.999d0)THEN
               mode_1d(k)=0                              ! Air-Water
               CALL single_phase_HTC(m,mode_1d(k))
            ELSEIF(tw.gt.cell%tl(i).and.cell%alphag(i).lt.0.999d0)THEN
               mode_1d(k)=2                              ! Liquid 1-phase
               CALL single_phase_HTC(m,mode_1d(k))
            ELSE
               mode_1d(k)=11                             ! Condensation    
               CALL condensation_HTC(m,mode_1d(k),ht_area_1d,qc)
            ENDIF
         ELSEIF(dt_sat.le.0.0d0)THEN
            mode_1d(k)=0                                  ! Air-Water
            CALL single_phase_HTC(m,mode_1d(k))
         ELSEIF(tw.lt.cell%tl(i))THEN
            mode_1d(k)=0                       
            CALL single_phase_HTC(m,mode_1d(k))
         ELSE
            IF(cell%alphag(i).ge.0.999d0)THEN
               mode_1d(k)=9                                ! Gas 1-phase
               CALL single_phase_HTC(m,mode_1d(k))
            ELSEIF(dt_sat.gt.600.0d0)THEN
               !mode_1d(k)=7~8                            ! Film boiling
               CALL CHF_calc(m)
               CALL trans_film_boiling_HTC(m,mode_1d(k))
               CALL subcooled_boiling(m,mul1d_o(m),ht_area_1d,qe)              ! Calculate Gamma_wall !!!
            ELSEIF(dt_sat.gt.100.0d0)THEN 
               !mode_1d(k)=5~8                              ! Transient or Film boiling
               CALL CHF_calc(m)
               CALL trans_film_boiling_HTC(m,mode_1d(k))
               CALL subcooled_boiling(m,mul1d_o(m),ht_area_1d,qe)              ! Calculate Gamma_wall !!!
            ELSE            
               CALL CHF_calc(m)
               CALL nucl_boiling_HTC(m,mode_1d(k))
               IF(qflux_t.ge.chf)THEN
                  !mode_1d(k)=5~8                            ! Transient or Film boiling
                  CALL trans_film_boiling_HTC(m,mode_1d(k))
               ELSEIF(qflux_t.gt.0.0d0)THEN
                  !mode_1d(k)=3,4
                  IF(cell%tl(i).lt.cell%ts(i)-0.1d0)THEN
                     mode_1d(k)=3                            ! Subcooled Boiling
                  ELSE
                     mode_1d(k)=4                            ! Nucleate Boling
                  ENDIF
               ELSE
                  mode_1d(k)=1                               ! Do nothing in MARS!!!
               ENDIF 
               CALL subcooled_boiling(m,mul1d_o(m),ht_area_1d,qe)              ! Calculate Gamma_wall !!!
            ENDIF
         ENDIF !End of mode selection
!
!........Mode correction according to the exsistance of NC gas
!
!         IF(cell%quala(i).ge.1.d-9)mode_1d(k)=mode_1d(k)+20         
!
!........Relaxation
!
!
!........Energy Partitioning !!!cyj
!
         qporous_liq(i)=qporous_liq(i)+(qflux_l-qe)*ht_area_1d
         qporous_gas(i)=qporous_gas(i)+(qflux_g+qc)*ht_area_1d
         IF(chf.ne.0)chfr(i)=qflux_t/chf   
!
!........Save calculated HTC in HS array
!
         IF(rfluid.eq.1)THEN
            hll_1d(k)=HTC_tl
            hgl_1d(k)=HTC_tg
            hstl_1d(k)=HTC_tst
            hspl_1d(k)=HTC_tgp
!
            tll_1d(k)=cell%tl(i)
            tgl_1d(k)=cell%tg(i)
            tstl_1d(k)=cell%tst(i)
            tspl_1d(k)=cell%ts(i)               
         ELSEIF(rfluid.eq.2)THEN
            hlr_1d(k)=HTC_tl
            hgr_1d(k)=HTC_tg
            hstr_1d(k)=HTC_tst
            hspr_1d(k)=HTC_tgp
!
            tlr_1d(k)=cell%tl(i)
            tgr_1d(k)=cell%tg(i)
            tstr_1d(k)=cell%tst(i)
            tspr_1d(k)=cell%ts(i)           
         ENDIF
!
      ENDDO
!
      END SUBROUTINE rv_wall_HT_1d
