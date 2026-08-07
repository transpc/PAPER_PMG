!
      SUBROUTINE non_rv_wall_HT_2d
!
!     This routine calculates wall heat transfer coefficient and heat flux at Fluid-solid interface.
!     No RV heat structures are needed.       
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
      USE SOLID_DATA   , ONLY: solid      
      USE Zparam       , ONLY: ndim
      USE STM_TBL_cupid, ONLY: st_tbl,nt,ndxstd,pcrit
      USE Zqvol        , ONLY: qporous_gas,qporous_liq
      USE Zmodel       , ONLY: qconden
      USE Znormal      , ONLY: sa_walll
      USE Znum_cell    , ONLY: istart_nf,nf_number_nb,nf_number_id,istart_nfs
      USE Zwall_HTC    , ONLY: tsat_t,dt_sat,sigma,hfg,sat_hfp,h_liq,hfg_p,sat_hfp,qual_eq,vfg, & !pik-halden-dbug mul_o, &
                               mflux_gasa,mflux_liqa,mflux_tota,qflux_t,qflux_l,qflux_g,chf,incnd,tw
      USE Zvector      , ONLY: vl_n,vg_n
      USE Zvec_index   , ONLY: left_nf,right_fsw      
      USE Zvec_geo     , ONLY: fac1_fsw,fac_fsw
      USE Zvec_index_solid , ONLY: flux_fsw
      USE Zio_unit     , ONLY: unit_log                             
      USE Zwall_HTC     , ONLY: dia_rod,dia_rod_cfd
      USE Zcore         , ONLY: myrank
      USE Zzone         , ONLY: ncell_fluid
!            
      IMPLICIT NONE
!
      INTEGER :: i,ii,kk,mode,f_direc
      INTEGER :: nv,nf_number,len,istart0,istart,i0,i1
      LOGICAL :: err
      LOGICAL,SAVE :: initial=.true.
      REAL(8) :: x_flow,h_gas,h_mix
      REAL(8) :: qe,qc,vx,vy,vz,partition,tl1,tg1,condw
      REAL(8) :: ali_tmp,agi_tmp,akli,akgi      
      REAL(8) :: s(36)
      REAL(8) :: dia_rod_save,qflux_liq,qflux_gas 
      REAL(8),ALLOCATABLE,SAVE ::mul_o(:)
! 
!.....use dia_rod_cfd instead of dia_rod, which is used in RV models.
!      
      dia_rod_save=dia_rod
      dia_rod=dia_rod_cfd
!
      chf=0.0d0
      qconden(:)=0.0d0
!
      IF(initial)THEN
         initial=.FALSE.
         IF(myrank.eq.0)WRITE(*,"(11x,a)")'--non_rv_wall_HT_2d is turned on for f-s wall.'
         ALLOCATE(mul_o(ncell_fluid))
         mul_o(:)=0.0d0
      ENDIF   
!
!.....Build summation info for non,inl,fsw,ctw,chw
!
      nf_number_nb=1
      nf_number_id(1)=5
      istart_nfs(1)=0
!
!.....Fluid-Solid interface
!
      nv=1
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)         ! Fluid cell number
         kk=right_fsw(i)        ! Solid cell number
         ali_tmp=cell%alphal(ii)
         agi_tmp=cell%alphag(ii)       
         akli=ali_tmp*cell%condl(ii)
         akgi=agi_tmp*cell%condg(ii)
         partition=1.0d0
         IF(agi_tmp.ge.0.9d0)partition=1.0d0-10.0d0*(0.1d0-ali_tmp)
!
!........Surface wall temp. estimation from heat structure
!
         tl1=cell%tl_o(ii)
         tg1=cell%tg_o(ii)
         tw=((partition*cell%condl(ii)*tl1+(1.0d0-partition)*cell%condg(ii)*tg1)*fac1_fsw(i)   &
             +fac_fsw(i)*solid%conds(kk)*solid%tsol_o(kk)) 
         condw=((partition*cell%condl(ii)          &
             +(1.0d0-partition)*cell%condg(ii))*fac1_fsw(i)+fac_fsw(i)*solid%conds(kk))
         tw=tw/condw
!
         tsat_t=cell%tst(ii)   
         dt_sat=tw-tsat_t
         qe=0.0d0
         qc=0.0d0
!
         hfg=cell%hgsat(ii)-cell%hlsat(ii)
         sat_hfp=cell%hlsat(ii)    
         h_liq=cell%el(ii)+cell%p(ii)/cell%rhol(ii)  
         sigma=cell%sigma(ii)              
!
!........Determine the main flow direction
!    
         IF(cell%alphag(ii).ge.0.999d0)THEN    ! gas 1-phase flow
            vx=DABS(vg_n(ii,1))
            vy=DABS(vg_n(ii,2))
            IF(ndim.eq.3)vz=DABS(vg_n(ii,3))
         ELSE                                  ! Liquid phase dominent flow
            vx=DABS(vl_n(ii,1))
            vy=DABS(vl_n(ii,2))
            IF(ndim.eq.3)vz=DABS(vl_n(ii,3))            
         ENDIF 
         IF(vx.ge.vy)THEN
            f_direc=1      
            IF(ndim.eq.3)THEN
               IF(vz.ge.vx)f_direc=3
            ENDIF   
         ELSE
            f_direc=2
            IF(ndim.eq.3)THEN
               IF(vz.ge.vy)f_direc=3
            ENDIF 
         ENDIF         
!
!........Calculate thermal equilibrium quality(qual_eq), liquid enthalpy  on total pressure (sat_hfp), heat of vaporization on total pressure (hfg_p)
!                
         mflux_gasa=DMAX1(DABS(vg_n(ii,f_direc)*cell%rhog(ii)*cell%alphag(ii)),0.001d0) 
         mflux_liqa=DMAX1(DABS(vl_n(ii,f_direc)*cell%rhol(ii)*cell%alphal(ii)),0.001d0) 
         IF(cell%quala(ii).gt.1.d-9)THEN
            s(2)=cell%p(ii)
            IF(s(2).gt.pcrit)THEN
               WRITE(*,*)        'Pressure exceeds saturation bound'
               WRITE(unit_log,*) 'Pressure exceeds saturation bound'
               s(2)=DMIN1(pcrit,s(2))
            ENDIF
            s(9)=0.0d0
            CALL sth2x2_cupid(st_tbl(ndxstd),st_tbl(ndxstd+nt),s,err)
            hfg_p=s(16)-s(15)
            sat_hfp=s(15)
            x_flow=(mflux_gasa+0.01*cell%alphag(ii)*cell%rhog(ii))/(mflux_liqa+mflux_gasa+0.01*cell%rhom(ii))
            h_gas=cell%eg(ii)+cell%p(ii)/cell%rhog(ii)
            IF(s(2).lt.pcrit)THEN
               h_mix=h_liq+x_flow*(h_gas-h_liq)
               qual_eq=(h_mix-sat_hfp)/hfg_p
            ELSE
               qual_eq=cell%quals(ii)
            ENDIF         
         ELSE
            hfg_p=hfg
            qual_eq=cell%quals(ii)
         ENDIF         
!
!........Set mass flux and relative velocity in major direction (1-D base)
!   
         mflux_tota=mflux_gasa+mflux_liqa 
         vfg=DABS(vg_n(ii,f_direc)-vl_n(ii,f_direc))          
!
!-----------------------------------------------------------------------------------------      
!      
!........Mode Selection
!
         IF(cell%quala(ii).gt.0.99999999d0)THEN
            mode=0                                 ! Air-Water
            CALL single_phase_HTC(ii,mode)
         ELSEIF(cell%p(ii).gt.pcrit)THEN
            mode=10                                ! Critical fluid
            CALL single_phase_HTC(ii,mode)
         ELSEIF(tw.lt.cell%ts(ii)-0.001d0)THEN
            IF(cell%alphag(ii).lt.0.1d0)THEN
               mode=2                              ! Liquid 1-phase
               CALL single_phase_HTC(ii,mode)
            ELSEIF(cell%quala(ii).gt.0.999d0)THEN
               mode=0                              ! Air-Water
               CALL single_phase_HTC(ii,mode)
            ELSEIF(tw.gt.cell%tl(ii).and.cell%alphag(ii).lt.0.999d0)THEN
               mode=2                              ! Liquid 1-phase
               CALL single_phase_HTC(ii,mode)
            ELSE
               mode=11                             ! Condensation    
               incnd=0
               CALL condensation_HTC(ii,mode,sa_walll(ii),qc)
            ENDIF
         ELSEIF(dt_sat.le.0.0d0)THEN
            mode=0                                  ! Air-Water
            CALL single_phase_HTC(ii,mode)
         ELSEIF(tw.lt.cell%tl(ii))THEN
            mode=0                                  ! Air-Water
            CALL single_phase_HTC(ii,mode)
         ELSE
            IF(cell%alphag(ii).ge.0.999d0)THEN
               mode=9                               ! Gas 1-phase
               CALL single_phase_HTC(ii,mode)
            ELSEIF(dt_sat.gt.600.0d0)THEN
               !mode=7~8                            ! Film boiling
               CALL CHF_calc(ii)
               CALL trans_film_boiling_HTC(ii,mode)
               CALL subcooled_boiling(ii,mul_o(ii),sa_walll(ii),qe)                 ! Calculate Gamma_wall !!!
            ELSEIF(dt_sat.gt.100.0d0)THEN 
               !mode(i)=5~8                         ! Transient or Film boiling
               CALL CHF_calc(ii)
               CALL trans_film_boiling_HTC(ii,mode)
               CALL subcooled_boiling(ii,mul_o(ii),sa_walll(ii),qe)                 ! Calculate Gamma_wall !!!
            ELSE            
               CALL CHF_calc(ii)
               CALL nucl_boiling_HTC(ii,mode)
               IF(qflux_t.ge.chf)THEN
                  !mode(i)=5~8                      ! Transient or Film boiling
                  CALL trans_film_boiling_HTC(ii,mode)
                  CALL subcooled_boiling(ii,mul_o(ii),sa_walll(ii),qe)              ! Calculate Gamma_wall !!!
               ELSEIF(qflux_t.gt.0.0d0)THEN
                  !mode(i)=3,4
                  IF(cell%tl(ii).lt.cell%tst(ii))THEN
                     mode=3                         ! Subcooled Boiling
                  ELSE
                     mode=4                         ! Nucleate Boling
                  ENDIF
                  CALL subcooled_boiling(ii,mul_o(ii),sa_walll(ii),qe)              ! Calculate Gamma_wall !!!
               ELSE
                  mode=1                            ! Do nothing in MARS!!!
               ENDIF 
            ENDIF
         ENDIF !End of mode selection         
!
!........Energy Partitioning 
!                      
         qflux_liq=(qflux_l-qe)*sa_walll(ii)
         qflux_gas=(qflux_g+qc)*sa_walll(ii)
         flux_fsw(i)=qflux_liq+qflux_gas  
!         IF(qconden(ii).ne.0)flux_fsw(i)=0.0d0       ! heat transfer is considered in qconden when the condensation occurs
         qporous_liq(ii)=qporous_liq(ii)+qflux_liq
         qporous_gas(ii)=qporous_gas(ii)+qflux_gas
! 
      ENDDO
      RETURN
!
      dia_rod=dia_rod_save
!      
      END SUBROUTINE non_rv_wall_HT_2d
      
      
