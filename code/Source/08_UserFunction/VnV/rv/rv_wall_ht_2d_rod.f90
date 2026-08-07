!
      SUBROUTINE rv_wall_ht_2d_rod
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
!     Add 20 to mode number if quala .gt. 1.d-9 
!     Add 40 to mode number if Reflood=1
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell
      USE STM_TBL_cupid, ONLY: st_tbl,    &
                               nt,ndxstd, &
                               pcrit
      USE Zqvol        , ONLY: qporous_gas,qporous_liq
      USE Zvector      , ONLY: vl_n,vg_n
      USE Zwall_HTC    , ONLY: mode,f_direc,c_direc1,c_direc2,                        &
                               reflood,                                               &
                               chfr,dt_sat_rv,chf_rv,k_grid,dis_grid,                 &
                               f_pp_axial,horiz_chf,angle_horiz,twall_rv,             &
                               qflux_t,qflux_l,qflux_g,HTC_tl,HTC_tst,HTC_tg,HTC_tgp, &
                               mflux_liqa,mflux_gasa,mflux_tota,tw, &
                               chf,tsat_t,dt_sat,hfg,sat_hfp,h_liq,sigma,hfg_p,qual_eq,vfg, &
                               mul_o
      USE Zwall_HTC    , ONLY: qflux_l0,qflux_g0
      USE Zrv_hts_2d   , ONLY: nr_2d 
      USE Zrv_hts_2d   , ONLY: hlr_2f,hgr_2f,hstr_2f,hspr_2f,tlr_2f,tgr_2f,tstr_2f,tspr_2f, &
                               t_fuel
      USE Zrv_ncell    , ONLY: ncell_fuel_rod,cupid_cell_hts2d,  &
                               channel_cell_hts2d
      !OPR1000 rod-scale                                
      USE Zmpi         , ONLY: ncell_fp,maxmt_fluid
      USE Zzone        , ONLY: ncell_fluid,nzone
      USE Zcore        , ONLY: np,myrank
      USE Zparam       , ONLY: pi
      USE Znum_cell    , ONLY: i_neigh,neigh
      USE Ztimecon     , ONLY: itim,time      
      USE zcoord2      , ONLY: cell_leng
      USE Zqvol        , ONLY: gamma_wall
      USE zcoord1      , ONLY: xloc
      USE Zvec_geo     , ONLY: xn_nf
      USE Zwall_HTC    , ONLY: gamma_wall_rod
      USE Zwall_HTC    , ONLY: qflux_l0,qflux_g0,qf0,qf1,qg0,qg1,gw0,gw1
      USE Zcoord3      , ONLY: volp
      USE Zrv_hts_2d   , ONLY: ri_2d
      USE Zporous  , ONLY: chn_type
      USE Zrv_hts_2d   , ONLY: power_2d   
      USE Zvector      , ONLY: ug_o
      USE Zconst2      , ONLY: hydraulicd   
      USE Zporous,only:qsum00
      USE Zio_unit     , ONLY: unit_log
      !
      IMPLICIT NONE
!
      INTEGER i,k,m,j0
!
      LOGICAL err
!
      REAL(8) s(36),qe,qc,ep
      REAL(8) x_flow,h_gas,h_mix
!.....Local arrays
      REAL(8),SAVE,ALLOCATABLE :: xn1(:),xn2(:)
      LOGICAL,SAVE :: init_xn=.true.
!
      !Additional connectivity - jrlee
      INTEGER m1,m2,m3,m4,j,sc
      REAL(8) ts,tst,p,quala,quals
      REAL(8) tg,hgsat,eg,rhog,alphag,vg
      REAL(8) tl,hlsat,el,rhol,alphal,vl
      REAL(8) rhom,hg
      REAL(8) ts0(ncell_fp),tst0(ncell_fp),p0(ncell_fp),quala0(ncell_fp),quals0(ncell_fp)
      REAL(8) tg0(ncell_fp),hgsat0(ncell_fp),rhog0(ncell_fp),alphag0(ncell_fp),vg0(ncell_fp)
      REAL(8) tl0(ncell_fp),hlsat0(ncell_fp),el0(ncell_fp),rhol0(ncell_fp),alphal0(ncell_fp),vl0(ncell_fp)
      !Single-phase HTC
      REAL(8) condg,visg,cpg,betag
      REAL(8) condl,visl,cpl,betal,vlf
      REAL(8) condg0(ncell_fp),visg0(ncell_fp),cpg0(ncell_fp),betag0(ncell_fp)
      REAL(8) condl0(ncell_fp),visl0(ncell_fp),cpl0(ncell_fp),betal0(ncell_fp)
      REAL(8) vl00(ncell_fp)
      !Nucleate boiling
      REAL(8) sigma0(ncell_fp)
      !Condensation
      REAL(8) hg0(ncell_fp)
      !film boiling refoold HTC
      REAL(8) hd,ug
      REAL(8) hd0(ncell_fp)
      REAL(8) ug0(ncell_fp)
!
      REAL(8) xn,yn

      real(8) delts
      REAL(8) qporous_liq_tot
      REAL(8) qporous_gas_tot
      REAL(8) mass_in2
      REAL(8) mass_in7
!      
      twall_rv(:,1)=0.d0
      twall_rv(:,2)=0.d0
      twall_rv(:,3)=0.d0
      DO i=1,ncell_fuel_rod        
         m=channel_cell_hts2d(i)
         twall_rv(m,1)=twall_rv(m,1)+t_fuel(i,1)
        !twall_rv(m,2)=twall_rv(m,2)+t_fuel(i,5)
         twall_rv(m,3)=twall_rv(m,3)+t_fuel(i,nr_2d)
      ENDDO

!.....pcrit is constant read in stread.f90 why change?
!     pcrit=22.4d6
      ep=1.d-10 
!
!.....Initialization of variables for Reflood calculation
!
      chf=0.d0
!
!.....Save wall tempertures calculated in Heat Structure module at the old time step. 
!     This variable is used only for Post-processing not for calculation.
!
!.....Main Loop for Wall HTC Calculation
!
!
      IF(np.gt.1) THEN
         CALL communicate_1d(cell%ts   , &
                             cell%tst  , &
                             cell%hgsat, &
                             cell%hlsat, &
                             cell%quals, &
                             cell%quala, &
                             cell%p    , &
                             cell%el      )
         CALL communicate_1d(cell%lviscosg, &
                             cell%lviscosl, &
                             cell%cpg     , &
                             cell%cpl     , &
                             cell%betag   , &
                             cell%betal   , &
                             cell%alphag  , &
                             cell%alphal     )
      !Single-phase HTC
         CALL communicate_1d(cell%lcondg  , &
                             cell%lcondl  , &
                             cell%rhog    , &
                             cell%rhol    , &
                             cell%sigma   , &  !sigma-> nucleat boiling
                             cell%hg      , &
                             cell%tl      , &
                             cell%tg         )  !hg   -> condensate
      !film boiling refoold HTC
         CALL communicate_1d(hydraulicd   , &
                             ug_o            )
      ENDIF

!
      DO m1=1,ncell_fp
         ts0(m1)    =0.d0
         tst0(m1)   =0.d0
         tg0(m1)    =0.d0
         tl0(m1)    =0.d0
         hgsat0(m1) =0.d0
         hlsat0(m1) =0.d0
         el0(m1)    =0.d0
         p0(m1)     =0.d0
         rhog0(m1)  =0.d0
         rhol0(m1)  =0.d0
         quala0(m1) =0.d0
         quals0(m1) =0.d0
         alphag0(m1)=0.d0
         alphal0(m1)=0.d0
         vg0(m1)    =0.d0
         vl0(m1)    =0.d0
         !Single-phase HTC
         condg0(m1) =0.d0
         condl0(m1) =0.d0
         visg0(m1)  =0.d0
         visl0(m1)  =0.d0
         cpg0(m1)   =0.d0
         cpl0(m1)   =0.d0
         betag0(m1) =0.d0
         betal0(m1) =0.d0
         vl00(m1)   =0.d0
         !Nucleate boiling
         sigma0(m1) =0.d0
         !Condensation
         hg0(m1)    =0.d0
         !Transition film refood
         hd0(m1)    =0.d0
         ug0(m1)    =0.d0
      ENDDO
!
!.....Get all  the xn scalar from xn_nf vector
!
      IF(init_xn)then
         init_xn=.false.
         ALLOCATE(xn1(maxmt_fluid),xn2(maxmt_fluid))
         xn1=0.d0
         xn2=0.d0
         CALL get_scalar_variable_n(xn_nf(1,1),xn1)
         CALL get_scalar_variable_n(xn_nf(1,2),xn2)
      ENDIF
!
!.....Vessel and/or core only
!     Averaging between i and k along x- direction
!
      DO m1=1,ncell_fluid
         IF(chn_type(m1).ne.0)then
!...........Get all the xn(j,1) for cell m1
!           CALL get_scalar_variable_n_i_ndim(xn_nf,xxn,m1,1)
            m2=0
            DO j=i_neigh(m1),i_neigh(m1+1)-1
               xn=xn1(j)
               IF(xn.gt.0.5d0) m2=neigh(j)
            ENDDO
            IF(m2.ne.0)then
               ts0(m1)    =(cell%ts(m1)     +cell%ts(m2)       )*0.5d0
               tst0(m1)   =(cell%tst(m1)    +cell%tst(m2)      )*0.5d0
               tg0(m1)    =(cell%tg_o(m1)   +cell%tg_o(m2)     )*0.5d0
               tl0(m1)    =(cell%tl_o(m1)   +cell%tl_o(m2)     )*0.5d0
               hgsat0(m1) =(cell%hgsat(m1)  +cell%hgsat(m2)    )*0.5d0
               hlsat0(m1) =(cell%hlsat(m1)  +cell%hlsat(m2)    )*0.5d0
               el0(m1)    =(cell%el(m1)     +cell%el(m2)       )*0.5d0
               p0(m1)     =(cell%p(m1)      +cell%p(m2)        )*0.5d0
               rhog0(m1)  =(cell%rhog(m1)   +cell%rhog(m2)     )*0.5d0
               rhol0(m1)  =(cell%rhol(m1)   +cell%rhol(m2)     )*0.5d0
               quala0(m1) =(cell%quala(m1)  +cell%quala(m2)    )*0.5d0
               quals0(m1) =(cell%quals(m1)  +cell%quals(m2)    )*0.5d0
               alphag0(m1)=(cell%alphag(m1) +cell%alphag(m2)   )*0.5d0
               alphal0(m1)=(cell%alphal(m1) +cell%alphal(m2)   )*0.5d0
               vg0(m1)    =(vg_n(m1,f_direc)+vg_n(m2,f_direc)  )*0.5d0
               vl0(m1)    =(vl_n(m1,f_direc)+vl_n(m2,f_direc)  )*0.5d0
               !Single-phase HTC
               condg0(m1) =(cell%lcondg(m1)  +cell%lcondg(m2)  )*0.5d0
               condl0(m1) =(cell%lcondl(m1)  +cell%lcondl(m2)  )*0.5d0
               visg0(m1)  =(cell%lviscosg(m1)+cell%lviscosg(m2))*0.5d0
               visl0(m1)  =(cell%lviscosl(m1)+cell%lviscosl(m2))*0.5d0
               cpg0(m1)   =(cell%cpg(m1)     +cell%cpg(m2)     )*0.5d0
               cpl0(m1)   =(cell%cpl(m1)     +cell%cpl(m2)     )*0.5d0
               betag0(m1) =(cell%betag(m1)   +cell%betag(m2)   )*0.5d0
               betal0(m1) =(cell%betal(m1)   +cell%betal(m2)   )*0.5d0
               vl00(m1)   =(vl_n(m1,3)       +vl_n(m2,3)       )*0.5d0
               !Nucleate boiling
               sigma0(m1) =(cell%sigma(m1)   +cell%sigma(m2)   )*0.5d0
               !Condensation 
               hg0(m1)    =(cell%hg(m1)      +cell%hg(m2)      )*0.5d0
               !Transition film refood
               hd0(m1)    =(hydraulicd(m1)   +hydraulicd(m2)   )*0.5d0
               ug0(m1)    =(ug_o(m1)         +ug_o(m2)         )*0.5d0
            ELSE
               ts0(m1)    =cell%ts(m1)     
               tst0(m1)   =cell%tst(m1)    
               tg0(m1)    =cell%tg_o(m1)   
               tl0(m1)    =cell%tl_o(m1)   
               hgsat0(m1) =cell%hgsat(m1)
               hlsat0(m1) =cell%hlsat(m1)  
               el0(m1)    =cell%el(m1)   
               p0(m1)     =cell%p(m1) 
               rhog0(m1)  =cell%rhog(m1) 
               rhol0(m1)  =cell%rhol(m1)
               quala0(m1) =cell%quala(m1) 
               quals0(m1) =cell%quals(m1) 
               alphag0(m1)=cell%alphag(m1)
               alphal0(m1)=cell%alphal(m1)
               vg0(m1)    =vg_n(m1,f_direc)
               vl0(m1)    =vl_n(m1,f_direc)
               !Single-phase HTC
               condg0(m1) =cell%lcondg(m1)
               condl0(m1) =cell%lcondl(m1) 
               visg0(m1)  =cell%lviscosg(m1)
               visl0(m1)  =cell%lviscosl(m1)
               cpg0(m1)   =cell%cpg(m1)     
               cpl0(m1)   =cell%cpl(m1) 
               betag0(m1) =cell%betag(m1)  
               betal0(m1) =cell%betal(m1) 
               vl00(m1)   =vl_n(m1,3)
               !Nucleate boiling
               sigma0(m1) =cell%sigma(m1)
               !Condensation 
               hg0(m1)    =cell%hg(m1)
               !Transition film refood
               hd0(m1)    =hydraulicd(m1)
               ug0(m1)    =ug_o(m1)
            ENDIF
         ELSE
            ts0(m1)    =cell%ts(m1)
            tst0(m1)   =cell%tst(m1)
            tg0(m1)    =cell%tg_o(m1)
            tl0(m1)    =cell%tl_o(m1)
            hgsat0(m1) =cell%hgsat(m1)
            hlsat0(m1) =cell%hlsat(m1)
            el0(m1)    =cell%el(m1)
            p0(m1)     =cell%p(m1)
            rhog0(m1)  =cell%rhog(m1)
            rhol0(m1)  =cell%rhol(m1)
            quala0(m1) =cell%quala(m1)
            quals0(m1) =cell%quals(m1)
            alphag0(m1)=cell%alphag(m1)
            alphal0(m1)=cell%alphal(m1)
            vg0(m1)    =vg_n(m1,f_direc)
            vl0(m1)    =vl_n(m1,f_direc)
            !Single-phase HTC
            condg0(m1) =cell%lcondg(m1)
            condl0(m1) =cell%lcondl(m1)
            visg0(m1)  =cell%lviscosg(m1)
            visl0(m1)  =cell%lviscosl(m1)
            cpg0(m1)   =cell%cpg(m1)
            cpl0(m1)   =cell%cpl(m1)
            betag0(m1) =cell%betag(m1)
            betal0(m1) =cell%betal(m1)
            vl00(m1)   =vl_n(m1,3)
            !Nucleate boiling
            sigma0(m1) =cell%sigma(m1)
            !Condensation
            hg0(m1)    =cell%hg(m1)
            !Transition film refood
            hd0(m1)    =hydraulicd(m1)
            ug0(m1)    =ug_o(m1)
         ENDIF
      ENDDO   

      IF(np.gt.1) THEN
         CALL communicate_1d(ts0   , &
                             tst0  , &
                             tg0   , &
                             tl0   , &
                             hgsat0, &
                             hlsat0, &
                             el0   , &
                             p0       )
         CALL communicate_1d(rhog0  , &
                             rhol0  , &
                             quala0 , &
                             quals0 , &
                             alphag0, &
                             alphal0, &
                             vg0    , &
                              vl0       )
      !Single-phase HTC
         CALL communicate_1d(condg0 , &
                             condl0 , &
                             visg0  , &
                             visl0  , &
                             cpg0   , &
                             cpl0   , &
                             betag0 , &
                             betal0    )
         CALL communicate_1d(vl00   , &
                             sigma0 , &
                             hg0    , &
                             hd0    , &
                             ug0       )
      ENDIF
                             
      do i=1,ncell_fuel_rod
         qflux_l0(i)=0.d0
         qflux_g0(i)=0.d0
         gamma_wall_rod(i)=0.d0
      enddo
      do i=1,ncell_fp
         qf0(i)=0.d0
         qf1(i)=0.d0
         qg0(i)=0.d0
         qg1(i)=0.d0
         gw0(i)=0.d0
         gw1(i)=0.d0
      enddo

      DO i=1,ncell_fuel_rod
         k=channel_cell_hts2d(i)
         m1=cupid_cell_hts2d(i)
         m=m1
         IF(chn_type(m1).eq.0)then
            write(*,*) 'm1 of fuel rod',i,'is',m1,'and nzone(m1)is',nzone(m1)
            stop 'err in rv_wall_HT_2d wow'
         ENDIF
         
!........skip guide tube region
         !IF(p3d_cupid(i).lt.1.d-3)then
         !   CYCLE
         !ENDIF


         ! Averaging i and k along y- direction
!........Get all the xn(j,2) for cell m1
!        CALL get_scalar_variable_n_i_ndim(xn_nf,yyn,m1,2)
         m4=0
         DO j=i_neigh(m1),i_neigh(m1+1)-1
            yn=xn2(j)
            IF(yn.lt.-0.5d0) then
              !j1=j-j0!
               m4=neigh(j)
            ENDIF
         ENDDO
         IF(m4.ne.0)then
            ts    =(ts0(m1)    +ts0(m4)    )*0.5d0
            tst   =(tst0(m1)   +tst0(m4)   )*0.5d0
            tg    =(tg0(m1)    +tg0(m4)    )*0.5d0
            tl    =(tl0(m1)    +tl0(m4)    )*0.5d0
            hgsat =(hgsat0(m1) +hgsat0(m4) )*0.5d0
            hlsat =(hlsat0(m1) +hlsat0(m4) )*0.5d0
            el    =(el0(m1)    +el0(m4)    )*0.5d0
            p     =(p0(m1)     +p0(m4)     )*0.5d0
            rhog  =(rhog0(m1)  +rhog0(m4)  )*0.5d0
            rhol  =(rhol0(m1)  +rhol0(m4)  )*0.5d0
            quala =(quala0(m1) +quala0(m4) )*0.5d0
            quals =(quals0(m1) +quals0(m4) )*0.5d0
            alphag=(alphag0(m1)+alphag0(m4))*0.5d0
            alphal=(alphal0(m1)+alphal0(m4))*0.5d0
            vg    =(vg0(m1)    +vg0(m4)    )*0.5d0
            vl    =(vl0(m1)    +vl0(m4)    )*0.5d0
            !Single-phase HTC
            condg =(condg0(m1) +condg0(m4) )*0.5d0
            condl =(condl0(m1) +condl0(m4) )*0.5d0
            visg  =(visg0(m1)  +visg0(m4)  )*0.5d0
            visl  =(visl0(m1)  +visl0(m4)  )*0.5d0
            cpg   =(cpg0(m1)   +cpg0(m4)   )*0.5d0
            cpl   =(cpl0(m1)   +cpl0(m4)   )*0.5d0
            betag =(betag0(m1) +betag0(m4) )*0.5d0
            betal =(betal0(m1) +betal0(m4) )*0.5d0
            vlf   =(vl00(m1)   +vl00(m4)   )*0.5d0
            !Nucleate boiling
            sigma =(sigma0(m1) +sigma0(m4) )*0.5d0
            !Condensation
            hg    =(hg0(m1)    +hg0(m4)    )*0.5d0
            !Transition film refood
            hd    =(hd0(m1)   +hd0(m4)     )*0.5d0
            ug    =(ug0(m1)   +ug0(m4)     )*0.5d0
         ELSE
            ts    =ts0(m1)    
            tst   =tst0(m1)   
            tg    =tg0(m1)    
            tl    =tl0(m1)    
            hgsat =hgsat0(m1) 
            hlsat =hlsat0(m1) 
            el    =el0(m1)    
            p     =p0(m1)     
            rhog  =rhog0(m1)  
            rhol  =rhol0(m1)  
            quala =quala0(m1) 
            quals =quals0(m1) 
            alphag=alphag0(m1)
            alphal=alphal0(m1)
            vg    =vg0(m1)    
            vl    =vl0(m1)    
            !Singlephase HTC
            condg =condg0(m1) 
            condl =condl0(m1) 
            visg  =visg0(m1)  
            visl  =visl0(m1)  
            cpg   =cpg0(m1)   
            cpl   =cpl0(m1)   
            betag =betag0(m1) 
            betal =betal0(m1) 
            vlf   =vl00(m1)   
            !Nucleate boiling
            sigma =sigma0(m1)
            !Condensation
            hg    =hg0(m1)
            !Transition film refood
            hd    =hd0(m1)
            ug    =ug0(m1)
         ENDIF
!
!........Load wall temp. from heat structure
!
         tw=t_fuel(i,nr_2d)    
         delts=max(0.d0,ts-tw)
! 
!........Set bagic properties
!             
         tsat_t=tst
         dt_sat=tw-tsat_t
         qe=0.d0
         qc=0.d0
!
         hfg=hgsat-hlsat
         sat_hfp=hlsat
         h_liq=el+p/rhol
!
!........Calculate thermal equilibrium quality(qual_eq), liquid enthalpy  on total pressure (sat_hfp), heat of vaporization on total pressure (hfg_p)
!
         IF(quala.gt.1.d-9)THEN
            s(2)=p
            IF(s(2).gt.pcrit)THEN
               WRITE(* ,*) 'Pressure exceeds saturation bound'
               WRITE(unit_log,*) 'Pressure exceeds saturation bound'
               s(2)=MIN(pcrit,s(2))
            ENDIF
            s(9)=0.d0
            CALL sth2x2_cupid(st_tbl(ndxstd),st_tbl(ndxstd+nt),s,err)
            hfg_p=s(16)-s(15)
            sat_hfp=s(15)
            rhom=alphag*rhog+alphal*rhol
            x_flow=(mflux_gasa+0.01*alphag*rhog)/(mflux_liqa+mflux_gasa+0.01*rhom)
            h_gas=eg+p/rhog
            IF(s(2).lt.pcrit)THEN
               h_mix=h_liq+x_flow*(h_gas-h_liq)
               qual_eq=(h_mix-sat_hfp)/hfg_p
            ELSE
               qual_eq=quals
            ENDIF         
         ELSE
            hfg_p=hfg
            qual_eq=quals
         ENDIF
!
!........Set mass flux and relative velocity in major direction (1-D base)
!   
         mflux_gasa=MAX(ABS(vg*rhog*alphag),0.001d0) 
         mflux_liqa=MAX(ABS(vl*rhol*alphal),0.001d0) 
         mflux_tota=mflux_gasa+mflux_liqa 
         vfg=max(0.001d0,ABS(vg-vl))
!
!........Calculate 'Reflood parameters'
!
!-----------------------------------------------------------------------------------------      
!      
!........Mode Selection
!
         IF(quala.gt.0.99999999d0)THEN
            IF(myrank.eq.0) write(*,*) 'quala.gt.0.999'
            mode(k)=0                                 ! Air-Water
            CALL single_phase_HTC_rod(i,mode(k),alphag,alphal,rhog,rhol,tg,tl,vg,vl,condg,condl,visg,visl,cpg,cpl,betag,betal)
         ELSEIF(p.gt.pcrit)THEN
            IF(myrank.eq.0) write(*,*) 'p.gt.pcrit'
            mode(k)=10                                ! Critical fluid
            CALL single_phase_HTC_rod(i,mode(k),alphag,alphal,rhog,rhol,tg,tl,vg,vl,condg,condl,visg,visl,cpg,cpl,betag,betal)
         ELSEIF(tw.lt.ts-0.001d0)THEN !original-cyj
            IF(alphag.lt.0.1d0)THEN
               mode(k)=2                              ! Liquid 1-phase
               CALL single_phase_HTC_rod(i,mode(k),alphag,alphal,rhog,rhol,tg,tl,vg,vl,condg,condl,visg,visl,cpg,cpl,betag,betal)
            ELSEIF(quala.gt.0.999d0)THEN
               IF(myrank.eq.0) write(*,*) 'quala.gt.0.999d0'
               mode(k)=0                              ! Air-Water
               CALL single_phase_HTC_rod(i,mode(k),alphag,alphal,rhog,rhol,tg,tl,vg,vl,condg,condl,visg,visl,cpg,cpl,betag,betal)
            ELSEIF(tw.gt.tl.and.alphag.lt.0.999d0)THEN
               mode(k)=2                              ! Liquid 1-phase
               CALL single_phase_HTC_rod(i,mode(k),alphag,alphal,rhog,rhol,tg,tl,vg,vl,condg,condl,visg,visl,cpg,cpl,betag,betal)
            ELSE
               IF(myrank.eq.0) write(*,*) 'condensation_HTC'
               mode(k)=11                             ! Condensation    
               CALL condensation_HTC_rod(i,mode(k),ts,quala,hg,hlsat,qc,                                                  &
                                         p,alphag,alphal,rhog,rhol,tg,tl,vg,vl,condg,condl,visg,visl,cpg,cpl,betag,betal)
            ENDIF
         ELSEIF(dt_sat.le.0.d0)THEN
            mode(k)=0                                  ! Air-Water
            CALL single_phase_HTC_rod(i,mode(k),alphag,alphal,rhog,rhol,tg,tl,vg,vl,condg,condl,visg,visl,cpg,cpl,betag,betal)

         ELSEIF(tw.lt.tl)THEN
            IF(myrank.eq.0) write(*,*) 'tw.lt.tl'
            mode(k)=0                                  ! Air-Water
            CALL single_phase_HTC_rod(i,mode(k),alphag,alphal,rhog,rhol,tg,tl,vg,vl,condg,condl,visg,visl,cpg,cpl,betag,betal)
         ELSE
            IF(alphag.ge.0.999d0)THEN
               IF(myrank.eq.0) write(*,*) 'Gas 1-phase'
               mode(k)=9                                ! Gas 1-phase
               CALL single_phase_HTC_rod(i,mode(k),alphag,alphal,rhog,rhol,tg,tl,vg,vl,condg,condl,visg,visl,cpg,cpl,betag,betal) 
            ELSEIF(dt_sat.gt.600.d0)THEN
               IF(myrank.eq.0) write(*,*) 'Film boiling'
               !mode(i)=7~8                            ! Film boiling
               CALL CHF_calc(m)
               IF(reflood.eq.1)THEN
                  CALL trans_film_reflood_HTC_rod(i,mode(k),quala,hd,ug,p,alphag,alphal,rhog,rhol,tg,tl,vg,vl,condg,condl,visg,visl,cpg,cpl,betag,betal)
               ELSE
                  CALL trans_film_boiling_HTC_rod(i,mode(k),quala,p,alphag,alphal,rhog,rhol,tg,tl,vg,vl,condg,condl,visg,visl,cpg,cpl,betag,betal)
               ENDIF 
              !CALL subcooled_boiling_rod(i,mul_o(k),qe,p,rhog,rhol,alphal)
               CALL subcooled_boiling_rod_sahazuber(i,mul_o(k),qe,p,rhog,rhol,alphal)
            ELSEIF(dt_sat.gt.100.d0)THEN 
               IF(myrank.eq.0) write(*,*) 'Transient or Film boiling'
               !mode(i)=5~8                              ! Transient or Film boiling
               CALL CHF_calc(m)
               IF(reflood.eq.1)THEN
                  CALL trans_film_reflood_HTC_rod(i,mode(k),quala,hd,ug,p,alphag,alphal,rhog,rhol,tg,tl,vg,vl,condg,condl,visg,visl,cpg,cpl,betag,betal)
               ELSE
                  CALL trans_film_boiling_HTC_rod(i,mode(k),quala,p,alphag,alphal,rhog,rhol,tg,tl,vg,vl,condg,condl,visg,visl,cpg,cpl,betag,betal)
               ENDIF
              !CALL subcooled_boiling_rod(i,mul_o(k),qe,p,rhog,rhol,alphal)
               CALL subcooled_boiling_rod_sahazuber(i,mul_o(k),qe,p,rhog,rhol,alphal)
            ELSE
               CALL CHF_calc(m)
               CALL nucl_boiling_HTC_rod(i,mode(k),p,alphag,alphal,rhog,rhol,tg,tl,vg,vl,condg,condl,visg,visl,cpg,cpl,betag,betal)
            !!!NEXT!!!
               chf=1.d20
               IF(qflux_t.ge.chf)THEN
                  IF(myrank.eq.0) then
                     write(*,*) 'qflux_t is ', qflux_t
                     write(*,*) 'mode(i)=5~8'
                  ENDIF
                  !mode(i)=5~8                            ! Transient or Film boiling
                  IF(reflood.eq.1)THEN
                     CALL trans_film_reflood_HTC_rod(i,mode(k),quala,hd,ug,p,alphag,alphal,rhog,rhol,tg,tl,vg,vl,condg,condl,visg,visl,cpg,cpl,betag,betal)
                  ELSE
                     CALL trans_film_boiling_HTC_rod(i,mode(k),quala,p,alphag,alphal,rhog,rhol,tg,tl,vg,vl,condg,condl,visg,visl,cpg,cpl,betag,betal)
                  ENDIF 
                 !CALL subcooled_boiling_rod(i,mul_o(k),qe,p,rhog,rhol,alphal)
                  CALL subcooled_boiling_rod_sahazuber(i,mul_o(k),qe,p,rhog,rhol,alphal)
               ELSEIF(qflux_t.gt.0.d0)THEN
                  !mode(i)=3,4
                  IF(tl.lt.tst)THEN
                     mode(k)=3                            ! Subcooled Boiling
                  ELSE
                     mode(k)=4                            ! Nucleate Boling
                  ENDIF
                 !CALL subcooled_boiling_rod(i,mul_o(k),qe,p,rhog,rhol,alphal)
                  CALL subcooled_boiling_rod_sahazuber(i,mul_o(k),qe,p,rhog,rhol,alphal)
               ELSE
                  mode(k)=1                               ! Do nothing in MARS!!!
               ENDIF 
!
            ENDIF
         ENDIF !End of mode selection
!
!........Energy Partitioning !!!cyj
!                 
         qflux_l0(i)=qflux_l-qe
         qflux_g0(i)=qflux_g+qc
         
         IF(chf.ne.0) chfr(k)=qflux_t/chf     
!       
!........Calculate 'Quenching Front,QF'
!
!
!........Save calculated HTC in HS array
!
         hlr_2f(i)=HTC_tl
         hgr_2f(i)=HTC_tg
         hstr_2f(i)=HTC_tst
         hspr_2f(i)=HTC_tgp
!
         tlr_2f(i)=tl !cell%tl(m)
         tgr_2f(i)=tg !cell%tg(m)
         tstr_2f(i)=tst !cell%tst(m)
         tspr_2f(i)=ts !cell%ts(m)
!
         dt_sat_rv(k)=dt_sat
         chf_rv(k)=chf                
!
         qf0(m1)=qflux_l0(i)*2.d0*pi*ri_2d(nr_2d)*cell_leng(m1,3)*0.25d0
         qg0(m1)=qflux_g0(i)*2.d0*pi*ri_2d(nr_2d)*cell_leng(m1,3)*0.25d0
         gw0(m1)=gamma_wall_rod(i)
         
      ENDDO !i=1,ncell_fuel_rod
      
      IF(np.gt.1) CALL communicate_1d(qf0, &
                                      qg0, &
                                      gw0)

      do m1=1,ncell_fluid
         IF(chn_type(m1).eq.0)cycle
!........Get all the xn(j,1) for cell m1
!        CALL get_scalar_variable_n_i_ndim(xn_nf,xxn,m1,1)
         m3=0
         DO j=i_neigh(m1),i_neigh(m1+1)-1
            xn=xn1(j)
            IF(xn .lt. -0.5d0) m3=neigh(j)
         enddo
         IF(m3.gt.0)then
            qf1(m1)=qf0(m3)
            qg1(m1)=qg0(m3)
            gw1(m1)=gw0(m3)
         ELSEIF(m3.eq.0)then
            qf1(m1)=0.d0
            qg1(m1)=0.d0
            gw1(m1)=0.d0
         ELSE
            STOP 'm3 is (-) in rv_wall_ht_rod'
         ENDIF
      enddo
      IF(np.gt.1) CALL communicate_1d(qf1, &
                                      qg1, &
                                      gw1)

      DO m=1,ncell_fluid
         IF(chn_type(m).eq.0)cycle
!........Get all the xn(j,1),xn(j,2) for cell m
!        CALL get_scalar_variable_n_i_ndim(xn_nf,xxn,m,1)
!        CALL get_scalar_variable_n_i_ndim(xn_nf,yyn,m,2)
!
!!!!!!!!!???? what if m1,m2,m3,m4 is never set in j loop???????
!
         m1=0
         m2=0
         m3=0
         m4=0
         DO j=i_neigh(m),i_neigh(m+1)-1
            xn=xn1(j)
            yn=xn2(j)
            IF(xn .gt. 0.5d0) m1=neigh(j)
            IF(yn .gt. 0.5d0) m2=neigh(j)
            IF(xn .lt.-0.5d0) m3=neigh(j)
            IF(yn .lt.-0.5d0) m4=neigh(j)
            IF(m1.lt.0 .or. m2.lt.0 .or. m3.lt.0 .or. m4.lt.0)then
               STOP 'm1~m4 are (-) in rv_wall_ht_rod'
            ENDIF
         ENDDO
         sc=chn_type(m)
         IF(sc.eq.3)then
            ! bottom wall
            IF(m4.eq.0)then
               ! bottom-left corner
               IF(m3.eq.0 .and. m1.ne.0 .and. m2.ne.0)then
                  qporous_liq(m)=qf0(m2)
                  qporous_gas(m)=qg0(m2)
                  gamma_wall(m) =gw0(m2)
               ! bottom-right corner
               ELSEIF(m3.ne.0 .and. m1.eq.0 .and. m2.ne.0)then
                  qporous_liq(m)=qf1(m2)
                  qporous_gas(m)=qg1(m2)
                  gamma_wall(m) =gw1(m2)
               ! bottom-middle plane   
               ELSEIF(m3.ne.0 .and. m1.ne.0 .and. m2.ne.0)then
                  qporous_liq(m)=qf0(m2)+qf1(m2)
                  qporous_gas(m)=qg0(m2)+qg1(m2)
                  gamma_wall(m) =gw0(m2)+gw1(m2)
               ENDIF
            ! top wall   
            ELSEIF(m2.eq.0)then
               ! top-left corner
               IF(m3.eq.0 .and. m1.ne.0 .and. m4.ne.0)then
                  qporous_liq(m)=qf0(m)
                  qporous_gas(m)=qg0(m)
                  gamma_wall(m) =gw0(m)
               ! top-right corner   
               ELSEIF(m1.eq.0 .and. m3.ne.0 .and. m4.ne.0)then
                  qporous_liq(m)=qf0(m3)
                  qporous_gas(m)=qg0(m3)
                  gamma_wall(m) =gw0(m3)
               ! top-middle plane   
               ELSEIF(m3.ne.0 .and. m1.ne.0 .and. m4.ne.0)then
                  qporous_liq(m)=qf0(m)+qf0(m3)
                  qporous_gas(m)=qg0(m)+qg0(m3)
                  gamma_wall(m) =gw0(m)+gw0(m3)
               ENDIF
            ! left wall   
            ELSEIF(m3.eq.0)then
               ! left-middel plane
               IF(m1.ne.0 .and. m2.ne.0 .and. m4.ne.0)then
                  qporous_liq(m)=qf0(m)+qf0(m2)
                  qporous_gas(m)=qg0(m)+qg0(m2)
                  gamma_wall(m) =gw0(m)+gw0(m2)
               ENDIF
            ! right wall   
            ELSEIF(m1.eq.0)then
               ! right-middle plane
               IF(m2.ne.0 .and. m3.ne.0 .and. m4.ne.0)then
                  qporous_liq(m)=qf0(m3)+qf1(m2)
                  qporous_gas(m)=qg0(m3)+qg1(m2)
                  gamma_wall(m) =gw0(m3)+gw1(m2)
               ENDIF
            ! inside subchannel   
            ELSE
               IF(m1.eq.0 .or. m2.eq.0 .or. m3.eq.0 .or. m4.eq.0)then
                  STOP 'check sc=3 in qvol_liq'
               ENDIF
               qporous_liq(m)=qf0(m)+qf0(m2)+qf0(m3)+qf1(m2)
               qporous_gas(m)=qg0(m)+qg0(m2)+qg0(m3)+qg1(m2)
               gamma_wall(m) =gw0(m)+gw0(m2)+gw0(m3)+gw1(m2)
            ENDIF
         ELSEIF(sc.eq.2)then
            ! bottom
            IF(m4.eq.0)then
               qporous_liq(m)=qf0(m2)+qf1(m2)
               qporous_gas(m)=qg0(m2)+qg1(m2)
               gamma_wall(m) =gw0(m2)+gw1(m2)
            ! top   
            ELSEIF(m2.eq.0)then
               qporous_liq(m)=qf0(m)+qf0(m3)
               qporous_gas(m)=qg0(m)+qg0(m3)
               gamma_wall(m) =gw0(m)+gw0(m3)
            ! right   
            ELSEIF(m1.eq.0)then
               qporous_liq(m)=qf0(m3)+qf1(m2)
               qporous_gas(m)=qg0(m3)+qg1(m2)
               gamma_wall(m) =gw0(m3)+gw1(m2)
            ! left   
            ELSEIF(m3.eq.0)then
               qporous_liq(m)=qf0(m)+qf0(m2)
               qporous_gas(m)=qg0(m)+qg0(m2)
               gamma_wall(m) =gw0(m)+gw0(m2)
            ! inside subchannel   
            ELSE
               IF(m1.eq.0 .or. m2.eq.0 .or. m3.eq.0 .or. m4.eq.0)then
                  STOP 'check sc=2 in qvol_liq'
               ENDIF
               qporous_liq(m)=qf0(m)+qf0(m2)+qf0(m3)+qf1(m2)
               qporous_gas(m)=qg0(m)+qg0(m2)+qg0(m3)+qg1(m2)
               gamma_wall(m) =gw0(m)+gw0(m2)+gw0(m3)+gw1(m2)
            ENDIF
         ELSE
            IF(sc.eq.0 .or. sc.gt.6)STOP 'check other sc in qvol_liq'
            qporous_liq(m)=qf0(m)+qf0(m2)+qf0(m3)+qf1(m2)
            qporous_gas(m)=qg0(m)+qg0(m2)+qg0(m3)+qg1(m2)
            gamma_wall(m) =gw0(m)+gw0(m2)+gw0(m3)+gw1(m2)
         ENDIF
         gamma_wall(m)=gamma_wall(m)/volp(m)         

      ENDDO 
!
!.....Total Energy check
      IF(MOD(itim,10).eq.0)then
         qporous_liq_tot=0.d0
         qporous_gas_tot=0.d0
         DO i=1,ncell_fluid
            IF(chn_type(i).ne.0)then
               qporous_liq_tot=qporous_liq_tot+qporous_liq(i)
               qporous_gas_tot=qporous_gas_tot+qporous_gas(i)
            ENDIF
         ENDDO
         IF(np.gt.1) THEN
            CALL allreducei_r1(qporous_liq_tot)
            CALL allreducei_r1(qporous_gas_tot)
         ENDIF
         IF(myrank.eq.0)write(5601,*)time,qporous_liq_tot/1.d6,power_2d/1.d6
         45 format(a,3(f10.2,1x))
         qsum00=qporous_liq_tot
      ENDIF
!
!.....Inlet Mass check
if(0)then
      IF(MOD(itim,10).eq.0)then
         mass_in2=0.d0
         mass_in7=0.d0
         DO i=1,ncell_fluid
            IF(chn_type(i).ne.0 .and. xloc(i,3).lt.0.3d0)then
               j0=i_neigh(i)-1
               DO j=i_neigh(i),i_neigh(i+1)-1
!                  IF(xn(j,i,3).lt.-0.5d0)then
!!                   mass_in2=mass_in2+cell%rhol(i)*sa(j,i)*dsqrt(vl_n(i,1)*vl_n(i,1)+vl_n(i,2)*vl_n(i,2)+vl_n(i,3)*vl_n(i,3))
!                  ENDIF
               ENDDO
            ENDIF
            IF(nzone(i).eq.7 .and. xloc(i,3).lt.0.3d0)then
               j0=i_neigh(i)-1
               DO j=i_neigh(i),i_neigh(i+1)-1
!                  IF(xn(j,i,3).lt.-0.5d0)then
!!                   mass_in7=mass_in7+cell%rhol(i)*sa(j,i)*dsqrt(vl_n(i,1)*vl_n(i,1)+vl_n(i,2)*vl_n(i,2)+vl_n(i,3)*vl_n(i,3))
!                  ENDIF
               ENDDO
            ENDIF
         ENDDO
         call allreducei_r1(mass_in2)
         call allreducei_r1(mass_in7)
         IF(myrank.eq.0)write(5700,45)'mass in zone2/zone7 is',mass_in2,mass_in7,time
         IF(myrank.eq.0)write(5701,*)time,mass_in2,mass_in7
         !call barrier_mpi
      ENDIF
endif
!
      END SUBROUTINE rv_wall_ht_2d_rod
