!
      SUBROUTINE choke_model
!
!     Henry-Fauske Critical Flow Model       
!
      USE Zcore           , ONLY: myrank 
      USE STM_TBL_cupid   , ONLY: st_tbl,       &
                                  nt,ndxstd,    &
                                  pcrit,        &
                                  cvao,uao,dcva
      USE Zconst2         , ONLY: dt
      USE Zrv_choke       , ONLY: choke,choke_pout,env_press_option,                                             &
                                  cpgas,cppf0,cppg0,gamma,ploss,pvzero,                                          &
                                  pzero,rgas,rnc,sgas,sliq,svap,tzero,vsubf0,vsubg0,                             &
                                  xnc,xzero,                                                                     &
                                  vl_choke,vg_choke,vl_choke_o,vg_choke_o,                                       &
                                  cell_leng_avg,ar_liq_avg,ar_gas_avg,p_avg,p_avg_out,pps_avg,vl_o_avg,vg_o_avg, &
                                  rhog_avg,eg_avg,el_avg,quala_avg,quals_avg,tl_avg,alphag_avg,alphal_avg,       &
                                  cvao_cell_avg,dcva_cell_avg,uao_cell_avg,ra_cell_avg,theta_avg,                &
                                  num_throatface,relax_choke,                                                    &
                                  choke_throat_area,throat_area,n_face_throat
      USE Zvec_major      , ONLY: mflux_l_nf,mflux_g_nf       
!
      IMPLICIT NONE
     
!
!.....Local variables for Henry-Fauske model.                              
!                              
      INTEGER :: i,i1       
      REAL(8) :: cpmix,cvgas,cvmix,cvvap,delt,dgcdp,               &
                 dvcdp,expn,gcrit,gmom,hgas,hfact,hmix,pbar,pgas,  &
                 treff,ugas,vcritt,                                &
                 xncx0,choke_throat_ratio
!                                                                       
!.....Local variables.                                                     
!
      REAL(8) :: denomhf 
      REAL(8) :: avrff,avrgg,avrho,ddx,signvc
      REAL(8) :: ttt,press,vbarr,ubarr,hbarr,betaa,kpa,cpp,quall,psatt,vsubff,vsubgg, &
                 usubff,usubgg,hsubf,hsubg,betf,betg,kpaf,kpag,cppf,cppg
      REAL(8) :: prop(36)
      REAL(8) :: gcrit_vap,gcrit_liq,dgcdp_vap,dgcdp_liq
      REAL(8) :: discharge
      REAL(8),SAVE :: gcrit0,dgcdp0
      REAL(8) :: choke_env
      REAL(8) :: mflux_l_avg,mflux_g_avg      
!      
      LOGICAL :: err
!
!     DATA pcrit   /22.4d6/
      DATA pbar    /1.d5/
      DATA treff   /273.15d0/
!                                                                       
!.....State properties                                                     
!
      EQUIVALENCE(prop( 1),ttt),     &
                 (prop( 2),press),   &
                 (prop( 3),vbarr),   &
                 (prop( 4),ubarr),   &
                 (prop( 5),hbarr),   &
                 (prop( 6),betaa),   &
                 (prop( 7),kpa),     &
                 (prop( 8),cpp),     &
                 (prop( 9),quall),   &
                 (prop(10),psatt),   &
                 (prop(11),vsubff),  &
                 (prop(12),vsubgg),  &
                 (prop(13),usubff),  &
                 (prop(14),usubgg),  &
                 (prop(15),hsubf),   &
                 (prop(16),hsubg),   &
                 (prop(17),betf),    &
                 (prop(18),betg),    &
                 (prop(19),kpaf),    &
                 (prop(20),kpag),    &
                 (prop(21),cppf),    &
                 (prop(22),cppg)
!
!.....No choke condition
!
      IF(vl_choke*vg_choke.lt.0.d0) RETURN
      IF(num_throatface.eq.0) RETURN
      IF(throat_area.eq.0.d0) RETURN    
!
!.....Set environment pressure
!      
      IF(env_press_option.eq.0) THEN
         choke_env=choke_pout
      ELSEIF(env_press_option.eq.1) THEN   
         choke_env=p_avg_out
      ENDIF
!
!.....Set choke_thrat_ratio
!            
      IF(choke_throat_area.ne.100.d0) THEN
         choke_throat_ratio=choke_throat_area/throat_area
      ELSE
         choke_throat_ratio=1.d0 
      ENDIF
!
!.....pcrit is constant read in stread.f90 why change?
!     pcrit=22.4d6
!
!.....Initialize 
!
      gcrit0=0.d0
      dgcdp0=0.d0
!         
      signvc=DSIGN(1.d0,vl_choke) 
      ddx=0.5d0*cell_leng_avg
!                                                                       
!.....Define mixture density & alfa*rho                               
!                                                                       
      avrff=ar_liq_avg
      avrgg=ar_gas_avg  
      avrho=avrff+avrgg
!                                                                       
!.....Set Stagnation Pressure to Cell-Center Value:                   
!                                                                       
      pzero=p_avg
!                                                                       
!.....Calculate "Pressure Losses" to Cell-Edge:                       
!     1) Add (or subtract) cell center momentum flux.              
!                                                                       
      ploss=0.5d0*signvc*(avrff*DABS(vl_o_avg)*vl_o_avg+avrgg*DABS(vg_o_avg)*vg_o_avg)
!                                                                       
!     2) Subtract pressure drop due to wall friction. 
!                                                                     
!     ploss=ploss-signvc*ddx*(fwalf(ick)*vl_choke+fwalg(ick)*vg_choke)    !fwalf=G/vol[kg/s/m3]
!                                                                       
!     3) subtract gravitational head for 1/2 cell.                 
!                                                                       
      ploss=ploss-9.81d0*ddx*avrho*DSIN(theta_avg)
!                                                                       
!.....No choke condition 2
!                    
      IF((pzero+MIN(0.d0,ploss)).lt.1.01d0*choke_env) THEN
         RETURN !no choking  !1.01*pzero needs be check later.
      ENDIF   
!                                                                       
!.....Saturation Properties at Vapor Partial Pressure                 
!
      pvzero=pps_avg
!         
!.....Check for supercritical pressure case                                
!                                                                     
      IF(pvzero.ge.pcrit)then 
         gamma=1.3d0 
         expn=(gamma+1.d0)/(gamma-1.d0) 
         gcrit=DSQRT(gamma*rhog_avg*(pzero+ploss)*(2.d0/(gamma+1.d0))**expn)
         dgcdp=gcrit/pzero 
         GOTO 160 
      ENDIF
!
!.....Sat. liquid condition (x2 stb)
!
      prop(2)=DMAX1(611.24d0,pvzero)
      prop(9)=0.d0
      CALL sth2x2_cupid(st_tbl(ndxstd),st_tbl(ndxstd+nt),prop,err)
!
      cppf0=cppf 
      cppg0=cppg 
      sliq=prop(25) 
      svap=prop(26) 
      vsubf0=vsubff 
      vsubg0=vsubgg 
      IF(avrho.lt.1.d-10) THEN
         RETURN  !no choking
      ENDIF   
      hmix=(avrgg*eg_avg+avrff*el_avg+pzero)/avrho
!
!.....Define Property Info With & Without NC Gas.                     
!               
      IF(quala_avg.lt.1.d-6) THEN
!                                                                     
!        1) Gas Phase is Pure Steam.                                     
!                                                                     
         xnc=0.d0
         gamma=1.3d0 
         cvao=715.d0 
         dcva=0.10d0 
         rnc=287.d0 
         cpgas=cppg0 
         rgas=0.d0
         sgas=svap 
!            
         xzero=(hmix-hsubf)/(hsubg-hsubf) 
!                                                                       
      ELSE 
!                                                                       
!        2) Gas Phase is Mixture of Steam and NC Gas.  
!                 
!........Set NC gas property constants:                            
!                                                                       
         cvao=cvao_cell_avg
         dcva=dcva_cell_avg
         uao=uao_cell_avg
         rnc=ra_cell_avg
!                                                                       
!........NC Gas Properties at Tsat(PVZERO).                        
!                                                                       
         delt=DMAX1(0.d0,ttt-250.d0) 
         ugas=cvao*ttt+0.5d0*dcva*delt**2+uao 
         pgas=DMAX1(612.d0,pzero-pvzero) 
         rgas=pgas/(rnc*ttt) 
         hgas=ugas+pgas/rgas 
         cvgas=cvao+dcva*delt 
         cpgas=cvgas+rnc 
         sgas=cpgas*DLOG(ttt/treff)-rnc*DLOG(pgas/pbar) 
!                                                                       
!........Define "stagnation" quality: use isenthalpic process: NC Mass Fraction remains constant.                           
!                                                                     
         xncx0=quals_avg*quala_avg
!                                                                       
!........Compute vapor mass fraction (store in XZERO for now).  
!                                                                       
         xzero=(hmix-xncx0*hgas-(1.d0-xncx0)*hsubf)/(hsubg-hsubf) 
         xzero=DMAX1(0.d0,xzero) 
!                                                                       
!........Correct Quality by adding NC mass fraction.            
!                                                                       
         xzero=xzero+xncx0 
         xzero=DMIN1(1.d0,xzero) 
!                                                                       
!........Backout NC Quality.                                    
!                                                                       
         xnc=xncx0/DMAX1(1.d-6,xzero) 
!                                                                       
!........Set Isentropic Expansion Coefficient:                     
!                                                                      
         cvvap=cppg/1.3d0 
         cvmix=quala_avg*cvgas+(1.d0-quala_avg)*cvvap 
         cpmix=quala_avg*cpgas+(1.d0-quala_avg)*cppg 
!                                                                       
         gamma=cpmix/cvmix 
            
      ENDIF
!                                                                       
!.....Set Non-Equilibrium Parameter:                                     
!                                                                       
      hfact=0.14d0
!
      IF(xzero.lt.1.d-6) THEN
!
!........Subcooled Liquid Critical Flow:                     
!
         tzero=DMIN1(tl_avg,ttt-0.001d0)
         CALL choke_gcsub(gcrit,dgcdp,hfact,err)
!
      ELSEIF(xzero.lt.0.998d0) THEN
         IF(quala_avg.ge.1.d-6)then 
!
!...........Single-Phase Vapor:          
!               
            tzero=DMIN1(tl_avg,ttt-0.001d0) 
            CALL choke_gcsub(gcrit_liq,dgcdp_liq,hfact,err)
!                                                          
            expn=(gamma+1.d0)/(gamma-1.d0)
            gcrit_vap=DSQRT(gamma*rhog_avg*(pzero+ploss)*(2.d0/(gamma+1.d0))**expn)
            dgcdp_vap=gcrit_vap/pzero 
!
            gcrit=gcrit_vap*alphag_avg+gcrit_liq*alphal_avg
            dgcdp=dgcdp_vap*alphag_avg+dgcdp_liq*alphal_avg
         ELSE
!
!...........Two-Phase Critical Flow:
!                           
            CALL choke_gctpm(gcrit,dgcdp,hfact,err) 
!                                                                       
            IF(err)then 
!               
!..............GCTPM did not converge, use default.                   
!
!               gcrit=avrho*v_da(kk)%sounde 
!               dgcdp=avrho * 0.15d0/(v_da(kk)%sounde*v_da(kk)%rho) 
               gcrit=gcrit0  !modified by LSJ due to no sounde calculation in H.F. model
               dgcdp=dgcdp0  !modified by LSJ due to no sounde calculation in H.F. model
!
               PRINT*,'warning: gctpm does not converged'                        
!                  
            ENDIF 
         ENDIF

      ELSE 
!                                                                       
!........Single-Phase Vapor:                                    
!                                                                       
         expn=(gamma+1.d0)/(gamma-1.d0) 
         gcrit=Dsqrt(gamma*rhog_avg*(pzero+ploss)*(2.d0/(gamma+1.d0))**expn)
         dgcdp=gcrit/pzero
!                                                                       
      ENDIF 
!
160   CONTINUE  !skipped to by supercritical pressure condition
!                                                                       
!.....Apply discharge coefficient: need to be modified later. 
!             
      discharge=1.d0            !no discharge coeff.
      gcrit=discharge*gcrit 
      dgcdp=discharge*dgcdp                                                           
!                                                                       
!.....Use throat area ratio to modify GCRIT.                       
!                     
      gcrit=gcrit*choke_throat_ratio
      dgcdp=dgcdp*choke_throat_ratio
!
!.....Update gcrit and dgcdp (LSJ)
!
      gcrit0=gcrit
      dgcdp0=dgcdp         
!                                                                       
!.....Mass Flux from Momentum Solution.                            
!new
!from cell face mass flux
       mflux_l_avg=0.d0
       mflux_g_avg=0.d0
       DO i=1,num_throatface  !number of faces for each throat tt
          i1=n_face_throat(i) !face index for each throat tt
          mflux_l_avg=mflux_l_avg+mflux_l_nf(i1)   
          mflux_g_avg=mflux_g_avg+mflux_g_nf(i1)   
       ENDDO
       mflux_l_avg=mflux_l_avg/DBLE(num_throatface)
       mflux_g_avg=mflux_g_avg/DBLE(num_throatface)
       gmom=ABS(mflux_g_avg+mflux_l_avg) !alphagf*rhogf*vgf/sf + alphalf*rholf*vlf/sf
!                                                                       
!.....Choking Test
!               
      IF(gmom.gt.gcrit)then 
!                                                                       
!........choke=pressure boundary critical condition on/off(T/F)    
!                                                                       
         choke=.true.   
!                                                                       
         vcritt=gcrit/avrho 
         dvcdp=dgcdp/avrho 
!
         denomhf=1.d0/dvcdp 
         denomhf=denomhf+((avrho*ddx)/dt) 
         dvcdp=1.d0/denomhf 
!               
         vl_choke=signvc*vcritt 
         vg_choke=vl_choke
!                                                                       
!........Under-Relax Choked Velocity.                                 
!        Damping factor is larger for liquid velocity decrease so "water hammer" problems avoided.                           
!                                                                       
         vl_choke=relax_choke*vl_choke+(1.d0-relax_choke)*vl_choke_o
         vg_choke=relax_choke*vg_choke+(1.d0-relax_choke)*vg_choke_o
         
         WRITE(*,100) myrank,gmom,gcrit,gmom-gcrit,vl_choke,vl_choke_o
      ELSE 
!                                                                       
!........Unchoked, use normal momentum solution.       
!                                                                       
         choke=.false.
         
!         WRITE(*,101) myrank,gmom,gcrit,gmom-gcrit,vl_choke,vl_choke_o
      ENDIF 
      
      vl_choke_o=vl_choke
      vg_choke_o=vg_choke  
      
100   FORMAT('>>choked<< ',i5,100(e15.7,4x))  
!101   FORMAT('>>not choked<< ',i5,100(e15.7,4x))
      

!            
      END SUBROUTINE choke_model
