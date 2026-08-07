!
      SUBROUTINE heat_partition_fluid
!
!     This routine defines heat partitioning model and calculate heat flux & wall
!     temperature using Newton-Rhapson method in fluid cells.
!
      USE VOL_DATA        , ONLY: cell
      USE Wall_DATA       , ONLY: face
      USE Zzone           , ONLY: ncell_fluid
      USE Zcore           , ONLY: myrank      
      USE Zparam          , ONLY: nin_max,pi      
      USE Zbc_index       , ONLY: nbcon,icell_type,iface_wall,iface_wall1
      USE Znum_cell       , ONLY: i_neigh
      USE Zb_condition    , ONLY: twall
      USE Ztimecon        , ONLY: time
      USE Zconst1         , ONLY: iheatpart,iturb,vv_prob
      USE Zconst2         , ONLY: grav      
      USE Zcoord3         , ONLY: volp
      USE Zface           , ONLY: Zequation,SST,Kepsilon,Kepsilon_RNG,Kepsilon_real,Laminar, &
                                  qqcell,qecell,qcgcell,qclcell,q1cell,ndensitycell
      USE Zheat_partition , ONLY: q1,qe,qq,qcg,qcl,kfactor,dlo,   &
                                  utaul,hconvg,hconvl
      USE Ziat            , ONLY: dbubble_init,iat_nucl
      USE Zmodel          , ONLY: dtl        
      USE Zpress          , ONLY: p
      USE Zqvol           , ONLY: gamma_wall,qwall_solid,t_bulk,t_plus,t_plus_bulk,nsiteden_o,dry_weight
      USE Zturb           , ONLY: utau,yplus,walln
      USE Zvector         , ONLY: ug_o,ul_o,vl_o
      USE Znormal         , ONLY: sa_walll
!
      IMPLICIT NONE
!.....Local variables
      INTEGER i,j,j0,nfcondition
!
      INTEGER n_iter
      INTEGER, SAVE:: fr_lift,fr_nden,fr_depfreq,fr_depdia      
!      
      LOGICAL const_temperature,const_heat_flux,solid_fluid_interface,l_iter
      LOGICAL, SAVE::initial      
!
      REAL(8) st,hi_gas,hi_liq
      REAL(8) delt_sup,delt_wl,delt_wg,delh_gsl,delr_lg
      REAL(8) bub_dept_dia,bub_site_dens,d_bub_site_dens,bub_dept_freq
      REAL(8) frac_al,d_frac_al,frac_ag,d_frac_ag,bub_wait_time,bub_dept_dia2,bub_dept_dia2_pi
      REAL(8) d_qcl,d_qcg,d_qq,d_qe
      REAL(8) h_quench,h_qc,evap_rate,d_evap_rate,evap_rc     
      REAL(8) ftw,d_ftw,twall_new,twall_old,delt,error
      REAL(8) tsup_1,tsup_2,frac_ag_1,frac_ag_2,d_frac_ag_1,d_frac_ag_2
      REAL(8) bub_site_dens_1,bub_site_dens_2,d_bub_site_dens_1,d_bub_site_dens_2          
      REAL(8) pr,gamma_plus,beta_plus,gamma_plus_bulk,yplus_const                                        ! variable for RPI t_bulk model
      REAL(8) hfg_unal,grav_unal,pr_unal,d_depart_min,d_depart_max,hconvl_unal                           ! variables for Unal departure diameter model
      REAL(8) c_unal,a_unal,b_unal,ul_unal,psi_unal,bpsi_unal,tsub,gamma_unal,alphal_unal,condw,rhow,cpw ! variables for Unal departure diameter model
      REAL(8) uslide,reb,gs,ja,xkt1,aa,csl,rstar,ss,vol_bubble,l_slide,re_coeff,ra_coeff                 ! variables for Lift-off diameter model
      REAL(8) rho_plus,rho_func,cavity_r,ndensity_tmp,d_fritz                                            ! variables for Kocamusta's Site density model          
      REAL(8) nd_q, nd_freq                                                                              ! variables for Situ departure frequency model
!
      DATA d_depart_min,d_depart_max/1.D-5,1.4D-3/ 
      DATA initial /.true./     
!
!.....Define surface conditions
!
      solid_fluid_interface=.false.
      const_temperature=.false.
      const_heat_flux=.false.
!
      DO i=1,ncell_fluid
         j0=i_neigh(i)-1
         IF(icell_type(i).ne.1) cycle
!
         j=iface_wall(i)
         nfcondition=nbcon(j+j0)
         IF(nfcondition.eq.-2)THEN
            solid_fluid_interface=.true.
         ELSEIF(nfcondition.ge.-nin_max)THEN
            const_temperature=.true.
         ELSE
            const_heat_flux=.true.
         ENDIF
         IF(nfcondition.eq.-1)THEN
            face%wall_fluxl_diff(i)=0.0d0
            face%wall_fluxg_diff(i)=0.0d0
            face%wall_fluxd_diff(i)=0.0d0
            face%ddepartw(i)=0.0d0
            face%ratio_evap(i)=0.0d0
            GOTO 100
         ENDIF
!
         j=iface_wall1(i)
         nfcondition=nbcon(j+j0)
!
!........Set the sub-models in HPM
!      
         IF(initial)THEN
            fr_lift=mod(iheatpart,10000)/1000        ! Lift off model;                  0=off,  1=On
            fr_nden=mod(iheatpart,1000)/100          ! Nucleation site density model:   0=Cole, 1=Lemmert and Chwala, 2=ocamustafaogullary, 3=Hibiki, 4=Modified 
            fr_depfreq=mod(iheatpart,100)/10         ! Departure frequency model:       0=Cole, 1=Situ
            fr_depdia=mod(iheatpart,10)              ! Depart diameter model;           1: Cole and Rosenhow (1968) 2: Fritz (1935) 3: Tolubinsky model 4: Unal  5: Constant  6: Kocamusta
            nsiteden_o(i)=0.0d0
            IF(fr_nden.ne.3)initial=.FALSE.         
         ENDIF
!
         IF(initial.and.fr_nden.eq.3)THEN 
            IF(myrank.eq.0)THEN
               WRITE(*,*) 'Hibiki site density model cannot be used for Newton-Rhapson method of heat partitioning model.'
               WRITE(*,*) 'Site density will be calculated using Cole model.'
            ENDIF
            fr_nden=0
            initial=.FALSE.
         ENDIF
!
         IF(const_temperature)THEN
!
!...........Set constant wall temperature
!
            face%twall_partition(i)=twall(-nfcondition)
            twall_old=twall(-nfcondition)
!
         ELSE
!
!...........Initial guess of wall temperature
!
            twall_old=face%twall_partition(i)
            IF(time.le.1.0d-4)twall_old=cell%ts(i)+0.5d0
!
         ENDIF
!
!........Initialize papameters
!
         delh_gsl=cell%hgsat(i)-cell%hl(i)
         delr_lg=cell%rhol(i)-cell%rhog(i)
         dry_weight(i)=dmax1(0.0d0,dmin1(1.0d0,(cell%alphag(i)-0.90d0)/(0.99d0-0.90d0)))   
!         dry_weight(i)=dmax1(0.0d0,dmin1(1.0d0,(alphag_cm-cell(i)%alphag)/(alphag_cm-alphag_bc)))            
!
!........Friction velocity for T+ calculation
!
         IF(iturb.eq.Zequation)THEN
            utaul=utau(i)
         ELSEIF(iturb.eq.SST)THEN
            utaul=utau(i)
         ELSEIF(iturb.eq.Kepsilon)THEN
            utaul=utau(i)
         ELSEIF(iturb.eq.Kepsilon_RNG)THEN
            utaul=utau(i)
         ELSEIF(iturb.eq.Kepsilon_real)THEN
            utaul=utau(i)
         ELSEIF(iturb.eq.Laminar)THEN
            utaul=ul_o(i)
         ELSE
            WRITE(*,*)'HP_by:Turbulence_Model should be (-1,0,2,3,4)!'
         ENDIF
!         
!........RPI Boiling model to find T_bulk at Y+=250        
!
         pr=cell%eviscosl(i)*cell%cpl(i)/cell%condl(i)
         beta_plus=(3.85d0*pr**(1.d0/3.d0)-1.3d0)**2.d0+2.12d0*log(pr)
!         
!........Calculate T+ near wall cell        
!      
         IF(yplus(i).eq.0.0d0) yplus(i)=cell%rhol(i)*utaul*walln(i)/cell%lviscosl(i)
         gamma_plus=1.d-2*(pr*yplus(i))**4.d0/(1.d0+5.d0*yplus(i)*pr**3.d0)
         t_plus(i)=pr*yplus(i)*exp(-gamma_plus)+(2.12d0*log(yplus(i))+beta_plus)*    &
                   exp(-1.0d0/gamma_plus)         
!         
!........Calculate T+ at Y+=250
!            
         yplus_const=250.d0
         gamma_plus_bulk=1.d-2*(pr*yplus_const)**4.d0/(1.d0+5.d0*yplus_const*pr**3.d0)
         t_plus_bulk(i)=pr*yplus_const*exp(-gamma_plus_bulk)+(2.12d0*log(yplus_const)+beta_plus)* &
                        exp(-1.0d0/gamma_plus_bulk)
!
!........Start of Iteration
!
         l_iter=.true.
         n_iter=0
         DO WHILE(l_iter)
!
            delt_sup=DMAX1(0.0d0,twall_old-cell%ts(i))
            delt_wl=twall_old-cell%tl(i)
            delt_wg=twall_old-cell%tg(i)      
!         
!...........Estimate T_bulk at Y+=250 
!                    
            t_bulk(i)=twall_old-t_plus_bulk(i)/t_plus(i)*(twall_old-cell%tl(i))
            IF(ISNAN(t_bulk(i))) t_bulk(i)=cell%tl(i)       
!      
!...........Bubble departure diameter,1: Cole and Rosenhow (1968) 2: Fritz (1935) 3: Tolubinsky model 4: Unal  5: Constant  6: Kocamusta
!
            SELECTCASE(fr_depdia)
            CASE(1)           
               bub_dept_dia=max(1.0d-5,1.5d-4*sqrt(cell%sigma(i)/9.806d0/delr_lg) &
                           *(cell%rhol(i)*cell%cpl(i)*cell%ts(i)/cell%rhog(i)/delh_gsl)**(5.0d0/4.0d0))
            CASE(2) 
               bub_dept_dia=0.0208d0*38.0d0*sqrt(cell%sigma(i)/9.806d0/delr_lg)
            CASE(3)     
               bub_dept_dia=0.0006d0*exp((t_bulk(i)-cell%ts(i))/45.d0) 
               IF(bub_dept_dia.ge.1.4d-3)bub_dept_dia=1.4d-3
            CASE(4)
!              IF(iturb.eq.Zequation)THEN
!                 utaul=utau(i)
!              ELSEIF(iturb.eq.Kepsilon_l.or.iturb.eq.Kepsilon_g.or.iturb.eq.Kepsilon_lg)THEN
!                 utaul=utau(i)
!              ELSEIF(iturb.eq.Laminar)THEN
!                 utaul=dsqrt(dot_product(vl_o(i,:),vl_o(i,:)))
!              ELSE
!                 WRITE(*,*)'HP_by:Turbulence_Model should be (-1,0,2)!'
!              ENDIF
               hconvl_unal=max(10.0d0,st*cell%rhol(i)*cell%cpl(i)*utaul)     
               IF(vv_prob.eq.'suboex'.or.vv_prob.eq.'subo')THEN 
                  tsub=cell%ts(i)-t_bulk(i) 
                  condw=9.8585d0+0.01565d0*twall_old-5.76527d-7*twall_old**2.0d0+7.59296d-10*twall_old**3.0d0  !J/mK
                  rhow=8525.9888d0-0.38477d0*twall_old-6.19869d-6*twall_old**2.0d0                      !kg/m3 
                  cpw=418.71007d0-0.63632d0*twall_old+0.00453d0*twall_old**2.0d0-9.71703d-6*twall_old**3.0d0+   &   !J/kg-K            
                      8.94271d-9*twall_old**4.0d0-2.8983d-12*twall_old**5.0d0
               ELSE
                  tsub=cell%ts(i)-t_bulk(i) 
                  condw=14.8d0    !J/mK
                  rhow=7850.0d0 !kg/m3 
                  cpw=486.0d0   !J/kg-K
               ENDIF    
               alphal_unal=cell%condl(i)/(cell%rhol(i)*cell%cpl(i))
               gamma_unal=condw*rhow*cpw/(cell%condl(i)*cell%rhol(i)*cell%cpl(i))
               hfg_unal=(cell%hgsat(i)-cell%hlsat(i))
               grav_unal=sqrt(dot_product(grav(:),grav(:)))
               pr_unal=cell%eviscosl(i)*cell%cpl(i)/cell%condl(i) !momentum diffusivity/heat diffusivity,nu/alpha,nu/(k/(rho*cp)
               c_unal=hfg_unal*cell%eviscosl(i)*(cell%cpl(i)/(0.013*hfg_unal*Pr_unal**1.7d0))**3.0d0 &
                     /(cell%sigma(i)/(cell%rhol(i)-cell%rhog(i))*grav_unal)**0.5d0
               a_unal=(qwall_solid(-nfcondition)-hconvl_unal*tsub)**(1.0d0/3.0d0)*cell%condl(i) &
                     /(2.0d0*c_unal**(1.d0/3.d0)*hfg_unal*sqrt(pi*alphal_unal)*cell%rhog(i)) &
                     *sqrt(gamma_unal)
               IF(a_unal.lt.0.d0)a_unal=0.0d0
               b_unal=tsub/(2.0d0*(1.0d0-cell%rhog(i)/cell%rhol(i)))
               ul_unal=sqrt(dot_product(vl_o(i,:),vl_o(i,:)))
               psi_unal=max1((ul_unal/0.61d0)**0.47d0, 1.0d0)   
               bpsi_unal=b_unal*psi_unal
               IF(bpsi_unal.gt.0.0d0)bub_dept_dia=2.42D-5*p(i)**0.709d0*a_unal/dsqrt(bpsi_unal)
            CASE(5) !Kocamustafaogullari(1983) core_catcher
               bub_dept_dia=0.0000264d0*38.0d0*dsqrt(cell%sigma(i)/9.806d0/delr_lg)*(delr_lg/cell%rhog(i))**0.9d0 
            CASE(6)
               bub_dept_dia=dbubble_init   
            ENDSELECT
            bub_dept_dia=dmin1(dmax1(d_depart_min,bub_dept_dia),d_depart_max)      
! 
!...........Bubble departure frequency: 0-Cole, 1-Situ(2007)         
!
            SELECTCASE(fr_depfreq)
            CASE(0)   
               bub_dept_freq=DSQRT(13.0747d0*delr_lg/(bub_dept_dia*cell%rhol(i)))    !4.0*9.806/3.0=13.0747
            CASE(1)   
               !non-dimensional nucleate boiling heat flux
               nd_q=(qqcell(i)+qecell(i))*bub_dept_dia/delh_gsl/cell%alphal(i)/cell%rhog(i)  
               !non-dimensional departure frequency
               nd_freq=4.06d0*nd_q**0.803d0                                              
               bub_dept_freq=nd_freq*cell%alphal(i)/bub_dept_dia**2.0d0       
            ENDSELECT
            IF(bub_dept_dia.le.0.0d0)bub_dept_freq=0.0d0
!
!...........Bubble wait time
!

            bub_wait_time=0.8d0/bub_dept_freq
            IF(bub_dept_freq.le.0.0d0)bub_wait_time=0.0d0           
!   
!...........Nucleation Site Density 
!     
            IF(delt_sup.le.1.d-8)THEN
               bub_site_dens=0.0d0
            ELSE
               SELECTCASE(fr_nden)
               CASE(0) 
                  ! Cole
                  bub_site_dens=12366.45d0*delt_sup**1.805d0       ! 185**1.805=12366.45d0
                  d_bub_site_dens=22321.44d0*delt_sup**0.805d0     ! (185**1.805)*1.805=22321.44d0         
               CASE(1) 
                  ! Lemmert and Chwala
                  bub_site_dens=15545.5405d0*delt_sup**1.805d0     ! 210**1.805=15545.5405d0
                  d_bub_site_dens=28059.7006d0*delt_sup**0.805d0   ! (210**1.805)*1.805=28059.7006d0     
               CASE(2) 
                  ! Kocamustafaogullary
                  rho_plus=delr_lg/cell%rhog(i)            
                  d_fritz=0.0208d0*38.0d0*dsqrt(cell%sigma(i)/9.806d0/delr_lg)
                  d_fritz=d_fritz*0.0012d0*rho_plus**0.9d0
                  cavity_r=2.0d0*cell%sigma(i)*cell%ts(i)/(cell%rhog(i)*delh_gsl)
                  cavity_r=cavity_r/(d_fritz/2)
                  rho_func=2.157d-7*rho_plus**(-3.2d0)*(1.0d0+0.0049d0*rho_plus)**4.13d0 
                  ndensity_tmp=rho_func*cavity_r**(-4.4d0)/d_fritz**2.0d0
                  bub_site_dens=ndensity_tmp*delt_sup**4.4d0 
                  d_bub_site_dens=4.4d0*ndensity_tmp*delt_sup**3.4d0           
               CASE(3) 
                  ! Hibiki model cannot be used in Newton-Rhapson Mothod
               CASE(4) 
                  ! Modified Cole        
                  bub_site_dens=12436.24d0*delt_sup**1.805d0        ! 12436.24=0.8d0*9.922D+5/10**1.805
                  d_bub_site_dens=22447.42d0*delt_sup**0.805d0      ! 22447.42=12436.24*1.805    
               ENDSELECT
            ENDIF        
!
!...........Bubble lift-off HKCHO-beg 
!
            IF(fr_lift.eq.1)THEN
               uslide=1.0d0*dsqrt(dot_product(vl_o(i,:),vl_o(i,:)))
               reb=uslide*cell%D1(i)*cell%rhol(i)/cell%lviscosl(i)
               gs=dabs(dsqrt(dot_product(vl_o(i,:),vl_o(i,:)))/walln(i))*(cell%D1(i)/2.)/uslide
               csl=3.877d0*dsqrt(gs)*(1./REB**2 + 0.014*gs**2)**0.25
               ja=cell%rhol(i)*cell%cpl(i)*(dmax1(delt_sup,0.0d0))/cell%rhog(i)/delh_gsl
               xkt1=cell%lcondl(i)/cell%rhol(i)/cell%cpl(i)
               aa=3.46d0*ja*dsqrt(xkt1/3.141592d0)
               IF(aa.gt.0.0d0)THEN
                  rstar=8.34d0*(csl*(bub_dept_dia*uslide/aa**2)**2)**(-0.7)
               ELSE
                  rstar=0.0d0
               ENDIF
               dlo=bub_dept_dia*(1.d0+rstar)
               dlo=dmin1(d_depart_max,25.0d0*bub_dept_dia,dlo)
               IF(bub_site_dens.ne.0.0) THEN
                  ss=1./dsqrt(bub_site_dens)
               ELSE
                  ss=-1.d0
               ENDIF       
               vol_bubble=pi*dlo**3/6.d0
               IF(aa.gt.0.0d0)THEN
                  l_slide=(1./32.)*(cell%rhol(i)/cell%rhog(i)-1.d0)*9.81d0*((dlo-bub_dept_dia)/aa)**4
               ELSE
                  l_slide=0.0d0        
               ENDIF         
               !re_coeff
               IF(ss.ge.0.0 .and. l_slide.gt.bub_dept_dia) THEN
                  IF(ss.lt.bub_dept_dia) THEN
                     re_coeff=(ss/dlo)**3
                  ELSEIF (ss.le.l_slide) THEN
                     re_coeff=(1.d0-(bub_dept_dia/dlo)**3)*(ss-bub_dept_dia)/(l_slide-bub_dept_dia)+(bub_dept_dia/dlo)**3
                  ELSE
                     re_coeff=1.d0
                  ENDIF
               ELSE
                  re_coeff=1.d0
               ENDIF
               !ra_coeff	        
               IF(ss.ge.0.0 .and. l_slide.gt.bub_dept_dia) THEN
                  IF(ss.lt.bub_dept_dia) THEN
                     ra_coeff=(ss/dlo)**2
                  ELSEIF (ss.le.l_slide) THEN
                     ra_coeff=(1.d0-(bub_dept_dia/dlo)**2)*(ss-bub_dept_dia)/(l_slide-bub_dept_dia)+(bub_dept_dia/dlo)**2
                  ELSE
                     ra_coeff=1.d0
                  ENDIF
               ELSE
                  ra_coeff=1.d0
               ENDIF
            ENDIF         
!
!...........Fraction of influential vapour area
!         
            bub_dept_dia2=bub_dept_dia*bub_dept_dia/4.0d0
            bub_dept_dia2_pi=pi*bub_dept_dia2
!
            frac_ag_1=0.8d0
            frac_ag_2=1.0d0
!
            SELECTCASE(fr_nden)
            CASE(0) 
               ! Cole
               tsup_1=frac_ag_1/(bub_dept_dia2_pi*kfactor)/12366.45d0   
               tsup_1=tsup_1**0.5540167d0                                ! 1.0/1.805=0.5540167d0
               d_frac_ag_1=bub_dept_dia2_pi*kfactor*22321.44d0*tsup_1**0.805d0       
               tsup_2=frac_ag_2/(bub_dept_dia2_pi*kfactor)/12366.45d0    
               tsup_2=tsup_2**0.5540167d0+dtl                            ! Add 3.0 to extend interpolation region
               d_frac_ag_2=0.0d0        
            CASE(1) 
               ! Lemmert and Chwala
               tsup_1=frac_ag_1/(bub_dept_dia2_pi*kfactor)/15545.5405d0         
               tsup_1=tsup_1**0.5540167d0                                ! 1.0/1.805=0.5540167d0
               d_frac_ag_1=bub_dept_dia2_pi*kfactor*28059.7006d0*tsup_1**0.805d0         
               tsup_2=frac_ag_2/(bub_dept_dia2_pi*kfactor)/15545.5405d0         
               tsup_2=tsup_2**0.5540167d0+dtl                            ! Add 3.0 to extend interpolation region
               d_frac_ag_2=0.0d0                
            CASE(2) 
               ! Kocamustafaogullary
               tsup_1=frac_ag_1/(bub_dept_dia2_pi*kfactor)/ndensity_tmp
               tsup_1=tsup_1**(1.0d0/4.4d0)                              
               d_frac_ag_1=bub_dept_dia2_pi*kfactor*4.4d0*ndensity_tmp*tsup_1**3.4d0       
               tsup_2=frac_ag_2/(bub_dept_dia2_pi*kfactor)/ndensity_tmp    
               tsup_2=tsup_2**(1.0d0/4.4d0)+dtl                          ! Add 3.0 to extend interpolation region
               d_frac_ag_2=0.0d0                    
            CASE(3) 
               ! Hibiki model cannot be used in Newton-Rhapson Mothod
            CASE(4) 
               ! Modified Cole  
               tsup_1=frac_ag_1/(bub_dept_dia2_pi*kfactor)/12436.24d0   
               tsup_1=tsup_1**0.5540167d0                                ! 1.0/1.805=0.5540167d0
               d_frac_ag_1=bub_dept_dia2_pi*kfactor*22447.42d0*tsup_1**0.805d0       
               tsup_2=frac_ag_2/(bub_dept_dia2_pi*kfactor)/12436.24d0    
               tsup_2=tsup_2**0.5540167d0+dtl                            ! Add 3.0 to extend interpolation region
               d_frac_ag_2=0.0d0   
            ENDSELECT
!
            IF(delt_sup.le.tsup_1)THEN
               frac_ag=bub_dept_dia2_pi*kfactor*bub_site_dens
               d_frac_ag=bub_dept_dia2_pi*kfactor*d_bub_site_dens
            ELSEIF(delt_sup.le.tsup_2)THEN
!
!..............Cubic interpolation between tsup_1 and tsup_2
!
               CALL cub_interp(delt_sup,tsup_1,tsup_2,frac_ag_1,frac_ag_2,d_frac_ag_1,d_frac_ag_2,frac_ag,d_frac_ag)
!
            ELSE
               frac_ag=frac_ag_2
               d_frac_ag=d_frac_ag_2
            ENDIF
!
!...........Fraction of influential liquid area
!
            frac_al=1.0d0-frac_ag
            d_frac_al=-d_frac_ag
            IF(frac_al.le.1.0d-4)THEN
               frac_al=1.0d-4
               d_frac_al=0.0d0
            ENDIF         
!
!...........Cubic interpolation of Bubble site density between frac_ag_1 and frac_ag_2
!
            frac_ag_1=0.8d0
            frac_ag_2=1.05d0   ! Add 0.05 to extend interpolation region
!         
            bub_site_dens_1=frac_ag_1/bub_dept_dia2_pi/kfactor 
            d_bub_site_dens_1=1.0d0/bub_dept_dia2_pi/kfactor          
            bub_site_dens_2=frac_ag_2/bub_dept_dia2_pi/kfactor
            d_bub_site_dens_2=0.0d0
!         
            IF(frac_ag.le.frac_ag_1)THEN
               bub_site_dens=frac_ag/bub_dept_dia2_pi/kfactor   
               d_bub_site_dens=1.0d0/bub_dept_dia2_pi/kfactor
            ELSEIF(frac_ag.le.frac_ag_2)THEN
!
!..............Cubic interpolation between frac_ag_1 andfrac_ag_2
!
               CALL cub_interp(frac_ag,frac_ag_1,frac_ag_2,bub_site_dens_1,bub_site_dens_2,d_bub_site_dens_1,d_bub_site_dens_2,bub_site_dens,d_bub_site_dens)
!
            ELSE
               bub_site_dens=bub_site_dens_2
               d_bub_site_dens=d_bub_site_dens_2
            ENDIF   
            IF(bub_dept_dia2_pi.le.0.0d0)THEN
               bub_site_dens=0.0d0 
               d_bub_site_dens=0.0d0
            ENDIF
!            bub_site_dens=(1.0d0-relax_na)*bub_site_dens+relax_na*nsiteden_o(i)         
!            nsiteden_o(i)=bub_site_dens            
!
!...........Liquid convective heat flux
!         
            st=0.0045d0
            IF(t_plus(i).eq.0.or.ISNAN(t_plus(i)))THEN
               hconvl=cell%rhol(i)*cell%cpl(i)*ul_o(i) 
            ELSE
               hconvl=cell%rhol(i)*cell%cpl(i)*utaul/t_plus(i)       
            ENDIF         
            qcl=frac_al*hconvl*delt_wl
            d_qcl=d_frac_al*hconvl*delt_wl+frac_al*hconvl
!
!...........Gas convective heat flux
!       
            hconvg=st*cell%rhog(i)*cell%cpg(i)*ug_o(i) 
            qcg=hconvg*delt_wg
            d_qcg=hconvg         
!
!...........Quenching heat flux
!
            h_qc=DSQRT(bub_wait_time*cell%lcondl(i)*cell%rhol(i)*cell%cpl(i)/pi)
            h_quench=2.0d0*bub_dept_freq*h_qc
            qq=frac_ag*h_quench*(twall_old-t_bulk(i))
            d_qq=d_frac_ag*h_quench*(twall_old-t_bulk(i))+frac_ag*h_quench
!
!...........Evaporation heat flux
!
            evap_rc=pi*bub_dept_dia**3.0d0*cell%rhog(i)/6.0d0
            evap_rate=evap_rc*bub_dept_freq*bub_site_dens
            d_evap_rate=evap_rc*bub_dept_freq*d_bub_site_dens
            qe=evap_rate*delh_gsl
            d_qe=d_evap_rate*delh_gsl
!
!...........Terminate the iteration for constant wall temperature
!
            IF(const_temperature)l_iter=.false.
!
!...........Find wall temperature when wall heat flux is given
!
            IF(solid_fluid_interface.or.const_heat_flux)THEN
!
!..............Set wall heat flux boundary conditions
!
               IF(solid_fluid_interface)THEN
!
!.................Solid-fluid interface
!.................Set wall boundary flux condition here
!
!.....AT fluid-solid interface (nbcon=-2)
!
                  IF(myrank.eq.0)THEN
                     PRINT*,'Heat partition using Newton-Raphson is not available in fluid-solid interface!'
                     STOP
                  ENDIF
!
               ELSE
!
!.................Constant wall heat flux
!
                  q1=qwall_solid(-nfcondition)
!
               ENDIF
!
               ftw=(1.0d0-dry_weight(i))*(qcl+qq+qe)+dry_weight(i)*qcg-q1
               d_ftw=(1.0d0-dry_weight(i))*(d_qcl+d_qq+d_qe)+dry_weight(i)*d_qcg
!
!..............Allow temperature change within 3K
!
               delt=ftw/d_ftw
               IF(delt.lt.-3.0d0)delt=-3.0d0
               IF(delt.gt.3.0d0)delt=3.0d0
!
               twall_new=twall_old-delt
!
!..............Turn off the heat partition model:
!................. if wall temperature falls below saturation temperature
!................ if liquid enthalpy is greater than saturated gas enthalpy
!
               IF((twall_new.lt.cell%ts(i)+1.0d-1.and.n_iter.gt.100).or.delh_gsl.le.0.0d0)THEN            
                  qq=0.0d0
                  qe=0.0d0
                  qcg=0.0d0               
                  qcl=q1/(1.0d0-dry_weight(i))
                  l_iter=.false.
                  IF(dry_weight(i).eq.1)THEN
                     qcl=0.0d0
                     qcg=q1
                  ENDIF
                  EXIT
               ENDIF
!
!..............Interation converges if error is less than 1.e-8
!
               error=DABS((twall_new-twall_old)/twall_old)
               IF(error.lt.1.0d-6)l_iter=.false.
!
!..............Physical properties are out of range
!
               IF(n_iter.ge.1000)THEN
                  twall_new=face%twall_partition(i)
                  IF(dry_weight(i).eq.1.0d0)THEN
                     qq=0.0d0
                     qe=0.0d0
                     qcl=0.0d0
                     qcg=qcgcell(i)
                     IF(time.lt.1.0d-4)qcg=q1    
                  ELSEIF(dry_weight(i).eq.0.0d0)THEN
                     qq=qqcell(i)
                     qe=qecell(i)
                     qcl=qclcell(i)
                     qcg=0.0d0
                     IF(time.lt.1.0d-4)qcl=q1             
                  ELSE
                     qq=qqcell(i)/(1.0d0-dry_weight(i))
                     qe=qecell(i)/(1.0d0-dry_weight(i))
                     qcl=qclcell(i)/(1.0d0-dry_weight(i))
                     qcg=qcgcell(i)/dry_weight(i)
                     IF(time.lt.1.0d-4)THEN
                        qcl=q1/(1.0d0-dry_weight(i))
                        qcg=q1/dry_weight(i)               
                     ENDIF             
                  ENDIF       
                  l_iter=.false.
                  EXIT
               ENDIF
!
!..............Update solution
!
               twall_old=twall_new
               n_iter=n_iter+1
!
            ENDIF
!
         ENDDO
!    
!      
!........Save parameters (Cell) to check the adequacy of heat partitioning model calculation
!
         qqcell(i)=qq*(1.0d0-dry_weight(i))
         qecell(i)=qe*(1.0d0-dry_weight(i))
         qclcell(i)=qcl*(1.0d0-dry_weight(i))
         qcgcell(i)=qcg*dry_weight(i)   
         q1cell(i)=q1  
         ndensitycell(i)=bub_site_dens
         IF(bub_dept_dia.ge.0)cell%Ddepart(i)=bub_dept_dia     
         cell%twall(i)=twall_new      
!      
!........Save parameters (Face) to check the adequacy of heat partitioning model calculation
!  
         face%twall_partition(i)=twall_new
         face%wall_fluxl_diff(i)=(qq+qcl)*(1.0d0-dry_weight(i))*sa_walll(i)
         face%wall_fluxg_diff(i)=qcg*dry_weight(i)*sa_walll(i)
         face%wall_fluxd_diff(i)=0.0d0    
         face%ddepartw(i)=bub_dept_dia
         face%bfreq(i)=bub_dept_freq      
         IF(q1.gt.0.0d0)THEN
            face%ratio_evap(i)=qe*(1.0d0-dry_weight(i))/q1
         ELSE
            face%ratio_evap(i)=0.0d0
         ENDIF
!
!........Calculate vapor generation rate(gamma_wall)
!
         IF((1.0d0-dry_weight(i))*qe.gt.0.0d0)THEN
            hi_gas=cell%hgsat(i)
            hi_liq=cell%hl(i)
            gamma_wall(i)=gamma_wall(i)+qe*(1.0d0-dry_weight(i))/(hi_gas-hi_liq)*sa_walll(i)/volp(i)
!
!...........IAT source by nuculate boiling
!
            IAT_nucl(i)=4.0d0*bub_dept_dia2_pi*bub_site_dens*bub_dept_freq*frac_ag*sa_walll(i)/volp(i)
            !IF(bub_dept_dia2_pi*bub_site_dens.ge.0.99d0)IAT_nucl(i)=IAT_nucl(i)*0.481d0
            IF(fr_lift.eq.1)IAT_nucl(i)=pi*dlo**2.d0*bub_site_dens*bub_dept_freq*frac_ag*sa_walll(i)/Volp(i)*ra_coeff            
!
         ENDIF
100      CONTINUE
      ENDDO
!
      RETURN 
      END SUBROUTINE heat_partition_fluid
