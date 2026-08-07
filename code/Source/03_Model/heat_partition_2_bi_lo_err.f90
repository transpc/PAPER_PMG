!
      Function heat_partition_2_bi_lo_err(deltaTs,deltaTl,deltaTg,deltarho,deltahlg &
                                    ,d_depart_init,fr_lift,fr_nden,fr_depfreq,fr_depdia,i,nfcondition,ra_coeff,tw)
!
!     This routine calculates partitioned heat with given wall temperature and 
!     return the heat error, err=q1-q2.
!     q1=assinged heat, q2=partitioned heat, twall=wall temperature.
!                                    
      USE VOL_DATA                    
      USE Zparam          , ONLY: pi,sqrt_pi
      USE Zconst1         , ONLY: iturb,vv_prob
      USE Zconst2         , ONLY: grav
      USE Zface
      USE Zheat_partition , ONLY: q1,qe,qq,qcg,qcl,bfreq,twait,kfactor,a_two,dlo,   &
                                   ndensity,d_depart,utaul,hconvg,hconvl
      USE Zpress          , ONLY: p
      USE Zqvol           , ONLY: qwall_solid,t_bulk,t_plus,t_plus_bulk,dry_weight
      USE Zturb           , ONLY: walln,utau,yplus
      USE Zvector         , ONLY: ug_o,ul_o,vl_o
      USE Zuserdefined    , ONLY: udfl_psbt_cfx_model
!
      IMPLICIT NONE
!      
      INTEGER i,nfcondition,fr_lift,fr_nden,fr_depfreq,fr_depdia
      REAL(8) deltarho,deltahlg
      REAL(8) deltaTs,deltaTl,deltaTg
      REAL(8) a_single, bub_area
      REAL(8) q2, st
      REAL(8) c_unal,a_unal,b_unal,ul_unal,psi_unal,bpsi_unal,tsub,gamma_unal,alphal_unal,condw,rhow,cpw
      REAL(8) hfg_unal,grav_unal,pr_unal,d_depart_min,d_depart_max,hconvl_unal
      REAL(8) heat_partition_2_bi_lo_err,d_depart_init
      REAL(8) uslide,reb,gs,ja,xkt1,aa,csl,rstar,ss,vol_bubble,l_slide,re_coeff,ra_coeff,tw 
      REAL(8) gamma_plus,beta_plus,gamma_plus_bulk,yplus_const    
      REAL(8) gas_const,angle_cont,mu,rho_plus,rho_func,cavity_r,d_fritz
      REAL(8) nd_q,nd_freq,unal_tmp      
!
      DATA d_depart_min,d_depart_max/1.D-5,1.4D-3/
!
!.....Friction velocity for T+ calculation
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
!.....RPI Boiling model to find T_bulk at Y+=250        
!
      pr_unal=cell%eviscosl(i)*cell%cpl(i)/cell%condl(i)
      beta_plus=(3.85d0*pr_unal**(1.d0/3.d0)-1.3d0)**2.d0+2.12d0*log(pr_unal)
      !Calculate T+ near wall cell        
      IF(yplus(i).eq.0.0d0) yplus(i)=cell%rhol(i)*utaul*walln(i)/cell%lviscosl(i)
      gamma_plus=1.d-2*(pr_unal*yplus(i))**4.d0/(1.d0+5.d0*yplus(i)*pr_unal**3.d0)
      t_plus(i)=pr_unal*yplus(i)*exp(-gamma_plus)+(2.12d0*log(yplus(i))+beta_plus)*    &
               exp(-1.0d0/gamma_plus)
      !Calculate T+ at Y+=250
      yplus_const=250.d0
      gamma_plus_bulk=1.d-2*(pr_unal*yplus_const)**4.d0/(1.d0+5.d0*yplus_const*pr_unal**3.d0)
      t_plus_bulk(i)=pr_unal*yplus_const*exp(-gamma_plus_bulk)+(2.12d0*log(yplus_const)+beta_plus)* &
                     exp(-1.0d0/gamma_plus_bulk)
      !Estimate T_bulk at Y+=250 
      t_bulk(i)=tw-t_plus_bulk(i)/t_plus(i)*(tw-cell%tl(i))   
      IF(ISNAN(t_bulk(i))) t_bulk(i)=cell%tl(i)                   
!      
!.....Convective heat transfer coefficient
!
       st=0.0045d0
       IF(t_plus(i).eq.0.or.ISNAN(t_plus(i)))THEN
          hconvl=cell%rhol(i)*cell%cpl(i)*ul_o(i) 
       ELSE
          hconvl=cell%rhol(i)*cell%cpl(i)*utaul/t_plus(i)       
       ENDIF
       hconvg=st*cell%rhog(i)*cell%cpg(i)*ug_o(i) 
!   
!.....Nucleation Site Density 
!     
      IF(deltaTs.le.1.d-8)THEN
         ndensity=0.0d0
      ELSE
         SELECTCASE(fr_nden)
         CASE(0) 
            ! Cole
            ndensity=(185.0d0*deltaTs)**1.805d0         
         CASE(1) 
            ! Lemmert and Chwala
            ndensity=(210.0d0*deltaTs)**1.805d0         
         CASE(2) 
            ! Kocamustafaogullary
            rho_plus=deltarho/cell%rhog(i)            
            d_fritz=0.0208d0*38.0d0*dsqrt(cell%sigma(i)/9.806d0/deltarho)
            d_fritz=d_fritz*0.0012d0*rho_plus**0.9d0
            cavity_r=2.0d0*cell%sigma(i)*cell%ts(i)/(cell%rhog(i)*deltahlg*(tw-cell%ts(i)))
            cavity_r=cavity_r/(d_fritz/2)
            rho_func=2.157d-7*rho_plus**(-3.2d0)*(1.0d0+0.0049d0*rho_plus)**4.13d0 
            ndensity=rho_func*cavity_r**(-4.4d0)/d_fritz**2.0d0           
         CASE(3) 
            ! Hibiki
            gas_const=462.0d0       !J/kgK for steam
            angle_cont=38.0d0/180.0d0*pi
            mu=0.722d0    !rad                        
            rho_plus=log10(deltarho/cell%rhog(i))
            rho_func=-0.01064d0+0.48246d0*rho_plus-0.22715d0*rho_plus**2.0d0+0.05468d0*rho_plus**3.0d0
            cavity_r=2.0d0*cell%sigma(i)*(1.0d0+cell%rhog(i)/cell%rhol(i))/cell%p(i)
            cavity_r=cavity_r/(exp(deltahlg*(tw-cell%ts(i))/(gas_const*cell%tg(i)*cell%ts(i)))-1.0d0)   !tg is replaced by tw
            cavity_r=DMAX1(0.0d0,cavity_r)
            ndensity=4.72d5*(1.0d0-exp(-angle_cont**2.0d0/8.0d0/mu**2.0d0))
            ndensity=ndensity*(exp(rho_func*2.5d-6/cavity_r)-1)      
            ndensity=DMAX1(0.0d0,ndensity)       
         CASE(4) 
         ! Modified Cole
            ndensity=(0.8d0*9.922D+5)*(deltaTs/10.0d0)**1.805d0          
         ENDSELECT
      ENDIF
!      
!.....Bubble departure diameter,1: Cole and Rosenhow (1968) 2: Fritz (1935) 3: Tolubinsky model 4: Unal  5: Constant  6: Kocamusta
!
      SELECTCASE(fr_depdia)
      CASE(1)           
         d_depart=dmax1(1.0d-5,1.5d-4*dsqrt(cell%sigma(i)/9.806d0/deltarho) &
                  *(cell%rhol(i)*cell%cpl(i)*cell%ts(i)/cell%rhog(i)/deltahlg)**(5.0d0/4.0d0))
      CASE(2) 
         d_depart=0.0208d0*38.0d0*dsqrt(cell%sigma(i)/9.806d0/deltarho)
      CASE(3)     
         d_depart=0.0006d0*exp((t_bulk(i)-cell%ts(i))/45.d0)                
      CASE(4)
         IF(iturb.eq.Laminar)  utaul=dsqrt(dot_product(vl_o(i,:),vl_o(i,:)))
         hconvl_unal=dmax1(10.0d0,st*cell%rhol(i)*cell%cpl(i)*utaul)     
         IF(udfl_psbt_cfx_model)THEN 
            tsub=cell%ts(i)-t_bulk(i) 
            condw=9.8585d0+0.01565d0*tw-5.76527d-7*tw**2.0d0+7.59296d-10*tw**3.0d0  !J/mK
            rhow=8525.9888d0-0.38477d0*tw-6.19869d-6*tw**2.0d0                      !kg/m3 
            cpw=418.71007d0-0.63632d0*tw+0.00453d0*tw**2.0d0-9.71703d-6*tw**3.0d0+   &   !J/kg-K            
                8.94271d-9*tw**4.0d0-2.8983d-12*tw**5.0d0
         ELSE
            tsub=cell%ts(i)-cell%tl(i)
            condw=14.8d0    !J/mK
            rhow=7850.0d0 !kg/m3 
            cpw=486.0d0   !J/kg-K
         ENDIF    
         alphal_unal=cell%condl(i)/(cell%rhol(i)*cell%cpl(i))
         gamma_unal=condw*rhow*cpw/(cell%condl(i)*cell%rhol(i)*cell%cpl(i))
         hfg_unal=(cell%hgsat(i)-cell%hlsat(i))
         grav_unal=dsqrt(dot_product(grav(:),grav(:)))
         pr_unal=cell%eviscosl(i)*cell%cpl(i)/cell%condl(i) !momentum diffusivity/heat diffusivity,nu/alpha,nu/(k/(rho*cp)
         c_unal=hfg_unal*cell%eviscosl(i)*(cell%cpl(i)/(0.013*hfg_unal*Pr_unal**1.7d0))**3.0d0 &
               /(cell%sigma(i)/(cell%rhol(i)-cell%rhog(i))*grav_unal)**0.5d0
         unal_tmp=qwall_solid(-nfcondition)-hconvl_unal*tsub
         IF(unal_tmp.le.0.d0)THEN
            a_unal=0.0d0
         ELSE
            a_unal=unal_tmp**(1.0d0/3.0d0)*cell%condl(i) &
               /(2.0d0*c_unal**(1.d0/3.d0)*hfg_unal*dsqrt(pi*alphal_unal)*cell%rhog(i)) &
               *dsqrt(gamma_unal)
         ENDIF   
         IF(a_unal.lt.0.d0)a_unal=0.0d0
         b_unal=tsub/(2.0d0*(1.0d0-cell%rhog(i)/cell%rhol(i)))
         ul_unal=dsqrt(dot_product(vl_o(i,:),vl_o(i,:)))
         psi_unal=dmax1((ul_unal/0.61d0)**0.47d0, 1.0d0)   
         bpsi_unal=b_unal*psi_unal
         IF(bpsi_unal.gt.0.0d0)d_depart=2.42D-5*p(i)**0.709d0*a_unal/dsqrt(bpsi_unal)   
      CASE(5) !Kocamustafaogullari(1983) core_catcher
         d_depart=0.0000264d0*38.0d0*dsqrt(cell%sigma(i)/9.806d0/deltarho)*(deltarho/cell%rhog(i))**0.9d0 
      CASE(6)
         d_depart=d_depart_init   
      ENDSELECT
!      
!.....Set minimum & maximum departure diameter      
!  
      IF(vv_prob.eq.'suboex')d_depart_max=5.0d-3
      d_depart=dmin1(dmax1(d_depart_min,d_depart),d_depart_max)
!
!.....Bubble lift-off HKCHO-beg 
!
      IF(fr_lift.eq.1)THEN
         !dlo
         uslide=1.0d0*dsqrt(dot_product(vl_o(i,:),vl_o(i,:)))
         reb=uslide*cell%D1(i)*cell%rhol(i)/cell%lviscosl(i)
         gs=dabs(dsqrt(dot_product(vl_o(i,:),vl_o(i,:)))/walln(i))*(cell%D1(i)/2.)/uslide
         csl=3.877d0*dsqrt(gs)*(1./REB**2 + 0.014*gs**2)**0.25
         ja=cell%rhol(i)*cell%cpl(i)*(dmax1(deltaTs,0.0d0))/cell%rhog(i)/deltahlg
         xkt1=cell%lcondl(i)/cell%rhol(i)/cell%cpl(i)
         aa=3.46d0*ja*dsqrt(xkt1/3.141592d0)
         IF(aa.gt.0.0d0.and.d_depart.ne.0.0d0)THEN
            rstar=8.34d0*(csl*(d_depart*uslide/aa**2)**2)**(-0.7)
         ELSE
            rstar=0.0d0
         ENDIF
         dlo=d_depart*(1.d0+rstar)
         dlo=dmin1(d_depart_max,25.0d0*d_depart,dlo)
         !d_depart=dlo
         IF(ndensity.ne.0.0) THEN
           ss=1./dsqrt(ndensity)
         ELSE
           ss=-1.d0
         ENDIF       
         vol_bubble=pi*dlo**3/6.d0
         IF(aa.gt.0.0d0)THEN
            l_slide=(1./32.)*(cell%rhol(i)/cell%rhog(i)-1.d0)*9.81d0*((dlo-d_depart)/aa)**4
         ELSE
            l_slide=0.0d0        
         ENDIF
         !re_coeff
         IF(ss.ge.0.0 .and. l_slide.gt.d_depart) THEN
            IF(ss.lt.d_depart) THEN
               re_coeff=(ss/dlo)**3
            ELSEIF (ss.le.l_slide) THEN
               re_coeff=(1.d0-(d_depart/dlo)**3)*(ss-d_depart)/(l_slide-d_depart)+(d_depart/dlo)**3
            ELSE
               re_coeff=1.d0
           ENDIF
         ELSE
            re_coeff=1.d0
         ENDIF
         !ra_coeff	        
         IF(ss.ge.0.0 .and. l_slide.gt.d_depart) THEN
            IF(ss.lt.d_depart) THEN
               ra_coeff=(ss/dlo)**2
            ELSEIF (ss.le.l_slide) THEN
               ra_coeff=(1.d0-(d_depart/dlo)**2)*(ss-d_depart)/(l_slide-d_depart)+(d_depart/dlo)**2
            ELSE
               ra_coeff=1.d0
            ENDIF
         ELSE
            ra_coeff=1.d0
         ENDIF
      ENDIF        
! 
!.....Bubble departure frequency: 0-Cole, 1-Situ(2007)         
!
      SELECTCASE(fr_depfreq)
      CASE(0)   
         !4.0*9.806/3.0=13.0747
         bfreq=DSQRT(13.0747d0*deltarho/(d_depart*cell%rhol(i)))  
      CASE(1)   
         !non-dimensional nucleate boiling heat flux
         nd_q=(qqcell(i)+qecell(i))*d_depart/deltahlg/cell%alphal(i)/cell%rhog(i)  
         !non-dimensional departure frequency
         nd_freq=4.06d0*nd_q**0.803d0                                              
         bfreq=nd_freq*cell%alphal(i)/d_depart**2.0d0       
      ENDSELECT
      IF(d_depart.le.0.0d0)bfreq=0.0d0
!
!.....Bubble wait time
!
      twait=0.8d0/bfreq
      IF(bfreq.le.0.0d0)twait=0.0d0      
!
!.....Two-phase area fraction
!      
      bub_area=pi*d_depart**2.0d0/4.0d0
      a_two=dmax1(0.0d0,dmin1(1.0d0,ndensity*bub_area*kfactor))   
      a_single=dmax1(1.0d-4,dmin1(1.0d0,1-a_two))              !Set minimum single phase area for stable calc.
!       
!.....Limit maximumum nucleation site density by area fraction
!
      IF(.false.)then
         IF(ndensity*bub_area.ge.1.0d0) ndensity=1.0d0/bub_area 
      ENDIF
!      
      IF(a_two.ge.1.d0)THEN
         ndensity=1.d0/kfactor/bub_area 
      ELSE
         dry_weight(i)=0.0d0
      ENDIF
      
!      ndensity=(1.0d0-relax_na)*ndensity+relax_na*nsiteden_o(i) 
!      nsiteden_o(i)=ndensity       
!      
!.....Quenching heat flux
!
      qq=2.0d0/sqrt_pi*dsqrt(twait*cell%lcondl(i)*cell%rhol(i)*cell%cpl(i))*bfreq*a_two*(tw-t_bulk(i))           
!
!.....Evaporation heat flux
!
      qe=ndensity*bfreq*pi/6.0d0*d_depart**3.0d0*cell%rhog(i)*deltahlg
      IF(fr_lift.eq.1) qe=ndensity*bfreq*pi/6.0d0*dlo**3.0d0*cell%rhog(i)*deltahlg *re_coeff               
      IF(udfl_psbt_cfx_model)THEN
         qcl=hconvl*a_single*deltaTl
! bug  qcg not calculater creates trouble in heat_partition_2_bi_lo
         qcg=0.d0
         q2=qq+qe+qcl 
      ELSE
         qcl=hconvl*a_single*deltaTl
         qcg=hconvg*1.0d0   *deltaTg
         q2=(1.0d0-dry_weight(i))*(qq+qe+qcl)+dry_weight(i)*qcg 
      ENDIF      
!    
      IF(nfcondition.eq.-2.or.nfcondition.le.-5)THEN
         heat_partition_2_bi_lo_err=q1-q2
      ELSE
         heat_partition_2_bi_lo_err=q2                  ! error = tatal heat flux for constant temp. BC
      ENDIF      
!
      END Function heat_partition_2_bi_lo_err
