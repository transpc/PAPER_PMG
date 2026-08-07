!
      SUBROUTINE choke_gctpm(gc,dgcdp,hfact,error) 
!
      USE STM_TBL_cupid  , ONLY: st_tbl,    &
                                 nt,ndxstd, &
                                 cvao,dcva
      USE Zrv_choke      , ONLY: cpgas,cppf0,cppg0,gamma,ploss,pvzero,        &
                                 pzero,rgas,rnc,sgas,sliq,svap,vsubf0,vsubg0, &
                                 xnc,xzero         
!
      IMPLICIT NONE
!
      INTEGER :: iter,n 
      LOGICAL :: error 
!
      REAL(8) :: coef,conv,cpmix,delp,delt,                                       &
                 dgcdp,dgmdp,dgtdp,dndp,dpmin,dsgdp,dsldp,dsvdp,dtdps,dxdp,dxdpl, &
                 dxdpv,dxdpv2,dxdp2,dxeqdp,gami,gamrat,gber,gc,gt,hfact,          &
                 pbar,pgas,ptmax,ptmin,pnew,poly,polyi,pterm,pthrot,pvap,         &
                 pvmin,ratio,rthrot,rthrtp,rvap,rzero,rzerop,                     &
                 sfgm,sfgm0,smix,szero,tref,vfg,vmixt,vmix0,                      &
                 vthrot,vthrtp,xeq,xgas,xmfvap,xvap,xvliq0
      REAL(8) :: ttt,press,vbarr,ubarr,hbarr,betaa,kpa,cpp,quall,psatt,vsubff,vsubgg, &
                 usubff,usubgg,hsubf,hsubg,betf,betg,kpaf,kpag,cppf,cppg
      REAL(8) argmnt 
      REAL(8) prop(36) 
!
      PARAMETER (conv=0.001d0) 
      PARAMETER (tref=273.15d0) 
      PARAMETER (pbar=1.0d5) 
      PARAMETER (pvmin=612.d0)       
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
      error=.false. 
!
!     ---------------------------------------      
!.....Properties at Upstream Conditions
!     --------------------------------------- 
!
!.....Vapor Mole Fraction & Partial Pressures.                        
!                                                                       
      xmfvap=pvzero/pzero 
      xmfvap=MAX(0.d0,MIN(1.d0,xmfvap)) 
      pvzero=MAX(pvmin,pvzero) 
!                                                                       
!.....Gas / Vapor Mixture Properties...                               
!                                                                       
      vmix0=1.d0/((1.d0/vsubg0)+rgas) 
      cpmix=xnc*cpgas+(1.d0-xnc)*cppg0 
      smix=xnc*sgas+(1.d0-xnc)*svap 
      sfgm0=svap-sliq 
!                                                                       
!.....Polytropic Expansion Coefficient...                             
!                                                                       
      gami=1.d0/gamma 
      poly=((1.d0-xzero)*cppf0/cpmix+1.d0)/((1.d0-xzero)*cppf0/cpmix+gami) 
!                                                                       
      polyi=1.d0/poly 
      gamrat=gamma/(gamma-1.d0) 
!                                                                       
!.....Two-Phase Mixture Properties...                                 
!                                                                       
      rzero=1.d0/(xzero*vmix0+(1.d0-xzero)*vsubf0) 
      szero=xzero*smix+(1.d0-xzero)*sliq 
!                                                                       
!.....Set terms that won't change.                                      
!                                                                       
      xvliq0=(1.d0-xzero)*vsubf0 
      rzerop=1.d0/(xvliq0+gamrat*xzero*vmix0) 
      pterm=pzero/rzerop+ploss/rzero 
      coef=xzero*vmix0*pzero**gami/gamma 
!                                                                       
      ptmax=MAX(pvmin,0.99d0*(pzero+MIN(0.d0,ploss)))       
!
!     --------------------------------------- 
!.....Set Minimum Value for Throat Pressure Iteration.                  
!     --------------------------------------- 
!                                                                       
      ptmin=max(pvmin,0.7d0*ptmax) 
      vmixt=vmix0*(pzero/ptmin)**gami 
      vthrot=xvliq0+xzero*vmixt 
      vthrtp=xvliq0+gamrat*xzero*vmixt 
      ptmin=pterm/(vthrtp+gamma*vthrot**2/(2.0d0*xzero*vmixt)) 
!                                                                       
      DO 100 n=1,20 
!                                                                       
!........mass flux from Bernoulli eq.                                    
!                                                                       
         vmixt=vmix0*(pzero/ptmin)**gami 
         rthrot=1.d0/(xvliq0+xzero*vmixt) 
         rthrtp=1.d0/(xvliq0+gamrat*xzero*vmixt) 
!                                                                       
         argmnt=pterm-ptmin/rthrtp 
         IF(argmnt.le.0.d0) THEN
            ptmax=ptmin 
            pnew=0.9d0*ptmin 
            GOTO 90 
         ENDIF 
!                                                                       
         gber=rthrot*SQRT(2.0d0*(pterm-ptmin/rthrtp)) 
!                                                                       
!........pressure at dG/dP = 0.                                          
!                                                                       
         pnew=(coef*gber**2)**(gamma/(gamma+1.d0)) 
!                                                                       
!........check convergence.                                              
!                                                                       
   90    IF(ABS(pnew-ptmin)/(pzero-pnew).gt.conv) THEN
            ptmin=MIN(pnew,ptmax) 
         ELSE 
            ptmin=pnew 
            GOTO 110 
         ENDIF 
!                                                                       
  100 ENDDO       
!    
!     ---------------------------------------  
!.....Iteration Loop for Throat Pressure.     
!     ---------------------------------------                           
!     
  110 CONTINUE       
      pthrot=MAX(ptmin,0.7d0*ptmax) 
      iter=1 
!
!.....Step 1:  Solve Bernoulli for Mass Flux.                           
!
   10 CONTINUE 
!                                                                       
      vmixt=vmix0*(pzero/pthrot)**gami 
      rthrot=1.d0/(xvliq0+xzero*vmixt) 
      rthrtp=1.d0/(xvliq0+gamrat*xzero*vmixt) 
!                                                                       
      argmnt=pterm-pthrot/rthrtp 
!
!.....check argument & decrease PTMAX if necessary.                  
!
      IF(argmnt.le.0.d0)then 
         ptmax=pthrot 
         pthrot=0.5d0*(ptmax+ptmin) 
         IF(iter.lt.20)then 
            iter=iter+1 
            GOTO 10 
         ELSE 
            error=.true. 
            RETURN 
         ENDIF 
      ENDIF 
!                                                                       
      gt=rthrot*sqrt(2.0d0*argmnt) 
!                                                                       
!.....Set Derivative of Mass Flux wrt Throat Pressure.                  
!                                                                       
      dgtdp=(xzero*vmixt*gt**2/(gamma*pthrot)-1.d0)*(rthrot/gt) 
!                                                                       
!.....Limit magnitude of derivative (don't let it go to zero).          
!                                                                       
      dgmdp=0.15d0*gt/MAX(1.0d+3,(pthrot-pzero)) 
      dgtdp=MIN(dgtdp,dgmdp)  
!      
!.....Step 2:  Compute Throat Equilibrium Quality.                      
!                                                                       
!.....Throat Saturation Properties: load pressure and zero quality
!                                                                       
      pvap=MAX(pvmin,xmfvap*pthrot) 
      pgas=MAX(pvmin,pthrot-pvap) 
      prop(2)=pvap 
      prop(9)=0.d0
!
      CALL sth2x2_cupid(st_tbl(ndxstd),st_tbl(ndxstd+nt),prop,error)
!
      rvap=1.d0/vsubgg 
      sliq=prop(25) 
      svap=prop(26) 
      dtdps=(vsubgg-vsubff)/(svap-sliq) 
      dsldp=-vsubff*betf+cppf*(dtdps/ttt) 
!                                                                       
!.....Non-Condensible Gas Properties...                               
!                                                                       
      rgas=pgas/(rnc*ttt) 
      delt=MAX(0.d0,ttt-250.0d0) 
      cpgas=cvao+dcva*delt+rnc 
      sgas=cpgas*LOG(ttt/tref)-rnc*LOG(pgas/pbar) 
!                                                                       
!.....Gas / Vapor Mixture Properties...                               
!                                                                       
      cpmix=xnc*cpgas+(1.d0-xnc)*cppg 
!	                                                                      
      sfgm=svap-sliq 
      vfg=vsubgg-vsubff 
!                                                                       
!.....Mass Fraction of NC Gas (constant):                             
!                                                                       
      xgas=xnc*xzero 
!                                                                       
!.....Mass Fraction of Vapor:                                         
!                                                                       
      xvap=(szero-xgas*sgas-(1.d0-xgas)*sliq)/sfgm 
      xvap=MAX(0.d0,xvap) 
!                                                                       
      xeq=xvap+xgas 
      xeq=MAX(xzero,MIN(1.d0,xeq)) 
      xvap=xeq-xgas 
!                                                                       
!.....Apply Non-Equilibrium Factor & set Derivative wrt Pressure.       
!                                                                       
      IF(xeq.gt.hfact)then 
!                                                                       
!........Use "Equilibrium" Option:                                       
!                                                                       
         dndp=0.d0
         ratio=1.d0 
!                                                                       
      ELSEIF(hfact.lt.1000.d0)then 
!                                                                       
!........Apply Non-Equilibrium Parameter:                                
!                                                                       
         dsgdp=(dtdps/ttt)*(cpgas-rnc*LOG(pgas/pbar))-(rgas/pthrot) 
         dsvdp=-vsubgg*betg+cppg*(dtdps/ttt) 
!                                                                       
         dxeqdp=-xmfvap*((1.d0-xeq)*dsldp+xvap*dsvdp+xgas*dsgdp)/sfgm 
         dxeqdp=MIN(0.d0,dxeqdp) 
!                                                                       
         dndp=dxeqdp/hfact 
!                                                                       
         ratio=MIN(1.d0,xeq/hfact) 
!                                                                       
      ELSE 
!                                                                       
!........Use "Frozen" Option:                                            
!                                                                       
         dndp=0.d0
         ratio=0.d0
!                                                                       
      ENDIF                  
!                                                                       
!.....Step 3:  Critical Mass Flux.                                      
!                                                                       
!.....Pressure Part of Quality Derivative.                              
!                                                                       
      dxdpl=-xmfvap*dsldp/sfgm 
      dxdpv=-cpmix*(polyi-gami)/(pthrot*sfgm0) 
!                                                                       
!.....Derivative of Quality wrt Pressure.                               
!                                                                       
      dxdp=(1.d0-xzero)*ratio*dxdpl+xzero*dxdpv 
      dxdp=MIN(0.d0,dxdp) 
!                                                                       
!.....Critical Mass Flux.                                               
!                                                                       
      gc=SQRT(1.d0/(xzero*vmixt/(poly*pthrot)-vfg*dxdp)) 
!                                                                       
!.....Approximate Derivative wrt Pressure.                              
!                                                                       
      dxdpv2=-dxdpv/pthrot 
      dxdp2=(1.d0-xzero)*(dxdpl*dndp)+xzero*dxdpv2 
!                                                                       
      dgcdp=0.5d0*gc**3*((vmixt/pthrot)*(xzero*(1.d0+gami)/(poly*pthrot)-gami*dxdp)+vfg*dxdp2)
!                                                                       
!.....Step 4:  Check Convergence & Update Throat Pressure.              
!                                                                       
!.....Check Convergence.                                                
!                                                                       
      IF(ABS(gt-gc)/gc.gt.0.001d0)then 
!                                                                       
         delp=(gt-gc)/(dgcdp-dgtdp) 
!                                                                       
!........Set iteration limits.                                             
!                                                                       
         IF(delp.gt.0.d0)then 
            ptmin=pthrot 
            pthrot=pthrot+MIN(delp,0.9d0*(ptmax-ptmin)) 
         ELSE 
            ptmax=pthrot 
            pthrot=pthrot+MAX(delp,0.9d0*(ptmin-ptmax)) 
         ENDIF 
!                                                                       
!........check pressure range.                                          
!                                                                       
         dpmin=0.5d0*conv*gc/MAX(1.0d-4,dgcdp) 
         IF((ptmax-ptmin).gt.dpmin)then 
            iter=iter+1 
            IF(iter.lt.20)then 
!            
!..............continue iterating on throat pressure.                    
!
               GOTO 10 
!               
            ELSE 
!            
!..............solution did not converge, set error flag.                
!
               error=.true. 
!               
            ENDIF 
         ENDIF 
!                                                                       
      ENDIF 
!                                                                       
!.....Use derivative from momentum eqn.                                 
!                                                                       
      dgcdp=ABS(dgtdp)       
!
      END SUBROUTINE choke_gctpm
