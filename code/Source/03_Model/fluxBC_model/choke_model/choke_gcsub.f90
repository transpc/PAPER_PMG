!
      SUBROUTINE choke_gcsub(gc,dgcdp,hfact,error) 
!
      USE STM_TBL_cupid , ONLY: st_tbl,             &
                                nt,np,ns,ns2,ndxstd
      USE Zrv_choke     , ONLY: ploss,pzero,sliq,svap,tzero
!
      IMPLICIT NONE
!            
!    Define Local Variables.                                             
!                                                                       
      INTEGER :: it,iter   
      LOGICAL :: error 
!      
      REAL(8) :: conv,delg,delp,delp1,                                              &
                 dgcdp,dgcdp1,dgcdp2,dgtdp,dndp,dpfi,dpmax,dpmin,dpsub,dsldp,dsvdp, &
                 dtdps,dvldp,dvvdp,dxdp,dxeqdp,gc,gt,hfact,ptmax,                   &
                 ptmin,psubc,pthrot,ratio,rzero,                                    &
                 sfg,szero,vfg,xeq
      REAL(8) :: ttt,press,vbarr,ubarr,hbarr,betaa,kpa,cpp,quall,psatt,vsubff,vsubgg, &
                 usubff,usubgg,hsubf,hsubg,betf,betg,kpaf,kpag,cppf,cppg
!
      REAL(8) :: prop(36) 
!
      PARAMETER (conv=0.001d0)                    
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
!.....Stagnation Properties:                                            
!
!.....Load in upstream pressure and temperature.                      
!                                                                       
!     initialize prop
      prop(:)=0.d0
      prop(2)=pzero 
      prop(1)=tzero 
      CALL sth2x3_cupid(prop,it,error,                     &  
                        st_tbl(ndxstd),                    &
                        st_tbl(ndxstd+nt),                 &
                        st_tbl(ndxstd+nt+np+13*ns+13*ns2))
      IF(error) RETURN
!
!.....Set liquid density, entropy, and saturation pressure at liquid temp.
!       
      rzero=1.d0/vbarr 
      szero=prop(24)       
!
!.....Since we now get tsat(P) instead of psat(T) for x3 calls             
!     Get the saturation pressure by a call to x2
!                                                                       
      prop(9)=1.d0 
      CALL sth2x1_cupid(st_tbl(ndxstd),st_tbl(ndxstd+nt),prop,error)
      psubc=prop(10)

!     
!.....Iteration Loop for Throat Pressure.                               
!
!.....Set initial bounds & guess.                                    
!                                                                       
      iter=1 
      ptmax=MIN(psubc,0.99d0*(pzero+MIN(0.d0,ploss))) 
      ptmin=6.12d+2 
!                                                                       
!.....Use Henry-Fauske at saturation to limit pressure undershoot.    
!                                                                       
      dpfi=2.571d-2*pzero**(1.1334d0) 
      dpsub=pzero-psubc 
      dpsub=(dpfi**4+dpsub**4)**0.25d0 
      pthrot=MAX(ptmin,(pzero+MIN(0.d0,ploss)-dpsub))       
!  
!----------------------------------------------------   
!.....Step 1:  Solve Bernoulli for Mass Flux.                           
!----------------------------------------------------
!     
   10 CONTINUE 
!
      gt=SQRT(2.d0*rzero*(pzero+ploss-pthrot)) 
!                                                                       
!.....Set Derivative of Mass Flux wrt Throat Pressure.                  
!                                                                       
      dgtdp=-rzero/gt 
!
!----------------------------------------------------
!.....Step 2:  Compute Throat Equilibrium Quality.  
!----------------------------------------------------
!  
!.....Throat Saturation Properties:                                     
!     Load in pressure and zero quality                            
!                                                                       
      prop(2)=pthrot 
      prop(9)=0.d0
      CALL sth2x2_cupid(st_tbl(ndxstd),st_tbl(ndxstd+nt),prop,error)
      IF(error) RETURN
!
      sliq=prop(25) 
      svap=prop(26) 
      vfg=vsubgg-vsubff 
      sfg=svap-sliq 
      dtdps=vfg/sfg 
      dsldp=-vsubff*betf+cppf*(dtdps/ttt) 
!                                                                       
      xeq=(szero-sliq)/sfg 
      xeq=MAX(1.d-8,MIN(1.d0,xeq))       
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
         dsvdp=-vsubgg*betg+cppg*(dtdps/ttt) 
         dxeqdp=-(xeq*dsvdp+(1.d0-xeq)*dsldp)/sfg 
         dxeqdp=MIN(0.d0,dxeqdp) 
         dndp=dxeqdp/hfact 
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
!----------------------------------------------------
!.....Step 3:  Critical Mass Flux.          
!----------------------------------------------------
!
      dvldp=-vsubff*kpaf 
!                                                                       
!.....Derivative of Quality wrt Pressure.                               
!                                                                       
      dxdp=-ratio*(dsldp/sfg) 
!                                                                       
!.....Critical Mass Flux.                                               
!                                                                       
      gc=sqrt(-1.d0/(dvldp+vfg*dxdp)) 
!                                                                       
!.....Approximate Derivative wrt Pressure.                              
!                                                                       
      dvvdp=vsubgg*(betg*dtdps-kpag) 
!                                                                       
      dgcdp=0.5d0*(vfg*dndp*(-dsldp/sfg)+dxdp*dvvdp)*gc**3 
!
!----------------------------------------------------
!     Step 4:  Check Convergence & Update Throat Pressure.  
!----------------------------------------------------
!      
!.....Check Convergence.  
!
      dpmax=ptmax-ptmin 
      delg=ABS(gt-gc) 
      delg=MIN(delg,dpmax*ABS(dgtdp)) 
!                                                                       
      IF(delg/gt.gt.conv)then 
!                                                                       
         delp1=(gt-gc)/(dgcdp-dgtdp) 
         dgcdp2=3.d0*dgcdp/gc 
         dgcdp1=MAX(0.25d0,MIN(1.5d0,(1.d0+0.5d0*dgcdp2*delp1)))*dgcdp
         delp=(gt-gc)/(dgcdp1-dgtdp)
!                                                                       
!........Set iteration limits.                                             
!                                                                       
         dpmin=ABS(conv*gt/dgtdp) 
!                                                                       
         IF(delp.gt.0.d0) THEN
            ptmin=pthrot 
            dpmax=ptmax-ptmin 
            pthrot=pthrot+MIN(delp,0.95d0*dpmax) 
         ELSE 
            ptmax=pthrot 
            dpmax=ptmax-ptmin 
            IF(ABS(gt-gc)/gt.gt.0.25d0) THEN
               pthrot=pthrot+MAX(MIN(-dpmin,delp),-0.95d0*dpmax) 
            ELSE 
               pthrot=pthrot+MAX(delp,-0.95d0*dpmax)
            ENDIF 
         ENDIF 
!                                                                       
         iter=iter+1 
         IF(iter.lt.20)then 
!
!...........Continue iterating on throat pressure.                       
!
            GOTO 10 
!            
         ELSE 
!         
!...........did not converge, set error flag and default value.          
!
            error=.true. 
            dpsub=MIN(dpsub,pzero+ploss) 
            gt=SQRT(2.d0*rzero*dpsub) 
            dgtdp=rzero/gt 
         ENDIF 
!                                                                       
      ENDIF 
!                                                                       
!.....Use Mass Flux from Momentum Solution:  the expression for         
!     the critical mass flux can be very sensitive to throat        
!     pressure and its derivative VERY LARGE.                       
!                                                                       
      gc=gt 
      dgcdp=ABS(dgtdp)       
!
      END SUBROUTINE choke_gcsub
