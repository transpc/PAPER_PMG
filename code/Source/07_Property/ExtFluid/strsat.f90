      SUBROUTINE strsat_cupid(a,icf,x,ysat,dpsdts,err) 
!deck strsat                                                            
!                                                                       
!                                                                       
!      strsat  - find saturation value of pressure or temperature which 
!                corresponds to given value of temperature or pressure; 
!                also find derivative of saturation pressure with       
!                respect to temperature                                 
!                                                                       
!      Calling sequence:                                                
!                                                                       
!                call  strsat (rp1,ip2,rp3,rp4,rp5,lp6)                 
!                                                                       
!      Parameters:                                                      
!                                                                       
!                rp1 = a      = steam tables (input)                    
!                                                                       
!                ip2 = icf    = control flag (input)                    
!                               1 = find saturation pressure for given  
!                                   temperature                         
!                               2 = find saturation temperature for     
!                                   given pressure                      
!                                                                       
!                rp3 = x      = temperature or pressure for which       
!                               saturation pressure or temperature is   
!                               desired (input)                         
!                                                                       
!                rp4 = ysat   = saturation pressure or temperature      
!                               (output)                                
!                                                                       
!                rp5 = dpsdts = derivative of saturation pressure with  
!                               respect to temperature (output)         
!                                                                       
!                lp6 = err    = error flag  (output)                    
!                                                                       
!      Note:                                                            
!                                                                       
!        The form of the interpolation function is                      
!                                                                       
!                         log (P) = -A/T + B                            
!                                                                       
!        where P is the saturation pressure at temperature T, and A and 
!        B are constants which must be determined at each interpolation.
!                                                                       
!        To derive the given functional form, the following simplifying 
!        assumptions are made:                                          
!                                                                       
!        (a)  Latent heat-of-vaporization is a linear function of       
!             pressure (only) which goes to zero at the critical        
!             point, i.e., htvap = c1 * ( Pcrit - P ).                  
!                                                                       
!        (b)  The vapor obeys the ideal gas law, i.e, the vapor         
!             volume is vvap = c2 * T / P.                              
!                                                                       
!        (c)  Liquid volume does not depend on P (incompressible) and   
!             is an increasing linear function of T which becomes       
!             equal to the vapor volume at the critical point,          
!             i.e, vliq = c2 * T / Pcrit.                               
!                                                                       
!        Applying these assumptions to the Clausius-Clapeyron equation  
!        dPbydT = htvap / ( T * ( vvap - vliq ) ) gives                 
!                                                                       
!        dPbydT = c1 * ( Pcrit - P ) / ( c2 * ( 1/P - 1/Pcrit ) * T**2 )
!                                                                       
!        or, since 1/P - 1/Pcrit = ( Pcrit - P ) / ( Pcrit * P ),       
!                                                                       
!                          dPbydT = A * P / T**2                        
!                                                                       
!        where A = c1 * Pcrit / c2.  Thus,                              
!                                                                       
!                         1/P * dPbydT = A / T**2                       
!                                                                       
!        or, integrating with respect to T,                             
!                                                                       
!                            log (P) = -A/T + B                         
!                                                                       
!        which is the given functional form.                            
!                                                                       
!        If (T1,P1) and (T2,P2) are the two points on the saturation    
!        line which define the lower and upper limits, respectively, of 
!        the interpolation interval, then we have                       
!                                                                       
!          (1)  log (P2) = -A/T2 + B                                    
!                                                                       
!        and                                                            
!                                                                       
!          (2)  log (P1) = -A/T1 + B                                    
!                                                                       
!        from which, subtracting (2) from (1),                          
!                                                                       
!          log (P2) - log (P1) = A * ( 1/T1 - 1/T2 )                    
!                                                                       
!        or                                                             
!                                                                       
!          (3)  A = ( log (P2) - log (P1) ) * T1 * T2 / ( T2 - T1 )     
!                                                                       
!        Adding (1) to (2) gives                                        
!                                                                       
!          log (P1) + log (P2) = -A * ( 1/T1 + 1/T2 ) + 2 * B           
!                                                                       
!        or                                                             
!                                                                       
!          (4)  B = 0.5 * ( log (P1) + log (P2)                         
!                                    + A * ( T1 + T2 ) / ( T1 * T2 ) )  
!                                                                       
!        which gives a better "average" value for B than just           
!        substituting (3) into either (1) or (2).                       
!                                                                       
      USE STM_TBL_cupid  , ONLY: nt,np,pcrit
!
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL(8) a(*),x,ysat,dpsdts 
      INTEGER icf 
      LOGICAL err 
      REAL(8) xmin,xmax,xtrip,xcrit,ystrip,yscrit,xlow,xhigh,x1,x2,ysx1,&
      ysx2,x1s,x2s,ysx1s,ysx2s,xa,ysata,xb,ysatb                        
      REAL(8) ac,bc,p,p1,p1log,p2,p2log,t,t1,t1t2,t2 
      INTEGER ix,iytbl,nx,kxs,kysat,nsat,l,i0,i1,i2,i 
!
!------------------------------------------------------------------------------
!
!      GAS Modification by Won-Jae Lee: MUST in case New GAS is added
!                                                                 
      INCLUDE 'stcom.h' 
      INCLUDE 'gastable.h'
!   D. LMR-K.S. Ha for liquid metal properties- lead-bismuth eutetic(nfluid=11)
      INCLUDE 'lmtable.h'
      INTEGER ntg,npg,nstg,nspg,it3bpg,it4bpg,it5bpg,nprpntg,it3p0g 
! ------------------------------------------------------------------------
      ntg=nt 
      npg=np
      nstg=nst
      nspg=nsp
      it3bpg=it3bp
      it4bpg=it4bp
      it5bpg=it5bp
      nprpntg=nprpnt
      it3p0g=it3p0
!      IF(nfluid.eq.3) then
!         ntg=ntgc
!         npg=npgc
!         nstg=nsc
!         nspg=ns2c
!         it3bpg=klpc
!         it4bpg=klp2c
!         it5bpg=llpc
!         nprpntg=nt5c
!         it3p0g=jplc
!      ELSEIF(nfluid.eq.4) then
!         ntg=ntgh
!         npg=npgh
!         nstg=nsh
!         nspg=ns2h
!         it3bpg=klph
!         it4bpg=klp2h
!         it5bpg=llph
!         nprpntg=nt5h
!         it3p0g=jplh
!      ELSEIF(nfluid.eq.5) then
!         ntg=ntgh2
!         npg=npgh2
!         nstg=nsh2
!         nspg=ns2h2
!         it3bpg=klph2
!         it4bpg=klp2h2
!         it5bpg=llph2
!         nprpntg=nt5h2
!         it3p0g=jplh2
!      ELSEIF(nfluid.eq.6) then
!         ntg=ntgo
!         npg=npgo
!         nstg=nso
!         nspg=ns2o
!         it3bpg=klpo
!         it4bpg=klp2o
!         it5bpg=llpo
!         nprpntg=nt5o
!         it3p0g=jplo
!      ELSEIF(nfluid.eq.7) then
!         ntg=ntgn
!         npg=npgn
!         nstg=nsn
!         nspg=ns2n
!         it3bpg=klpn
!         it4bpg=klp2n
!         it5bpg=llpn
!         nprpntg=nt5n
!         it3p0g=jpln
!      ELSEIF(nfluid.eq.8) then
!         ntg=ntgna
!         npg=npgna
!         nstg=nsna
!         nspg=ns2na
!         it3bpg=klpna
!         it4bpg=klp2na
!         it5bpg=llpna
!         nprpntg=nt5na
!         it3p0g=jplna
!      ELSEIF(nfluid.eq.11) then
!         ntg=ntlbe
!         npg=nplbe
!         nstg=nslbe
!         nspg=ns2lbe
!         it3bpg=klplbe
!         it4bpg=klp2lbe
!         it5bpg=llplbe
!         nprpntg=nt5lbe
!         it3p0g=jpllbe
!      ENDIF
! ------------------------------------------------------------------------
!                                                                       
!                                                                       
!--initialize variables and pointers according to control flag          
!                                                                       
!      xmin   = minimum allowed value of x                              
!      xmax   = maximum allowed value of x                              
!      xtrip  = triple point value of x                                 
!      xcrit  = critical point value of x                               
!      ystrip = triple point value of ysat                              
!      yscrit = critical point value of ysat                            
!      ix     = base pointer to table of x values                       
!      iytbl  = base pointer to table containing ysat values            
!      nx     = number of values of x which lie within the saturation   
!               region                                                  
!      kxs    = base pointer to saturation table of x values            
!      kysat  = base pointer to saturation table of ysat values         
!      nsat   = number of values of ysat variable which lie within the  
!               saturation region                                       
!                                                                       
      IF(icf.eq.1)then 
!                                                                       
         xmin=tmin 
         xmax=tmax 
         xtrip=ttrip 
         xcrit=tcrit 
         ystrip=ptrip 
         yscrit=pcrit 
         ix=0 
         iytbl=ntg 
         nx=nstg 
         kxs=it4bpg 
         kysat=it3bpg 
         nsat=nspg 
!                                                                       
      ELSEIF(icf.eq.2)then 
!                                                                       
         xmin=pmin 
         xmax=pmax 
         xtrip=ptrip 
         xcrit=pcrit 
         ystrip=ttrip 
         yscrit=tcrit 
         ix=nt 
         iytbl=0 
         nx=nspg 
         kxs=it3bpg 
         kysat=it4bpg 
         nsat=nstg 
!                                                                       
      ELSE 
!                                                                       
         GOTO 50 
!                                                                       
      ENDIF 
!                                                                       
!--check for valid x value                                              
!                                                                       
      xlow=max(xmin,xtrip) 
      xhigh=min(xmax,xcrit) 
      IF(x.lt.xlow.or.x.gt.xhigh) GOTO 50 
!                                                                       
!--search x table to find interval which contains x                     
!                                                                       
!--initialize pointers                                                  
!                                                                       
      l=kysat 
      i0=ix 
      i1=i0+1 
!                                                                       
!--if x lies below the first value in the x table, use the triple point 
!--values of x and ysat for the low point of the interpolation interval,
!--and the first table values of x and ysat for the high point of the   
!--interpolation interval                                               
!                                                                       
      IF(x.lt.a(i1))then 
         x1=xtrip 
         x2=a(i1) 
         ysx1=ystrip 
         ysx2=a(l+13) 
         GOTO 20 
      ENDIF 
!                                                                       
!--x is .ge. the lowest value in the x table                            
!                                                                       
      i1=i1+1 
      i2=i0+nx 
      DO 10 i=i1,i2 
         IF(a(i).le.x) GOTO 10 
         x1=a(i-1) 
         x2=a(i) 
         l=l+13*(i-i1+1) 
         ysx1=a(l) 
         ysx2=a(l+13) 
         GOTO 20 
   10 END DO 
!                                                                       
!--x is .ge. to the highest value in the x table;  if x is equal to the 
!--highest x table value, use the last two values in the x table        
!--for x1 and x2;  otherwise, try to use the critical value for x2;     
!--if the critical value cannot be used, set x2 to an illegally high    
!--value for possible use as an error flag later on                     
!                                                                       
      IF(x.eq.a(i2))then 
         x1=a(i2-1) 
         x2=a(i2) 
         l=l+13*(nx-1) 
         ysx1=a(l) 
         ysx2=a(l+13) 
      ELSE 
         x1=a(i2) 
         ysx1=a(l+13*nx) 
         IF(xcrit.le.xmax)then 
            x2=xcrit 
            ysx2=yscrit 
         ELSE 
            x2=xmax+1.0d0 
         ENDIF 
      ENDIF 
!                                                                       
!--search appropriate saturation (a.s.) table to find interval which    
!--contains x                                                           
!                                                                       
!--initialize pointers                                                  
!                                                                       
   20 l=iytbl 
      i0=kxs 
      i1=i0+13 
!                                                                       
!--if x lies below the first value in the a.s. table, use the triple    
!--point values of x and ysat for the low point of the interpolation    
!--interval, and the first a.s. table values of x and ysat for the high 
!--point of the interpolation interval                                  
!                                                                       
      IF(x.lt.a(i1))then 
         x1s=xtrip 
         x2s=a(i1) 
         ysx1s=ystrip 
         ysx2s=a(l+1) 
         GOTO 40 
      ENDIF 
!                                                                       
!--x is .ge. the lowest value in the a.s. table                         
!                                                                       
      i1=i1+13 
      i2=i0+13*nsat 
      DO 30 i=i1,i2,13 
         IF(a(i).le.x) GOTO 30 
         x1s=a(i-13) 
         x2s=a(i) 
         l=l+(i-i1)/13+1 
         ysx1s=a(l) 
         ysx2s=a(l+1) 
         GOTO 40 
   30 END DO 
!                                                                       
!--x is .ge. to the highest value in the a.s. table;  if x is equal to  
!--the highest a.s. table value, use the last two values in the a.s.    
!--table for x1s and x2s;  otherwise, try to use the critical value for 
!--x2s;  if the critical value cannot be used, set x2s to an illegally  
!--high value for possible use as an error flag later on                
!                                                                       
      IF(x.eq.a(i2))then 
         x1s=a(i2-13) 
         x2s=a(i2) 
         l=l+nsat-1 
         ysx1s=a(l) 
         ysx2s=a(l+1) 
      ELSE 
         x1s=a(i2) 
         ysx1s=a(l+nsat) 
         IF(xcrit.le.xmax)then 
            x2s=xcrit 
            ysx2s=yscrit 
         ELSE 
            x2s=xmax+1.0d0 
         ENDIF 
      ENDIF 
!                                                                       
!--verify that a valid x2 or x2s value is available for interpolation   
!                                                                       
   40 IF(min(x2,x2s).gt.xmax) GOTO 50 
!                                                                       
!--find best set of x and ysat values for interpolation                 
!                                                                       
      IF(x1.ge.x1s)then 
         xa=x1 
         ysata=ysx1 
      ELSE 
         xa=x1s 
         ysata=ysx1s 
      ENDIF 
!                                                                       
      IF(x2.le.x2s)then 
         xb=x2 
         ysatb=ysx2 
      ELSE 
         xb=x2s 
         ysatb=ysx2s 
      ENDIF 
!                                                                       
!--interpolate to find saturation value                                 
!                                                                       
      IF(icf.eq.1)then 
         t=x 
         t1=xa 
         t2=xb 
         p1=ysata 
         p2=ysatb 
         t1t2=t1*t2 
         p1log=log(p1) 
         p2log=log(p2) 
         ac=(p2log-p1log)*t1t2/(t2-t1) 
         bc=0.5d0*(p1log+p2log+ac*(t1+t2)/t1t2) 
         p=exp(bc-ac/t) 
         ysat=p 
      ELSE 
         p=x 
         p1=xa 
         p2=xb 
         t1=ysata 
         t2=ysatb 
         t1t2=t1*t2 
         p1log=log(p1) 
         p2log=log(p2) 
         ac=(p2log-p1log)*t1t2/(t2-t1) 
         bc=0.5d0*(p1log+p2log+ac*(t1+t2)/t1t2) 
         t=ac/(bc-log(p)) 
!  mod by Won-Jae Lee to set min temperature as tripple temperature
		 t=max(t,ttrip)
         ysat=t 
      ENDIF 
!                                                                       
!--get derivative of saturation pressure with respect to temperature    
!                                                                       
      dpsdts=ac*p/t**2 
!                                                                       
      err=.false. 
      GOTO 60 
!                                                                       
!--error                                                                
!                                                                       
   50 err=.true. 
!                                                                       
!--done                                                                 
!                                                                       
   60 RETURN 
!                                                                       
      END SUBROUTINE strsat_cupid                         
