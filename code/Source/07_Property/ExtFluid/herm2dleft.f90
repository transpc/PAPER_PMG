      SUBROUTINE herm2dleft_cupid(tables,ptr,a,arg,g,t,p,err) 
!define win32dvf                                                        
!define erf                                                             
!define fourbyt                                                         
!define hconden                                                         
!define impnon                                                          
!define in32                                                            
!define newnrc                                                          
!define ploc                                                            
!define sphaccm                                                         
!define unix                                                            
!define noselap                                                         
!define noextvol                                                        
!define noextv20                                                        
!define noextsys                                                        
!define noextjun                                                        
!define noextj20                                                        
!define noparcs                                                         
!define nonpa                                                           
!define nomap                                                           
!define logp                                                            
!deck herm2dleft                                                        
!                                                                       
!  $Id: herm2dleft.ff,v 1.1 2001/02/01 23:15:49 r5qa Exp dbarber $      
!                                                                       
!  Evaluates the 2D Hermite interpolating polynomial and several        
!  derivatives                                                          
!  This is a special version for extrapolating to the left assuming     
!  gpppp, gtttt, and cross products and higher derivatives are 0.0      
!  It is equivalent to setting the fifth and sixth rows and columns     
!  of the coefficient matrix to 0.0.  gam 1/10/2001                     
!                                                                       
!  Cognizant engineer:  rwt 5/25/1999                                   
!                                                                       
      IMPLICIT none 
!                                                                       
!***********************************************************************
!                                                                       
!  Declarations                                                         
!                                                                       
      INTEGER ptr 
      REAL(8) xx(6),yy(6),tables(*) 
!gam       real x,x2,y,y1,gfunct,gtfunct,gttfunct,                      
!gam     &      gpfunct,gptfunct,gppfunct                               
      REAL(8) x,x2,y,y1,gfuncte_cupid,gtfuncte_cupid,gttfuncte_cupid,gpfuncte_cupid,gptfuncte_cupid,  &
      gppfuncte_cupid                                                         
!if def,old                                                             
!      real diff1,diff2,diff3,diff4,diff5,diff6                         
!      real randy1,randy2,randy3,randy4,randy5,randy6                   
!endif                                                                  
      REAL(8) a,g,p,t,pii,pij,pji,pjj,tii,tij,tji,tjj 
      REAL(8) deltt,deltp,gtilde,dgtildedt,dgtildedp,d2gtildedt2,       &
      d2gtildedp2,d2gtildedpdt                                          
!if def,swesty                                                          
!      real psi0,psi1,psi2,z                                            
!      real dgdtii,dgdtij,dgdtji,                                       
!    &      dgdtjj,dgdpii,dgdpij,dgdpji,dgdpjj,d2gdp2ii,                
!    &      d2gdp2ij,d2gdp2ji,d2gdp2jj,d2gdt2ii,                        
!    &      d2gdt2ij,d2gdt2ji,d2gdt2jj,d2gdtdpii,                       
!    &      d2gdtdpij,d2gdtdpji,d2gdtdpjj,d3gdt2dpii,                   
!    &      d3gdt2dpij,d3gdt2dpji,d3gdt2dpjj,d3gdtdp2ii,                
!    &      d3gdtdp2ij,d3gdtdp2ji,d3gdtdp2jj,d4gdt2dp2ii,               
!    &      d4gdt2dp2ij,d4gdt2dp2ji,d4gdt2dp2jj,dpsi0,dpsi1,            
!    &      dpsi2,d2psi0,d2psi1,d2psi2,d3psi0,d3psi1,d3psi2,gii,gij,    
!    &      gji,gjj,giit,gijt,gjit,gjjt,giitt,gijtt,gjitt,              
!    &      gjjtt,giip,gijp,gjip,gjjp,giipp,gijpp,gjipp,gjjpp,          
!    &      giitp,gijtp,gjitp,gjjtp,giittp,gijttp,gjittp,               
!    &      gjjttp,giitpp,gijtpp,gjitpp,gjjtpp,giittpp,gijttpp,         
!    &      gjittpp,gjjttpp                                             
!      real psi0x,psi1x,psi2x,dpsi0x,dpsi1x,dpsi2x,d2psi0x,             
!    &      d2psi1x,d2psi2x,                                            
!    &      psi0x2,psi1x2,psi2x2,dpsi0x2,dpsi1x2,dpsi2x2,d2psi0x2,      
!    &      d2psi1x2,d2psi2x2,                                          
!    &      psi0y,psi1y,psi2y,dpsi0y,dpsi1y,dpsi2y,d2psi0y,             
!    &      d2psi1y,d2psi2y,                                            
!    &      psi0y1,psi1y1,psi2y1,dpsi0y1,dpsi1y1,dpsi2y1,d2psi0y1,      
!    &      d2psi1y1,d2psi2y1                                           
!      real d3gtildedt3,d3gtildedpdt2,d3gtildedp2dt,d4gtildedp2dt2      
!endif                                                                  
      DIMENSION a(11,4),g(10) 
      LOGICAL err 
      CHARACTER*(*)arg 
!if def,old                                                             
!      data diff1 / 0 /                                                 
!      data diff2 / 0 /                                                 
!      data diff3 / 0 /                                                 
!      data diff4 / 0 /                                                 
!      data diff5 / 0 /                                                 
!      data diff6 / 0 /                                                 
!      save diff1                                                       
!      save diff2                                                       
!      save diff3                                                       
!      save diff4                                                       
!      save diff5                                                       
!      save diff6                                                       
!endif                                                                  
!***********************************************************************
!                                                                       
!  Statement functions to define Hermite polynomials                    
!                                                                       
!if def,swesty                                                          
!      psi0(z) = 1.0 + (z**3.0)*(-10.0 + z*(15.0 - z*(6.0)))            
!      psi1(z) = z*(1.0 + z**2*(-6.0 + z*(8.0 + z*(-3.0))))             
!      psi2(z) = (0.5*z**2)*(1.0 + z*(-3.0 + z*(3.0 + z*(-1.0))))       
!                                                                       
!  First derivatives                                                    
!                                                                       
!      dpsi0(z) = z**2*(-30.0 + z*(60.0 + z*(-30.0)))                   
!      dpsi1(z) = 1.0 + z**2*(-18.0 + z*(32.0 + z*(-15.0)))             
!      dpsi2(z) = .5*z*(2.0 + z*(-9.0 + z*(12.0 + z*(-5.0))))           
!                                                                       
!  Second derivatives                                                   
!                                                                       
!      d2psi0(z) = z*(-60.0 + z*(180.0 + z*(-120.0)))                   
!      d2psi1(z) = z*(-36.0 + z*(96.0 + z*(-60.0)))                     
!      d2psi2(z) = 0.5*(2.0 + z*(-18.0 + z*(36.0 + z*(-20.0))))         
!                                                                       
!  Third derivatives                                                    
!                                                                       
!      d3psi0(z) = -60.0 + z*(360.0 + z*(-360.0))                       
!      d3psi1(z) = -36.0 + z*(192.0 + z*(-180.0))                       
!      d3psi2(z) = 0.5*(-18.0 + z*(72.0 + z*(-60.0)))                   
!endif                                                                  
!***********************************************************************
!                                                                       
!  Definitions                                                          
!                                                                       
!      call timstart ('herm2dleft')                                     
      err=.false. 
!                                                                       
      tii=a(1,1) 
      tij=a(1,2) 
      tji=a(1,3) 
      tjj=a(1,4) 
!if def,swesty                                                          
!      gii = a(3,1)                                                     
!      gij = a(3,2)                                                     
!      gji = a(3,3)                                                     
!      gjj = a(3,4)                                                     
!endif                                                                  
      pii=a(2,1) 
      pij=a(2,2) 
      pji=a(2,3) 
      pjj=a(2,4) 
!                                                                       
      deltt=(tij-tii) 
      deltp=(pji-pii) 
      x=(t-tii)/deltt 
      x2=1.0d0-x 
      y=(p-pii)/deltp 
      y1=1.0d0-y 
!if def,logp                                                            
      IF(arg.eq.'vapor')then 
         deltp=(log(pji)-log(pii)) 
         y=(log(p)-log(pii))/deltp 
         y1=1.0d0-y 
      ENDIF 
!endif                                                                  
!                                                                       
!if def,swesty                                                          
!      dgdtii      = a(6,1)                                             
!      dgdtij      = a(6,2)                                             
!      dgdtji      = a(6,3)                                             
!      dgdtjj      = a(6,4)                                             
!      dgdpii      = a(4,1)                                             
!      dgdpij      = a(4,2)                                             
!      dgdpji      = a(4,3)                                             
!      dgdpjj      = a(4,4)                                             
!      d2gdp2ii    = a(5,1)                                             
!      d2gdp2ij    = a(5,2)                                             
!      d2gdp2ji    = a(5,3)                                             
!      d2gdp2jj    = a(5,4)                                             
!      d2gdt2ii    = a(9,1)                                             
!      d2gdt2ij    = a(9,2)                                             
!      d2gdt2ji    = a(9,3)                                             
!      d2gdt2jj    = a(9,4)                                             
!      d2gdtdpii   = a(7,1)                                             
!      d2gdtdpij   = a(7,2)                                             
!      d2gdtdpji   = a(7,3)                                             
!      d2gdtdpjj   = a(7,4)                                             
!      d3gdt2dpii  = a(10,1)                                            
!      d3gdt2dpij  = a(10,2)                                            
!      d3gdt2dpji  = a(10,3)                                            
!      d3gdt2dpjj  = a(10,4)                                            
!      d3gdtdp2ii  = a(8,1)                                             
!      d3gdtdp2ij  = a(8,2)                                             
!      d3gdtdp2ji  = a(8,3)                                             
!      d3gdtdp2jj  = a(8,4)                                             
!      d4gdt2dp2ii = a(11,1)                                            
!      d4gdt2dp2ij = a(11,2)                                            
!      d4gdt2dp2ji = a(11,3)                                            
!      d4gdt2dp2jj = a(11,4)                                            
!endif                                                                  
!if def,logp                                                            
      IF(arg.eq.'vapor')then 
!  gyy                                                                  
!if def,swesty                                                          
!        d2gdp2ii = d2gdp2ii*pii**2 + dgdpii*pii                        
!        d2gdp2ij = d2gdp2ij*pij**2 + dgdpij*pij                        
!        d2gdp2ji = d2gdp2ji*pji**2 + dgdpji*pji                        
!        d2gdp2jj = d2gdp2jj*pjj**2 + dgdpjj*pjj                        
!  gyyt                                                                 
!        d3gdtdp2ii = d3gdtdp2ii*pii**2 + d2gdtdpii*pii                 
!        d3gdtdp2ij = d3gdtdp2ij*pij**2 + d2gdtdpij*pij                 
!        d3gdtdp2ji = d3gdtdp2ji*pji**2 + d2gdtdpji*pji                 
!        d3gdtdp2jj = d3gdtdp2jj*pjj**2 + d2gdtdpjj*pjj                 
!  gyytt                                                                
!        d4gdt2dp2ii = d4gdt2dp2ii*pii**2 + d3gdt2dpii*pii              
!        d4gdt2dp2ij = d4gdt2dp2ij*pij**2 + d3gdt2dpij*pij              
!        d4gdt2dp2ji = d4gdt2dp2ji*pji**2 + d3gdt2dpji*pji              
!        d4gdt2dp2jj = d4gdt2dp2jj*pjj**2 + d3gdt2dpjj*pjj              
!  gytt                                                                 
!        d3gdt2dpii = d3gdt2dpii*pii                                    
!        d3gdt2dpij = d3gdt2dpij*pij                                    
!        d3gdt2dpji = d3gdt2dpji*pji                                    
!        d3gdt2dpjj = d3gdt2dpjj*pjj                                    
!  gyt                                                                  
!        d2gdtdpii = d2gdtdpii*pii                                      
!        d2gdtdpij = d2gdtdpij*pij                                      
!        d2gdtdpji = d2gdtdpji*pji                                      
!        d2gdtdpjj = d2gdtdpjj*pjj                                      
!  gy                                                                   
!        dgdpii = dgdpii*pii                                            
!        dgdpij = dgdpij*pij                                            
!        dgdpji = dgdpji*pji                                            
!        dgdpjj = dgdpjj*pjj                                            
!endif                                                                  
!                                                                       
         p=log(p) 
         pii=log(pii) 
         pij=log(pij) 
         pji=log(pji) 
         pjj=log(pjj) 
      ENDIF 
!endif                                                                  
!                                                                       
!if def,swesty                                                          
!      giit    = dgdtii*deltt                                           
!      gijt    = dgdtij*deltt                                           
!      gjit    = dgdtji*deltt                                           
!      gjjt    = dgdtjj*deltt                                           
!      giitt   = d2gdt2ii*deltt**2                                      
!      gijtt   = d2gdt2ij*deltt**2                                      
!      gjitt   = d2gdt2ji*deltt**2                                      
!      gjjtt   = d2gdt2jj*deltt**2                                      
!      giip    = dgdpii*deltp                                           
!      gijp    = dgdpij*deltp                                           
!      gjip    = dgdpji*deltp                                           
!      gjjp    = dgdpjj*deltp                                           
!      giipp   = d2gdp2ii*deltp**2                                      
!      gijpp   = d2gdp2ij*deltp**2                                      
!      gjipp   = d2gdp2ji*deltp**2                                      
!      gjjpp   = d2gdp2jj*deltp**2                                      
!      giitp   = d2gdtdpii*deltp*deltt                                  
!      gijtp   = d2gdtdpij*deltp*deltt                                  
!      gjitp   = d2gdtdpji*deltp*deltt                                  
!      gjjtp   = d2gdtdpjj*deltp*deltt                                  
!      giittp  = d3gdt2dpii*deltp*deltt**2                              
!      gijttp  = d3gdt2dpij*deltp*deltt**2                              
!      gjittp  = d3gdt2dpji*deltp*deltt**2                              
!      gjjttp  = d3gdt2dpjj*deltp*deltt**2                              
!      giitpp  = d3gdtdp2ii*deltt*deltp**2                              
!      gijtpp  = d3gdtdp2ij*deltt*deltp**2                              
!      gjitpp  = d3gdtdp2ji*deltt*deltp**2                              
!      gjjtpp  = d3gdtdp2jj*deltt*deltp**2                              
!      giittpp = d4gdt2dp2ii*(deltt**2)*(deltp**2)                      
!      gijttpp = d4gdt2dp2ij*(deltt**2)*(deltp**2)                      
!      gjittpp = d4gdt2dp2ji*(deltt**2)*(deltp**2)                      
!      gjjttpp = d4gdt2dp2jj*(deltt**2)*(deltp**2)                      
!                                                                       
!      psi0x  = psi0(x)                                                 
!      psi1x  = psi1(x)                                                 
!      psi2x  = psi2(x)                                                 
!      psi0x2 = psi0(x2)                                                
!      psi1x2 = psi1(x2)                                                
!      psi2x2 = psi2(x2)                                                
!      psi0y  = psi0(y)                                                 
!      psi1y  = psi1(y)                                                 
!      psi2y  = psi2(y)                                                 
!      psi0y1 = psi0(y1)                                                
!      psi1y1 = psi1(y1)                                                
!      psi2y1 = psi2(y1)                                                
!                                                                       
!      dpsi0x  = dpsi0(x)                                               
!      dpsi1x  = dpsi1(x)                                               
!      dpsi2x  = dpsi2(x)                                               
!      dpsi0x2 = dpsi0(x2)                                              
!      dpsi1x2 = dpsi1(x2)                                              
!      dpsi2x2 = dpsi2(x2)                                              
!      dpsi0y  = dpsi0(y)                                               
!      dpsi1y  = dpsi1(y)                                               
!      dpsi2y  = dpsi2(y)                                               
!      dpsi0y1 = dpsi0(y1)                                              
!      dpsi1y1 = dpsi1(y1)                                              
!      dpsi2y1 = dpsi2(y1)                                              
!                                                                       
!      d2psi0x  = d2psi0(x)                                             
!      d2psi1x  = d2psi1(x)                                             
!      d2psi2x  = d2psi2(x)                                             
!      d2psi0x2 = d2psi0(x2)                                            
!      d2psi1x2 = d2psi1(x2)                                            
!      d2psi2x2 = d2psi2(x2)                                            
!      d2psi0y  = d2psi0(y)                                             
!      d2psi1y  = d2psi1(y)                                             
!      d2psi2y  = d2psi2(y)                                             
!      d2psi0y1 = d2psi0(y1)                                            
!      d2psi1y1 = d2psi1(y1)                                            
!      d2psi2y1 = d2psi2(y1)                                            
!endif                                                                  
!                                                                       
!***********************************************************************
! calculate gtilde                                                      
!                                                                       
!      if (arg .eq. 'vapor') then                                       
!if def,swesty                                                          
!      gtilde = gii*psi0x*psi0y                                         
!    &        + gij*psi0x2*psi0y                                        
!    &        + gji*psi0x*psi0y1                                        
!    &        + gjj*psi0x2*psi0y1                                       
!    &        + giit*psi1x*psi0y                                        
!    &        - gijt*psi1x2*psi0y                                       
!    &        + gjit*psi1x*psi0y1                                       
!    &        - gjjt*psi1x2*psi0y1                                      
!    &        + giitt*psi2x*psi0y                                       
!    &        + gijtt*psi2x2*psi0y                                      
!    &        + gjitt*psi2x*psi0y1                                      
!    &        + gjjtt*psi2x2*psi0y1                                     
!    &        + giip*psi0x*psi1y                                        
!    &        + gijp*psi0x2*psi1y                                       
!    &        - gjip*psi0x*psi1y1                                       
!    &        - gjjp*psi0x2*psi1y1                                      
!    &        + giipp*psi0x*psi2y                                       
!    &        + gijpp*psi0x2*psi2y                                      
!    &        + gjipp*psi0x*psi2y1                                      
!    &        + gjjpp*psi0x2*psi2y1                                     
!    &        + giitp*psi1x*psi1y                                       
!    &        - gijtp*psi1x2*psi1y                                      
!    &        - gjitp*psi1x*psi1y1                                      
!    &        + gjjtp*psi1x2*psi1y1                                     
!    &        + giittp*psi2x*psi1y                                      
!    &        + gijttp*psi2x2*psi1y                                     
!    &        - gjittp*psi2x*psi1y1                                     
!    &        - gjjttp*psi2x2*psi1y1                                    
!    &        + giitpp*psi1x*psi2y                                      
!    &        - gijtpp*psi1x2*psi2y                                     
!    &        + gjitpp*psi1x*psi2y1                                     
!    &        - gjjtpp*psi1x2*psi2y1                                    
!    &        + giittpp*psi2x*psi2y                                     
!    &        + gijttpp*psi2x2*psi2y                                    
!    &        + gjittpp*psi2x*psi2y1                                    
!    &        + gjjttpp*psi2x2*psi2y1                                   
!endif                                                                  
!if -def,swesty                                                         
      xx(1)=1.0d0 
      xx(2)=x 
      xx(3)=x**2 
      xx(4)=x**3 
      xx(5)=x**4 
      xx(6)=x**5 
      yy(1)=1.0d0 
      yy(2)=y 
      yy(3)=y**2 
      yy(4)=y**3 
      yy(5)=y**4 
      yy(6)=y**5 
!gam       gtilde = gfunct(yy,xx,tables,ptr)                            
      gtilde=gfuncte_cupid(yy,xx,tables,ptr) 
!endif                                                                  
!         if (abs ((gtilde-randy1)/gtilde) .gt. diff1) then             
!           diff1 = abs ((gtilde-randy1)/gtilde)                        
!           IF(myrank.eq.0)write(90,*) gtilde,randy1,diff1,x,y,p,t                    
!         endif                                                         
!       endif                                                           
!***********************************************************************
! calculate d(gtilde)/dt                                                
!                                                                       
!if def,swesty                                                          
!      dgtildedt = gii*psi0y*dpsi0x                                     
!    &           - gij*psi0y*dpsi0x2                                    
!    &           + gji*psi0y1*dpsi0x                                    
!    &           - gjj*psi0y1*dpsi0x2                                   
!    &           + giit*psi0y*dpsi1x                                    
!    &           + gijt*psi0y*dpsi1x2                                   
!    &           + gjit*psi0y1*dpsi1x                                   
!    &           + gjjt*psi0y1*dpsi1x2                                  
!    &           + giitt*psi0y*dpsi2x                                   
!    &           - gijtt*psi0y*dpsi2x2                                  
!    &           + gjitt*psi0y1*dpsi2x                                  
!    &           - gjjtt*psi0y1*dpsi2x2                                 
!    &           + giip*psi1y*dpsi0x                                    
!    &           - gijp*psi1y*dpsi0x2                                   
!    &           - gjip*psi1y1*dpsi0x                                   
!    &           + gjjp*psi1y1*dpsi0x2                                  
!    &           + giipp*psi2y*dpsi0x                                   
!    &           - gijpp*psi2y*dpsi0x2                                  
!    &           + gjipp*psi2y1*dpsi0x                                  
!    &           - gjjpp*psi2y1*dpsi0x2                                 
!    &           + giitp*psi1y*dpsi1x                                   
!    &           + gijtp*psi1y*dpsi1x2                                  
!    &           - gjitp*psi1y1*dpsi1x                                  
!    &           - gjjtp*psi1y1*dpsi1x2                                 
!    &           + giittp*psi1y*dpsi2x                                  
!    &           - gijttp*psi1y*dpsi2x2                                 
!    &           - gjittp*psi1y1*dpsi2x                                 
!    &           + gjjttp*psi1y1*dpsi2x2                                
!    &           + giitpp*psi2y*dpsi1x                                  
!    &           + gijtpp*psi2y*dpsi1x2                                 
!    &           + gjitpp*psi2y1*dpsi1x                                 
!    &           + gjjtpp*psi2y1*dpsi1x2                                
!    &           + giittpp*psi2y*dpsi2x                                 
!    &           - gijttpp*psi2y*dpsi2x2                                
!    &           + gjittpp*psi2y1*dpsi2x                                
!    &           - gjjttpp*psi2y1*dpsi2x2                               
!      dgtildedt = dgtildedt/deltt                                      
!endif                                                                  
!       if (arg .eq. 'vapor') then                                      
!         randy2 = gtfunct(y,x,c)                                       
!         randy2 = randy2/deltt                                         
!if -def, swesty                                                        
!gam       dgtildedt = gtfunct(yy,xx,tables,ptr)                        
      dgtildedt=gtfuncte_cupid(yy,xx,tables,ptr) 
      dgtildedt=dgtildedt/deltt 
!endif                                                                  
!         if (abs ((dgtildedt-randy2)/dgtildedt) .gt. diff2) then       
!           diff2 = abs ((dgtildedt-randy2)/dgtildedt)                  
!           IF(myrank.eq.0)write(91,*) dgtildedt,randy2,diff2,x,y,p,t                 
!         endif                                                         
!       endif                                                           
!***********************************************************************
! calculate d(gtilde)/dp                                                
!                                                                       
!if def,swesty                                                          
!      dgtildedp = gii*dpsi0y*psi0x                                     
!    &           + gij*dpsi0y*psi0x2                                    
!    &           - gji*dpsi0y1*psi0x                                    
!    &           - gjj*dpsi0y1*psi0x2                                   
!    &           + giit*dpsi0y*psi1x                                    
!    &           - gijt*dpsi0y*psi1x2                                   
!    &           - gjit*dpsi0y1*psi1x                                   
!    &           + gjjt*dpsi0y1*psi1x2                                  
!    &           + giitt*dpsi0y*psi2x                                   
!    &           + gijtt*dpsi0y*psi2x2                                  
!    &           - gjitt*dpsi0y1*psi2x                                  
!    &           - gjjtt*dpsi0y1*psi2x2                                 
!    &           + giip*dpsi1y*psi0x                                    
!    &           + gijp*dpsi1y*psi0x2                                   
!    &           + gjip*dpsi1y1*psi0x                                   
!    &           + gjjp*dpsi1y1*psi0x2                                  
!    &           + giipp*dpsi2y*psi0x                                   
!    &           + gijpp*dpsi2y*psi0x2                                  
!    &           - gjipp*dpsi2y1*psi0x                                  
!    &           - gjjpp*dpsi2y1*psi0x2                                 
!    &           + giitp*dpsi1y*psi1x                                   
!    &           - gijtp*dpsi1y*psi1x2                                  
!    &           + gjitp*dpsi1y1*psi1x                                  
!    &           - gjjtp*dpsi1y1*psi1x2                                 
!    &           + giittp*dpsi1y*psi2x                                  
!    &           + gijttp*dpsi1y*psi2x2                                 
!    &           + gjittp*dpsi1y1*psi2x                                 
!    &           + gjjttp*dpsi1y1*psi2x2                                
!    &           + giitpp*dpsi2y*psi1x                                  
!    &           - gijtpp*dpsi2y*psi1x2                                 
!    &           - gjitpp*dpsi2y1*psi1x                                 
!    &           + gjjtpp*dpsi2y1*psi1x2                                
!    &           + giittpp*dpsi2y*psi2x                                 
!    &           + gijttpp*dpsi2y*psi2x2                                
!    &           - gjittpp*dpsi2y1*psi2x                                
!    &           - gjjttpp*dpsi2y1*psi2x2                               
!      dgtildedp = dgtildedp/deltp                                      
!endif                                                                  
!       if (arg .eq. 'liquid') then                                     
!         randy3 = gpfunct(y,x,c)                                       
!         randy3 = randy3/deltp                                         
!if -def, swesty                                                        
!gam       dgtildedp = gpfunct(yy,xx,tables,ptr)                        
      dgtildedp=gpfuncte_cupid(yy,xx,tables,ptr) 
      dgtildedp=dgtildedp/deltp 
!endif                                                                  
!         if (abs ((dgtildedp-randy3)/dgtildedp) .gt. diff3) then       
!           diff3 = abs ((dgtildedp-randy3)/dgtildedp)                  
!           IF(myrank.eq.0)write(92,*) dgtildedp,randy3,diff3,x,y,p,t                 
!         endif                                                         
!       endif                                                           
!***********************************************************************
! calculate d2(gtilde)/dt2                                              
!if def,swesty                                                          
!      d2gtildedt2 = gii*psi0y*d2psi0x                                  
!    &             + gij*psi0y*d2psi0x2                                 
!    &             + gji*psi0y1*d2psi0x                                 
!    &             + gjj*psi0y1*d2psi0x2                                
!    &             + giit*psi0y*d2psi1x                                 
!    &             - gijt*psi0y*d2psi1x2                                
!    &             + gjit*psi0y1*d2psi1x                                
!    &             - gjjt*psi0y1*d2psi1x2                               
!    &             + giitt*psi0y*d2psi2x                                
!    &             + gijtt*psi0y*d2psi2x2                               
!    &             + gjitt*psi0y1*d2psi2x                               
!    &             + gjjtt*psi0y1*d2psi2x2                              
!    &             + giip*psi1y*d2psi0x                                 
!    &             + gijp*psi1y*d2psi0x2                                
!    &             - gjip*psi1y1*d2psi0x                                
!    &             - gjjp*psi1y1*d2psi0x2                               
!    &             + giipp*psi2y*d2psi0x                                
!    &             + gijpp*psi2y*d2psi0x2                               
!    &             + gjipp*psi2y1*d2psi0x                               
!    &             + gjjpp*psi2y1*d2psi0x2                              
!    &             + giitp*psi1y*d2psi1x                                
!    &             - gijtp*psi1y*d2psi1x2                               
!    &             - gjitp*psi1y1*d2psi1x                               
!    &             + gjjtp*psi1y1*d2psi1x2                              
!    &             + giittp*psi1y*d2psi2x                               
!    &             + gijttp*psi1y*d2psi2x2                              
!    &             - gjittp*psi1y1*d2psi2x                              
!    &             - gjjttp*psi1y1*d2psi2x2                             
!    &             + giitpp*psi2y*d2psi1x                               
!    &             - gijtpp*psi2y*d2psi1x2                              
!    &             + gjitpp*psi2y1*d2psi1x                              
!    &             - gjjtpp*psi2y1*d2psi1x2                             
!    &             + giittpp*psi2y*d2psi2x                              
!    &             + gijttpp*psi2y*d2psi2x2                             
!    &             + gjittpp*psi2y1*d2psi2x                             
!    &             + gjjttpp*psi2y1*d2psi2x2                            
!      d2gtildedt2 = d2gtildedt2/deltt**2                               
!endif                                                                  
!       if (arg .eq. 'vapor') then                                      
!         randy4 = gttfunct(y,x,c)                                      
!         randy4 = randy4/deltt**2                                      
!if -def, swesty                                                        
!gam       d2gtildedt2 = gttfunct(yy,xx,tables,ptr)                     
      d2gtildedt2=gttfuncte_cupid(yy,xx,tables,ptr) 
      d2gtildedt2=d2gtildedt2/deltt**2 
!endif                                                                  
!         if (abs ((d2gtildedt2-randy4)/d2gtildedt2) .gt. diff4) then   
!           diff4 = abs ((d2gtildedt2-randy4)/d2gtildedt2)              
!           IF(myrank.eq.0)write(93,*) d2gtildedt2,randy4,diff4,x,y,p,t               
!         endif                                                         
!       endif                                                           
!***********************************************************************
! calculate d2(gtilde)/dp2                                              
!                                                                       
!if def,swesty                                                          
!      d2gtildedp2 = gii*d2psi0y*psi0x                                  
!    &             + gij*d2psi0y*psi0x2                                 
!    &             + gji*d2psi0y1*psi0x                                 
!    &             + gjj*d2psi0y1*psi0x2                                
!    &             + giit*d2psi0y*psi1x                                 
!    &             - gijt*d2psi0y*psi1x2                                
!    &             + gjit*d2psi0y1*psi1x                                
!    &             - gjjt*d2psi0y1*psi1x2                               
!    &             + giitt*d2psi0y*psi2x                                
!    &             + gijtt*d2psi0y*psi2x2                               
!    &             + gjitt*d2psi0y1*psi2x                               
!    &             + gjjtt*d2psi0y1*psi2x2                              
!    &             + giip*d2psi1y*psi0x                                 
!    &             + gijp*d2psi1y*psi0x2                                
!    &             - gjip*d2psi1y1*psi0x                                
!    &             - gjjp*d2psi1y1*psi0x2                               
!    &             + giipp*d2psi2y*psi0x                                
!    &             + gijpp*d2psi2y*psi0x2                               
!    &             + gjipp*d2psi2y1*psi0x                               
!    &             + gjjpp*d2psi2y1*psi0x2                              
!    &             + giitp*d2psi1y*psi1x                                
!    &             - gijtp*d2psi1y*psi1x2                               
!    &             - gjitp*d2psi1y1*psi1x                               
!    &             + gjjtp*d2psi1y1*psi1x2                              
!    &             + giittp*d2psi1y*psi2x                               
!    &             + gijttp*d2psi1y*psi2x2                              
!    &             - gjittp*d2psi1y1*psi2x                              
!    &             - gjjttp*d2psi1y1*psi2x2                             
!    &             + giitpp*d2psi2y*psi1x                               
!    &             - gijtpp*d2psi2y*psi1x2                              
!    &             + gjitpp*d2psi2y1*psi1x                              
!    &             - gjjtpp*d2psi2y1*psi1x2                             
!    &             + giittpp*d2psi2y*psi2x                              
!    &             + gijttpp*d2psi2y*psi2x2                             
!    &             + gjittpp*d2psi2y1*psi2x                             
!    &             + gjjttpp*d2psi2y1*psi2x2                            
!      d2gtildedp2 = d2gtildedp2/deltp**2                               
!endif                                                                  
!       if (arg .eq. 'vapor') then                                      
!         randy5 = gppfunct(y,x,c)                                      
!         randy5 = randy5/deltp**2                                      
!if -def, swesty                                                        
!gam       d2gtildedp2 = gppfunct(yy,xx,tables,ptr)                     
      d2gtildedp2=gppfuncte_cupid(yy,xx,tables,ptr) 
      d2gtildedp2=d2gtildedp2/deltp**2 
!endif                                                                  
!         if (abs ((d2gtildedp2-randy5)/d2gtildedp2) .gt. diff5) then   
!           diff5 = abs ((d2gtildedp2-randy5)/d2gtildedp2)              
!           IF(myrank.eq.0)write(94,*) d2gtildedp2,randy5,diff5,x,y,p,t               
!         endif                                                         
!       endif                                                           
!***********************************************************************
! calculate d2(gtilde)/dtdp                                             
!if def,swesty                                                          
!      d2gtildedpdt = gii*dpsi0y*dpsi0x                                 
!    &              - gij*dpsi0y*dpsi0x2                                
!    &              - gji*dpsi0y1*dpsi0x                                
!    &              + gjj*dpsi0y1*dpsi0x2                               
!    &              + giit*dpsi0y*dpsi1x                                
!    &              + gijt*dpsi0y*dpsi1x2                               
!    &              - gjit*dpsi0y1*dpsi1x                               
!    &              - gjjt*dpsi0y1*dpsi1x2                              
!    &              + giitt*dpsi0y*dpsi2x                               
!    &              - gijtt*dpsi0y*dpsi2x2                              
!    &              - gjitt*dpsi0y1*dpsi2x                              
!    &              + gjjtt*dpsi0y1*dpsi2x2                             
!    &              + giip*dpsi1y*dpsi0x                                
!    &              - gijp*dpsi1y*dpsi0x2                               
!    &              + gjip*dpsi1y1*dpsi0x                               
!    &              - gjjp*dpsi1y1*dpsi0x2                              
!    &              + giipp*dpsi2y*dpsi0x                               
!    &              - gijpp*dpsi2y*dpsi0x2                              
!    &              - gjipp*dpsi2y1*dpsi0x                              
!    &              + gjjpp*dpsi2y1*dpsi0x2                             
!    &              + giitp*dpsi1y*dpsi1x                               
!    &              + gijtp*dpsi1y*dpsi1x2                              
!    &              + gjitp*dpsi1y1*dpsi1x                              
!    &              + gjjtp*dpsi1y1*dpsi1x2                             
!    &              + giittp*dpsi1y*dpsi2x                              
!    &              - gijttp*dpsi1y*dpsi2x2                             
!    &              + gjittp*dpsi1y1*dpsi2x                             
!    &              - gjjttp*dpsi1y1*dpsi2x2                            
!    &              + giitpp*dpsi2y*dpsi1x                              
!    &              + gijtpp*dpsi2y*dpsi1x2                             
!    &              - gjitpp*dpsi2y1*dpsi1x                             
!    &              - gjjtpp*dpsi2y1*dpsi1x2                            
!    &              + giittpp*dpsi2y*dpsi2x                             
!    &              - gijttpp*dpsi2y*dpsi2x2                            
!    &              - gjittpp*dpsi2y1*dpsi2x                            
!    &              + gjjttpp*dpsi2y1*dpsi2x2                           
!      d2gtildedpdt = d2gtildedpdt/(deltp*deltt)                        
!endif                                                                  
!       if (arg .eq. 'vapor') then                                      
!         randy6 = gptfunct(y,x,c)                                      
!         randy6 = randy6/(deltp*deltt)                                 
!if -def, swesty                                                        
!gam       d2gtildedpdt = gptfunct(yy,xx,tables,ptr)                    
      d2gtildedpdt=gptfuncte_cupid(yy,xx,tables,ptr) 
      d2gtildedpdt=d2gtildedpdt/(deltp*deltt) 
!endif                                                                  
!         if (abs ((d2gtildedpdt-randy6)/d2gtildedpdt) .gt. diff6) then 
!           diff6 = abs ((d2gtildedpdt-randy6)/d2gtildedpdt)            
!           IF(myrank.eq.0)write(95,*) d2gtildedpdt,randy6,diff6,x,y,p,t              
!         endif                                                         
!       endif                                                           
!***********************************************************************
! calculate d3(gtilde)/dt3                                              
!      d3gtildedt3 = gii*psi0y*d3psi0x                                  
!    &             - gij*psi0y*d3psi0x2                                 
!    &             + gji*psi0y1*d3psi0x                                 
!    &             - gjj*psi0y1*d3psi0x2                                
!    &             + giit*psi0y*d3psi1x                                 
!    &             + gijt*psi0y*d3psi1x2                                
!    &             + gjit*psi0y1*d3psi1x                                
!    &             + gjjt*psi0y1*d3psi1x2                               
!    &             + giitt*psi0y*d3psi2x                                
!    &             - gijtt*psi0y*d3psi2x2                               
!    &             + gjitt*psi0y1*d3psi2x                               
!    &             - gjjtt*psi0y1*d3psi2x2                              
!    &             + giip*psi1y*d3psi0x                                 
!    &             - gijp*psi1y*d3psi0x2                                
!    &             - gjip*psi1y1*d3psi0x                                
!    &             + gjjp*psi1y1*d3psi0x2                               
!    &             + giipp*psi2y*d3psi0x                                
!    &             - gijpp*psi2y*d3psi0x2                               
!    &             + gjipp*psi2y1*d3psi0x                               
!    &             - gjjpp*psi2y1*d3psi0x2                              
!    &             + giitp*psi1y*d3psi1x                                
!    &             + gijtp*psi1y*d3psi1x2                               
!    &             - gjitp*psi1y1*d3psi1x                               
!    &             - gjjtp*psi1y1*d3psi1x2                              
!    &             + giittp*psi1y*d3psi2x                               
!    &             - gijttp*psi1y*d3psi2x2                              
!    &             - gjittp*psi1y1*d3psi2x                              
!    &             + gjjttp*psi1y1*d3psi2x2                             
!    &             + giitpp*psi2y*d3psi1x                               
!    &             + gijtpp*psi2y*d3psi1x2                              
!    &             + gjitpp*psi2y1*d3psi1x                              
!    &             + gjjtpp*psi2y1*d3psi1x2                             
!    &             + giittpp*psi2y*d3psi2x                              
!    &             - gijttpp*psi2y*d3psi2x2                             
!    &             + gjittpp*psi2y1*d3psi2x                             
!    &             - gjjttpp*psi2y1*d3psi2x2                            
!      d3gtildedt3 = d3gtildedt3/deltt**3                               
!***********************************************************************
! calculate d3(gtilde)/dt2dp                                            
!      d3gtildedpdt2 = gii*dpsi0y*d2psi0x                               
!    &               + gij*dpsi0y*d2psi0x2                              
!    &               - gji*dpsi0y1*d2psi0x                              
!    &               - gjj*dpsi0y1*d2psi0x2                             
!    &               + giit*dpsi0y*d2psi1x                              
!    &               - gijt*dpsi0y*d2psi1x2                             
!    &               - gjit*dpsi0y1*d2psi1x                             
!    &               + gjjt*dpsi0y1*d2psi1x2                            
!    &               + giitt*dpsi0y*d2psi2x                             
!    &               + gijtt*dpsi0y*d2psi2x2                            
!    &               - gjitt*dpsi0y1*d2psi2x                            
!    &               - gjjtt*dpsi0y1*d2psi2x2                           
!    &               + giip*dpsi1y*d2psi0x                              
!    &               + gijp*dpsi1y*d2psi0x2                             
!    &               + gjip*dpsi1y1*d2psi0x                             
!    &               + gjjp*dpsi1y1*d2psi0x2                            
!    &               + giipp*dpsi2y*d2psi0x                             
!    &               + gijpp*dpsi2y*d2psi0x2                            
!    &               - gjipp*dpsi2y1*d2psi0x                            
!    &               - gjjpp*dpsi2y1*d2psi0x2                           
!    &               + giitp*dpsi1y*d2psi1x                             
!    &               - gijtp*dpsi1y*d2psi1x2                            
!    &               + gjitp*dpsi1y1*d2psi1x                            
!    &               - gjjtp*dpsi1y1*d2psi1x2                           
!    &               + giittp*dpsi1y*d2psi2x                            
!    &               + gijttp*dpsi1y*d2psi2x2                           
!    &               + gjittp*dpsi1y1*d2psi2x                           
!    &               + gjjttp*dpsi1y1*d2psi2x2                          
!    &               + giitpp*dpsi2y*d2psi1x                            
!    &               - gijtpp*dpsi2y*d2psi1x2                           
!    &               - gjitpp*dpsi2y1*d2psi1x                           
!    &               + gjjtpp*dpsi2y1*d2psi1x2                          
!    &               + giittpp*dpsi2y*d2psi2x                           
!    &               + gijttpp*dpsi2y*d2psi2x2                          
!    &               - gjittpp*dpsi2y1*d2psi2x                          
!    &               - gjjttpp*dpsi2y1*d2psi2x2                         
!      d3gtildedpdt2 = d3gtildedpdt2/(deltt**2*deltp)                   
!***********************************************************************
! calculate d3(gtilde)/dp2dt                                            
!      d3gtildedp2dt = gii*d2psi0y*dpsi0x                               
!    &               - gij*d2psi0y*dpsi0x2                              
!    &               + gji*d2psi0y1*dpsi0x                              
!    &               - gjj*d2psi0y1*dpsi0x2                             
!    &               + giit*d2psi0y*dpsi1x                              
!    &               + gijt*d2psi0y*dpsi1x2                             
!    &               + gjit*d2psi0y1*dpsi1x                             
!    &               + gjjt*d2psi0y1*dpsi1x2                            
!    &               + giitt*d2psi0y*dpsi2x                             
!    &               - gijtt*d2psi0y*dpsi2x2                            
!    &               + gjitt*d2psi0y1*dpsi2x                            
!    &               - gjjtt*d2psi0y1*dpsi2x2                           
!    &               + giip*d2psi1y*dpsi0x                              
!    &               - gijp*d2psi1y*dpsi0x2                             
!    &               - gjip*d2psi1y1*dpsi0x                             
!    &               + gjjp*d2psi1y1*dpsi0x2                            
!    &               + giipp*d2psi2y*dpsi0x                             
!    &               - gijpp*d2psi2y*dpsi0x2                            
!    &               + gjipp*d2psi2y1*dpsi0x                            
!    &               - gjjpp*d2psi2y1*dpsi0x2                           
!    &               + giitp*d2psi1y*dpsi1x                             
!    &               + gijtp*d2psi1y*dpsi1x2                            
!    &               - gjitp*d2psi1y1*dpsi1x                            
!    &               - gjjtp*d2psi1y1*dpsi1x2                           
!    &               + giittp*d2psi1y*dpsi2x                            
!    &               - gijttp*d2psi1y*dpsi2x2                           
!    &               - gjittp*d2psi1y1*dpsi2x                           
!    &               + gjjttp*d2psi1y1*dpsi2x2                          
!    &               + giitpp*d2psi2y*dpsi1x                            
!    &               + gijtpp*d2psi2y*dpsi1x2                           
!    &               + gjitpp*d2psi2y1*dpsi1x                           
!    &               + gjjtpp*d2psi2y1*dpsi1x2                          
!    &               + giittpp*d2psi2y*dpsi2x                           
!    &               - gijttpp*d2psi2y*dpsi2x2                          
!    &               + gjittpp*d2psi2y1*dpsi2x                          
!    &               - gjjttpp*d2psi2y1*dpsi2x2                         
!      d3gtildedp2dt = d3gtildedp2dt/(deltt*deltp**2)                   
!***********************************************************************
! calculate d4(gtilde)/dp2dt2                                           
!      d4gtildedp2dt2 = gii*d2psi0y*d2psi0x                             
!    &                + gij*d2psi0y*d2psi0x2                            
!    &                + gji*d2psi0y1*d2psi0x                            
!    &                + gjj*d2psi0y1*d2psi0x2                           
!    &                + giit*d2psi0y*d2psi1x                            
!    &                - gijt*d2psi0y*d2psi1x2                           
!    &                + gjit*d2psi0y1*d2psi1x                           
!    &                - gjjt*d2psi0y1*d2psi1x2                          
!    &                + giitt*d2psi0y*d2psi2x                           
!    &                + gijtt*d2psi0y*d2psi2x2                          
!    &                + gjitt*d2psi0y1*d2psi2x                          
!    &                + gjjtt*d2psi0y1*d2psi2x2                         
!    &                + giip*d2psi1y*d2psi0x                            
!    &                + gijp*d2psi1y*d2psi0x2                           
!    &                - gjip*d2psi1y1*d2psi0x                           
!    &                - gjjp*d2psi1y1*d2psi0x2                          
!    &                + giipp*d2psi2y*d2psi0x                           
!    &                + gijpp*d2psi2y*d2psi0x2                          
!    &                + gjipp*d2psi2y1*d2psi0x                          
!    &                + gjjpp*d2psi2y1*d2psi0x2                         
!    &                + giitp*d2psi1y*d2psi1x                           
!    &                - gijtp*d2psi1y*d2psi1x2                          
!    &                - gjitp*d2psi1y1*d2psi1x                          
!    &                + gjjtp*d2psi1y1*d2psi1x2                         
!    &                + giittp*d2psi1y*d2psi2x                          
!    &                + gijttp*d2psi1y*d2psi2x2                         
!    &                - gjittp*d2psi1y1*d2psi2x                         
!    &                - gjjttp*d2psi1y1*d2psi2x2                        
!    &                + giitpp*d2psi2y*d2psi1x                          
!    &                - gijtpp*d2psi2y*d2psi1x2                         
!    &                + gjitpp*d2psi2y1*d2psi1x                         
!    &                - gjjtpp*d2psi2y1*d2psi1x2                        
!    &                + giittpp*d2psi2y*d2psi2x                         
!    &                + gijttpp*d2psi2y*d2psi2x2                        
!    &                + gjittpp*d2psi2y1*d2psi2x                        
!    &                + gjjttpp*d2psi2y1*d2psi2x2                       
!      d4gtildedp2dt2 = d4gtildedp2dt2/((deltt**2)*(deltp**2))          
!***********************************************************************
!if def,logp                                                            
      IF(arg.eq.'vapor')then 
         p=exp(p) 
         d2gtildedp2=(d2gtildedp2-dgtildedp)/p**2 
         dgtildedp=dgtildedp/p 
         d2gtildedpdt=d2gtildedpdt/p 
      ENDIF 
!endif                                                                  
!  fill in g array                                                      
      g(1)=gtilde 
      g(2)=dgtildedp 
      g(3)=d2gtildedp2 
      g(4)=dgtildedt 
      g(5)=d2gtildedpdt 
!      g(6)  = d3gtildedp2dt                                            
      g(7)=d2gtildedt2 
!      g(8)  = d3gtildedpdt2                                            
!      g(9)  = d4gtildedp2dt2                                           
!      g(10) = d3gtildedt3                                              
!***********************************************************************
!      call timstop ('herm2dleft')                                      
      RETURN 
      END SUBROUTINE herm2dleft_cupid                     
