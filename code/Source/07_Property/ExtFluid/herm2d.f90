      SUBROUTINE herm2d_cupid(tables,ptr,a,arg,g,t,p,err) 
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
!deck herm2d                                                            
!                                                                       
!  $Id: herm2d.ff,v 1.2 2001/04/02 17:41:41 dbarber Exp dbarber $       
!                                                                       
!  Evaluates the 2D Hermite interpolating polynomial and several        
!  derivatives                                                          
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
      REAL(8) x,y,gfunct_cupid,gtfunct_cupid,gttfunct_cupid,gpfunct_cupid,gptfunct_cupid,gppfunct_cupid 
      REAL(8) a,g,p,t,pii,pij,pji,pjj,tii,tij,tji,tjj 
      REAL(8) deltt,deltp,gtilde,dgtildedt,dgtildedp,d2gtildedt2,       &
      d2gtildedp2,d2gtildedpdt                                          
      DIMENSION a(11,4),g(10) 
      LOGICAL err 
      CHARACTER*(*)arg 
!                                                                       
!  Statement functions to define Hermite polynomials                    
!                                                                       
!***********************************************************************
!                                                                       
!  Definitions                                                          
!                                                                       
!      call timstart ('herm2d')                                         
      err=.false. 
!                                                                       
      tii=a(1,1) 
      tij=a(1,2) 
      tji=a(1,3) 
      tjj=a(1,4) 
!                                                                       
      pii=a(2,1) 
      pij=a(2,2) 
      pji=a(2,3) 
      pjj=a(2,4) 
!                                                                       
      deltt=(tij-tii) 
      deltp=(pji-pii) 
      x=(t-tii)/deltt 
      y=(p-pii)/deltp 
!if def,logp                                                            
      IF(arg.eq.'vapor')then 
         deltp=(log(pji)-log(pii)) 
         y=(log(p)-log(pii))/deltp 
      ENDIF 
!endif                                                                  
!                                                                       
!                                                                       
!***********************************************************************
! calculate gtilde                                                      
!                                                                       
!      if (arg .eq. 'vapor') then                                       
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
      gtilde=gfunct_cupid(yy,xx,tables,ptr) 
!***********************************************************************
! calculate d(gtilde)/dt                                                
!                                                                       
      dgtildedt=gtfunct_cupid(yy,xx,tables,ptr) 
      dgtildedt=dgtildedt/deltt 
!***********************************************************************
! calculate d(gtilde)/dp                                                
!                                                                       
      dgtildedp=gpfunct_cupid(yy,xx,tables,ptr) 
      dgtildedp=dgtildedp/deltp 
!***********************************************************************
! calculate d2(gtilde)/dt2                                              
      d2gtildedt2=gttfunct_cupid(yy,xx,tables,ptr) 
      d2gtildedt2=d2gtildedt2/deltt**2 
!***********************************************************************
! calculate d2(gtilde)/dp2                                              
!                                                                       
      d2gtildedp2=gppfunct_cupid(yy,xx,tables,ptr) 
      d2gtildedp2=d2gtildedp2/deltp**2 
!***********************************************************************
! calculate d2(gtilde)/dtdp                                             
      d2gtildedpdt=gptfunct_cupid(yy,xx,tables,ptr) 
      d2gtildedpdt=d2gtildedpdt/(deltp*deltt) 
!***********************************************************************
!if def,logp                                                            
      IF(arg.eq.'vapor')then 
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
!      call timstop ('herm2d')                                          
      RETURN 
      END SUBROUTINE herm2d_cupid                         
