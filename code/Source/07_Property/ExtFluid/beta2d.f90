      SUBROUTINE beta2d_cupid(tables,ptr,a,arg,g,t,p,err) 
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
!deck beta2d                                                            
!                                                                       
!  $Id: beta2d.ff,v 1.1 2001/04/19 18:50:57 dbarber Exp dbarber $       
!                                                                       
!  Evaluates the 2D Hermite interpolating polynomial and several        
!  derivatives                                                          
!                                                                       
!  Cognizant engineer:  rex 3/28/2001                                   
!                                                                       
      IMPLICIT none 
!                                                                       
!***********************************************************************
!                                                                       
!  Declarations                                                         
!                                                                       
      INTEGER ptr 
      REAL(8) xx(6),yy(6),tables(*) 
      REAL(8) x,y,gptfunct_cupid 
      REAL(8) a,g,p,t,pii,pij,pji,pjj,tii,tij,tji,tjj 
      REAL(8) deltt,deltp,d2gtildedpdt 
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
!      call timstart ('beta2d')                                         
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
!                                                                       
!***********************************************************************
! calculate d2(gtilde)/dtdp                                             
      d2gtildedpdt=gptfunct_cupid(yy,xx,tables,ptr) 
      d2gtildedpdt=d2gtildedpdt/(deltp*deltt) 
!***********************************************************************
!if def,logp                                                            
      IF(arg.eq.'vapor')then 
         d2gtildedpdt=d2gtildedpdt/p 
      ENDIF 
!endif                                                                  
!  fill in g array                                                      
      g(5)=d2gtildedpdt 
!                                                                       
!***********************************************************************
!      call timstop ('beta2d')                                          
      RETURN 
      END SUBROUTINE beta2d_cupid                         
