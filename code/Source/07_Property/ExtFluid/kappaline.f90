      SUBROUTINE kappaline_cupid(tables,ptr,a,arg,g,t,p,err) 
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
!deck kappaline                                                         
!                                                                       
!  $Id: kappaline.ff,v 1.1 2001/04/02 17:41:41 dbarber Exp dbarber $    
!                                                                       
!  Evaluates the g(3) in a linear mode. kappa = -g(3)/g(2)              
!                                                                       
!  Cognizant engineer:  rex 3/20/2001                                   
!                                                                       
      IMPLICIT none 
!                                                                       
!***********************************************************************
!                                                                       
!  Declarations                                                         
!                                                                       
      INTEGER ptr 
      REAL(8) xx(6),yy(6),tables(*) 
      REAL(8) x,y,gppfunctline_cupid,gpfunctline_cupid 
!                                                                       
      REAL(8) a,g,p,t,pii,pij,pji,pjj,tii,tij,tji,tjj 
      REAL(8) deltt,deltp,d2gtildedp2,dgtildedp 
!                                                                       
      DIMENSION a(11,4),g(10) 
      LOGICAL err 
      CHARACTER*(*)arg 
!  Definitions                                                          
!                                                                       
!      call timstart ('kappaline')                                      
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
!***********************************************************************
! calculate d(gtilde)/dp                                                
!                                                                       
      dgtildedp=gpfunctline_cupid(yy,xx,tables,ptr) 
      dgtildedp=dgtildedp/deltp 
!***********************************************************************
! calculate d2(gtilde)/dp2                                              
!                                                                       
      d2gtildedp2=gppfunctline_cupid(yy,xx,tables,ptr) 
      d2gtildedp2=d2gtildedp2/deltp**2 
!***********************************************************************
!if def,logp                                                            
      IF(arg.eq.'vapor')then 
         d2gtildedp2=(d2gtildedp2-dgtildedp)/p**2 
      ENDIF 
!endif                                                                  
!                                                                       
      g(3)=d2gtildedp2 
!                                                                       
!***********************************************************************
!      call timstop ('kappaline')                                       
      RETURN 
      END SUBROUTINE kappaline_cupid                      
