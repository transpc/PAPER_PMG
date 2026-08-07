      SUBROUTINE cubic_cupid(a,x,y,err) 
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
!deck cubic                                                             
!                                                                       
!  $Id: cubic.ff,v 1.1 2001/02/01 23:15:49 r5qa Exp dbarber $           
!                                                                       
!  Evaluates the 1D Hermite interpolating polynomial y(x)               
!                                                                       
!  Cognizant engineer:  rwt 5/25/1999                                   
!                                                                       
      IMPLICIT none 
!                                                                       
!***********************************************************************
!                                                                       
!  Local variables                                                      
!                                                                       
      REAL(8) a(8),delt,psi0,psi1,dydxi,dydxj,psi0x1,psi1x1,psi0x2,     &
      psi1x2,xi,xj,x1,x2,x,yi,yj,yit,yjt,y,z                            
      LOGICAL err 
!***********************************************************************
!                                                                       
!  Statement functions to define Hermite basis functions                
!                                                                       
      psi0(z)=1.0d0+(z**2*(-3.0d0+z*(2.0d0))) 
      psi1(z)=z*(1.0d0+z*(-2.0d0+z*(1.0d0))) 
!***********************************************************************
!                                                                       
!  Executable statements                                                
!                                                                       
!      call timstart ('cubic')                                          
      err=.false. 
      xi=a(1) 
      xj=a(5) 
      delt=(xj-xi) 
      x1=(x-xi)/delt 
      x2=1.0d0-x1 
!***********************************************************************
!                                                                       
!  Definitions                                                          
!                                                                       
      yi=a(2) 
      yj=a(6) 
      dydxi=a(3) 
      dydxj=a(7) 
!                                                                       
      yit=dydxi*delt 
      yjt=dydxj*delt 
!                                                                       
      psi0x1=psi0(x1) 
      psi1x1=psi1(x1) 
      psi0x2=psi0(x2) 
      psi1x2=psi1(x2) 
!***********************************************************************
!                                                                       
! calculate y at x                                                      
!                                                                       
      y=yi*psi0x1+yj*psi0x2+yit*psi1x1-yjt*psi1x2 
!***********************************************************************
!      call timstop ('cubic')                                           
      RETURN 
      END SUBROUTINE cubic_cupid                          
