      SUBROUTINE satprops_cupid(xy,p1,p2,tables,x,y,err) 
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
!deck satprops                                                          
!                                                                       
!  $Id: satprops.ff,v 1.2 2001/04/25 14:55:49 dbarber Exp dbarber $     
!                                                                       
!  Evaluates the saturation properties using Hermite polonomials        
!                                                                       
!  Cognizant engineer:  rwt 5/25/1999                                   
!  modification of herm1d; gam 1/25/1000                                
!                                                                       
      IMPLICIT none 
!                                                                       
!***********************************************************************
!                                                                       
!  Local variables                                                      
!                                                                       
      INTEGER g 
      REAL(8) xy(*),tables(*),x,y(*) 
      REAL(8) x1 
      REAL(8) delt,delt2 
      INTEGER p1,p2,p1p1,p2p1,p1p2,p2p2 
      LOGICAL err 
      INTEGER nprop,i 
      DATA  nprop /12/ 
!                                                                       
      err=.false. 
      delt=xy(2)-xy(1) 
      delt2=delt**2 
      x1=(x-xy(1))/delt 
!                                                                       
      DO 10 i=1,nprop 
         p1p1=p1+1 
         p2p1=p2+1 
         p1p2=p1+2 
         p2p2=p2+2 
!                                                                       
         g=p1+2 
         y(i)=tables(g+1)+x1*(tables(g+2)+x1*(tables(g+3)+x1*(tables(g+ &
         4)+x1*(tables(g+5)+x1*tables(g+6)))))                          
!                                                                       
!  bump the pointers to the next guy in the table                       
!                                                                       
         p1=p1+9 
         p2=p2+9 
   10 END DO 
!                                                                       
      RETURN 
      END SUBROUTINE satprops_cupid                       
