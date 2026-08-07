      FUNCTION gpfunctline_cupid(deltapres,deltatemp,tables,ptr) 
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
!deck gpfunctline                                                       
!                                                                       
!  $Id: gpfunctline.ff,v 1.2 2001/04/19 18:50:57 dbarber Exp dbarber $  
!                                                                       
!  Evaluate gp (dg/dp) as a function of dimensionless                   
!    deltapres and deltatemp given the                                  
!    coefficients of the bi-quintic polynomial.                         
!  coefficients of the bi-quintic polynomial.                           
!  You can pass in a 1-D vector coeff as long as the first index is for 
!  The row index (first index) for c(ip,jt) corresponds to the deltapres
!  The column index (second index) for c(ip,jt) corresponds to the delta
!                                                                       
!  Cognizant engineer:  rex 3/20/2001                                   
!                                                                       
      IMPLICIT none 
      REAL(8) tables(*) 
!                                                                       
      REAL(8) gpfunctline_cupid 
      REAL(8) deltapres(*),deltatemp(*) 
!                                                                       
      INTEGER ptr 
!                                                                       
!                                                                       
      gpfunctline_cupid=1.0d0*deltapres(1)*(tables(ptr+7)*deltatemp(1)+tables(&
      ptr+8)*deltatemp(2)+tables(ptr+9)*deltatemp(3))+2.0d0*deltapres(2)&
      *(tables(ptr+13)*deltatemp(1)+tables(ptr+14)*deltatemp(2)+tables( &
      ptr+15)*deltatemp(3))                                             
      RETURN 
      END FUNCTION gpfunctline_cupid                      
