      FUNCTION gtfuncte_cupid(deltapres,deltatemp,tables,ptr) 
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
!deck gtfuncte                                                          
!                                                                       
!  Evaluate gt (dg/dt) as a function of dimensionless                   
!    deltapres and deltatemp given the                                  
!    coefficients of the bi-quintic polynomial.                         
!  coefficients of the bi-quintic polynomial.                           
!  You can pass in a 1-D vector coeff as long as the first index is for 
!  The row index (first index) for c(ip,jt) corresponds to the deltapres
!  The column index (second index) for c(ip,jt) corresponds to the delta
!                                                                       
      IMPLICIT none 
      REAL(8) tables(*) 
!      real dg                                                          
      REAL(8) gtfuncte_cupid 
      REAL(8) deltapres(*),deltatemp(*) 
!      real coeff(6,6)                                                  
!      integer ip, jt, ptr                                              
      INTEGER ptr 
!                                                                       
!      gtfunct = 0.0                                                    
!      do 20 ip = 2, 6                                                  
!        do 10 jt = 1, 6                                                
!          gtfunct = gtfunct                                            
!    &             + coeff(ip,jt)                                       
!    &               *(ip-1)*(deltatemp(ip-1))                          
!    &               *(deltapres(jt))                                   
!          dg = coeff(2,jt)*1.0*deltatemp(1)                            
!    &        + coeff(3,jt)*2.0*deltatemp(2)                            
!    &        + coeff(4,jt)*3.0*deltatemp(3)                            
!    &        + coeff(5,jt)*4.0*deltatemp(4)                            
!    &        + coeff(6,jt)*5.0*deltatemp(5)                            
!        gtfunct = gtfunct + dg*deltapres(jt)                           
! 10     continue                                                       
! 20   continue                                                         
!                                                                       
      gtfuncte_cupid=deltapres(1)*(tables(ptr+2) * 1.0d0*deltatemp(1)+tables( &
      ptr+3) * 2.0d0*deltatemp(2)+tables(ptr+4) * 3.0d0*deltatemp(3))+  &
      deltapres(2)*(tables(ptr+8) * 1.0d0*deltatemp(1)+tables(ptr+9) *  &
      2.0d0*deltatemp(2)+tables(ptr+10) * 3.0d0*deltatemp(3))+deltapres(&
      3)*(tables(ptr+14) * 1.0d0*deltatemp(1)+tables(ptr+15) * 2.0d0*   &
      deltatemp(2)+tables(ptr+16) * 3.0d0*deltatemp(3))                 
      gtfuncte_cupid=gtfuncte_cupid+deltapres(4)*(tables(ptr+20) * 1.0d0*deltatemp( &
      1)+tables(ptr+21) * 2.0d0*deltatemp(2)+tables(ptr+22) * 3.0d0*    &
      deltatemp(3))+deltapres(5)*(tables(ptr+26) * 1.0d0*deltatemp(1)+  &
      tables(ptr+27) * 2.0d0*deltatemp(2)+tables(ptr+28) * 3.0d0*       &
      deltatemp(3))+deltapres(6)*(tables(ptr+32) * 1.0d0*deltatemp(1)+  &
      tables(ptr+33) * 2.0d0*deltatemp(2)+tables(ptr+34) * 3.0d0*       &
      deltatemp(3))                                                     
      RETURN 
      END FUNCTION gtfuncte_cupid                         
