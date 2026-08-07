      FUNCTION gttfuncte_cupid(deltapres,deltatemp,tables,ptr) 
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
!deck gttfuncte                                                         
!                                                                       
!  Evaluate gtt (d2g/dt2) as a function of dimensionless                
!    deltapres and deltatemp given the                                  
!    coefficients of the bi-quintic polynomial.                         
!  coefficients of the bi-quintic polynomial.                           
!  You can pass in a 1-D vector coeff as long as the first index is for 
!  The row index (first index) for c(ip,jt) corresponds to the deltapres
!  The column index (second index) for c(ip,jt) corresponds to the delta
!                                                                       
      IMPLICIT none 
      REAL(8) gttfuncte_cupid,tables(*) 
      REAL(8) deltapres(*),deltatemp(*) 
!      real coeff(6,6)                                                  
!      integer ip, jt                                                   
      INTEGER ptr 
!                                                                       
!      dg = 0.0                                                         
!      gttfunct = 0.0                                                   
!      do 20 ip = 3, 6                                                  
!        do 10 jt = 1, 6                                                
!          gttfunct = gttfunct                                          
!    &              + coeff(ip,jt)                                      
!    &                *(ip-1)*(ip-2)*(deltatemp(ip-2))                  
!    &                *(deltapres(jt))                                  
!          dg = coeff(3,jt)*2.0*1.0*deltatemp(1)                        
!    &        + coeff(4,jt)*3.0*2.0*deltatemp(2)                        
!    &        + coeff(5,jt)*4.0*3.0*deltatemp(3)                        
!    &        + coeff(6,jt)*5.0*4.0*deltatemp(4)                        
!          gttfunct = gttfunct + dg*deltapres(jt)                       
! 10     continue                                                       
! 20   continue                                                         
!                                                                       
      gttfuncte_cupid=deltapres(1)*(tables(ptr+3) * 2.0d0 * 1.0d0*deltatemp(1)&
      +tables(ptr+4) * 3.0d0 * 2.0d0*deltatemp(2))+deltapres(2)*(tables(&
      ptr+9) * 2.0d0 * 1.0d0*deltatemp(1)+tables(ptr+10) * 3.0d0 *      &
      2.0d0*deltatemp(2))+deltapres(3)*(tables(ptr+15) * 2.0d0 * 1.0d0* &
      deltatemp(1)+tables(ptr+16) * 3.0d0 * 2.0d0*deltatemp(2))         
      gttfuncte_cupid=gttfuncte_cupid+deltapres(4)*(tables(ptr+21) * 2.0d0 * 1.0d0* &
      deltatemp(1)+tables(ptr+22) * 3.0d0 * 2.0d0*deltatemp(2))+        &
      deltapres(5)*(tables(ptr+27) * 2.0d0 * 1.0d0*deltatemp(1)+tables( &
      ptr+28) * 3.0d0 * 2.0d0*deltatemp(2))+deltapres(6)*(tables(ptr+33)&
      * 2.0d0 * 1.0d0*deltatemp(1)+tables(ptr+34) * 3.0d0 * 2.0d0*      &
      deltatemp(2))                                                     
      RETURN 
      END FUNCTION gttfuncte_cupid                        
