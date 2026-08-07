      FUNCTION gppfuncte_cupid(deltapres,deltatemp,tables,ptr) 
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
!deck gppfuncte                                                         
!                                                                       
!  Evaluate gpp (d2g/dp2) as a function of dimensionless                
!    deltapres and deltatemp given the                                  
!    coefficients of the bi-quintic polynomial.                         
!  You can pass in a 1-D vector coeff as long as the first index is for 
!  The row index (first index) for c(ip,jt) corresponds to the deltapres
!  The column index (second index) for c(ip,jt) corresponds to the delta
!                                                                       
      IMPLICIT none 
      REAL(8) tables(*) 
!      real dg                                                          
      REAL(8) gppfuncte_cupid 
      REAL(8) deltapres(*),deltatemp(*) 
!      real coeff(6,6)                                                  
!      integer ip, jt                                                   
      INTEGER ptr 
!                                                                       
!      gppfunct = 0.0                                                   
!      do 20 ip = 1, 6                                                  
!        do 10 jt = 3, 6                                                
!          gppfunct = gppfunct                                          
!    &              + coeff(ip,jt)                                      
!    &                *(deltatemp(ip))                                  
!    &                *(jt-1)*(jt-2)*(deltapres(jt-2))                  
!          dg =                                                         
!    &        + coeff(ip,3)*2.0*1.0*deltapres(1)                        
!    &        + coeff(ip,4)*3.0*2.0*deltapres(2)                        
!    &        + coeff(ip,5)*4.0*3.0*deltapres(3)                        
!    &        + coeff(ip,6)*5.0*4.0*deltapres(4)                        
! 10     continue                                                       
!        gppfunct = gppfunct + dg*deltatemp(ip)                         
! 20   continue                                                         
      gppfuncte_cupid=2.0d0 * 1.0d0*deltapres(1)*(tables(ptr+13)*deltatemp(1)+&
      tables(ptr+14)*deltatemp(2)+tables(ptr+15)*deltatemp(3)+tables(   &
      ptr+16)*deltatemp(4))+3.0d0 * 2.0d0*deltapres(2)*(tables(ptr+19)* &
      deltatemp(1)+tables(ptr+20)*deltatemp(2)+tables(ptr+21)*deltatemp(&
      3)+tables(ptr+22)*deltatemp(4))                                   
      gppfuncte_cupid=gppfuncte_cupid+4.0d0 * 3.0d0*deltapres(3)*(tables(ptr+25)*   &
      deltatemp(1)+tables(ptr+26)*deltatemp(2)+tables(ptr+27)*deltatemp(&
      3)+tables(ptr+28)*deltatemp(4))+5.0d0 * 4.0d0*deltapres(4)*(      &
      tables(ptr+31)*deltatemp(1)+tables(ptr+32)*deltatemp(2)+tables(   &
      ptr+33)*deltatemp(3)+tables(ptr+34)*deltatemp(4))                 
!                                                                       
      RETURN 
      END FUNCTION gppfuncte_cupid                        
