      FUNCTION gpfuncte_cupid(deltapres,deltatemp,tables,ptr) 
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
!deck gpfuncte                                                          
!                                                                       
!  Evaluate gp (dg/dp) as a function of dimensionless                   
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
      REAL(8) gpfuncte_cupid 
      REAL(8) deltapres(*),deltatemp(*) 
!      real coeff(6,6)                                                  
!      integer ip, jt                                                   
      INTEGER ptr 
!                                                                       
!      gpfunct = 0.0                                                    
!      do 20 ip = 1, 6                                                  
!        do 10 jt = 2, 6                                                
!          gpfunct = gpfunct                                            
!    &             + coeff(ip,jt)                                       
!    &               *(deltatemp(ip))                                   
!    &               *(jt-1)*(deltapres(jt-1))                          
!          dg = coeff(ip,2)*1.0*deltapres(1)                            
!    &        + coeff(ip,3)*2.0*deltapres(2)                            
!    &        + coeff(ip,4)*3.0*deltapres(3)                            
!    &        + coeff(ip,5)*4.0*deltapres(4)                            
!    &        + coeff(ip,6)*5.0*deltapres(5)                            
! 10     continue                                                       
!        gpfunct = gpfunct + dg*deltatemp(ip)                           
! 20   continue                                                         
!                                                                       
      gpfuncte_cupid=1.0d0*deltapres(1)*(tables(ptr+7)*deltatemp(1)+tables(   &
      ptr+8)*deltatemp(2)+tables(ptr+9)*deltatemp(3)+tables(ptr+10)*    &
      deltatemp(4))+2.0d0*deltapres(2)*(tables(ptr+13)*deltatemp(1)+    &
      tables(ptr+14)*deltatemp(2)+tables(ptr+15)*deltatemp(3)+tables(   &
      ptr+16)*deltatemp(4))+3.0d0*deltapres(3)*(tables(ptr+19)*         &
      deltatemp(1)+tables(ptr+20)*deltatemp(2)+tables(ptr+21)*deltatemp(&
      3)+tables(ptr+22)*deltatemp(4))                                   
      gpfuncte_cupid=gpfuncte_cupid+4.0d0*deltapres(4)*(tables(ptr+25)*deltatemp(1)+&
      tables(ptr+26)*deltatemp(2)+tables(ptr+27)*deltatemp(3)+tables(   &
      ptr+28)*deltatemp(4))+5.0d0*deltapres(5)*(tables(ptr+31)*         &
      deltatemp(1)+tables(ptr+32)*deltatemp(2)+tables(ptr+33)*deltatemp(&
      3)+tables(ptr+34)*deltatemp(4))                                   
      RETURN 
      END FUNCTION gpfuncte_cupid                         
