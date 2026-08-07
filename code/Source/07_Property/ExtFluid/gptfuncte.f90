      FUNCTION gptfuncte_cupid(deltapres,deltatemp,tables,ptr) 
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
!deck gptfuncte                                                         
!                                                                       
!  Evaluate gpt (d2g/dpdt) as a function of dimensionless               
!    deltapres and deltatemp given the                                  
!    coefficients of the bi-quintic polynomial.                         
!  coefficients of the bi-quintic polynomial.                           
!  You can pass in a 1-D vector coeff as long as the first index is for 
!  The row index (first index) for c(ip,jt) corresponds to the deltapres
!  The column index (second index) for c(ip,jt) corresponds to the delta
!                                                                       
      IMPLICIT none 
      REAL(8) gptfuncte_cupid,tables(*) 
      REAL(8) deltapres(*),deltatemp(*) 
!      real coeff(6,6)                                                  
!      integer ip, jt                                                   
      INTEGER ptr 
!                                                                       
!      gptfunct = 0.0                                                   
!      do 20 ip = 2, 6                                                  
!        do 10 jt = 2, 6                                                
!          gptfunct = gptfunct                                          
!    &              + coeff(ip,jt)                                      
!    &                *(ip-1)*(deltatemp(ip-1))                         
!    &                *(jt-1)*(deltapres(jt-1))                         
! 10     continue                                                       
! 20   continue                                                         
      gptfuncte_cupid=1.0d0*deltapres(1)*(tables(ptr+8) * 1.0d0*deltatemp(1)+ &
      tables(ptr+9) * 2.0d0*deltatemp(2)+tables(ptr+10) * 3.0d0*        &
      deltatemp(3))+2.0d0*deltapres(2)*(tables(ptr+14) * 1.0d0*         &
      deltatemp(1)+tables(ptr+15) * 2.0d0*deltatemp(2)+tables(ptr+16) * &
      3.0d0*deltatemp(3))+3.0d0*deltapres(3)*(tables(ptr+20) * 1.0d0*   &
      deltatemp(1)+tables(ptr+21) * 2.0d0*deltatemp(2)+tables(ptr+22) * &
      3.0d0*deltatemp(3))+4.0d0*deltapres(4)*(tables(ptr+26) * 1.0d0*   &
      deltatemp(1)+tables(ptr+27) * 2.0d0*deltatemp(2)+tables(ptr+28) * &
      3.0d0*deltatemp(3))+5.0d0*deltapres(5)*(tables(ptr+32) * 1.0d0*   &
      deltatemp(1)+tables(ptr+33) * 2.0d0*deltatemp(2)+tables(ptr+34) * &
      3.0d0*deltatemp(3))                                               
!                                                                       
      RETURN 
      END FUNCTION gptfuncte_cupid                        
