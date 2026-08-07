      FUNCTION gfuncte_cupid(deltapres,deltatemp,tables,ptr) 
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
!deck gfuncte                                                           
!                                                                       
!  Evaluate g as a function of dimensionless deltapres and deltatemp giv
!  coefficients of the bi-quintic polynomial.                           
!  You can pass in a 1-D vector coeff as long as the first index is for 
!  The row index (first index) for c(ip,jt) corresponds to the deltatemp
!  The column index (second index) for c(ip,jt) corresponds to the delta
!                                                                       
      IMPLICIT none 
      INTEGER ptr 
      REAL(8) gfuncte_cupid,tables(*) 
!      real dg                                                          
      REAL(8) deltapres(*),deltatemp(*) 
!      real coeff(6,6)                                                  
!      integer ip, jt                                                   
!                                                                       
!      gfunct = 0.0                                                     
!      do 20 ip = 1, 6                                                  
!      dg = 0.0                                                         
!        do 10 jt = 1, 6                                                
!          gfunct = gfunct                                              
!    &            + coeff(ip,jt)                                        
!    &              *(deltatemp(ip))                                    
!    &              *(deltapres(jt))                                    
!          dg = dg + coeff(ip,jt)*deltapres(jt)                         
! 10     continue                                                       
!        gfunct = gfunct + dg*deltatemp(ip)                             
! 20   continue                                                         
      gfuncte_cupid=deltapres(1)*(tables(ptr+1)*deltatemp(1)+tables(ptr+2)*   &
      deltatemp(2)+tables(ptr+3)*deltatemp(3)+tables(ptr+4)*deltatemp(4)&
      )+deltapres(2)*(tables(ptr+7)*deltatemp(1)+tables(ptr+8)*         &
      deltatemp(2)+tables(ptr+9)*deltatemp(3)+tables(ptr+10)*deltatemp( &
      4))+deltapres(3)*(tables(ptr+13)*deltatemp(1)+tables(ptr+14)*     &
      deltatemp(2)+tables(ptr+15)*deltatemp(3)+tables(ptr+16)*deltatemp(&
      4))                                                               
      gfuncte_cupid=gfuncte_cupid+deltapres(4)*(tables(ptr+19)*deltatemp(1)+tables( &
      ptr+20)*deltatemp(2)+tables(ptr+21)*deltatemp(3)+tables(ptr+22)*  &
      deltatemp(4))+deltapres(5)*(tables(ptr+25)*deltatemp(1)+tables(   &
      ptr+26)*deltatemp(2)+tables(ptr+27)*deltatemp(3)+tables(ptr+28)*  &
      deltatemp(4))+deltapres(6)*(tables(ptr+31)*deltatemp(1)+tables(   &
      ptr+32)*deltatemp(2)+tables(ptr+33)*deltatemp(3)+tables(ptr+34)*  &
      deltatemp(4))                                                     
!                                                                       
      RETURN 
      END FUNCTION gfuncte_cupid                          
