      SUBROUTINE setsngprops_cupid(getprops) 
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
!deck setsngprops                                                       
!                                                                       
!  $Id: setsngprops.ff,v 1.1 2001/02/01 23:15:49 r5qa Exp dbarber $     
!                                                                       
!  Cognizant engineer:  rwt                                             
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER i 
      LOGICAL getprops(26) 
!                                                                       
      CALL clearprops_cupid(getprops) 
!                                                                       
      DO 10 i=1,9 
         getprops(i)=.true. 
   10 END DO 
!                                                                       
!  Set flag for entropy                                                 
!                                                                       
      getprops(24)=.true. 
!                                                                       
      RETURN 
      END SUBROUTINE setsngprops_cupid                    
