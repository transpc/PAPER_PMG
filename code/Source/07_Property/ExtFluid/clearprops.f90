      SUBROUTINE clearprops_cupid(getprops) 
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
!deck clearprops                                                        
!                                                                       
!  $Id: clearprops.ff,v 1.1 2001/02/01 23:15:49 r5qa Exp dbarber $      
!                                                                       
!  Cognizant engineer:  rwt                                             
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER i 
      LOGICAL getprops(26) 
!                                                                       
      DO 10 i=1,26 
         getprops(i)=.false. 
   10 END DO 
!                                                                       
      RETURN 
      END SUBROUTINE clearprops_cupid                     
