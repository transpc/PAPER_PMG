      SUBROUTINE movehiright_cupid(itleftphi,itleftplo) 
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
!deck movehiright                                                       
!                                                                       
!  $Id: movehiright.ff,v 1.1 2001/02/01 23:15:49 r5qa Exp dbarber $     
!                                                                       
!  Cognizant engineer:  rwt                                             
!                                                                       
      IMPLICIT none 
!                                                                       
      INCLUDE 'stcom.h' 
      INCLUDE 'newstcom.h' 
      INCLUDE 'gibbpnt.h' 
!                                                                       
      INTEGER itleftplo,itleftphi 
!                                                                       
      itleftphi=itleftplo 
!                                                                       
      RETURN 
      END SUBROUTINE movehiright_cupid                    
