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
!comdeck magic                                                          
!                                                                       
!  $Id: magic.hh,v 1.3 2001/04/02 17:41:41 dbarber Exp dbarber $        
!                                                                       
!  magic constant                                                       
!                                                                       
      INTEGER sizetables1 
      PARAMETER (sizetables1=1070000) 
!                                                                       
      COMMON/magic_cupid/magic,tables1 
      REAL(8) magic,tables1(sizetables1) 
!                                                                       
!***********************************************************************
!                                                                       
! Data dictionary for local variables                                   
!                                                                       
! Number of common variables =   1                                      
!                                                                       
! i=integer r=real l=logical c=character                                
!***********************************************************************
! Type Name   Definition                                                
!-----------------------------------------------------------------------
!  r  magic     = magic number, signals a point for which there are no  
!                 data values because of:                               
!                 thermal stability limit (Cv<0),                       
!                 mechanical stablity limit (v*kappa<0),                
!                 negative internal energy (u<0), or                    
!                 negative specific volume (v<0)                        
!********************************************************************** 
