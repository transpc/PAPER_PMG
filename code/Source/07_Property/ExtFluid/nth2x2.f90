      SUBROUTINE nth2x2_cupid(a,s,err) 
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
!deck nth2x2                                                            
!                                                                       
!  Compute water thermodynamic properties as a function                 
!  of pressure and quality                                              
!                                                                       
!  $Id: nth2x2.ff,v 1.1 2001/02/01 23:15:49 r5qa Exp dbarber $          
!                                                                       
      USE STM_TBL_cupid  , ONLY: pcrit
!
      IMPLICIT none 
!                                                                       
      REAL(8) a(*),s(*) 
      LOGICAL err,getprops(26) 
!                                                                       
      INCLUDE 'stcom.h' 
      INCLUDE 'newstcom.h' 
      INCLUDE 'sparms.h' 
!                                                                       
!  check for valid input                                                
!                                                                       
      err=s(pres).lt.ptrip.or.s(pres).gt.pcrit.or.s(qual).lt.0.0d0.or.s(&
      qual).gt.1.0d0                                                    
      IF(err)return 
!                                                                       
      CALL setallprops_cupid(getprops) 
      getprops(pres)=.false. 
      getprops(qual)=.false. 
!                                                                       
      CALL getsatatp_cupid(a,getprops,s(pres),s,err) 
!                                                                       
      s(temp)=s(tsat) 
!                                                                       
      RETURN 
      END SUBROUTINE nth2x2_cupid                         
