      SUBROUTINE nth2x1_cupid(a,s,err) 
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
!deck nth2x1                                                            
!                                                                       
!  $Id: nth2x1.ff,v 1.1 2001/02/01 23:15:49 r5qa Exp dbarber $          
!                                                                       
!  compute water thermodynamic properties as a function                 
!  of temperature and quality                                           
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
!  make sure temperature and quality are within range                   
!                                                                       
      err=s(temp).lt.ttrip.or.s(temp).gt.tcrit.or.s(qual).lt.0.0d0.or.s(&
      qual).gt.1.0d0                                                    
      IF(err)return 
!                                                                       
!  get all saturation properties at this temperature and quality        
!                                                                       
      CALL setallprops_cupid(getprops) 
      getprops(temp)=.false. 
      getprops(qual)=.false. 
!                                                                       
      CALL getsatatt_cupid(a,getprops,s(temp),s,err) 
!                                                                       
      s(pres)=s(psat) 
!                                                                       
      RETURN 
      END SUBROUTINE nth2x1_cupid                         
