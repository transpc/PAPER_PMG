      SUBROUTINE nth2x3_cupid(a,s,it,err) 
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
!deck nth2x3                                                            
!                                                                       
!  $Id: nth2x3.ff,v 1.2 2001/04/25 14:55:49 dbarber Exp dbarber $       
!                                                                       
!  compute water thermodynamic properties as a function of temperature  
!  and pressure                                                         
!                                                                       
      USE STM_TBL_cupid  , ONLY: pcrit
!
      IMPLICIT none 
!                                                                       
      REAL(8) a(*),s(*) 
      INTEGER it 
      LOGICAL err,getprops(26) 
!                                                                       
      INCLUDE 'stcom.h' 
      INCLUDE 'sparms.h' 
      INCLUDE 'newstcom.h' 
!                                                                       
!   check for valid input                                               
      err=s(temp).lt.tmin.or.s(temp).gt.tmax.or.s(pres).le.pmin.or.s(   &
      pres).gt.pmax                                                     
      IF(err)return 
!                                                                       
      IF(s(pres).ge.pcrit)then 
         it=4 
      ELSEIF(s(pres).lt.ptrip.or.s(temp).gt.tcrit)then 
         it=3 
         s(qual)=1.0d0 
      ELSE 
!                                                                       
!  get saturation temperature at this pressure                          
!                                                                       
         CALL gettsat_cupid(a,s(pres),s,err) 
         IF(err)return 
!                                                                       
!  determine if we are in liquid or vapor region                        
!                                                                       
         IF(s(temp).le.s(tsat))then 
            it=1 
            s(qual)=0.0d0 
         ELSE 
            it=3 
            s(qual)=1.0d0 
         ENDIF 
      ENDIF 
!                                                                       
!  now that the region is determined, get the properties                
!  first get the saturation properties at this pressure                 
!  if not in supercritical pressure region (it=4)                       
!  second get the single phase properties for the region we are in      
!                                                                       
      IF(it.ne.4)then 
         CALL setallprops_cupid(getprops) 
         CALL getsatatp_cupid(a,getprops,s(pres),s,err) 
      ENDIF 
!                                                                       
      CALL setsngprops_cupid(getprops) 
      IF(it.eq.1)then 
         CALL getsnglatpt_cupid(a,'liquid',getprops,s(temp),s(pres),s,err) 
      ELSEIF(it.eq.3)then 
         CALL getsnglatpt_cupid(a,'vapor',getprops,s(temp),s(pres),s,err) 
      ELSEIF(it.eq.4)then 
         CALL getsnglatpt_cupid(a,'supercritical',getprops,s(temp),s(pres),s, &
         err)                                                           
      ENDIF 
!                                                                       
      RETURN 
      END SUBROUTINE nth2x3_cupid                         
