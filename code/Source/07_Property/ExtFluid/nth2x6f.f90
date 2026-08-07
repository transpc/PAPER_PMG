      SUBROUTINE nth2x6f_cupid(a,s,it,err,satflag) 
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
!deck nth2x6f                                                           
!                                                                       
!  $Id: nth2x6f.ff,v 1.2 2001/04/25 14:55:49 dbarber Exp dbarber $      
!                                                                       
!  Compute water thermodynamic properties as a function of pressure     
!  and internal energy.                                                 
!  Coomputes the saturation properties also.                            
!                                                                       
      USE STM_TBL_cupid  , ONLY: pcrit
!
      IMPLICIT none 
!                                                                       
      REAL(8) a(*),s(*),satt 
      INTEGER it
      LOGICAL err,getprops(26) 
      CHARACTER*(*)satflag 
!                                                                       
      INCLUDE 'stcom.h' 
      INCLUDE 'sparms.h' 
      INCLUDE 'newstcom.h' 
!     INCLUDE 'efiless.h' 
!                                                                       
!     INCLUDE 'efilesd.h' 
!                                                                       
!  check for valid input                                                
!                                                                       
      err=s(pres).le.pmin.or.s(pres).gt.pmax 
      IF(err)return 
!                                                                       
!  check if below the triple point pressure                             
!  if true, use the vapor table                                         
!                                                                       
      IF(s(pres).lt.ptrip)then 
         IF(it.eq.5)then 
            err=.true. 
            RETURN 
         ENDIF 
         it=3 
!                                                                       
!  check if above the critical point pressure                           
!  if true, us the supercritical pressure table                         
!                                                                       
      ELSEIF(s(pres).ge.pcrit)then 
         it=4 
!                                                                       
!  compute the saturation temperature                                   
!                                                                       
      ELSE 
         IF(satflag.eq.'6')then 
            CALL gettsat_cupid(a,s(pres),s,err) 
            satt=s(tsat) 
            IF(err)return 
         ELSE 
            satt=s(tsat) 
         ENDIF 
!                                                                       
!  compute the remaining saturation properties                          
!                                                                       
         CALL setallprops_cupid(getprops) 
         getprops(ubar)=.false. 
         getprops(pres)=.false. 
         IF(it.eq.5)then 
            s(qual)=0.0d0 
         ELSEIF(it.eq.6)then 
            s(qual)=1.0d0 
         ENDIF 
!                                                                       
         CALL getsatatp_cupid(a,getprops,s(pres),s,err) 
!                                                                       
         s(tsat)=satt 
         IF(err)return 
!                                                                       
!  check if subcooled liquid                                            
!                                                                       
         IF(it.eq.5)then 
            it=1 
         ELSEIF(it.eq.6)then 
            it=3 
         ENDIF 
      ENDIF 
!                                                                       
!  get the single phase properties from the appripriate table           
!                                                                       
      IF(it.eq.1)then 
         CALL setsngprops_cupid(getprops) 
         CALL getsnglatpu_cupid(a,'liquid',getprops,s(ubar),s(pres),s,err) 
      ELSEIF(it.eq.3)then 
         CALL setsngprops_cupid(getprops) 
         CALL getsnglatpu_cupid(a,'vapor',getprops,s(ubar),s(pres),s,err) 
      ELSEIF(it.eq.4)then 
         CALL setsngprops_cupid(getprops) 
         CALL getsnglatpu_cupid(a,'supercritical',getprops,s(ubar),s(pres),s, &
         err)                                                           
      ENDIF 
!                                                                       
      RETURN 
      END SUBROUTINE nth2x6f_cupid                        
