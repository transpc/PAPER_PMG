      SUBROUTINE mixprops_cupid(s,getprops) 
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
!deck mixprops                                                          
!                                                                       
!  $Id: mixprops.ff,v 1.1 2001/02/01 23:15:49 r5qa Exp dbarber $        
!                                                                       
!  Returns the two-phase mixture properties                             
!  Does not compute s(ubar) if getprops(ubar) is false                  
!                                                                       
!  Cognizant engineer:  rwt                                             
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL(8) s(*) 
      LOGICAL getprops(*) 
      REAL(8) sq,sq1 
      INCLUDE 'sparms.h' 
!                                                                       
!  liquid side                                                          
!                                                                       
      IF(s(qual).eq.0.0d0)then 
         s(vbar)=s(vsubf) 
         IF(getprops(hbar))s(hbar)=s(hsubf) 
         s(beta)=s(betaf) 
         s(kapa)=s(kapaf) 
         s(cp)=s(cpf) 
         s(entpy)=s(entpyf) 
         IF(getprops(ubar))s(ubar)=s(usubf) 
!                                                                       
!  vapor side                                                           
!                                                                       
      ELSEIF(s(qual).eq.1.0d0)then 
         s(vbar)=s(vsubg) 
         IF(getprops(hbar))s(hbar)=s(hsubf) 
         s(beta)=s(betag) 
         s(kapa)=s(kapag) 
         s(cp)=s(cpg) 
         s(entpy)=s(entpyg) 
         IF(getprops(ubar))s(ubar)=s(usubg) 
!                                                                       
!  two-phase mixture                                                    
!                                                                       
      ELSE 
         sq=s(qual) 
         sq1=1.0d0-sq 
         s(vbar)=sq*s(vsubg)+sq1*s(vsubf) 
         IF(getprops(hbar))then 
            s(hbar)=sq*s(hsubg)+sq1*s(hsubf) 
         ENDIF
         s(beta)=sq*s(betag)+sq1*s(betaf) 
         s(kapa)=sq*s(kapag)+sq1*s(kapaf) 
         s(cp)=sq*s(cpg)+sq1*s(cpf) 
         s(entpy)=sq*s(entpyg)+sq1*s(entpyf) 
         IF(getprops(ubar))then 
            s(ubar)=sq*s(usubg)+sq1*s(usubf) 
         ENDIF 
      ENDIF 
!                                                                       
      RETURN 
      END SUBROUTINE mixprops_cupid                       
