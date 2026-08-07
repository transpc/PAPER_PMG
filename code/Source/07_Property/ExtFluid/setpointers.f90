      SUBROUTINE setpointers_cupid(arg,ofirstpres,olastpres,ptableprop,       &
      stableprop,ntableprop,ptablelimits,ntablelimits)                  
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
!deck setpointers                                                       
!                                                                       
!  $Id: setpointers.ff,v 1.1 2001/02/01 23:15:49 r5qa Exp dbarber $     
!                                                                       
!  Cognizant engineer:  rwt                                             
!                                                                       
      IMPLICIT none 
!                                                                       
      INCLUDE 'stcom.h' 
      INCLUDE 'newstcom.h' 
      INCLUDE 'gibbpnt.h' 
!                                                                       
      INTEGER ofirstpres,olastpres,ntablelimits,ntableprop,ptablelimits,&
      ptableprop,stableprop                                             
      CHARACTER*(*)arg 
      LOGICAL err 
!                                                                       
      err=.false. 
!                                                                       
      IF(arg.eq.'liquid')then 
         ofirstpres=ofirstsatpres 
         olastpres=olastsatpres 
         ptableprop=ptable5 
         stableprop=stable5 
         ntableprop=ntable5 
         ptablelimits=ptable8 
         ntablelimits=ntable8 
      ELSEIF(arg.eq.'vapor')then 
         ofirstpres=1 
         olastpres=olastsatpres 
         ptableprop=ptable6 
         stableprop=stable6 
         ntableprop=ntable6 
         ptablelimits=ptable9 
         ntablelimits=ntable9 
      ELSEIF(arg.eq.'supercritical')then 
         ofirstpres=ofirstsupcritpres 
         olastpres=olastsupcritpres 
         ptableprop=ptable7 
         stableprop=stable7 
         ntableprop=ntable7 
         ptablelimits=ptable10 
         ntablelimits=ntable10 
      ELSE 
         err=.true. 
      ENDIF 
!                                                                       
      RETURN 
      END SUBROUTINE setpointers_cupid                    
