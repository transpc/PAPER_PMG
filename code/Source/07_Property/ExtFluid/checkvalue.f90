      SUBROUTINE checkvalue_cupid(arg,p,g,s,getprops,tleft,tright,deltt,err) 
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
!deck checkvalue                                                        
!                                                                       
!  $Id: checkvalue.ff,v 1.2 2001/04/02 17:41:41 dbarber Exp dbarber $   
!                                                                       
!  Cognizant engineer:  rwt                                             
!                                                                       
      USE STM_TBL_cupid  , ONLY: pcrit
!
      IMPLICIT none 
!                                                                       
      INCLUDE 'stcom.h' 
      INCLUDE 'newstcom.h' 
      INCLUDE 'gibbpnt.h' 
      INCLUDE 'sparms.h' 
!                                                                       
      REAL(8) cvg,deltt,g(11),p,s(*),slope,tleft,tright,tnaught,ttest 
      LOGICAL getprops(*),err 
      CHARACTER*(*)arg 
!                                                                       
      err=.false. 
!                                                                       
      IF(arg.eq.'liquid')then 
!                                                                       
!  Liquid must stay above ice line                                      
!                                                                       
         IF(s(temp).lt.tmin)then 
            err=.true. 
         ENDIF 
!                                                                       
!  Spinodal line test for superheated liquid                            
!  595.0K is the spinodal line temperature                              
!  at the triple point pressure                                         
!                                                                       
         tnaught=595.0d0 
         slope=(tcrit-tnaught)/(log(pcrit)-log(ptrip)) 
         ttest=tnaught+slope*(log(p)-log(ptrip)) 
         IF(s(temp).gt.ttest)then 
            err=.true. 
         ENDIF 
!                                                                       
      ELSEIF(arg.eq.'vapor')then 
!                                                                       
!  Vapor cannot exceed maximum table temperature                        
!                                                                       
         IF(s(temp).gt.tmax)then 
            err=.true. 
         ENDIF 
!                                                                       
!  Spinodal line test for subcooled vapor                               
!                                                                       
!  Test against 3 criteria for determining a good fit.                  
!  (1) specific volume > 0.0                                            
!  (2) specific heat at constant volume > 0.0                           
!  (3) specific volume*kappa > 0.0                                      
!  These are the same criteria used in the generator.                   
!                                                                       
         cvg=-s(temp)*(g(7)-(g(5)**2/g(3))) 
!rex+ 23 Feb 2001                                                       
         IF(cvg.lt.0.0d0.and.g(3).lt.0.0d0.and.s(cp).gt.0.0d0)then 
!          cvg = cp + t*g(5)**2/g(3)                                    
!          use cp to get cv, assume cp/cv=1.1, solve for g5             
            g(5)=sqrt(-0.1d0*s(cp)*g(3)/s(temp)) 
            s(beta)=g(5)/g(2) 
            cvg=-s(temp)*(g(7)-(g(5)**2/g(3))) 
         ENDIF 
!rex-                                                                   
         IF((getprops(vbar).and.(s(vbar).le.0.0d0)).or.(cvg.le.0.0d0)   &
         .or.(-g(3).le.0.0d0))then                                      
            err=.true. 
         ENDIF 
      ELSE 
!                                                                       
!  Supercritical fluid test                                             
!                                                                       
         IF((s(temp).lt.ttrip).or.(s(temp).gt.tmax))err=.true. 
!                                                                       
      ENDIF 
!                                                                       
!  Test to see if interpolated/extrapolated temperature is within       
!  an acceptable interval for all phases                                
!                                                                       
      tleft=tright-deltt 
      IF(s(temp).lt.(tleft-deltt))then 
!        err = .true.                                                   
      ENDIF 
      IF(s(temp).gt.(tright+deltt))then 
!        err = .true.                                                   
      ENDIF 
!                                                                       
      RETURN 
      END SUBROUTINE checkvalue_cupid                     
