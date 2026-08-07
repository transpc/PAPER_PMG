      INTEGER temp,pres,vbar,ubar,hbar,beta,kapa,cp,qual,tsat,psat,     &
      vsubf,vsubg,usubf,usubg,hsubf,hsubg,betaf,betag,kapaf,kapag,cpf,  &
      cpg,is23,entpy,entpyf,entpyg                                      
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
!comdeck sparms                                                         
!                                                                       
!  $Id: sparms.hh,v 1.1 2001/02/01 23:15:49 r5qa Exp dbarber $          
!                                                                       
!  These are the parameters that can be used                            
!  for subscripts in the s array                                        
!                                                                       
      PARAMETER (temp=1,pres=2,vbar=3,ubar=4,hbar=5,beta=6,kapa=7,cp=8, &
      qual=9,tsat=10,psat=10,vsubf=11,vsubg=12,usubf=13,usubg=14,hsubf= &
      15,hsubg=16,betaf=17,betag=18,kapaf=19,kapag=20,cpf=21,cpg=22,    &
      is23=23,entpy=24,entpyf=25,entpyg=26)                             
!                                                                       
