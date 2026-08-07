      SUBROUTINE getsatatt_cupid(tables,getprops,t,s,err) 
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
!deck getsatatt                                                         
!                                                                       
!  $Id: getsatatt.ff,v 1.2 2001/04/25 14:55:49 dbarber Exp dbarber $    
!                                                                       
!  Returns the saturation properties specified saturation               
!  temperature. This routine replaces sth2x0.F and sth2x1.F             
!  and sth2x2.F                                                         
!                                                                       
!  Cognizant engineer:  rwt                                             
!                                                                       
      IMPLICIT none 
!                                                                       
      INCLUDE 'stcom.h' 
      INCLUDE 'sparms.h' 
      INCLUDE 'newstcom.h' 
      INCLUDE 'gibbpnt.h' 
!                                                                       
      REAL(8) a(2),delt,delt2,s(*),tleft,tright,t,tables(*),xy(2),y(  &
      12),z                                                             
      INTEGER g,i,itleft,itleftlimit,itright,itrightlimit,ol,opl,opr,   &
      opsattleft,opsattright,or                                         
      LOGICAL err,getprops(26) 
!                                                                       
!      call timstart ('getsatatt')                                      
!                                                                       
      err=.false. 
!                                                                       
      itleftlimit=ofirstsattemp 
      itrightlimit=olastsattemp 
!                                                                       
      CALL getindex_cupid(tables,t,ptable1,itleftlimit,itrightlimit,itleft,   &
      err)                                                              
!                                                                       
      IF(err)return 
!                                                                       
      itright=itleft+1 
!                                                                       
      tleft=tables(ptable1+itleft) 
      tright=tables(ptable1+itright) 
!                                                                       
      opsattleft=ptable3+(itleft-itleftlimit)*ntable3 
      opsattright=ptable3+(itright-itleftlimit)*ntable3 
!                                                                       
!  Load a, b, and c arrays for 2-point formula                          
!                                                                       
      a(1)=tleft 
      a(2)=tright 
!                                                                       
      i=1 
      opl=opsattleft+9*(i-1)+1 
      opr=opsattright+9*(i-1)+1 
      delt=tright-tleft 
      delt2=(tright-tleft)**2 
      z=(t-tleft)/(tright-tleft) 
      g=opl+2 
      s(psat)=tables(g+1)+z*(tables(g+2)+z*(tables(g+3)+z*(tables(g+4)+ &
      z*(tables(g+5)+z*tables(g+6)))))                                  
!if def,logp                                                            
      s(psat)=exp(s(psat)) 
!endif                                                                  
!                                                                       
      i=2 
      ol=opsattleft+9*(i-1)+1 
      or=opsattright+9*(i-1)+1 
!                                                                       
      xy(1)=tleft 
      xy(2)=tright 
!                                                                       
      CALL satprops_cupid(xy,ol,or,tables,t,y,err) 
!                                                                       
      s(vsubf)=y(1) 
      s(usubf)=y(2) 
      s(betaf)=y(3) 
      s(kapaf)=y(4) 
      s(cpf)=y(5) 
      s(entpyf)=y(6) 
!                                                                       
      s(vsubg)=y(7) 
      s(usubg)=y(8) 
      s(betag)=y(9) 
      s(kapag)=y(10) 
      s(cpg)=y(11) 
      s(entpyg)=y(12) 
!                                                                       
!  Compute enthalpy saturation values                                   
!                                                                       
      s(hsubf)=s(usubf)+s(psat)*s(vsubf) 
      s(hsubg)=s(usubg)+s(psat)*s(vsubg) 
!                                                                       
!  Compute the mixture properties                                       
!                                                                       
      CALL mixprops_cupid(s,getprops) 
!                                                                       
!      call timstop ('getsatatt')                                       
      RETURN 
      END SUBROUTINE getsatatt_cupid                      
