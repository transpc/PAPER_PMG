      SUBROUTINE gettsat_cupid(tables,p,s,err) 
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
!deck gettsat                                                           
!                                                                       
!  $Id: gettsat.ff,v 1.2 2001/04/25 14:55:49 dbarber Exp dbarber $      
!                                                                       
!  Returns the saturation temperature given pressure.                   
!  This routine is similiar to getsatatp except that it only            
!  returns the saturation temperature.                                  
!                                                                       
!  Cognizant engineer:  rwt,gam                                         
!                                                                       
      IMPLICIT none 
!                                                                       
      INCLUDE 'stcom.h' 
      INCLUDE 'sparms.h' 
      INCLUDE 'newstcom.h' 
      INCLUDE 'gibbpnt.h' 
!                                                                       
      REAL(8) delt,p,phi,plo,s(*),t,tables(*),z 
      INTEGER iphi,iplo,ofirstpres,olastpres,otsatplo,otsatphi 
      LOGICAL err 
!                                                                       
!      call timstart ('gettsat')                                        
!                                                                       
      err=.false. 
!                                                                       
      ofirstpres=ofirstsatpres 
      olastpres=olastsatpres 
!                                                                       
      CALL getindex_cupid(tables,p,ptable2,ofirstpres,olastpres,iplo,err) 
!                                                                       
      IF(err)return 
!                                                                       
      iphi=iplo+1 
      plo=tables(ptable2+iplo) 
      phi=tables(ptable2+iphi) 
!                                                                       
      otsatplo=ptable4+(iplo-ofirstpres)*ntable4 
      otsatphi=ptable4+(iphi-ofirstpres)*ntable4 
!                                                                       
!  compute saturation temperature using the 2-point formula             
!                                                                       
      delt=log(phi)-log(plo) 
      z=(log(p)-log(plo))/delt 
      t=tables(otsatplo+4)+z*(tables(otsatplo+5)+z*(tables(otsatplo+6)+ &
      z*(tables(otsatplo+7)+z*(tables(otsatplo+8)+z*tables(otsatplo+9)))&
      ))                                                                
!                                                                       
      s(tsat)=t 
!                                                                       
!      call timstop ('gettsat')                                         
!                                                                       
      RETURN 
      END SUBROUTINE gettsat_cupid                        
