      SUBROUTINE strpu2_cupid(a,s,it,err) 
!deck strpu2                                                            
!                                                                       
!                                                                       
!      strpu2  - compute thermodynamic properties as a function of      
!                pressure and internal energy, when the saturation      
!                temperature is already known                           
!                                                                       
!      Calling sequence:                                                
!                                                                       
!                call  strpu2 (rp1,rp2,ip3,lp4)                         
!                                                                       
!      Parameters:                                                      
!                                                                       
!                rp1 = a   = steam tables (input)                       
!                                                                       
!                rp2 = s   = array into which the computed              
!                            thermodynamic properties are stored        
!                            (input,output)                             
!                                                                       
!                ip3 = it  = flag indicating physical state of fluid,   
!                            i.e., liquid, vapor, superheated vapor     
!                            (output)                                   
!                                                                       
!                lp4 = err = error flag (output)                        
!                                                                       
!                                                                       
!      This routine adapted from sth2xf (entry point in sth2x6) written 
!      by R. J. Wagner for light water steam tables                     
!                                                                       
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL(8) a(*),s(*) 
      INTEGER it 
      LOGICAL err 
!                                                                       
!                                                                       
      EXTERNAL strpu_cupid 
!                                                                       
!                                                                       
!--get thermodynamic properties                                         
!                                                                       
      it=1 
      CALL strpu_cupid(a,s,it,err) 
!                                                                       
!--done                                                                 
!                                                                       
      RETURN 
!                                                                       
!                                                                       
      END SUBROUTINE strpu2_cupid                         
