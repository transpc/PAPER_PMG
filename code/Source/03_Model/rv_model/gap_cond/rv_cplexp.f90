      SUBROUTINE rv_cplexp(cltave,dtdt,hoop,irupt,eplas,block) 
!
!  Calculates the plastic strain using data in nureg-0630.              
!  Assumes that the gap conductance model is implemented.               
!  Skip the calculations if the hoop stress is negative.                
!                                                                       
! ******** input parameters *********                                   
!  cltave = average clad temperature (k)                                
!  dtdt   = clad heatup rate (k/s)                                      
!  hoop   = clad hoop stress (pa)                                       
!  irupt  = rupture flag = 1 if rod has ruptured, 0 otherwise           
!           set in ht1tdp after return from madata if block > 0.0       
! ******** output parameters *********                                  
!  eplas  = plastic strain in the clad (m/m)                            
!  block  = flow channel blockage (%) (calc when failure occurs)        
!         = 0.0 if no fuel rod rupture has occurred                     
!                            calculate rupture temperature              
!                            and plastic deformation temperature        
!                            tplask = plastic deformation temp (k)      
!                            trupk  = rupture temp (k)                  
!
      IMPLICIT NONE
!                   
      REAL(8) cltave,dtdt,hoop,eplas,block 
      INTEGER irupt                                                     
!                                                                       
!.....Local variables.                                                     
!
      REAL(8) tplask,trupk 
!                                                                       
      IF(hoop.le.0.0d0.or.irupt.ne.0)then 
         eplas=0.0d0 
         block=0.0d0 
      ELSE 
!
!........Compute tplask, trupk
!      
         CALL rv_ruplas(hoop,dtdt,tplask,trupk) 
!
!........Compute the deformation of the claddig: eplas, block
!
         CALL rv_plstrn(cltave,tplask,trupk,dtdt,eplas,block) 
!         
      ENDIF 
!      
      RETURN 
      END SUBROUTINE rv_cplexp
