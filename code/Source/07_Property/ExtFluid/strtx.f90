      SUBROUTINE strtx_cupid(a,s,err) 
!deck strtx                                                             
!                                                                       
!                                                                       
!      strtx   - compute thermodynamic properties as a function of      
!                temperature and quality                                
!                                                                       
!      Calling sequence:                                                
!                                                                       
!                call  strtx (rp1,rp2,lp3)                              
!                                                                       
!      Parameters:                                                      
!                                                                       
!                rp1 = a   = steam tables (input)                       
!                                                                       
!                rp2 = s   = array into which the computed              
!                            thermodynamic properties are stored        
!                            (input,output)                             
!                                                                       
!                lp3 = err = error flag (output)                        
!                                                                       
!                                                                       
!      This routine adapted from sth2x1 routine written by R. J. Wagner 
!      for light water steam tables                                     
!                                                                       
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL(8) a(*),s(*) 
      LOGICAL err 
!                                                                       
!                                                                       
      REAL(8) dpsdts 
      EXTERNAL strsat_cupid,strx_cupid 
!                                                                       
!                                                                       
!--get saturation pressure                                              
!                                                                       
      CALL strsat_cupid(a,1,s(1),s(10),dpsdts,err) 
      IF(err) GOTO 10 
      s(2)=s(10) 
!                                                                       
!--compute thermodynamic properties as a function of quality            
!                                                                       
      CALL strx_cupid(a,s,err) 
!                                                                       
!--done                                                                 
!                                                                       
   10 RETURN 
!                                                                       
!                                                                       
      END SUBROUTINE strtx_cupid                          
