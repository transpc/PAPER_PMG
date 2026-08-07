      SUBROUTINE rv_celmdr(ctemp,ey,xnu) 
!
!     calculates cladding young's modulus, and poisson's ratio     
!
      IMPLICIT NONE
!                                                                       
!     input                                                           
!        ctemp   - cladding temperature                              
!     output                                                          
!        ey      - young's modulus                                   
!        xnu     - poisson's ratio                                   
!                                                                       
!     taken from matpro-11, revision 1       
!
      REAL(8) ctemp,ey,xnu      
!
      IF(ctemp.lt.1090.d0) THEN
         ey=1.088d11-5.475d7*ctemp 
      ELSEIF(ctemp.le.1240.d0) THEN
         ey=4.912d10-4.827d7*(ctemp-1090.d0) 
      ELSE
         ey=DMAX1(1.0d10,9.21d10-4.05d7*ctemp) 
      ENDIF   
! 
      xnu=0.30d0 
!  
      RETURN 
      END SUBROUTINE rv_celmdr                         
