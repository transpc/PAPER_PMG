      SUBROUTINE surftn_cupid(n,ifluid,temp,sigma) 
!
!     jjj  10/31/1997
!     bsw-transp : copy from MARS. 2006.4.5 
!                                                                       
!     surftn_cupid  - compute surface tension for given fluid                
!                                                                       
!     Author:   J. E. Tolli, EG&G Idaho, Inc.                          
!     Date:     7/89                                                   
!     Modified: 12/97 (IAPS '75 for light water) by WJL
!     Modified: 03/03 (CO2, He) by WJL
!                                                                 
!                                                                       
      IMPLICIT NONE 
!
      REAL(8) sigma(*),temp(*) 
      INTEGER ifluid
!                                                                       
!  Local variables:                                                     
!                                                                       
      REAL(8) term,term1
      INTEGER i,n 
!     INCLUDE 'contrl.h'
!                                                                       
!  Data statements:                                                     
!                                                                       
!  Execution:                                                           
!                                                                       
!                                                                       
!                                                                       
!--branch to correlation for given fluid                                
!                                                                       
      IF(ifluid.eq.1.or.ifluid.eq.15)then 
!
!   Light Water
!
!   Correlation from IAPS '75 Release
!   Reference:  ASME '93 Steam Table, Appendix 7
!   Modified by Won-Jae Lee on Dec. '97
!
!
         DO  i=1,n
               term=max(647.15d0-temp(i),0.0d0) 
               term1=term/647.15d0
               sigma(i)=0.2358d0*(1.d0-0.625d0*term1)*term1**1.256d0                               
!  Set minimum value to 1.d-9 (for supercritical fluid): Won-Jae Lee
               sigma(i)=max(sigma(i),1.d-9)
         END DO 
!                                                                       
      ELSEIF(ifluid.eq.2)then 
!                                                                       
!                                                                       
!--heavy water                                                          
!                                                                       
!--correlation from Flowtran program (Savannah River);                  
!--FORTRAN coding from C. S. Miller and R. J. Wagner, EG&G Idaho, Inc.  
!                                                                       
         DO  i=1,n
            term=min(max(temp(i)-273.15d0,1.0d-6),366.33d0) 
            sigma(i)=80.61755d-3-0.22006382d-3*term 
!  Set minimum value to 1.d-9 (for supercritical fluid): Won-Jae Lee
            sigma(i)=max(sigma(i),1.d-9)
         END DO 
!---------------------------------------------------------------------------
      ELSEIF(ifluid.eq.3)then 
!                                                                       
!   CO2 Surface Tension - Modified by Won-Jae Lee on Oct. '02
!	Rathjen, W. and Straub, J.,
!	 "Temperature dependence of surface tension, coexistence curve, and vapor
!	 pressure of CO2, CClF3, CBrF3, and SF6,"
!	 Chapter 18 in:  Heat Transfer in Boiling, New York, Academic Press,
!	 pp. 425-451, 1977.
!
!	sigma in N/m
!
!                                                                       
         DO  i=1,n
            term1=min(max(temp(i),216.58d0),304.1282d0) 
            term=1.d0-term1/304.1282d0   ! Diff from Critical Temperature by R&S (304.17 K)
            sigma(i)=0.084497d0*term**1.28d0 
!
!  Set minimum value to 1.d-9 (for supercritical fluid): Won-Jae Lee
               sigma(i)=max(sigma(i),1.d-9)
!
         END DO 
!---------------------------------------------------------------------------
      ELSEIF(ifluid.eq.4)then 
!                                                                       
!  He-4 Surface Tension by Won-Jae Lee 
!  Ref: Chemical Properties Handbook 
!      sigma = A * (1 - T/Tc)**n
!  http://www.knovel.com/knovel2/Toc.jsp?SpaceID=10093&BookID=49
!  This model is Not EFFECTIVE in the calculation but NEEDED for completeness
!      Application Ranges: 1.76 K < T(K) < 5.2 K
!    units: t (K), surftn (N/m))                                             
!
!                                                                       
         DO  i=1,n
            term=min(max(temp(i),1.76d0),5.2d0) 
            term=1.d0-term/5.2d0 
            sigma(i)=0.511d0*term**1.003d0 
!  Set minimum value to 1.d-9 (for supercritical fluid): Won-Jae Lee
            sigma(i)=max(sigma(i),1.d-9)
         END DO 
!
!---------------------------------------------------------------------------
      ELSEIF(ifluid.eq.5)then 
!                                                                       
!  H2 Surface Tension by Won-Jae Lee 
!  NOT IMPLEMENTED YET - 
!
!  Input error for ifluid
!
!
!---------------------------------------------------------------------------
      ELSEIF(ifluid.eq.6)then 
!                                                                       
!  O2 Surface Tension by Won-Jae Lee 
!  Lemmon E.W. & Penoncello, S.G. "The surface tension of air and air component Mixture"
!  Adv. Cryo. Eng. v39, 1994 
!
         DO  i=1,n
            term1=min(max(temp(i),54.361d0),154.581d0) 
            term=1.d0-term1/154.581d0   ! Diff from Critical Temperature by R&S (304.17 K)
            sigma(i)=0.038612652d0*term**1.228d0 
!
!  Set minimum value to 1.d-9 (for supercritical fluid): Won-Jae Lee
            sigma(i)=max(sigma(i),1.d-9)
         END DO 
!
!---------------------------------------------------------------------------
      ELSEIF(ifluid.eq.7)then 
!                                                                       
!  N2 Surface Tension by Won-Jae Lee 
!  Lemmon E.W. & Penoncello, S.G. "The surface tension of air and air component Mixture"
!  Adv. Cryo. Eng. v39, 1994 
!
         DO  i=1,n
            term1=min(max(temp(i),63.151d0),126.192d0) 
            term=1.d0-term1/126.192d0   ! Diff from Critical Temperature by R&S (304.17 K)
            sigma(i)=0.029324108d0*term**1.259d0 
!
!  Set minimum value to 1.d-9 (for supercritical fluid): Won-Jae Lee
            sigma(i)=max(sigma(i),1.d-9)
         END DO 
!
!---------------------------------------------------------------------------
      ELSEIF(ifluid.eq.8)then 
!>>
!  D. LMR-K.S. Ha for liquid metal properties - Sodium
!  Na Surface Tension
!
!  melting temperature : 93 C
         DO  i=1,n
            term=max(temp(i)-273.15d0,93.d0)
            sigma(i)=0.2067d0-1.d-4*term
         END DO 
!<<
!---------------------------------------------------------------------------
      ELSEIF(ifluid.eq.9)then 
!                                                                       
!  Argon Surface Tension by Won-Jae Lee 
!  Lemmon E.W. & Penoncello, S.G. "The surface tension of air and air component Mixture"
!  Adv. Cryo. Eng. v39, 1994 
!
         DO  i=1,n
            term1=min(max(temp(i),83.8058d0),150.6633d0) 
            term=1.d0-term1/150.6633d0   ! Diff from Critical Temperature by R&S (304.17 K)
            sigma(i)=0.037898063d0*term**1.278d0 
!
!  Set minimum value to 1.d-9 (for supercritical fluid): Won-Jae Lee
            sigma(i)=max(sigma(i),1.d-9)
         END DO 
!
!---------------------------------------------------------------------------
      ELSEIF(ifluid.eq.10)then 
!                                                                       
!  Air Surface Tension by Won-Jae Lee 
!  NOT IMPLEMENTED YET - 
!
!  Input error for ifluid
!
!
!---------------------------------------------------------------------------
      ELSEIF(ifluid.eq.11)then 
!>>
!  D. LMR-K.S. Ha for liquid metal properties - LBE
!  Pb-Bi Surface Tension
!
         DO  i=1,n
            term=max(temp(i)-273.15d0,93.d0)
            sigma(i)=0.437d0-6.6d-5*term
         END DO 
!<<
!---------------------------------------------------------------------------
      ELSEIF(ifluid.eq.13)then 
!                                                                       
!  FLIBE Surface Tension by Won-Jae Lee 
!  NOT IMPLEMENTED YET - 
!
!  Input error for ifluid
!
!
!
!---------------------------------------------------------------------------
      ELSEIF(ifluid.eq.14)then 
!                                                                       
!  FLINABE Surface Tension by Won-Jae Lee 
!  NOT IMPLEMENTED YET - 
!
!  Input error for ifluid
!
!
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
      ELSEIF(ifluid.eq.16)then 
!>>
!   R12
!
         DO  i=1,n
            term=max(temp(i),120.0d0) 
            sigma(i)=0.0565d0*(1.d0-0.0026d0*term)**1.26403d0       !!!CYJ 200107 [N/m], R2=1
            sigma(i)=max(sigma(i),1.d-9)
         END DO 
!<<
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
      ELSEIF(ifluid.eq.17)then 
!>>
!  r134a
!
         DO  i=1,n
            term=max(temp(i),170.0d0) 
            sigma(i)=0.05801d0*(1.d0-0.00267d0*term)**1.24099d0       !!!CYJ 200107 [N/m], R2=1
            sigma(i)=max(sigma(i),1.d-9)
         END DO 
!<<
!---------------------------------------------------------------------------
      ELSE 
!                                                                       
!  Input error for ifluid                                               
!                                                                       
      ENDIF  
!
      END SUBROUTINE surftn_cupid
