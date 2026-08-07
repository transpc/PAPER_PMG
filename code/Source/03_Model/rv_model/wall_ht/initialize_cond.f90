!     
      SUBROUTINE initialize_cond
!
!     This routine conductivity of solid
!              
      USE SOLID_DATA   , ONLY: solid
      USE Zzone        , ONLY: ncell_cond,num_max_zone,nmaterial_c
!
      IMPLICIT NONE
!
      INTEGER i
      INTEGER iokr,iokk
      REAL(8) CpVol,Condu         
!
!.....Call & save material properties
!    
      DO i=1,ncell_cond
         IF(IABS(nmaterial_c(i)).lt.50) THEN
            CALL mat_prop(IABS(nmaterial_c(i)),solid%tsol(i),CpVol,Condu,iokr,iokk)
            solid%conds(i)=Condu
         ENDIF
      ENDDO      
      END SUBROUTINE initialize_cond
