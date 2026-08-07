!
      SUBROUTINE initialize_specific_variables_solid
!
!     Define initial user-defined solid variables 
!       
      USE SOLID_DATA   , ONLY: solid
      USE Zparam       , ONLY: pi
      USE Zconst1      , ONLY: vv_prob
      USE Zcoord3      , ONLY: aporous
      USE Zqvol        , ONLY: qvol_ice_solid
      USE Zzone        , ONLY: ncell_cond,nmaterial_c
      USE Znum_cell    , ONLY: n_fluid
!
      IMPLICIT NONE
!
      INTEGER i,ii
      INTEGER iOKr,iOKk
!      
      REAL(8) CpVol,Condu         
!  
!.....PAFS-POOL, SMALL-POOL, SMALL-POOL-3D
!  
      IF (vv_prob.eq.'PAFS-POOL') THEN
         DO i=1,ncell_cond
            IF(abs(nmaterial_c(i)).lt.50)THEN
               CALL mat_prop(ABS(nmaterial_c(i)),solid%tsol(i),CpVol,Condu,iOKr,iOKk)
               solid%rhocps(i)=CpVol/100.d0
               solid%conds(i)=Condu*1.d0
               IF(qvol_ice_solid(i).eq.0.d0) solid%conds(i)=Condu
             ENDIF
             solid%matnum(i)=ABS(nmaterial_c(i))
          ENDDO 
      ENDIF
!   
!.....sgbundle
!  
      IF (vv_prob.eq.'sgbundle') THEN
         DO i=1,ncell_cond
            ii=n_fluid(i)
            aporous(ii)=0.03d0
            IF(iabs(nmaterial_c(i)).lt.50) THEN
               CALL mat_prop(ABS(nmaterial_c(i)),solid%tsol(i),CpVol,Condu,iOKr,iOKk)
               solid%rhocps(i)=CpVol
               solid%conds(i)=Condu
             ENDIF
             solid%matnum(i)=ABS(nmaterial_c(i))         
          ENDDO
      ENDIF     
!                      
      END SUBROUTINE initialize_specific_variables_solid
