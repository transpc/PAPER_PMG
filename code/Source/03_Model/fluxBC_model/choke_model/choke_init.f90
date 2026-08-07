!
      SUBROUTINE choke_init

      USE Ztimecon      , ONLY: time
      USE Zrv_choke     , ONLY: time_cflow_on,choke,choke,choke_update
!
      IMPLICIT NONE

      LOGICAL, SAVE::initial_throatfinding
      DATA initial_throatfinding /.true./
!
!......Find throat
!
      IF(initial_throatfinding) THEN
         CALL choke_define_fluxBC
         initial_throatfinding=.false.
      ENDIF  
!
!......Initialize choke_updae
!
!      choke=Local update for velocity field (LOGICAL)
!      choke_update=Global upodate for pressure field (INTEGER)        
!      
      IF(time.ge.time_cflow_on) THEN
         choke=.false.
         choke_update=0
      ENDIF
!            
      END SUBROUTINE choke_init
