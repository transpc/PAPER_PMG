!
      SUBROUTINE choke_main
      
      USE Zrv_choke    , ONLY: time_cflow_on,vl_choke,vg_choke
      USE Ztimecon     , ONLY: time
      USE Zcore        , ONLY: np
!
      IMPLICIT NONE
!
!......Check choke model On/off
!
      vl_choke=0.d0
      vg_choke=0.d0
!      
      IF(time.lt.time_cflow_on) RETURN
!
!......Averaged parameter
!
      CALL choke_cell_avg
!
!......Critical flow model: choke, vl_choke calculate          
!
      CALL choke_model

      IF(np.gt.1) THEN
          call allreducei_r1(vl_choke)
          call allreducei_r1(vg_choke)
      ENDIF       
!         
      RETURN
      END SUBROUTINE choke_main