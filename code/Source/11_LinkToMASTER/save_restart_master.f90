!
      SUBROUTINE save_restart_master(nout)
!
      USE Zcore        , ONLY: myrank
      USE Zconst1      , ONLY: cplmaster
      USE Zconst2      , ONLY: dt
      USE Zio_unit     , ONLY: unit_saveout
      USE Ztimecon     , ONLY: time,itim
!
      IMPLICIT NONE
!
      INTEGER nout
!
!.....Save restart file, in case of MASTER/CUPID
!
      CALL restart_write(nout)
      IF(myrank.eq.0)WRITE(unit_saveout,*)dt,time,nout,itim         
      IF(cplmaster.eq.1)THEN
          CALL t_masterC_assembly(5)
      ELSEIF(cplmaster.eq.2)THEN    
          CALL write_trod
          CALL t_masterC_rod(5)
      ENDIF    
!
      END SUBROUTINE save_restart_master
!
!
      SUBROUTINE save_restart_master_only
!
      USE Zconst1      , ONLY: cplmaster
!
      IMPLICIT NONE
!
!.....Save restart file, in case of MASTER/CUPID
!
      IF(cplmaster.eq.1)THEN
          CALL t_masterC_assembly(5)
      ELSEIF(cplmaster.eq.2)THEN    
          CALL write_trod
          CALL t_masterC_rod(5)
      ENDIF    
!
      END SUBROUTINE save_restart_master_only
