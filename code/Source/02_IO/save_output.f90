!
      SUBROUTINE save_output(nout)
!
!     This routine write time step, problem time, iteration number on the screen
!
      USE Zcore        , ONLY: myrank
      USE Zconst1      , ONLY: save_option,noutput,restart,cplmaster
      USE Zconst2      , ONLY: dt
      USE Zio_unit     , ONLY: unit_saveout
      USE Ztimecon     , ONLY: time,treststep,t_end_ctrl,treststep_ctrl,nctrl,ctrl_opt,itim
!
      IMPLICIT NONE
!
      LOGICAL, SAVE::initial, move_step
!
      INTEGER nout,restart_nout,i,ictrl2,ictrl2_opt
      INTEGER, SAVE::nout2,ictrl2_o
!      
      REAL(8) restart_time,time_out,time_out_next
      REAL(8), SAVE:: time_out_start
!
      SAVE restart_time, restart_nout
!
      DATA initial,move_step /.true.,.true./      
!      
      IF(nout.lt.30000000)THEN
         IF(itim.eq.1)THEN
            nout=nout+1 
         ENDIF
!
!........1: step based output, 0: time based output
!
         IF(save_option.eq.1)THEN   
            IF(MOD(itim,noutput).eq.0)THEN
               CALL restart_write(nout)
               IF(myrank.eq.0)WRITE(unit_saveout,*)dt,time,nout,itim
               nout=nout+1               
               IF(cplmaster.gt.0) CALL save_restart_master_only 
            ENDIF
         ELSE 
            IF(restart.ne.0)THEN
               IF(initial)THEN
                  restart_time=time-dt
               ENDIF
            ELSE
               restart_time=0
               restart_nout=0
               IF(initial)THEN               
                 time_out_start=time-dt
               ENDIF
            ENDIF                      
 !           
            ictrl2_opt=0      
            IF(ctrl_opt.eq.0)THEN
                time_out=time_out_start+treststep*nout
                time_out_next=time_out_start+treststep*(nout+1)              
                ictrl2=0             
            ELSE   
               ictrl2=1
               time_out=0.0d0
               time_out_next=0.0d0               
               t_end_ctrl(0)=0.0d0 
               treststep_ctrl(0)=1.0d0 
!
!..............Find ictrl and save as ictrl2
!
               DO i=1,nctrl
                  ictrl2=i
                  IF (time.gt.t_end_ctrl(i-1) .and. time.le.t_end_ctrl(i)) EXIT
               ENDDO  
               IF(initial)ictrl2_o=0  
               IF (ictrl2_o.ne.ictrl2)THEN
                  nout2=nout-1               
                  ictrl2_opt=1   
!                  IF(initial.eq..false.) move_step=.false.
                  IF(.not.(initial)) move_step=.false.
               ENDIF
               ictrl2_o=ictrl2
!
!..............Find time to write restart file
!
               IF(restart.eq.0) then
                  IF(ictrl2.eq.1) nout2=0 
               ENDIF
!               
               IF(restart.ne.0.and.move_step)THEN               
                  time_out=restart_time+treststep_ctrl(ictrl2)*(nout-nout2)
                  time_out_next=restart_time+treststep_ctrl(ictrl2)*(nout+1-nout2)                  
               ELSE             
                  time_out=t_end_ctrl(ictrl2-1)+treststep_ctrl(ictrl2)*(nout-nout2)
                  time_out_next=t_end_ctrl(ictrl2-1)+treststep_ctrl(ictrl2)*(nout+1-nout2)
               ENDIF
!
            ENDIF 
!
            IF((time.gt.time_out.and.time.le.time_out_next).or.(ictrl2.ne.1.and.ictrl2_opt.eq.1))THEN 
               IF(restart.ne.0) nout=nout+1
               CALL restart_write(nout) 
               IF(myrank.eq.0) WRITE(unit_saveout,*)dt,time,nout,itim               
!               IF(udfl_gid_output) CALL GID_out_cell_lsj(time,nout)                           
               IF(restart.eq.0) nout=nout+1
               IF((time.le.time_out.or.time.gt.time_out_next).and.(ictrl2.ne.1.and.ictrl2_opt.eq.1))nout2=nout2+1
               IF(cplmaster.gt.0) CALL save_restart_master_only 
            ENDIF           
!            
         ENDIF
!
      ENDIF
!
      initial=.false.               
!
      RETURN
      END SUBROUTINE save_output
