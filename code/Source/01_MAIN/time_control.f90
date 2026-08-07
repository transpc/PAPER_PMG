!
      SUBROUTINE time_control_start
!
!     This routine defines initial time step and CFL number
!      
      USE Zparam          , ONLY: ndim
      USE Zcore           , ONLY: np,myrank
      USE Zconst1         , ONLY: vv_prob      
      USE Znum_cell       , ONLY: ncell      
      USE Ztimecon        , ONLY: ctrl_opt,ictrl,nctrl,time,t_end,t_end_ctrl,              &
                                   toutstep,toutstep_ctrl,cfl_ratio_max,cfl_ratio_max_ctrl, &
                                   smac,smac_ctrl,dt_opt,dt_opt_ctrl,dt_max,dt_max_ctrl,ictrl
      USE Ztplot          , ONLY: time2view,tplot_dt                                   
      USE viewData_common , ONLY: nframe,time2plot   
      USE Zio_unit        , ONLY: unit_log
!       
      IMPLICIT NONE
!
      LOGICAL, SAVE::INITIAL
!
      DATA INITIAL /.TRUE./      
!      
      IF(ctrl_opt.eq.0.and.INITIAL)then
         IF(myrank.eq.0)THEN
            WRITE(*,*)'          *********************************************************'
            WRITE(*,*)'          Time control off'
            WRITE(*,*)'          t_end,toutstep,cfl_ratio_max,smac,dt_opt,dt_max,   np'
            WRITE(*,31)t_end,toutstep,cfl_ratio_max,smac,dt_opt,dt_max,np
            WRITE(*,*)'          *********************************************************'
            WRITE(unit_log,*)'          *********************************************************'
            WRITE(unit_log,*)'          Time control off'
            WRITE(unit_log,*)'          t_end,toutstep,cfl_ratio_max,smac,dt_opt,dt_max,   np'
            WRITE(unit_log,31)t_end,toutstep,cfl_ratio_max,smac,dt_opt,dt_max,np
            WRITE(unit_log,*)'          *********************************************************'            
         ENDIF      
      ENDIF
   31 FORMAT(10x, 3(1pe10.3),1x,2i3,1x,1pe10.3,2x,1i3)          
!   
      IF(ctrl_opt.eq.0) GOTO 1
!
      IF(INITIAL) THEN
         INITIAL=.FALSE.
      ELSE
         GOTO 1
      ENDIF
!     
      t_end_ctrl(0)=0.d0
      DO ictrl=1,nctrl
         IF (time.ge.t_end_ctrl(ictrl-1) .and. time.lt.t_end_ctrl(ictrl))THEN
            t_end=t_end_ctrl(ictrl)
            toutstep=toutstep_ctrl(ictrl)
            cfl_ratio_max=cfl_ratio_max_ctrl(ictrl)
            smac=smac_ctrl(ictrl)
            dt_opt=dt_opt_ctrl(ictrl)
            dt_max=dt_max_ctrl(ictrl)
            IF(myrank.eq.0)THEN
               WRITE(*,*)'          *********************************************************'
               WRITE(*,*)'          Time control start:',ictrl
               WRITE(*,*)'          t_end,toutstep,cfl_ratio_max,smac,dt_opt,dt_max,ictrl,np'
               WRITE(*,30)t_end,toutstep,cfl_ratio_max,smac,dt_opt,dt_max,ictrl,np
               WRITE(*,*)'          *********************************************************'
               WRITE(unit_log,*)'          *********************************************************'
               WRITE(unit_log,*)'          Time control start:',ictrl
               WRITE(unit_log,*)'          t_end,toutstep,cfl_ratio_max,smac,dt_opt,dt_max,ictrl,np'
               WRITE(unit_log,30)t_end,toutstep,cfl_ratio_max,smac,dt_opt,dt_max,ictrl,np
               WRITE(unit_log,*)'          *********************************************************'               
            ENDIF
            EXIT
         ENDIF
      ENDDO
      ictrl=ictrl+1
!
   30 FORMAT(10x, 3(1pe10.3),1x,2i3,1x,1pe10.3,2x,2i3) 
!
    1 CONTINUE
!
      IF(myrank.eq.0)THEN
         vv_prob=trim(vv_prob)
         WRITE(*,40)vv_prob,ndim,ncell,np
         WRITE(unit_log,40)vv_prob,ndim,ncell,np         
         WRITE(*,*)'          ITIM-ITER     TIME        DP_MAX    TOTAL_MASS        DT'
      ENDIF
   40 FORMAT(11x,'Problem: ',1A20,1I1,'D',1x,1I7,1x,'cells',1x,1i3,1x,'cores')
!
!.....Print data for ParaView
!
      nframe = 0
      time2plot =toutstep
      IF (time2plot <= time) time2plot = time + toutstep
      CALL write_fieldview  
!
!.....Print data for real time view
!
      time2view=tplot_dt 
      IF (time2view <= time) time2view = time + tplot_dt
      CALL write_fieldview_tplot      
!
      RETURN
      END SUBROUTINE time_control_start
!
!------------------------------------------------------------------------------
!
      SUBROUTINE time_control
!
!     This routine assigns multiple time step and multiple CFL number
!      
      USE Zconst1         , ONLY: cplmaster
      USE Zcore           , ONLY: np,myrank
      USE Ztimecon        , ONLY: ctrl_opt,ictrl,nctrl,time,t_end,t_end_ctrl,                    &
                                   toutstep,toutstep_ctrl,cfl_ratio_max,cfl_ratio_max_ctrl, &
                                   smac,smac_ctrl,dt_opt,dt_opt_ctrl,dt_max,dt_max_ctrl
      USE Ztplot          , ONLY: time2view,tplot_dt                                    
      USE viewData_common , ONLY: nframe,time2plot,viwunit    
      USE Zio_unit        , ONLY: unit_log
!      
      IMPLICIT NONE
      INTEGER writefopt
      DATA writefopt/1/
!
      IF(ctrl_opt.eq.0)goto 1
!
! bug ictrl is greater than nctrl for last time control intervall
      IF(ictrl.gt.nctrl) goto 1
      IF (time.ge.t_end_ctrl(ictrl-1) .and. time.lt.t_end_ctrl(ictrl))THEN
         t_end=t_end_ctrl(ictrl)
         toutstep=toutstep_ctrl(ictrl)
         cfl_ratio_max=cfl_ratio_max_ctrl(ictrl)
         smac=smac_ctrl(ictrl)
         dt_opt=dt_opt_ctrl(ictrl)
         dt_max=dt_max_ctrl(ictrl)
         IF(myrank.eq.0)THEN
            WRITE(*,*)'          *********************************************************'
            WRITE(*,*)'          Time_control:',ictrl
            WRITE(*,*)'          t_end,toutstep,cfl_ratio_max,smac,dt_opt,dt_max,ictrl,np'
            WRITE(*,30)t_end,toutstep,cfl_ratio_max,smac,dt_opt,dt_max,ictrl,np
            WRITE(*,*)'          *********************************************************'
            WRITE(unit_log,*)'          *********************************************************'
            WRITE(unit_log,*)'          Time_control:',ictrl
            WRITE(unit_log,*)'          t_end,toutstep,cfl_ratio_max,smac,dt_opt,dt_max,ictrl,np'
            WRITE(unit_log,30)t_end,toutstep,cfl_ratio_max,smac,dt_opt,dt_max,ictrl,np
            WRITE(unit_log,*)'          *********************************************************'            
         ENDIF
         ictrl=ictrl+1
      ENDIF
   30 FORMAT(10x, 3(1pe10.3),1x,2i3,1x,1pe10.3,2x,2i3) 
  
    1 CONTINUE
!
!.... Print data for Paraview
!
      IF(time >= time2plot)THEN
         nframe = nframe + 1 
         CALL write_fieldview          
         time2plot = time2plot+toutstep
         writefopt=0
      ENDIF   
      IF(time >= time2view)THEN
         nframe = nframe + 1 
         CALL write_fieldview_tplot
         time2view = time2view+tplot_dt
         writefopt=0
      ENDIF            
      IF(time.ge.t_end) then
         IF(writefopt)then
            nframe=nframe+1
            CALL write_fieldview          
         ENDIF   
         CLOSE (viwUnit) 
      ENDIF 
!
      IF(cplmaster.gt.0)CALL change_tend_cupid        
!
      RETURN
      END SUBROUTINE time_control
!----------------------------------------------------------------------
      SUBROUTINE change_tend_cupid
!
      USE Ztimecon        , ONLY: ictrl,nctrl,time,t_end,t_end_ctrl
      USE Zconst1         , ONLY: restart
      USE Zcore           , ONLY: myrank
      USE MASTER4         , ONLY: dmaster_pass,mas_delay
      USE Zio_unit        , ONLY: unit_log
!     
      IMPLICIT NONE
!      
      LOGICAL,SAVE::initial
!
      INTEGER i
!      
      REAL(8) time_end     
!      
      DATA initial/.TRUE./
!      
      IF(dmaster_pass.eq.0)RETURN
!
      IF(restart.eq.0)THEN      
!      
         IF(initial)THEN
            initial=.FALSE.
            time_end=time+mas_delay
            t_end=time_end
            IF(myrank.eq.0)THEN
               WRITE(*,"(a,11f12.6)")'t_end=',t_end
               WRITE(unit_log,"(a,11f12.6)")'t_end=',t_end
            ENDIF   
            DO i=ictrl,nctrl
               t_end_ctrl(i)=time_end
               IF(myrank.eq.0)THEN
                  WRITE(*,"(a,1i5,1f12.6)")'i,t_end_ctrl(i)=',i,t_end_ctrl(i)
                  WRITE(unit_log,"(a,1i5,1f12.6)")'i,t_end_ctrl(i)=',i,t_end_ctrl(i)
               ENDIF   
            ENDDO 
         ENDIF
!      
      ENDIF   
!
      RETURN    
      ENDSUBROUTINE change_tend_cupid
