!
     SUBROUTINE restart_write_cupid_master   
!     
     USE Zconst1   ,only:restart,cplmars,cplmaster
     USE Zconst2   ,only:dt
     USE Zio_unit  ,ONLY:unit_saveout,unit_log
     USE Zcore     ,only:myrank
     USE Ztimecon  ,only:time,itim,t_end,itim_restart
     USE Zmars     ,ONLY:time_mars                  
     USE MASTER4   ,ONLY:i_flag,TTIME,DTTR,TFC,TFS,TCOO,DCOO,BCOO,NXY_TH,NZ_TH
!     
     IMPLICIT NONE
!     
!DEC$IF defined (MCC)      
!DEC$ELSEIF defined (MCC_DLL)      
      !dec$ attributes dllexport :: restart_write_cupid_master      
!DEC$ELSEIF defined (SPACE)          
      !dec$ attributes dllexport :: restart_write_cupid_master
!DEC$ENDIF
!
     INTEGER nout,i,j
!     
     IF(cplmars.lt.1)RETURN
!     
     nout=1111
     IF(myrank.eq.0)WRITE(*,"(a,2e12.5,1x,2i8)")'Restart_write_cupid: dt,time,nout,itim=',dt,time,nout,itim 
     IF(myrank.eq.0)WRITE(unit_log,"(a,2e12.5,1x,2i8)")'Restart_write_cupid: dt,time,nout,itim=',dt,time,nout,itim 
     CALL restart_write(nout,1)
     IF(myrank.eq.0)WRITE(unit_saveout,*)dt,time,nout,itim         
!     
     IF(cplmaster.eq.1)THEN
        CALL t_masterC_assembly(5)
     ELSEIF(cplmaster.eq.2)THEN
        CALL t_masterC_rod(5)
     ENDIF    
!
     RETURN
     ENDSUBROUTINE restart_write_cupid_master                
