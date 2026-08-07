!DEC$IF defined (SPACE)
    SUBROUTINE trans_cupid(time_c)
!
      USE Zbicg     ,only:pbcgind,pbcgind_max,pbcgsig 
      USE Zconst1   ,only:iat,mboron,restart
      USE Zconst2   ,only:i_repeat,dt,dtr
      USE Zcore     ,only:np,myrank
      USE Ztimecon  ,only:ctrl_opt,smac,time,itim,t_end
      USE Zzone     ,only:ncell_cond_all
!      
      IMPLICIT NONE
!
      INCLUDE 'c3com.h'
      INCLUDE 'c3com_space.h'       
!dec$ attributes dllexport :: trans_cupid
      !INCLUDE 'c3trans.h'
      !INCLUDE 'c3mserr.h'
      !INCLUDE 'contrller.h'
!DEC$IF defined (mpi_flag)
      INCLUDE 'mpif.h'
!DEC$ENDIF      
      INTEGER:: iv,ip,is
      logical:: s_flg(2),repeat
      REAL(8):: tarray(2), tresult, DTIME, time1
      COMMON/CUPID_MAIN/iv,ip,is, &
                        s_flg,repeat, &
                        tarray, tresult, DTIME, time1 
!        
      DATA iv,ip,is/1,2,3/
!
      REAL(8) time_c
      INTEGER nout,restart_out
      COMMON/CUPID_LOCAL/nout,restart_out      
!
      INTEGER cupid__ekd__,cupid__ekd__release,cupid__ekd__debug
      COMMON/Zdist/cupid__ekd__,cupid__ekd__release,cupid__ekd__debug   
!
!     SAPCE coupling instead of c3trans.h, c3mserr.h, contrller.h
      INTEGER iofail, jflag, i_mass,nmfail ! nstep_c !by pik
      REAL(8) timet, dtvmx, varer
!     end space 
      cupid__ekd__=0
      cupid__ekd__release=0
      cupid__ekd__debug=0
!
      CALL print_logo
!
      CALL c3com_copy_S2C      
!
      CALL c3com_print_space
!
!.....Prepare for a time step advancement
!
      IF(i_where.eq.1)THEN
!      나중에 필요 변수를 SPACE 변수로 setup 해야 함
         WRITE(*,*)'i_where==1,starts'
         iofail=0 
         jflag=0 
         i_mass=0 
         nstep_c=0 
         IF(restart.eq.0)itim=0 
         nout=0
         !timet=0.0d0  !only for new problem
         !timet=timehy !from mars
         timet=time    !from cupid
         restart_out=0
         i_repeat=0
         flag_cobra=.true.!pik-back-ins 
         
!         
!........Controlling cupid time constants
!
         CALL time_control_start
!
!.....Prepare to begin a time-step advancement (100 continue)                      
!
         WRITE(*,*)'i_where==1,ends'
      ELSEIF(i_where.eq.2)THEN
!  
         WRITE(*,*)'i_where==2,starts'

         IF(jflag.eq.0)then 
            i_mass=i_mass+1 
            nmfail=0 
            itim=itim+1 !back-pik-ins
         ENDIF 
         IF(.not.flag_cobra)then !back-pik-ins 
            CALL scalar_reset  
            IF(myrank.eq.0)WRITE(*,"(11x,a,l,l,1i3)")'##i_where2, scalar_reset: R,C,j flags=',flag_relap,flag_cobra,jflag          
         ENDIF           
!                                                                       
!........return here for smaller delt if Outer iteration fails             
!
         125 CONTINUE 
         nstep_c=nstep_c+1 
         dtvmx=0.0d0 
         flag_cobra=.true. 
         IF(jflag.ne.0)then
            flag_cobra=.false. 
         ENDIF   
!
!.....Prepare for a time-step advancement: determine dt
!
         WRITE(*,*)'i_where==2,ends'

      ELSEIF(i_where.eq.3)THEN
               WRITE(*,*)'i_where==3,starts'
!      
         IF(.not.flag_relap.and.flag_cobra)then 
            jflag=3 
            timet=timet-dt_super 
            CALL scalar_reset            
            IF(myrank.eq.0)WRITE(*,"(11x,a,l,l,1i3)")'## i_where3, scalar_reset: R,C,j flags=',flag_relap,flag_cobra,jflag    
         ENDIF 
!
!........Calculate dt in CUPID
!         
         !atlas-delete CALL set_dt 
         IF(.not.flag_cobra)then !pbcgind.gt.0.and.pbcgind.le.pbcgind_max
            dt=dt*0.5d0
            dtr=1.0d0/dt            
            IF(myrank.eq.0)WRITE(*,"(11x,a,1pe20.10)")'## i_where3, reduce dt=',dt
         ENDIF   
         CALL set_dt !atlas-ins
         dt_cobra=dt 
!                                                                       
!.....Set fluid boundary conditions                          
!      
         s3dt_cobra=dt_cobra !cupid_space_debug
         WRITE(*,*)'i_where==4,ends'

      ELSEIF(i_where.eq.4)THEN
         WRITE(*,*)'i_where==4,starts'

         iofail=0 !back-pik-ins
         dt=dt_super 
         dtr=1.0d0/dt         
!
!........Transfer informations from CUPID to MARS
!           
         CALL cupid2mars   !MCC-jjj
!
!........Reflect user condition
!         
         CALL user_def_inp(itim)
!
!........Store old value 
!         
         CALL shift_solutions
!
!.....Advance a time step in CUPID
!
         WRITE(*,*)'i_where==4,ends'

      ELSEIF(i_where.eq.5)THEN
         WRITE(*,*)'i_where==5,starts'
      
!
         varer=0.0d0 
!
!........Comunicate major parameters for parallel computing
!
         IF(np.gt.1) CALL communicate_allb  
!
!........Calculcate physical models
!
         CALL calc_models
!
!........Set the donor properties at the CUPID-MARS interface
!
         CALL donor_CMI(timet)
!
         98 pcnv=.true.
!
!........Explicit momentum calculation
!
         CALL calc_momentum
!
!........Coupled Scalar and Pressure calculation
!
         CALL calc_scalar(s_flg,pcnv) 
         
         IF(.not.pcnv) THEN
            iter_p=iter_p+1
            CALL prn_iter
            GOTO 98
          ENDIF         
!
!........Calculate the 1D/3D interface
!
         CALL v1d3d_cupid   
!      
!.....Calculation failure in MARS, Reset  !MCC-jjj-NEXT
!          
         WRITE(*,*)'i_where==5,ends'

      ELSEIF(i_where.ge.6.and.i_where.le.8)THEN
!
!.....Air-repeat in MARS, Reset           !MCC-jjj-NEXT
!
      ELSEIF(i_where.eq.10)THEN
!
!.....Now, finish CUPID part and update time
!
      ELSEIF(i_where.eq.9)THEN
         WRITE(*,*)'i_where==9,starts'
      
!      
         !atlas IF(ncell_cond_all.gt.0) CALL calc_solid  !because it can be failed!     
!
!........Check scalar errors & goto 125 
!
         IF(pbcgind.gt.0.and.pbcgind.le.pbcgind_max)then ! back-pik-ins
            IF(dt.lt.1.d-6)THEN
               PRINT*,'          *** Pressure iteration failed in CUPID! ***'
               WRITE(*,"(a,1i5,1e12.5)")'          pbcgind,dt=',pbcgind,dt
               PAUSE
               STOP
            ENDIF
            pbcgsig=1
            jflag=2
            goto 125 ! => cobra_flase=false => reduce dt_cobra => not increase dt_cobra during 15 run
         ENDIF
!         IF(pbcgind.gt.0.and.pbcgind.le.pbcgind_max)THEN
!             CALL check_iteration
         !  jflag=2
         !  goto 125
         !ENDIF  
!
         jflag=0 !back-pik-ins 
         flag_cobra=.true. 
!      
!........Prepare to solve the boron transport equation
!........Update Material properties
!........Update Boron concentration
!
!........Conduction calculation
!
         IF(ncell_cond_all.gt.0) CALL calc_solid      
!
!........Controlling time constants
!
         CALL time_control 
!
!........Save output 
!
         CALL user_def_output(nout) 
!
!........Calculation Finishes at t_end
!
!        IF(time.gt.t_end)exit
!
!........Advance time
!         
         timet=timet+dt
 	     time_c=timet
         time=time_c
!
         WRITE(*,*)'i_where==9,ends'

      ENDIF
!
      CALL c3com_print_space
!
      IF(flag_cobra)then
         jflag_cobra=1
      ELSE
         jflag_cobra=0
      ENDIF
!
      RETURN
      END SUBROUTINE trans_cupid 
!-----------------------------------------------------------------------------------
      SUBROUTINE print_logo
!
      USE Zcore    ,ONLY: myrank
!            
      IMPLICIT NONE
!      
      LOGICAL,SAVE::initial
!      
      DATA initial/.TRUE./
!      
      IF(initial.eq..FALSE.)RETURN
!      
      IF(myrank.eq.0)THEN
          WRITE(*,*)''                               
          WRITE(*,*)''                               
          WRITE(*,*)' SSSSS  PPPPP    A    CCCCC  EEEEE     //  CCC  U   U  PPPPP   IIII  DDDD  '
          WRITE(*,*)' S      P   P   A A   C      E        // CC     U   U  P    P   II   D   D '
          WRITE(*,*)' SSSSS  PPPPP  AAAAA  C      EEEEE   //  C      U   U  PPPPP    II   D   D '
          WRITE(*,*)'     S  P      A   A  C      E      //   CC     U   U  P        II   D   D '
          WRITE(*,*)' SSSSS  P      A   A  CCCCC  EEEEE //      CCC   UUU   P       IIII  DDDD  ' 
          WRITE(*,*)
          WRITE(*,*)''     
          WRITE(*,*)'                                      v1.0                    '     
          WRITE(*,*)''     
          WRITE(*,*)'                                Copyright (c) 2018            '     
          WRITE(*,*)'                                       KAERI                   '     
          WRITE(*,*)''     
          WRITE(*,*)'                                All Rights Reserved           '      
      ENDIF
!
     

!          
      initial=.FALSE.
!      
      RETURN
      ENDSUBROUTINE print_logo           
!DEC$ENDIF      
