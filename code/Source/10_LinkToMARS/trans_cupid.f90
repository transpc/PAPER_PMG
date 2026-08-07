!
    SUBROUTINE trans_cupid
!
      USE Zinterface
      USE Zbicg     ,only:pbcgind,pbcgind_max,pbcgsig 
      USE Zconst1   ,only:restart
      USE Zconst2   ,only:i_repeat,dt,dtr
      USE Zcore     ,only:np,myrank
      USE Ztimecon  ,only:dt_opt,time,itim,iter_p,itim_restart
      USE Zzone     ,only:ncell_cond_all
!      
      IMPLICIT NONE
!
!DEC$IF defined (MCC)
      INCLUDE 'c3com.h'     !i_where,dt_super,dt_cobra,flag_cobra,flag_relap,nstep_c
!DEC$ELSEIF defined (MCC_DLL)
      INCLUDE 'c3com.h'     !i_where,dt_super,dt_cobra,flag_cobra,flag_relap,nstep_c
      !dec$ attributes dllexport :: trans_cupid
!DEC$ELSEIF defined (SPACE)
      INCLUDE 'c3com_space.h'
      !dec$ attributes dllexport :: trans_cupid
!DEC$ENDIF
      LOGICAL s_flg(3),repeat,pcnv
      INTEGER nout,restart_out
      INTEGER iofail, jflag, i_mass,nmfail  !Instead of c3trans.h, c3mserr.h, contrller.h
      REAL(8) timet, dtvmx, varer           !Instead of c3trans.h, c3mserr.h, contrller.h
      COMMON/cupid_main_/s_flg,repeat,pcnv,nout,restart_out
      COMMON/trans_cupid_/iofail,jflag,i_mass,nmfail,timet,dtvmx,varer                      
!
      repeat=.FALSE. 
!
      CALL print_logo
!
!DEC$IF defined (SPACE)
      CALL c3com_copy_S2C      
!DEC$ENDIF
!
!.....Prepare for a time step advancement
!
      IF(i_where.eq.1)THEN
!      
         iofail=0 
         jflag=0 
         i_mass=0 
         !timet=c3time_sys    !from mars
         timet=time           !from cupid
         flag_cobra=.true.    !pik-back-ins 
!
!........from cupid_main
         i_repeat=0
         restart_out=0
         nout=0
         IF(restart.eq.0)itim=0 
!         
!........Controlling cupid time constants
!
         CALL time_control_start
!
!.....Prepare to begin a time-step advancement (100 continue)                      
!
      ELSEIF(i_where.eq.2)THEN
!         
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
         dtvmx=0.0d0 
         flag_cobra=.true. 
         IF(jflag.ne.0)then
            flag_cobra=.false. 
         ENDIF   
         IF(.not.flag_cobra)then !back-pik-ins 
            IF(myrank.eq.0)WRITE(*,"(11x,a,l,l,1i3)")'##i_where2,125, R,C,j flags=',flag_relap,flag_cobra,jflag          
         ENDIF  
!
!.....Prepare for a time-step advancement: determine dt
!
      ELSEIF(i_where.eq.3)THEN
!      
         IF(.not.flag_relap.and.flag_cobra)then 
            jflag=3 
            timet=timet-dt_super 
            CALL scalar_reset            
            IF(myrank.eq.0)WRITE(*,"(11x,a,l,l,1i3,1e12.4)")'## i_where3, scalar_reset: R,C,j flags,dt_super=',flag_relap,flag_cobra,jflag,dt_super    
            IF(myrank.eq.0)WRITE(97,"(11x,a,l,l,1i3,1e12.4)")'## i_where3, scalar_reset: R,C,j flags,dt_super=',flag_relap,flag_cobra,jflag,dt_super    
         ENDIF 
!
!........Calculate dt in CUPID
!         
         !atlas-delete CALL set_dt 
         IF(.not.flag_cobra)then !pbcgind.gt.0.and.pbcgind.le.pbcgind_max
            dt=dt*0.5d0
            dtr=1.0d0/dt            
            IF(myrank.eq.0)WRITE(*,"(11x,a,1pe20.10)")'## i_where3, reduce dt=',dt
            IF(myrank.eq.0)WRITE(97,"(11x,a,1pe20.10)")'## i_where3, reduce dt=',dt
         ENDIF  
!
         iter_p=0
         itim_restart=itim          
         CALL set_dt !atlas-ins
         dt_cobra=dt 
!DEC$IF defined (SPACE)
         s3dt_cobra=dt_cobra 
!DEC$ENDIF  
!                                                                       
!.....Set fluid boundary conditions                          
!      
      ELSEIF(i_where.eq.4)THEN
         iofail=0 !back-pik-ins
         dt=dt_super 
         dtr=1.0d0/dt 
         CALL broadcast_i(iofail,1)
         CALL broadcast_r(dt,1)        
         CALL broadcast_r(dtr,1)        
!  
!........Change end times of MARS and CUPID, Check trip
!   
         CALL change_tend_mars_cupid ! see cupid_master.in
!
!........Transfer informations from CUPID to MARS
!           
         CALL cupid2mars   !MCC-jjj
         IF(.not.flag_relap.or..not.flag_cobra)CALL scalar_reset_cupvols !mcc-mpi
!
!........Reflect user condition
!         
         CALL user_def_inp(1)
!
!........Store old value 
!         
         CALL shift_solutions
         CALL shift_solutions_cupvols !mcc-mpi
!
!.....Advance a time step in CUPID
!
      ELSEIF(i_where.eq.5)THEN
!
      99 continue  
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
!         
         IF(.not.pcnv) THEN
            iter_p=iter_p+1
            CALL prn_iter
            GOTO 98
          ENDIF  
!
!........Check scalar errors
!
         IF(dt_opt.ge.2) THEN
            CALL check_scalar(s_flg,repeat)
            !IF(repeat) GOTO 99
         ENDIF                 
!
!........Calculate the 1D/3D interface
!
         CALL v1d3d_cupid                         
!      
!.....Calculation failure in MARS, Reset  !MCC-jjj-NEXT
!          
      ELSEIF(i_where.ge.6.and.i_where.le.8)THEN
!
!.....Air-repeat in MARS, Reset           !MCC-jjj-NEXT
!
      ELSEIF(i_where.eq.10)THEN
!
!.....Now, finish CUPID part and update time
!
      ELSEIF(i_where.eq.9)THEN
!      
         !atlas IF(ncell_cond_all.gt.0) CALL calc_solid  !because it can be failed!     
!
!........Check scalar errors & goto 125 
!
         IF(repeat)THEN
            jflag=1
            GOTO 125 ! => flag_cobra=.false. => reduce dt_cobra => not increase dt_cobra during 15 run
         ENDIF
         IF(pbcgind.gt.0.and.pbcgind.le.pbcgind_max)then ! back-pik-ins
            IF(dt.lt.1.d-8)THEN
               PRINT*,'          *** Pressure iteration failed in CUPID! ***'
               WRITE(*,"(a,1i5,1e12.5)")'          pbcgind,dt=',pbcgind,dt
               PAUSE
               STOP
            ENDIF
            pbcgsig=1
            jflag=2
            GOTO 125 ! => flag_cobra=.false. => reduce dt_cobra => not increase dt_cobra during 15 run
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
         time=timet
!
      ENDIF
!      
!DEC$IF defined (SPACE)
      IF(flag_cobra)then
         jflag_cobra=1
      ELSE
         jflag_cobra=0
      ENDIF
!DEC$ENDIF
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
!DEC$IF defined (MCC) 
!      
      IF(myrank.eq.0)THEN
          WRITE(*,*)''                               
          WRITE(*,*)''                               
          WRITE(*,*)'       CCCC U    U  PPPPP  IIIIII DDDD       // M   M    A    RRRRR  SSSSS'     
          WRITE(*,*)'     CC     U    U  P    P   II   D  DD     //  MM MM   A A   R   R  S    '     
          WRITE(*,*)'     C      U    U  PPPPP    II   D    D   //   M M M  AAAAA  RRRRR  SSSSS'     
          WRITE(*,*)'     CC     U    U  P        II   D  DD   //    M   M  A   A  R R        S'     
          WRITE(*,*)'       CCCC  UUUU   P      IIIIII DDDD   //     M   M  A   A  R  RR  SSSSS'     
          WRITE(*,*)
          WRITE(*,*)'                    M   M    A    SSSSS  TTTTTT EEEEE  RRRRR   '     
          WRITE(*,*)'                    MM MM   A A   S        TT   E      R    R   '     
          WRITE(*,*)'                    M M M  AAAAA  SSSSS    TT   EEEEE  RRRRR   '     
          WRITE(*,*)'                    M   M  A   A      S    TT   E      R R      '     
          WRITE(*,*)'                    M   M  A   A  SSSSS    TT   EEEEE  R  RR   '  
          WRITE(*,*)''     
          WRITE(*,*)'                                      v1.0                    '     
          WRITE(*,*)''     
          WRITE(*,*)'                                Copyright (c) 2019            '     
          WRITE(*,*)'                                       KAERI                   '     
          WRITE(*,*)''     
          WRITE(*,*)'                                All Rights Reserved           '      
      ENDIF
!
!DEC$ELSEIF defined (MCC_DLL) 
!      
      IF(myrank.eq.0)THEN
          WRITE(*,*)''                               
          WRITE(*,*)''                             
          WRITE(*,*)' M   M   A   SSSS TTTTTT EEEE RRRR     //  CCCC U    U  PPPPP  IIIIII DDDD      // M   M    A    RRRR   SSSSS'     
          WRITE(*,*)' MM MM  A A  S      TT   E    R   R   // CC     U    U  P    P   II   D  DD    //  MM MM   A A   R   R  S    '     
          WRITE(*,*)' M M M AAAAA SSSS   TT   EEEE RRRR   //  C      U    U  PPPPP    II   D    D  //   M M M  AAAAA  RRRR   SSSSS'     
          WRITE(*,*)' M   M A   A    S   TT   E    R R   //   CC     U    U  P        II   D  DD  //    M   M  A   A  R R        S'     
          WRITE(*,*)' M   M A   A SSSS   TT   EEEE R  RR//      CCCC  UUUU   P      IIIIII DDDD  //     M   M  A   A  R  RR  SSSSS'     
          WRITE(*,*)
          WRITE(*,*)'               '     
          WRITE(*,*)'                '     
          WRITE(*,*)'               '     
          WRITE(*,*)'                '     
          WRITE(*,*)'               '  
          WRITE(*,*)''     
          WRITE(*,*)'                                      v1.0                    '     
          WRITE(*,*)''     
          WRITE(*,*)'                                Copyright (c) 2019            '     
          WRITE(*,*)'                                       KAERI                  '     
          WRITE(*,*)''     
          WRITE(*,*)'                                All Rights Reserved           '      
      ENDIF      
!                
!          
!DEC$ELSEIF defined (SPACE)
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
          WRITE(*,*)'                                Copyright (c) 2019            '     
          WRITE(*,*)'                                       KAERI                   '     
          WRITE(*,*)''     
          WRITE(*,*)'                                All Rights Reserved           '      
      ENDIF 
!DEC$ENDIF      
!
      initial=.FALSE.
!      
      RETURN
      ENDSUBROUTINE print_logo           
