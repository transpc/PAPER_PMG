!
      SUBROUTINE cupid_main
!
!     CUPID 1.9 main routine
!
      USE Vol_DATA
      USE Zcore        , ONLY: np,myrank,cupid_alone,cupid_mars
      USE Zconst1      , ONLY: cplmaster
      USE Zconst2      , ONLY: dt,i_repeat
      USE Ztimecon     , ONLY: dt_opt,time,itim_restart,cfl_ratio,t_end,itim,iter_p,repeat_smac,smac
      USE Zzone        , ONLY: ncell_cond_all
      USE Zio_unit     , ONLY: unit_log
      USE unitManager  , ONLY: createUnit,printAllUnits
      USE Zuserdefined, ONLY : MG_solver
!
      IMPLICIT NONE
!
!     INCLUDE '../MASTER/master_c.h'
!
!DEC$IF defined (mpi_flag)
      INCLUDE 'mpif.h'
!DEC$ELSE
      INTEGER clock_count,clock_max,clock_rate
!DEC$ENDIF
      INTEGER itim_s,itim_e,nout,restart_out,wunit
!
      LOGICAL s_flg(3),repeat,pcnv
!
      REAL(8) time_begin,time_end
      REAL(8) :: t1,t2
! 
      cupid_alone=1
      cupid_mars=0
      MG_solver = .true.
      ! MG_solver = .false.
!
      time=0.d0
      itim=0
      nout=0
      restart_out=0
      i_repeat=0
      s_flg=.false.
      repeat_smac=.false.
!
!.....OPEN input and output files
!
      CALL open_files
!
!.....Initialize some user-defined variables after read problem name
!
      nout=-1
      CALL read_flow(nout)     
      nout=0
!
!.....input DATA and initialization
!
      CALL read_grid
      
      CALL vectorize_index
      IF(ncell_cond_all.gt.0) CALL vectorize_index_solid
      CALL vectorize_index_connectivity
      CALL vectorize_geo_var
      CALL vectorize_major_flux
      CALL vectorize_scalar_upwind
!!to enable vector code
!
      CALL read_flow(nout)
!            
      CALL vectorize_deallocate_face
!
!============================start time loop===================================
!
!.....Controlling time constants
!
      CALL time_control_start
!
      itim_s=itim+1
      itim_e=itim+1e8
!
!.....Check Wallclock time
!
!DEC$IF defined (mpi_flag)
      time_begin=MPI_Wtime()
!DEC$ELSE
      CALL SYSTEM_CLOCK(clock_count,clock_rate,clock_max)
      time_begin=REAL(clock_count,kind=8)/REAL(clock_rate,kind=8)
!DEC$ENDIF
!
!.....Time-marching Loop
!
      DO 100 itim=itim_s,itim_e

!
         iter_p=0
!
!........Set variable time step
!
         itim_restart=itim
         IF(cfl_ratio.gt.0.d0) CALL set_dt
         time=time+dt
!
!........Problem specific user defined input
!
         CALL user_def_inp(1)
!
!........shift solutions in time
!
         CALL shift_solutions
!
   99    CONTINUE
!
!........Communicate major parameters for parallel computing
!
         IF(np.gt.1) CALL communicate_allb
!
!........Calculcate physical models
!  
         CALL calc_models
!
   98    CONTINUE
         pcnv=.true.
!         
!........Explicit momentum calculation
!
! test-CPU

       IF(myrank.EQ.0) THEN
       !DEC$IF defined (mpi_flag)
         t1 = MPI_Wtime()
       !DEC$ENDIF
       ENDIF

         CALL calc_momentum

! test-CPU
       IF(myrank.EQ.0) THEN
       !DEC$IF defined (mpi_flag)
         t2 = MPI_Wtime()

!         write(*,*)'test here'
       write(101,*) itim,t2-t1  
       !DEC$ENDIF
       ENDIF

!
!........Coupled Scalar and Pressure calculation

!
       IF(myrank.EQ.0) THEN
       !DEC$IF defined (mpi_flag)
         t1 = MPI_Wtime()
       !DEC$ENDIF
       ENDIF

         CALL calc_scalar(s_flg,pcnv)

! test-CPU
       IF(myrank.EQ.0) THEN
       !DEC$IF defined (mpi_flag)
         t2 = MPI_Wtime()

!         write(*,*) 'test here'
       write(101,*) t2-t1
       !DEC$ENDIF
       ENDIF
!
         IF(smac.eq.3.and.repeat_smac) GOTO 99
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
            IF(repeat) GOTO 99
         ENDIF
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
         IF(time.gt.t_end)exit
!
  100 END DO
!
!.....Check wallclock time
!
!DEC$IF defined (mpi_flag)
      time_end=MPI_Wtime()
!DEC$ELSE
      CALL SYSTEM_CLOCK(clock_count,clock_rate,clock_max)
      time_end=REAL(clock_count,kind=8)/REAL(clock_rate,kind=8)
!DEC$ENDIF
!
!.....Print out the total wallclock time
!
      WRITE(*,*) '          Wallclock Time=',time_end-time_begin
      IF(myrank.eq.0) THEN
         WRITE(unit_log,*) '   '
         WRITE(unit_log,*) '## End of Calculation. CPU Time=',time_end-time_begin       
         WRITE(unit_log,*) '   '  
         WRITE(*,*) '          *** End of Calculation ***'
         CLOSE(unit_log)
         wunit=createUnit('zsuccess')
         wunit=97
         OPEN(wunit,file='zsuccess.out',status='replace')
         WRITE(wunit,*)'Success!'
         CLOSE(wunit)
      ENDIF
!
      IF(cplmaster.gt.0) CALL save_restart_master(nout)
!
      RETURN
      END SUBROUTINE cupid_main

