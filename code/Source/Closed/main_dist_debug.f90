!
      PROGRAM main
!
      USE Zcore , ONLY: np,myrank
!
      IMPLICIT NONE
!      
      INTEGER:: argc, iargc
      INTEGER cupid__ekd__,cupid__ekd__release,cupid__ekd__debug 
      INTEGER unit_log
      
!
      LOGICAL, SAVE:: initial_expire
!      
      CHARACTER(30):: argv 
!      
      DATA initial_expire /.TRUE./
!        
      COMMON/Zdist/cupid__ekd__,cupid__ekd__release,cupid__ekd__debug         
!
!.....Activate MPI COMM WORLD
!
!DEC$IF defined (mpi_flag)
      INCLUDE 'mpif.h'
      INTEGER status(MPI_STATUS_SIZE),tag,ierr,rc
!
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)
!DEC$ELSE
      np=1
      myrank=0
!DEC$ENDIF
!
!.....Set the distribution version
!
      cupid__ekd__=1
      cupid__ekd__release=0
      IF(cupid__ekd__release.eq.1)THEN
         cupid__ekd__debug=0
      ELSE
         cupid__ekd__debug=1
      ENDIF      
      IF(cupid__ekd__.eq.0)THEN
         cupid__ekd__release=0
         cupid__ekd__debug=0
      ENDIF   
!
!.....Check the serial number and expire date
!      
      IF(cupid__ekd__release.eq.1)THEN 
         argc=iargc()
         IF(argc==0) THEN
            WRITE(*,*)'          !!! This CUPID is for a workbench. !!!'
            WRITE(unit_log,*)'          !!! This CUPID is for a workbench. !!!'
            STOP
         ELSE  
            CALL getarg(argc,argv)
            IF(argv.eq.'QKRDLRRB305755353ekd')THEN 
            ELSE
                WRITE(*,*)'          !!! This CUPID is for a workbench. !!!'
                WRITE(unit_log,*)'          !!! This CUPID is for a workbench. !!!'
                STOP              
            ENDIF
         ENDIF   
         IF(initial_expire) THEN
             CALL check_expire
             initial_expire=.false.
         ENDIF
      ENDIF
!
      IF(myrank.eq.0) THEN 
         WRITE(*,*)
         WRITE(*,*)   
         WRITE(*,*)'                                                             '     
         WRITE(*,*)'                    CCCC  U       U  PPPPPP   IIIIII  DDDDD    '     
         WRITE(*,*)'                  CC      U       U  P     P    II    D    DD  '     
         WRITE(*,*)'                 C        U       U  P     P    II    D      D '     
         WRITE(*,*)'                 C        U       U  PPPPPP     II    D      D '     
         WRITE(*,*)'                  CC       U     U   P          II    D    DD  '     
         WRITE(*,*)'                    CCCC    UUUUU    P        IIIIII  DDDDD    '     
         WRITE(*,*)'                                                             '     
         WRITE(*,*)'                                     v2.20                   '     
         WRITE(*,*)'                                                             '     
         WRITE(*,*)'                                                             '     
         WRITE(*,*)'                               Copyright (c) 2018            '     
         WRITE(*,*)'                                     KAERI                   '     
         WRITE(*,*)'                                                             '     
         WRITE(*,*)'                               All Rights Reserved           '      
         WRITE(*,*)'                                                             '     
         WRITE(*,*)'                                                             '
      ENDIF
!
!.....CUPID v2.20
!         
      CALL cupid_main
!
!DEC$IF defined (mpi_flag)
      CALL MPI_FINALIZE(rc)
!DEC$ENDIF
!
      IF(cupid__ekd__.eq.1.and.myrank.eq.0)PAUSE !distribution
      STOP
      END PROGRAM main
!
!------------------------------------------------------------------------------
!
      SUBROUTINE check_expire
!
      CHARACTER(8) :: date
      CHARACTER(10):: time
      CHARACTER(5) :: zone
      INTEGER,DIMENSION(8) :: values
      INTEGER(8) :: expire_date
!
      CALL date_and_time(date,time,zone,values)
      CALL date_and_time(date=date,zone=zone)
      CALL date_and_time(time=time)
      CALL date_and_time(values=values)
      !print '(a,2x,a,2x,a)',date,time,zone
      !print '(8i5)',values
!
      expire_date=values(1)*10000+values(2)*100+values(3)
      
      IF(expire_date.gt.20201231)THEN
         WRITE(*,*)'          !!! License for CUPID 2.0 is expired. !!!'
         WRITE(unit_log,*)'          !!! License for CUPID 2.0 is expired. !!!'
         STOP
      ENDIF
!
      END SUBROUTINE
    