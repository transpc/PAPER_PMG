!
      SUBROUTINE stread_cupid(n,nuse,record,a) 
! 
!      stread_cupid  - read and initialize steam tables from thermodynamic 
!                properties file DATA 
! 
!      Calling sequence: 
! 
!                CALL  stread_cupid (ip1,ip2,cp3,rp4)
! 
!      Parameters: 
! 
!                ip1 = n      = FORTRAN unit number from which 
!                               thermodynamic properties file DATA is 
!                               read (input) 
! 
!                ip2 = nUSE   = number of words available in rp4 for 
!                               storage of steam tables (input) 
!                             = number of words of rp4 actually needed f
!                               steam tables (output) 
!                             = -1 IF error detected during steam table 
!                                read or initialization (output) 
! 
!                cp3 = record = CHARACTER array into which information 
!                               about the steam tables and generating 
!                               program is placed (output) 
! 
!                rp4 = a      = array into which steam tables are read 
!                               (output) 
! 
!      I/O units: 
! 
!                ip1 (input);  see above 
! 
!                * (default output) 
! 
!      This routine adapted from sth2xi routine written by R. J. Wagner 
!      for light water steam tables 
! 
      USE STM_TBL_cupid  , ONLY: pcrit
!
      IMPLICIT NONE 
! 
      INCLUDE 'stcom.h' 
!
      INTEGER n,nuse
      INTEGER i,ios,ntot 
!
      REAL a(*) 
!
      CHARACTER(80) record(*) 
! 
!.....rewind thermodynamic properties file 
! 
      REWIND n 
! 
!.....get thermodynamic properties file title, and information about the generating program 
! 
      READ(n,end=10,err=20,iostat=ios)record(1) 
      READ(n,end=10,err=20,iostat=ios)record(2) 
! 
!.....get triple point and critical point DATA, minimum and maximum temperatures and pressures, table statistics, and table pointers 
! 
      READ(n,end=10,err=20,iostat=ios) ttrip,ptrip,vtrip,tcrit,pcritt,                                     &
                                       vcrit,tmin,pmin,tmax,pmax,ntt,npp,nst,nsp,it3bp,it4bp,it5bp,nprpnt, &
                                       it3p0
      pcrit=pcritt
! 
!.....get number of words in steam tables 
! 
      READ(n,end=10,err=20,iostat=ios)ntot 
! 
!.....check number of words in steam tables against number of words available for steam tables storage 
! 
      IF(ntot.gt.nuse) GOTO 30 
      nuse=ntot 
! 
!.....get steam tables 
! 
      READ(n,end=10,err=20,iostat=ios) (a(i),i=1,ntot) 
      GOTO 50 
! 
!.....premature END of DATA encountered 
! 
   10 WRITE(*,1001) 
      GOTO 40 
! 
!.....error reading steam table DATA  
! 
   20 WRITE(*,1002)ios 
      GOTO 40 
! 
!.....insufficient space 
!                                                                       
   30 WRITE(*,1003) 
! 
   40 nuse=-1 
! 
!.....done 
! 
   50 RETURN 
! 
 1001 FORMAT  ('0***** END of DATA encountered reading thermodynamic ', &
     &           'property file') 
 1002 FORMAT  ('0***** read error encountered reading thermodynamic ',  &
     &           'property file, iostat =',i4) 
 1003 FORMAT  ('0***** insufficient space furnished for thermodynamic ',&
     &           'property file') 
! 
      END SUBROUTINE stread_cupid 
