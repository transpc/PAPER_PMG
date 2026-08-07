      SUBROUTINE newstread_cupid(n,nuse,record,a) 
!define win32dvf                                                        
!define erf                                                             
!define fourbyt                                                         
!define hconden                                                         
!define impnon                                                          
!define in32                                                            
!define newnrc                                                          
!define ploc                                                            
!define sphaccm                                                         
!define unix                                                            
!define noselap                                                         
!define noextvol                                                        
!define noextv20                                                        
!define noextsys                                                        
!define noextjun                                                        
!define noextj20                                                        
!define noparcs                                                         
!define nonpa                                                           
!define nomap                                                           
!define logp                                                            
!deck newstread                                                         
!                                                                       
!  $Id: newstread.ff,v 1.2 2001/04/02 17:41:41 dbarber Exp dbarber $    
!                                                                       
!  newstread reads and initializes steam tables from new                
!  thermodynamic properties file data for light water                   
!                                                                       
!  Cognizant engineer: rwt.                                             
!                                                                       
      USE STM_TBL_cupid  , ONLY: pcrit
!
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL(8) a(*) 
      INTEGER n,nuse 
      CHARACTER(80) record(*) 
!                                                                       
!                                                                       
      INCLUDE 'stcom.h' 
      INCLUDE 'newstcom.h' 
      INCLUDE 'gibbpnt.h' 
!                                                                       
      INTEGER i,ios,ntot 
!                                                                       
!--rewind thermodynamic properties file                                 
!                                                                       
      REWIND n 
!                                                                       
!--get thermodynamic properties file title, and information about the   
!--generating program                                                   
!                                                                       
      READ(n,end=10,err=20,iostat=ios)record(1) 
      READ(n,end=10,err=20,iostat=ios)record(2) 
!                                                                       
!--get triple point and critical point data, minimum and maximum        
!--temperatures and pressures, table statistics, and table pointers     
!                                                                       
      READ(n,end=10,err=20,iostat=ios)ttrip,ptrip,vtrip,tcrit,pcrit,    &
      vcrit,tmin,tmax,pmin,pmax                                         
!                                                                       
      READ(n,end=10,err=20,iostat=ios)ntemp,npres,nsubcrittemp,         &
      nsupcrittemp,nsattemp,ofirstsattemp,olastsattemp,otemptrip1,      &
      otempcrit1,nsubcritpres,nsupcritpres,nsatpres,ofirstsatpres,      &
      olastsatpres,oprestrip2,oprescrit2,ofirstsupcritpres,             &
      olastsupcritpres                                                  
!                                                                       
      READ(n,end=10,err=20,iostat=ios)ptable1,ltable1,ntable1,stable1,  &
      ptable2,ltable2,ntable2,stable2,ptable3,ltable3,ntable3,stable3,  &
      ptable4,ltable4,ntable4,stable4,ptable5,ltable5,ntable5,stable5,  &
      ptable6,ltable6,ntable6,stable6,ptable7,ltable7,ntable7,stable7,  &
      ptable8,ltable8,ntable8,stable8,ptable9,ltable9,ntable9,stable9,  &
      ptable10,ltable10,ntable10,stable10                               
!                                                                       
!--get number of words in steam tables                                  
!                                                                       
      READ(n,end=10,err=20,iostat=ios)ntot 
!                                                                       
!--check number of words in steam tables against number of words        
!--available for steam tables storage                                   
!                                                                       
!      IF(myrank.eq.0)write(8,1101)nuse,ntot 
! 1101 FORMAT  (' from newstread: nuse=',i8,', ntot=',i8) 
!                                                                       
      IF(ntot.gt.nuse) GOTO 30 
      nuse=ntot 
!                                                                       
!--get steam tables                                                     
!                                                                       
      READ(n,end=10,err=20,iostat=ios)(a(i),i=1,ntot) 
      GOTO 50 
!                                                                       
!--premature end of data encountered                                    
!                                                                       
   10 WRITE(*,1001) 
!      IF(myrank.eq.0)write(8,1001) 
      GOTO 40 
!                                                                       
!--error reading steam table data                                       
!                                                                       
   20 WRITE(*,1002)ios 
!      IF(myrank.eq.0)write(8,1002)ios 
      GOTO 40 
!                                                                       
!--insufficient space                                                   
!                                                                       
   30 WRITE(*,1003) 
!      IF(myrank.eq.0)write(8,1003) 
!                                                                       
   40 nuse=-1 
!                                                                       
!--done                                                                 
!                                                                       
   50 RETURN 
!                                                                       
!                                                                       
 1001 FORMAT  ('0***** end of data encountered reading thermodynamic ', &
     &           'property file')                                       
 1002 FORMAT  ('0***** read error encountered reading thermodynamic ',  &
     &           'property file, iostat =',i4)                          
 1003 FORMAT  ('0***** insufficient space furnished for thermodynamic ',&
     &           'property file, increase ntotal in stgh2onew1.')       
!                                                                       
!                                                                       
      END SUBROUTINE newstread_cupid                     
