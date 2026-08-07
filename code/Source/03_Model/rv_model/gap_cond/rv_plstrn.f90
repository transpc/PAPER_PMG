      SUBROUTINE rv_plstrn(tempk,tplask,trupk,h,eplas,block) 
!
!deck plstrn                                                            
!                                                                       
!     calculates zircaloy plastic strain as function of temperature and 
!     stress.                                                           
!     calculates flow blockage when rupture occurs                      
!in32 itab1                                                             
!in32 itab2                                                             
!in32 itab3                                                             
!in32 itab4                                                             
!in32end                                                                
!                                                                       
! ********** input ***************                                      
!    tempk  = cladding temperature (deg k)                              
!    tplask = cladding temp at which plastic strain begins (deg k)      
!    trupk  = cladding rupture temp (deg k)                             
!    h      = cladding heatup rate (deg k/sec)                          
! ********* output ***************                                      
!    eplas  = plastic strain  --  unitless ratio of lengths             
! output only at clad failure                                           
!    block  = unit cell flow blockage (percent)                         
!                                                                       
!      this subroutine is similar to the one in frapt6 ver 11           
!      this subroutine includes the strain and blockage tables from     
!      nureg-0630.                                                      
!      two sets of tables are input as data                             
!      statements, slow ramp and fast ramp correlations.                
!      the slow ramp corresponds to less than or equal to 10 c/sec.     
!      the fast ramp is for greater than or equal to 25 c/sec.          
!      interpolation is used for ramp rates between 10 and 25 c/sec.    
!                                                                       
!if def,in32,1                                                          
!      IMPLICIT real (8)(a-h,o-z) 
      IMPLICIT NONE 
!                                                                       
!      INCLUDE 'ufiles.h' 
!                                                                       
      INTEGER itab1(2,3),itab2(2,3),itab3(2,3),itab4(2,3) 
      REAL(8) tab1(50),tab2(50),tab3(50),tab4(50) 
      REAL(8) eplas
      REAL(8) degkc,tempc,tempk,tplasc,tplask,trupc,trupk,exrs,exrf,h,exr,block,blocs,blocf
      REAL(8) exrs1(1),exrf1(1),blocs1(1),blocf1(1)
      EQUIVALENCE(exrs,exrs1)
      EQUIVALENCE(exrf,exrf1)
      EQUIVALENCE(blocs,blocs1)
      EQUIVALENCE(blocf,blocf1)
      LOGICAL ierr 
!                                                                       
!           10 c/sec                                                    
!    tab1-  slow ramp rupture strain vs rupture temperature (deg c)     
!if -def,in32,1                                                         
!      data itab1/-2,50,1/                                              
!if def,in32,1                                                          
      DATA itab1/0,-2,0,50,0,1/ 
      DATA tab1/600.0d0,0.10d0,625.0d0,0.11d0,650.0d0,0.13d0,675.0d0,   &
      0.20d0,700.0d0,0.45d0,725.0d0,0.67d0,750.0d0,0.82d0,775.0d0,      &
      0.89d0,800.0d0,0.90d0,825.0d0,0.89d0,850.0d0,0.82d0,875.0d0,      &
      0.67d0,900.0d0,0.48d0,925.0d0,0.28d0,950.0d0,0.25d0,975.0d0,      &
      0.28d0,1000.0d0,0.33d0,1025.0d0,0.35d0,1050.0d0,0.33d0,1075.0d0,  &
      0.25d0,1100.0d0,0.14d0,1125.0d0,0.11d0,1150.0d0,0.10d0,1175.0d0,  &
      0.10d0,1200.0d0,0.10d0/                                           
!           25 c/sec                                                    
!    tab2-  fast ramp rupture strain vs rupture temperature (deg c)     
!if -def,in32,1                                                         
!      data itab2/-2,50,1/                                              
!if def,in32,1                                                          
      DATA itab2/0,-2,0,50,0,1/ 
      DATA tab2/600.0d0,0.10d0,625.0d0,0.10d0,650.0d0,0.12d0,675.0d0,   &
      0.15d0,700.0d0,0.20d0,725.0d0,0.28d0,750.0d0,0.38d0,775.0d0,      &
      0.48d0,800.0d0,0.57d0,825.0d0,0.60d0,850.0d0,0.60d0,875.0d0,      &
      0.57d0,900.0d0,0.45d0,925.0d0,0.28d0,950.0d0,0.25d0,975.0d0,      &
      0.28d0,1000.0d0,0.35d0,1025.0d0,0.48d0,1050.0d0,0.77d0,1075.0d0,  &
      0.80d0,1100.0d0,0.77d0,1125.0d0,0.39d0,1150.0d0,0.26d0,1175.0d0,  &
      0.26d0,1200.0d0,0.36d0/                                           
!           10 c/sec                                                    
!    tab3-  slow ramp flow blockage (%) vs rupture temp (deg c)         
!if -def,in32,1                                                         
!      data itab3/-2,50,1/                                              
!if def,in32,1                                                          
      DATA itab3/0,-2,0,50,0,1/ 
      DATA tab3/600.0d0,6.5d0,625.0d0,7.0d0,650.0d0,8.4d0,675.0d0,      &
      13.8d0,700.0d0,33.5d0,725.0d0,52.5d0,750.0d0,65.8d0,775.0d0,      &
      71.0d0,800.0d0,71.5d0,825.0d0,71.0d0,850.0d0,65.8d0,875.0d0,      &
      52.5d0,900.0d0,35.7d0,925.0d0,20.0d0,950.0d0,18.0d0,975.0d0,      &
      20.0d0,1000.0d0,24.1d0,1025.0d0,25.7d0,1050.0d0,24.1d0,1075.0d0,  &
      18.0d0,1100.0d0,9.2d0,1125.0d0,7.0d0,1150.0d0,6.5d0,1175.0d0,     &
      6.5d0,1200.0d0,6.5d0/                                             
!           25 c/sec                                                    
!    tab4-  fast ramp flow blockage (%) vs rupture temp (deg c)         
!if -def,in32,1                                                         
!      data itab4/-2,50,1/                                              
!if def,in32,1                                                          
      DATA itab4/0,-2,0,50,0,1/ 
      DATA tab4/600.0d0,6.5d0,625.0d0,6.5d0,650.0d0,7.5d0,675.0d0,      &
      10.0d0,700.0d0,13.8d0,725.0d0,20.0d0,750.0d0,27.5d0,775.0d0,      &
      35.7d0,800.0d0,43.3d0,825.0d0,46.0d0,850.0d0,46.0d0,875.0d0,      &
      43.3d0,900.0d0,33.5d0,925.0d0,20.0d0,950.0d0,18.0d0,975.0d0,      &
      20.0d0,1000.0d0,25.7d0,1025.0d0,35.7d0,1050.0d0,61.6d0,1075.0d0,  &
      64.5d0,1100.0d0,61.6d0,1125.0d0,28.5d0,1150.0d0,18.3d0,1175.0d0,  &
      18.3d0,1200.0d0,26.2d0/                                           
!                                                                       
!                            conversion from degress c to k             
      DATA degkc/273.15d0/ 
!                            convert temperatures to degree c      
      tempc=tempk-degkc 
      tplasc=tplask-degkc 
      trupc=trupk-degkc 
!
!     initialize

      eplas=0.d0
      
!      print*,tempc,tplasc,trupc
!                                                                       
!  Determine if there is any pre-rupture plastic deformation; if none,  
!  then return.                                                         
      IF(tempc.lt.tplasc) GOTO 20 
!  Interpolate slow ramp rate rupture strain.                           
      call polat_cupid(itab1,tab1,trupc,exrs1,ierr) 
      IF(ierr)exrs=0.0d0 
      call polat_cupid(itab2,tab2,trupc,exrf1,ierr) 
      IF(ierr)exrf=0.0d0 
!  Interpolate rupture strain for given heatup rate.                    
      IF(h.gt.10.0d0.and.h.lt.25.0d0)exr=exrs+(exrf-exrs)*(h-10.0d0)/   &
      15.0d0                                                            
      IF(h.le.10.0d0)exr=exrs 
      IF(h.ge.25.0d0)exr=exrf 
!  Check for clad temp above rupture temp.                              
      IF(tempc.ge.trupc) GOTO 30 
!  Clad temp is above rupture temp, calculate plastic strain.           
!  Plastic strain = 25% of the rupture strain.                          
      eplas=exr * 0.25d0*DEXP(-0.0153d0*(trupc-tempc)) 
   20 block=0.0d0 
!      
      GOTO 999 
!  Rupture indicated.                                                   
   30 eplas=exr 
!  Note- the 0.25 factor is removed at rupture.                         
!  Interpolate slow ramp rate blockage.                                 
      call polat_cupid(itab3,tab3,trupc,blocs1,ierr) 
      IF(ierr)blocs=0.0d0 
!  Interpolate fast ramp rate blockage.                                 
      call polat_cupid(itab4,tab4,trupc,blocf1,ierr) 
      IF(ierr)blocf=0.0d0 
!  Interpolate blockage for given heatup rate.                          
      IF(h.gt.10.0d0.and.h.lt.25.0d0)block=blocs+(blocf-blocs)*(h-      &
      10.0d0)/15.0d0                                                    
      IF(h.le.10.0d0)block=blocs 
      IF(h.ge.25.0d0)block=blocf 
!
!      print*,'222',eplas,block
!      
  999 CONTINUE 
!
!      print*,'lsj'
!      print*,eplas,block  
! 
      RETURN 
      END SUBROUTINE rv_plstrn           
!
!==============================================================================
!
      SUBROUTINE polat_cupid(itble,ftble,arg,val,err) 
!jjj
!                                                                       
!  Performs linear interpolation.                                       
!                                                                       
      IMPLICIT NONE
!                                                                       
!  Performs linear interpolation on data in ftble using arg as search   
!  arguments.  Interpolated values returned in val.                     
!  itble has three words:                                               
!     itble(1) - the number of items per set.  If negative, error if    
!                arg is not bounded by independent variable of ftble.   
!                If positive, arg not bounded is not an error and       
!                returned values are appropriate end values in ftble.   
!     itble(2) - the total number of items.                             
!     itble(3) - the last used subscript.                               
!  ftble has table values.  The first word of each set in the table is  
!  the search variable.                                                 
!                                                                       
      INTEGER itble(2,*) 
      REAL(8) ftble(*),arg,val(*) 
      LOGICAL err 
      INTEGER i,ie,in,ig,k 
      REAL(8) r 
!                                                                       
      err=.false. 
      in=abs(itble(2,1)) 
      ie=itble(2,2) 
      i=itble(2,3) 
      IF(ie.eq.in) GOTO 20 
   10 IF(arg-ftble(i))13,25,11 
   13 IF(i.eq.1) GOTO 20 
      i=i-in 
      GOTO 10 
!
   11 ig=i+in 
   14 IF(ig.gt.ie) GOTO 20 
      IF(ftble(ig)-arg)17,24,12 
   17 i=ig 
      GOTO 11 
!
   12 r=(arg-ftble(i))/(ftble(ig)-ftble(i)) 
      DO k=2,in 
         val(k-1)=r*(ftble(ig+k-1)-ftble(i+k-1))+ftble(i+k-1) 
      END DO 
      GOTO 1000 
!
   20 IF(itble(2,1).ge.0) GOTO 25 
      err=.true. 
      RETURN 
!
   24 i=ig 
   25 DO k=2,in 
         val(k-1)=ftble(i+k-1) 
      END DO 
 1000 itble(2,3)=i 
! 
      RETURN 
      END SUBROUTINE polat_cupid    
