!comdeck stcom                                                          
!                                                                       
!if -def,in32,1                                                         
!      parameter  ( lstcom = 19 )                                       
!if def,in32,1                                                          
      INTEGER lstcom 
      PARAMETER (lstcom=19)   !PARAMETER (lstcom=15) 
!      Note:  if the length of /stcom/ is changed, the value of lstcom  
!             must be changed accordingly (lstcom = real word length)   
!                                                                       
      COMMON/stcom_cupid/ttrip,ptrip,vtrip,tcrit,pcritt,vcrit,tmin,pmin,tmax,  &
      pmax,ntt,npp,nst,nsp,it3bp,it4bp,it5bp,nprpnt,it3p0                 
!  Note:  This common block should agree with the common block of the   
!  same name on the Environmental Library.                              
      REAL (8) ttrip,ptrip,vtrip,tcrit,pcritt,vcrit,tmin,pmin,tmax,pmax 
      INTEGER ntt,npp,nst,nsp,it3bp,it4bp,it5bp,nprpnt,it3p0 
!
!     REAL (8) epsilon                                                    
!     DATA epsilon/1.0d-8/
