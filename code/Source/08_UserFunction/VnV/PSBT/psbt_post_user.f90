!
      SUBROUTINE psbt_post
!
!     This subroutine write calculation results for PSBT V&V problem 
!
      USE VOL_DATA        , ONLY: cell
      USE Zzone           , ONLY: ncell_fluid,ncell_fluid_all
      USE Zcore           , ONLY: myrank 
      USE Ztimecon        , ONLY: time
      USE Zcoord1         , ONLY: xloc_tmp
      USE Zcoord3         , ONLY: vol    
!
      IMPLICIT NONE
!      
!.....Local variables
      INTEGER :: i,na
      LOGICAL,SAVE :: INITIAL=.TRUE.
!.....Local arrays
      REAL(8) :: psbt2(13,3)
      REAL(8),DIMENSION(:),ALLOCATABLE :: alphag_all,vol_all,alphag_x_vol_all
!      
      na=ncell_fluid_all 
      IF(myrank.eq.0)THEN
         ALLOCATE(alphag_all(na),vol_all(na))
      ELSE
         ALLOCATE(alphag_all(1),vol_all(1))
      ENDIF
!
      psbt2(:,:)=0.0
      CALL gatherv_r(vol        ,ncell_fluid,vol_all   ,ncell_fluid_all,0)
      CALL gatherv_r(cell%alphag,ncell_fluid,alphag_all,ncell_fluid_all,0)
!
!.....Open output file for writing
!
      IF(INITIAL) THEN
         IF(myrank.eq.0)THEN 
            OPEN(111, file='PSBT_ref.dat')              
            WRITE(111,5021) 
         ENDIF
         INITIAL=.FALSE.
      ENDIF
 5021    FORMAT(4x,'Time(s)         h=0.3m        h=0.4m         h=0.5m         h=0.6m          h=0.7m   &
          h=0.8m         h=0.9m         h=1.0m         h=1.1m         h=1.2m         h=1.3m         h=1.4m         h=1.5m')                           
      IF(myrank.eq.0)THEN 
!     
!........Calculate the area averaged void fraction at 13 different heights
!
         ALLOCATE(alphag_x_vol_all(na))
         DO i=1,ncell_fluid_all
            alphag_x_vol_all(i)=alphag_all(i)*vol_all(i)           
         ENDDO   
         DO i=1,ncell_fluid_all
            IF    (xloc_tmp(i,3).ge.0.29d0.and.xloc_tmp(i,3).lt.0.31d0)THEN
               psbt2(1,1)=psbt2(1,1)+alphag_x_vol_all(i)           
               psbt2(1,2)=psbt2(1,2)+vol_all(i)        
            ELSEIF(xloc_tmp(i,3).ge.0.39d0.and.xloc_tmp(i,3).lt.0.41d0)THEN
               psbt2(2,1)=psbt2(2,1)+alphag_x_vol_all(i)           
               psbt2(2,2)=psbt2(2,2)+vol_all(i)          
            ELSEIF(xloc_tmp(i,3).ge.0.49d0.and.xloc_tmp(i,3).lt.0.51d0)THEN
               psbt2(3,1)=psbt2(3,1)+alphag_x_vol_all(i)           
               psbt2(3,2)=psbt2(3,2)+vol_all(i)          
            ELSEIF(xloc_tmp(i,3).ge.0.59d0.and.xloc_tmp(i,3).lt.0.61d0)THEN
               psbt2(4,1)=psbt2(4,1)+alphag_x_vol_all(i)           
               psbt2(4,2)=psbt2(4,2)+vol_all(i)          
            ELSEIF(xloc_tmp(i,3).ge.0.69d0.and.xloc_tmp(i,3).lt.0.71d0)THEN
               psbt2(5,1)=psbt2(5,1)+alphag_x_vol_all(i)           
               psbt2(5,2)=psbt2(5,2)+vol_all(i)          
            ELSEIF(xloc_tmp(i,3).ge.0.79d0.and.xloc_tmp(i,3).lt.0.81d0)THEN
               psbt2(6,1)=psbt2(6,1)+alphag_x_vol_all(i)           
               psbt2(6,2)=psbt2(6,2)+vol_all(i)          
            ELSEIF(xloc_tmp(i,3).ge.0.89d0.and.xloc_tmp(i,3).lt.0.91d0)THEN
               psbt2(7,1)=psbt2(7,1)+alphag_x_vol_all(i)           
               psbt2(7,2)=psbt2(7,2)+vol_all(i)          
            ELSEIF(xloc_tmp(i,3).ge.0.99d0.and.xloc_tmp(i,3).lt.1.01d0)THEN
               psbt2(8,1)=psbt2(8,1)+alphag_x_vol_all(i)           
               psbt2(8,2)=psbt2(8,2)+vol_all(i)          
            ELSEIF(xloc_tmp(i,3).ge.1.09d0.and.xloc_tmp(i,3).lt.1.11d0)THEN
               psbt2(9,1)=psbt2(9,1)+alphag_x_vol_all(i)           
               psbt2(9,2)=psbt2(9,2)+vol_all(i)          
            ELSEIF(xloc_tmp(i,3).ge.1.19d0.and.xloc_tmp(i,3).lt.1.21d0)THEN
               psbt2(10,1)=psbt2(10,1)+alphag_x_vol_all(i)
               psbt2(10,2)=psbt2(10,2)+vol_all(i)          
            ELSEIF(xloc_tmp(i,3).ge.1.29d0.and.xloc_tmp(i,3).lt.1.31d0)THEN
               psbt2(11,1)=psbt2(11,1)+alphag_x_vol_all(i)
               psbt2(11,2)=psbt2(11,2)+vol_all(i)       
            ELSEIF(xloc_tmp(i,3).ge.1.39d0.and.xloc_tmp(i,3).lt.1.41d0)THEN
               psbt2(12,1)=psbt2(12,1)+alphag_x_vol_all(i)
               psbt2(12,2)=psbt2(12,2)+vol_all(i) 
            ELSEIF(xloc_tmp(i,3).ge.1.49d0.and.xloc_tmp(i,3).lt.1.51d0)THEN
               psbt2(13,1)=psbt2(13,1)+alphag_x_vol_all(i)
               psbt2(13,2)=psbt2(13,2)+vol_all(i)   
            ENDIF
         ENDDO   
         DO i=1,13 
            psbt2(i,3)=psbt2(i,1)/psbt2(i,2)
         ENDDO                          
!
         WRITE(111,411) time, (psbt2(i,3),i=1,13)
         DEALLOCATE(alphag_x_vol_all)
      ENDIF 
!     
      DEALLOCATE(alphag_all,vol_all)
!
  411 FORMAT(14(e14.7,1x))    
!  
      END SUBROUTINE psbt_post
