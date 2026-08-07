!
      SUBROUTINE apr1400_lbloca_out
!
      USE Zmpi            , ONLY: jperm
      USE Zcore           , ONLY: myrank  
      USE Zparam          , ONLY: ndim 
      USE Zrv_ncell       , ONLY: ncell_fluid_core,ncell_fluid_core_all,num_ch,cupid_cell_channel
      USE Ztimecon        , ONLY: time
      USE Zwall_HTC       , ONLY: twall_rv
      USE Zzone           , ONLY: ncell_fluid_all
      USE Zrv_hts_2d      , ONLY: nz0_2d   
      USE Zcoord1         , ONLY: xloc_tmp
!      
      IMPLICIT NONE
!
!....Local variables
      INTEGER i,j,na,nrv1,nrv2,ii
      LOGICAL,SAVE :: initial=.true.
      REAL(8) temp_an
      REAL(8) dat_rv1_rodx
!.....Local arrays
      INTEGER :: cupid_cell_channel_gl(ncell_fluid_core),icell_rod(num_ch)
      REAL(8) :: dat_rv1_rod(num_ch)
!.....Allocatable local arrays
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: cupid_cell_channel_tmp
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: xloc_rod,dat_rv1
!
      na=ncell_fluid_all
      nrv1=ncell_fluid_core
      nrv2=ncell_fluid_core_all
!
      IF(initial) THEN
         IF(myrank.eq.0) OPEN(335,file='VD15_apr1400_lbloca_PCT.dat')  
         IF(myrank.eq.0) THEN
            ALLOCATE(cupid_cell_channel_tmp(nz0_2d,num_ch))
         ELSE
            ALLOCATE(cupid_cell_channel_tmp(1,1))
         ENDIF
         DO i=1,nrv1
            cupid_cell_channel_gl(i)=jperm(cupid_cell_channel(i))
         ENDDO
         CALL gatherv_i(cupid_cell_channel_gl,nrv1,cupid_cell_channel_tmp,nrv2,2)
!........PCT for each channel
         IF(myrank.eq.0) THEN
            ALLOCATE(xloc_rod(num_ch,ndim))
            DO i=1,num_ch
               ii=cupid_cell_channel_tmp(1,i)
               icell_rod(i)=ii
               xloc_rod(i,:)=xloc_tmp(ii,:)
            ENDDO          
            WRITE(335,1000) num_ch,(icell_rod(i),i=1,num_ch)
            WRITE(335,1001) num_ch,(xloc_rod(i,1),i=1,num_ch)
            WRITE(335,1001) num_ch,(xloc_rod(i,2),i=1,num_ch)
            DEALLOCATE(xloc_rod)
         ENDIF    
         initial=.false.
         DEALLOCATE(cupid_cell_channel_tmp)
      ENDIF
!         
      temp_an=time             
!
      IF(myrank.eq.0) THEN
         ALLOCATE(dat_rv1(nz0_2d,num_ch))
      ELSE
         ALLOCATE(dat_rv1(1,1))
      ENDIF
      CALL gatherv_r(twall_rv,nrv1,dat_rv1,nrv2,2)
!
      IF(myrank.eq.0) THEN
         DO i=1,num_ch
            dat_rv1_rodx=dat_rv1(1,i)
            DO j=2,nz0_2d         
               dat_rv1_rodx=MAX(dat_rv1_rodx,dat_rv1(j,i))
            ENDDO
            dat_rv1_rod(i)=dat_rv1_rodx
         ENDDO   
         WRITE(335,1002) temp_an,(dat_rv1_rod(i),i=1,num_ch)
      ENDIF
      DEALLOCATE(dat_rv1)
!
 1000 FORMAT(1x,1i20,1000(2x,10i14))          
 1001 FORMAT(1x,1i20,1000(2x,10e14.5))          
 1002 FORMAT(1x,1e20.10,1000(2x,10e14.5))          
!
      END SUBROUTINE apr1400_lbloca_out   
