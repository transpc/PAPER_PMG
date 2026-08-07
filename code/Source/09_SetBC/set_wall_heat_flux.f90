!
      SUBROUTINE wall_heat_flux
!
!.....This routine defines the wall heat flux
!
      USE Zconst2  , ONLY: stime_hup,stime_hflat
      USE Zqvol    , ONLY: qwall_solid ,qwall_origin
      USE Ztimecon , ONLY: time
!
      IMPLICIT NONE
!
      LOGICAL, SAVE::initial
!
      DATA INITIAL /.TRUE./ 
!      
      IF(initial)THEN
         initial=.FALSE.
         qwall_origin=qwall_solid(5)
      ENDIF
!
!.....Linearly increase the heat flux in time
!
      IF(time.ge.stime_hup.and.time.lt.stime_hflat)THEN
         qwall_solid(5)=qwall_origin*(time-stime_hup)/(stime_hflat-stime_hup)
      ELSEIF(time.lt.stime_hup)THEN
         qwall_solid(5)=0.d0
      ELSE
         qwall_solid(5)=qwall_origin
      ENDIF
!
      END SUBROUTINE wall_heat_flux
!     
!-----------------------------------------------------------------------------------------      
!
      SUBROUTINE udfn_hflux_bc_profile(icell,h_profile)
!
!     This routine defines the wall heat flux disribution
!
      USE Zparam          , ONLY: pi 
      USE Zcoord1         , ONLY: xloc_c        
      USE Zio_unit        , ONLY: unit_log
!
      IMPLICIT NONE
!
      INTEGER icell,profile_direc
!
      LOGICAL, SAVE::initial
!
      REAL(8) x,h_profile
      REAL(8) loc_min,loc_max,length_tot
!      
      DATA INITIAL /.TRUE./       
      DATA profile_direc /3/                                ! select x,y,z dimension 
      DATA length_tot /1.0d0/                               ! totla length to normalize 
!  
      loc_min=-pi*0.5d0
      loc_max=pi*0.5d0
!                                                     
      x=xloc_c(icell,profile_direc)/length_tot                 ! Normalize length
      x=(xloc_c(icell,profile_direc)-0.5d0)*pi                 ! make x to become -pi/2 <= x < pi/2                                                     
      x=DMAX1(loc_min,DMIN1(loc_max,x))
      h_profile=COS(x)
!      
      !WRITE(123,*)xloc_c(icell,profile_direc),x,h_profile
!             
      IF(h_profile.lt.0.0d0.or.h_profile.gt.1.0d0)THEN
         WRITE(*,*) 'h_profile=',h_profile,'x=',x
         WRITE(unit_log,*) 'h_profile=',h_profile,'x=',x
         STOP
      ENDIF     
!      
      IF(INITIAL)WRITE(*,*)'              Wall Heat Flux Profile Option is ON'
      INITIAL=.FALSE.
!
      END SUBROUTINE udfn_hflux_bc_profile          
!     
!-----------------------------------------------------------------------------------------      
!
      SUBROUTINE udfn_hflux_bc_profile_chw
!
!     This routine defines the wall heat flux disribution
!
      USE Zcore        , ONLY: np,myrank
      USE Znum_cell    , ONLY: istart_nf
      USE Zvec_index   , ONLY: left_nf
      USE Zparam       , ONLY: pi,pio2
      USE Zcoord1      , ONLY: xloc           
      USE Zuserdefined , ONLY: udfl_hflux_bc_profile,hflux_bc_profile_chw
      USE Zio_unit     , ONLY: unit_log
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,icell
      INTEGER :: nf_number,istart,len,i1
      INTEGER :: kill
      LOGICAL,SAVE :: initial=.TRUE.
      INTEGER :: profile_direc=3                     ! select x,y,z dimension 
      REAL(8) :: length_tot=1.d0                     ! totla length to normalize
      REAL(8) :: x,h_profile
      REAL(8) :: loc_min,loc_max
!      
!  
      loc_min=-pio2
      loc_max= pio2
!                                                     
      IF(INITIAL) THEN
         if(myrank.eq.0) WRITE(*,*)'              Wall Heat Flux Profile Option is ON'
         nf_number=7
         istart=istart_nf(1,nf_number)
         len    =istart_nf(2,nf_number)
         IF(len.gt.0) ALLOCATE(hflux_bc_profile_chw(len))
         IF(udfl_hflux_bc_profile) THEN
            kill=0
            DO i=1,len
               i1=istart+i
               icell=left_nf(i1)
               x=xloc(icell,profile_direc)/length_tot                 ! Normalize length
               x=(xloc(icell,profile_direc)-0.5d0)*pi                 ! make x to become -pi/2 <= x < pi/2 
               x=MAX(loc_min,MIN(loc_max,x))
               h_profile=COS(x)
!      
              !WRITE(123,*)xloc(icell,profile_direc),x,h_profile
!             
               IF(h_profile.lt.0.d0.or.h_profile.gt.1.d0)THEN
                  WRITE(*,*) myrank,'h_profile=',h_profile,'x=',x
                  WRITE(unit_log,*) 'h_profile=',h_profile,'x=',x
                  kill=1
                  exit
               ENDIF     
               hflux_bc_profile_chw(i)=h_profile
            ENDDO
            IF(np.gt.1) CALL allreducei_i1(kill)
            IF(kill.gt.0) THEN
               CALL barrier_mpi
               CALL finalize_mpi
               STOP
            ENDIF
         ELSE
            DO i=1,len
               hflux_bc_profile_chw(i)=1.d0
            ENDDO
         ENDIF
         INITIAL=.FALSE.
      ENDIF     
!
      END SUBROUTINE udfn_hflux_bc_profile_chw
!     
!-----------------------------------------------------------------------------------------      
!
      SUBROUTINE udfn_hflux_bc_profile_chw_c
!
!     This routine defines the wall heat flux disribution
!
      USE Zcore        , ONLY: np,myrank
      USE Znum_cell    , ONLY: istartc_nf
      USE Zvec_index_solid , ONLY: left_solid_nf
      USE Zparam       , ONLY: pi,pio2
      USE Zcoord1      , ONLY: xloc_c           
      USE Zuserdefined , ONLY: udfl_hflux_bc_profile,hflux_bc_profile_chw_c
      USE Zio_unit     , ONLY: unit_log
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,icell
      INTEGER :: nf_number,istart,len,i1
      INTEGER :: kill
      LOGICAL,SAVE :: initial=.TRUE.
      INTEGER :: profile_direc=3                     ! select x,y,z dimension 
      REAL(8) :: length_tot=1.d0                     ! totla length to normalize
      REAL(8) :: x,h_profile
      REAL(8) :: loc_min,loc_max
!      
      loc_min=-pio2
      loc_max= pio2
!                                                     
      IF(INITIAL) THEN
         if(myrank.eq.0) WRITE(*,*)'              Wall Heat Flux Profile Option is ON'
         nf_number=3
         istart=istartc_nf(1,nf_number)
         len    =istartc_nf(2,nf_number)
         IF(len.gt.0) ALLOCATE(hflux_bc_profile_chw_c(len))
         IF(udfl_hflux_bc_profile) THEN
            kill=0
            DO i=1,len
               i1=istart+i
               icell=left_solid_nf(i1)
               x=xloc_c(icell,profile_direc)/length_tot                 ! Normalize length
               x=(xloc_c(icell,profile_direc)-0.5d0)*pi                 ! make x to become -pi/2 <= x < pi/2 
               x=MAX(loc_min,MIN(loc_max,x))
               h_profile=COS(x)
!      
               !WRITE(123,*)xloc_c(icell,profile_direc),x,h_profile
!             
               IF(h_profile.lt.0.d0.or.h_profile.gt.1.d0)THEN
                  WRITE(*,*) myrank,'h_profile=',h_profile,'x=',x
                  WRITE(unit_log,*) 'h_profile=',h_profile,'x=',x
                  kill=1
                  exit
               ENDIF     
               hflux_bc_profile_chw_c(i)=h_profile
            ENDDO
            IF(np.gt.1) CALL allreducei_i1(kill)
            IF(kill.gt.0) THEN
               CALL barrier_mpi
               CALL finalize_mpi
               STOP
            ENDIF
         ELSE
            DO i=1,len
               hflux_bc_profile_chw_c(i)=1.d0
            ENDDO
         ENDIF
         INITIAL=.FALSE.
      ENDIF     
!
      END SUBROUTINE udfn_hflux_bc_profile_chw_c
