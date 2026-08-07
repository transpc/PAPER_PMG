      MODULE Zconst2
!
      IMPLICIT NONE
      SAVE
!
      INTEGER i_repeat,iprn
      REAL(8) dt,dtr,dt_old,stime_hup,stime_hflat,imomwall
      REAL(8) l_old
      REAL(8) stime_vup,stime_vflat
      REAL(8) ggc ! magnitude of gravity vector, grav(:).
      REAL(8),Allocatable::grav(:),gfactor(:),vl_origin(:,:),vg_origin(:,:),vd_origin(:,:)
      REAL(8),Allocatable::hydraulicd(:),sl(:,:)
      INTEGER::hydraulicd_init=0 
      logical lgravity
!
      END MODULE Zconst2