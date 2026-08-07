!
      SUBROUTINE udfn_subchannel_info
!
!     Define user-defined geometry data for subchannel analysis
!   
      USE Zinterface
      USE Zmpi        , ONLY: jperm,ncell_fp,maxmt_fluid
      USE Zzone       , ONLY: ncell_fluid,ncell_fluid_all
      USE Zcore       , ONLY: np,myrank
      USE Znum_cell   , ONLY: i_neigh
      USE Zconst2     , ONLY: hydraulicd
      USE Zcoord2     , ONLY: cell_leng
      USE Zcoord3     , ONLY: volp        
      USE Znormal     , ONLY: xn
      USE Zrv_ncell   , ONLY: chn_nx,chn_ny
      USE Zporous     , ONLY: chn_type,chn_type_tmp,sgap
      USE Zporous     , ONLY: mv_fac
      USE Zporous     , ONLY: n_mv,h_mv,mv_loc
      USE Zio_unit    , ONLY: unit_log
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: mv
      INTEGER :: mv1,mv2,mv3,mv4
      INTEGER :: ix,iy,i0
      INTEGER :: i,ii,j
      
!.....Local allocatable variables
      INTEGER,DIMENSION(:,:,:),ALLOCATABLE :: mv_fac_tmp
!      
!.....Other subchannel model V&V (ex.RPI2x2)      
!
      IF(.not.ALLOCATED(chn_type_tmp)) then
         ALLOCATE(chn_type_tmp(ncell_fluid_all))
         chn_type_tmp=0
         ALLOCATE(chn_type(ncell_fp))
         chn_type=0
         IF(myrank.eq.0)then 
            WRITE(* ,*)'           NO defined subchannel type (RPI2x2) '
            WRITE(97,*)'           NO defined subchannel type (RPI2x2) '
         ENDIF
         RETURN
      ENDIF
      
      ALLOCATE(chn_type(ncell_fp))
      chn_type=0

!........chn_type_tmp NOT allocated check bug correction
      DO ii=1,ncell_fluid
         i=jperm(ii)
         chn_type(ii)=chn_type_tmp(i)
      ENDDO

      IF(np.gt.1) THEN
         CALL communicate_1d_int(chn_type)
         CALL communicate_1d(hydraulicd)
         CALL communicate_1d(volp)
         CALL communicate_2d(cell_leng)
        !CALL communicate_2d(sl)
         CALL communicate_2d(sgap)
      ENDIF
!
!.....Mixing Vane and Spacer Grid 
!
      OPEN(1020, file='mixing_vane.in', status='old', iostat=mv)
      IF(mv.ne.0)then
         IF(myrank.eq.0)then
            WRITE(* ,*)'NO mixing vane information'
            WRITE(unit_log,*)'NO mixing vane information'
         ENDIF
      ELSEIF(mv.eq.0)then
        !ALLOCATE(mv_fac(ns,nn))
         ALLOCATE(mv_fac(maxmt_fluid))
         mv_fac=0
         ix=maxval(chn_nx(:))
         iy=maxval(chn_ny(:))
         i0=ix*iy
         
         ALLOCATE(mv_fac_tmp(ix,iy,4))
         mv_fac_tmp=0
         DO i=1,i0
            READ(1020,*)ii,ix,iy,mv1,mv2,mv3,mv4
            mv_fac_tmp(ix,iy,1)=mv1
            mv_fac_tmp(ix,iy,2)=mv2
            mv_fac_tmp(ix,iy,3)=mv3
            mv_fac_tmp(ix,iy,4)=mv4
         ENDDO
         
         mv_fac=0
         DO i=1,ncell_fluid
            ii=jperm(i)
            ix=chn_nx(ii)
            iy=chn_ny(ii)

            DO j=i_neigh(i),i_neigh(i+1)-1
               IF(xn(j,1).gt. 0.5d0)mv_fac(j)=mv_fac_tmp(ix,iy,1)
               IF(xn(j,1).lt.-0.5d0)mv_fac(j)=mv_fac_tmp(ix,iy,2)
               IF(xn(j,2).gt. 0.5d0)mv_fac(j)=mv_fac_tmp(ix,iy,3)
               IF(xn(j,2).lt.-0.5d0)mv_fac(j)=mv_fac_tmp(ix,iy,4)
            ENDDO
         ENDDO
         DEALLOCATE(mv_fac_tmp)
!
!........Mixing Vane Location
         READ(1020,*) n_mv
         READ(1020,*) h_mv
         ALLOCATE(mv_loc(n_mv))
         DO i=1,n_mv
            READ(1020,*)mv_loc(i)
         ENDDO

      ENDIF
!
      END SUBROUTINE udfn_subchannel_info
