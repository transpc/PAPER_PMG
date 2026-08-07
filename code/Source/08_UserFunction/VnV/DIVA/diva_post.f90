!
      SUBROUTINE diva_post
!
!     This subroutine write calculation results for DIVA V&V problem 
!
      USE Zzone           , ONLY: ncell_fluid,ncell_fluid_all
      USE Zcore           , ONLY: myrank 
      USE Ztimecon        , ONLY: time
      USE Zcoord1         , ONLY: xloc_tmp
      USE Zvector         , ONLY: vl_n
!
      IMPLICIT NONE 
!      
!.....Local variables
      INTEGER :: i,na
      LOGICAL,SAVE :: INITIAL=.TRUE.
!.....Local arrays
      REAL(8),DIMENSION(:),ALLOCATABLE :: dat
!
!.....Open output file for writing
!
      IF (INITIAL.and.myrank.eq.0) THEN
         OPEN(111, file='diva_post.dat')
         INITIAL=.FALSE.
      ENDIF
!     
      na=ncell_fluid_all
      IF(myrank.eq.0) THEN
         ALLOCATE(dat(na))
      ELSE
         ALLOCATE(dat(1))
      ENDIF
      CALL gatherv_r(vl_n(1,1),ncell_fluid,dat,na,0)
!
      IF(time.gt.10.d0)THEN         
         IF(myrank.eq.0)THEN         
            DO i=1,ncell_fluid_all
               IF(xloc_tmp(i,1).lt.0.1d0) WRITE(111,411)xloc_tmp(i,3),dat(i)
            ENDDO
         ENDIF
      ENDIF
      DEALLOCATE(dat)
!
  411 FORMAT(2(e14.7,1x))        
!      
      END SUBROUTINE diva_post





!!
!      SUBROUTINE diva_post(nplot)
!!
!!     This subroutine write calculation results for DIVA V&V problem 
!!
!      USE VOL_DATA                    
!      USE Zparam          , ONLY: ndim
!      USE Zcoord1         , ONLY: xloc      
!      USE Zcoord3         , ONLY: vol      
!      USE Zcore           , ONLY: np,myrank 
!      USE Ztimecon        , ONLY: time
!      USE Zvector         , ONLY: vl_n
!      USE Zzone           , ONLY: ncell_fluid,ncell_fluid_all
!!
!      IMPLICIT NONE 
!!
!!DEC$IF defined (mpi_flag)
!      INCLUDE 'mpif.h'
!!DEC$ENDIF        
!!      
!!      
!      INTEGER i,nplot,m,n,k,num_grid_x,num_grid_z
!!      
!      LOGICAL, SAVE::INITIAL
!!            
!      REAL(8) point_x,point_z,vl_mag,vl_ang
!      REAL(8), ALLOCATABLE::temp(:),tempall(:),alphal_all(:)
!      REAL(8), ALLOCATABLE::xlocall(:,:),vol_all(:),vlp(:,:),sum_vol(:,:),sum_vlx(:,:),sum_vlz(:,:),sum_ml(:,:)       
!!
!      DATA INITIAL /.TRUE./
!      DATA num_grid_x,num_grid_z/19,7/
!!     
!      ALLOCATE(sum_vol(num_grid_x,num_grid_z),sum_vlx(num_grid_x,num_grid_z),sum_vlz(num_grid_x,num_grid_z),sum_ml(num_grid_x,num_grid_z))
!      ALLOCATE(tempall(ncell_fluid_all),vol_all(ncell_fluid_all),alphal_all(ncell_fluid_all))
!      ALLOCATE(temp(ncell_fluid))     
!      ALLOCATE(xlocall(ndim,ncell_fluid_all),vlp(ndim,ncell_fluid_all))   
!!
!      vlp(:,:)=0.0d0
!      xlocall(:,:)=0.0d0
!      vol_all(:)=0.0d0 
!      sum_vol(:,:)=0.0d0
!      sum_vlx(:,:)=0.0d0
!      sum_vlz(:,:)=0.0d0
!      sum_ml(:,:)=0.0d0
!      alphal_all(:)=0.0d0
!!      
!      IF(np.gt.1)THEN
!         DO k=1,ndim
!            temp(:)=xloc(k,:)
!            CALL allgatherv_r(temp,tempall,ncell_fluid,ncell_fluid_all,0)
!            xlocall(k,:)=tempall(:)
!         ENDDO
!         DO k=1,ndim
!            temp(:)=vl_n(k,:)
!            CALL allgatherv_r(temp,tempall,ncell_fluid,ncell_fluid_all,0)
!            vlp(k,:)=tempall(:)
!         ENDDO
!         CALL allgatherv_r(vol,vol_all,ncell_fluid,ncell_fluid_all,0)
!         CALL allgatherv_r(cell%alphal,alphal_all,ncell_fluid,ncell_fluid_all,0)
!      ELSE
!         vlp(:,:)=vl_n(:,:)
!         xlocall(:,:)=xloc(:,:)
!         vol_all(:)=vol(:)
!         alphal_all(:)=cell%alphal(:)
!      ENDIF        
!!
!!.....Open output file for writing
!!
!      IF (INITIAL.and.myrank.eq.0) THEN
!         OPEN(111, file='diva_vector.dat')                                  
!         INITIAL=.FALSE.
!      ENDIF
!!
!      DO i=1,ncell_fluid_all
!         IF(xlocall(2,i).ge.0.0d0.and.xlocall(2,i).lt.0.2d0.and.alphal_all(i).ge.0.05)THEN         ! Define the calculating area usign the depth of channel and void fraction
!            point_z=0.28d0                                                                          ! z coordinate of 1st measuring grid
!            DO m=1,num_grid_z
!               point_x=0.54d0                                                                       ! x coordinate of 1st measuring grid
!               IF(xlocall(3,i).ge.point_z.and.xlocall(3,i).lt.point_z+0.03d0)THEN                  ! Measuring grid size of z axis is 0.03m
!                  DO n=1,num_grid_x
!                     IF(xlocall(1,i).ge.point_x.and.xlocall(1,i).lt.point_x+0.03d0)THEN            ! Measuring grid size of x axis is 0.03m
!                        sum_vlx(n,m)=sum_vlx(n,m)+vlp(1,i)*vol_all(i)*alphal_all(i)
!                        sum_vlz(n,m)=sum_vlz(n,m)+vlp(3,i)*vol_all(i)*alphal_all(i)
!                        sum_ml(n,m)=sum_ml(n,m)+alphal_all(i)*vol_all(i)
!                        sum_vol(n,m)=sum_vol(n,m)+vol_all(i)*alphal_all(i)
!                     ENDIF
!                     point_x=point_x+0.03d0
!                  ENDDO         
!               ENDIF   
!               point_z=point_z+0.03d0
!            ENDDO           
!         ENDIF
!      ENDDO
!
!!    
!!.....Write interesting parameters
!!           
!      IF(myrank.eq.0)THEN 
!         WRITE(111,411) time
!         point_z=0.28d0+0.015d0           
!         DO m=1,num_grid_z
!            point_x=0.54d0+0.015d0 
!            DO n=1,num_grid_x
!               IF(sum_vol(n,m).eq.0)THEN
!                  vl_mag=0.0d0
!               ELSE
!                  vl_mag=dSQRT((sum_vlx(n,m)/sum_vol(n,m))**2.d0+(sum_vlz(n,m)/sum_vol(n,m))**2.d0)  ! Volume averaged vector magnitude
!               ENDIF
!!               
!               IF(vl_mag.le.0.1) vl_mag=0.0d0    
!!               
!               IF(sum_vlz(n,m).eq.0.or.sum_vlx(n,m).eq.0)THEN
!                  !vl_ang=0.0d0
!               ELSE
!                  vl_ang=DATAN2(sum_vlz(n,m),sum_vlx(n,m))                                          
!               ENDIF                         
!               IF(sum_ml(n,m)/9.0d-4.le.0.2d-3)sum_ml(n,m)=0.0d0                                    !Cut-out thin film compared to Exp.          
!               WRITE(111,412) point_x,point_z,vl_ang,vl_mag,sum_ml(n,m)/9.0d-4
!               point_x=point_x+0.03d0 
!            ENDDO
!            point_z=point_z+0.03d0
!         ENDDO
!      ENDIF 
!!
!  411 FORMAT('time=',1(e14.7,1x))    
!  412 FORMAT(5(e14.7,1x))    
!!      
!      DEALLOCATE(sum_vol,sum_vlx,sum_vlz,sum_ml)
!      DEALLOCATE(tempall,vol_all,alphal_all)
!      DEALLOCATE(temp)
!      DEALLOCATE(xlocall,vlp)
!!  
!      RETURN 
!      END SUBROUTINE diva_post
!
!!
!      SUBROUTINE diva_post_zero(nplot)
!!
!!     This subroutine write calculation results for DIVA V&V problem when the liquid inlet flow is zero
!!
!      USE VOL_DATA                    
!      USE Zparam          , ONLY: ndim
!      USE Zcoord1         , ONLY: xloc      
!      USE Zcoord3         , ONLY: vol      
!      USE Zcore           , ONLY: np,myrank 
!      USE Ztimecon        , ONLY: time
!      USE Zvector         , ONLY: vl_n
!      USE Zzone           , ONLY: ncell_fluid,ncell_fluid_all
!!
!      IMPLICIT NONE 
!!
!!DEC$IF defined (mpi_flag)
!      INCLUDE 'mpif.h'
!!DEC$ENDIF        
!!      
!!      
!      INTEGER i,nplot,m,n,k,num_grid_x,num_grid_z
!      REAL(8) point_x,point_z,vl_mag,vl_ang
!!      
!      LOGICAL, SAVE::INITIAL
!!
!      DATA INITIAL /.TRUE./
!      DATA num_grid_x,num_grid_z/12,20/
!!            
!      REAL(8), ALLOCATABLE::temp(:),tempall(:),alphal_all(:)
!      REAL(8), ALLOCATABLE::xlocall(:,:),vol_all(:),vlp(:,:),sum_vol(:,:),sum_vlx(:,:),sum_vlz(:,:),sum_ml(:,:)       
!     
!      ALLOCATE(sum_vol(num_grid_x,num_grid_z),sum_vlx(num_grid_x,num_grid_z),sum_vlz(num_grid_x,num_grid_z),sum_ml(num_grid_x,num_grid_z))
!      ALLOCATE(tempall(ncell_fluid_all),vol_all(ncell_fluid_all),alphal_all(ncell_fluid_all))
!      ALLOCATE(temp(ncell_fluid))     
!      ALLOCATE(xlocall(ndim,ncell_fluid_all),vlp(ndim,ncell_fluid_all))   
!!
!      vlp(:,:)=0.0d0
!      xlocall(:,:)=0.0d0
!      vol_all(:)=0.0d0 
!      sum_vol(:,:)=0.0d0
!      sum_vlx(:,:)=0.0d0
!      sum_vlz(:,:)=0.0d0
!      sum_ml(:,:)=0.0d0
!      alphal_all(:)=0.0d0
!      IF(np.gt.1)THEN
!         DO k=1,ndim
!            temp(:)=xloc(k,:)
!            CALL allgatherv_r(temp,tempall,ncell_fluid,ncell_fluid_all,0)
!            xlocall(k,:)=tempall(:)
!         ENDDO
!         DO k=1,ndim
!            temp(:)=vl_n(k,:)
!            CALL allgatherv_r(temp,tempall,ncell_fluid,ncell_fluid_all,0)
!            vlp(k,:)=tempall(:)
!         ENDDO
!         CALL allgatherv_r(vol,vol_all,ncell_fluid,ncell_fluid_all,0)
!         CALL allgatherv_r(cell%alphal,alphal_all,ncell_fluid,ncell_fluid_all,0)
!      ELSE
!         vlp(:,:)=vl_n(:,:)
!         xlocall(:,:)=xloc(:,:)
!         vol_all(:)=vol(:)
!         alphal_all(:)=cell%alphal(:)
!      ENDIF        
!!
!!.....Open output file for writing
!!
!      IF (INITIAL.and.myrank.eq.0) THEN
!         OPEN(111, file='diva_vector_zero.dat')                                  
!         INITIAL=.FALSE.
!      ENDIF
!!
!      DO i=1,ncell_fluid_all
!         IF(xlocall(2,i).ge.0.0d0.and.xlocall(2,i).lt.0.2d0.and.alphal_all(i).ge.0.05)THEN
!            point_z=0.20d0          
!            DO m=1,num_grid_z
!               point_x=0.60d0 
!               IF(xlocall(3,i).ge.point_z.and.xlocall(3,i).lt.point_z+0.02d0)THEN                   ! Measuring grid size of z axis is 0.02m
!                  DO n=1,num_grid_x
!                     IF(xlocall(1,i).ge.point_x.and.xlocall(1,i).lt.point_x+0.02d0)THEN             ! Measuring grid size of x axis is 0.02m
!                        sum_vlx(n,m)=sum_vlx(n,m)+vlp(1,i)*vol_all(i)*alphal_all(i)
!                        sum_vlz(n,m)=sum_vlz(n,m)+vlp(3,i)*vol_all(i)*alphal_all(i)
!                        sum_ml(n,m)=sum_ml(n,m)+alphal_all(i)*vol_all(i)
!                        sum_vol(n,m)=sum_vol(n,m)+vol_all(i)*alphal_all(i)
!                     ENDIF
!                     point_x=point_x+0.02d0
!                  ENDDO         
!               ENDIF   
!               point_z=point_z+0.02d0
!            ENDDO           
!         ENDIF
!      ENDDO
!!    
!!.....Write interesting parameters
!!           
!      IF(myrank.eq.0)THEN 
!         WRITE(111,411) time
!         point_z=0.20d0+0.01d0           
!         DO m=1,num_grid_z
!            point_x=0.60d0+0.01d0 
!            DO n=1,num_grid_x
!               IF(sum_vol(n,m).eq.0)THEN
!                  vl_mag=0.0d0
!               ELSE
!                  vl_mag=dSQRT((sum_vlx(n,m)/sum_vol(n,m))**2.d0+(sum_vlz(n,m)/sum_vol(n,m))**2.d0)
!               ENDIF
!!               
!               IF(vl_mag.le.0.1) vl_mag=0.0d0    
!!               
!               IF(sum_vlz(n,m).eq.0.or.sum_vlx(n,m).eq.0)THEN
!                  !vl_ang=0.0d0
!               ELSE
!                  vl_ang=DATAN2(sum_vlz(n,m),sum_vlx(n,m))
!               ENDIF                           
!               IF(sum_ml(n,m)/4.0d-4.le.0.2d-3)sum_ml(n,m)=0.0d0                                     !Cut-out thin film compared to Exp.
!               WRITE(111,412) point_x,point_z,vl_ang,vl_mag,sum_ml(n,m)/4.0d-4
!               point_x=point_x+0.02d0 
!            ENDDO
!            point_z=point_z+0.02d0
!         ENDDO
!      ENDIF 
!!
!  411 FORMAT('time=',1(e14.7,1x))    
!  412 FORMAT(5(e14.7,1x))    
!!      
!      DEALLOCATE(sum_vol,sum_vlx,sum_vlz,sum_ml)
!      DEALLOCATE(tempall,vol_all,alphal_all)
!      DEALLOCATE(temp)
!      DEALLOCATE(xlocall,vlp)
!!  
!      RETURN 
!      END SUBROUTINE diva_post_zero
