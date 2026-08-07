
      SUBROUTINE rocom_output_cylinder_user(time,stavg,ftavg)
!
!     These subroutine is very special print routine for secific problems
!     Don't USE for another problem
!     Save lateral output for Origin
!
      USE VOL_DATA        , ONLY: cell
      USE Zmpi            , ONLY: ncell_fp
      USE Zzone           , ONLY: ncell_fluid,ncell_fluid_all   
      USE Zcore           , ONLY: myrank
      USE Zparam          , ONLY: ndim
      USE Zcoord1         , ONLY: xloc_tmp
      USE Zcoord2         , ONLY: xloc_xfc_radius_min,xloc_xfc_radius_max
      USE Zcoord3         , ONLY: vol
      USE Zrocom_specific , ONLY: tl_space_avg,tl_space_min
      USE Zvector         , ONLY: vl_n
!
      IMPLICIT NONE 
!
!.....Input
      REAL(8) :: time
      REAL(8) :: stavg,ftavg
!.....Local variables
      INTEGER :: i,level,na
      INTEGER :: Nradius      
!     INTEGER :: time_avg_opt  
      LOGICAL,SAVE :: INITIAL_p=.TRUE.
      REAL(8),SAVE :: tli_sum_num  
      REAL(8) :: radius(2),h1,h2
      REAL(8) :: radiusmax,radiusmin,height,hzero
      REAL(8) :: tl_vol_sum(2),tl_vol_avg(2),vol_sum(2),tl_min(2)
!.....Local arrays
      REAL(8),DIMENSION(:),ALLOCATABLE :: pr,tlall,cboron,volall
      REAL(8),DIMENSION(:),ALLOCATABLE :: radiusmax_all,radiusmin_all
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: vlp
      REAL(8),DIMENSION(:),SAVE,ALLOCATABLE :: tli_sum,pi_sum,cboroni_sum
      REAL(8),DIMENSION(:),SAVE,ALLOCATABLE :: tli_ave,pi_ave,cboroni_ave
      REAL(8),DIMENSION(:,:),SAVE,ALLOCATABLE :: vli_sum,vli_ave
!            
      DATA Nradius,hzero,h1,h2,radius(1),radius(2)/2,1.8d0,-1.1905,-0.3225,0.438,0.499/
!         
      na=ncell_fluid_all
      IF(myrank.eq.0) THEN
         ALLOCATE(vlp(na,ndim))
         ALLOCATE(pr(na),tlall(na),cboron(na),volall(na))
         ALLOCATE(radiusmax_all(na),radiusmin_all(na))
      ELSE
         ALLOCATE(vlp(1,ndim))
         ALLOCATE(pr(1),tlall(1),cboron(1),volall(1))
         ALLOCATE(radiusmax_all(1),radiusmin_all(1))
      ENDIF
!
      CALL gatherv_r_2d(vl_n,ncell_fp,vlp    ,ncell_fluid,na,0)
      CALL gatherv_r(cell%p     ,ncell_fluid,pr    ,na,0)
      CALL gatherv_r(cell%tl    ,ncell_fluid,tlall ,na,0)
      CALL gatherv_r(cell%cboron,ncell_fluid,cboron,na,0)
      CALL gatherv_r(vol        ,ncell_fluid,volall,na,0)
!            
      CALL gatherv_r(xloc_xfc_radius_max,ncell_fluid,radiusmax_all,na,0)
      CALL gatherv_r(xloc_xfc_radius_min,ncell_fluid,radiusmin_all,na,0)
!
      IF(myrank.eq.0) THEN
!      
!........DO just one time     
!
         IF(ndim.ne.3) RETURN 
!      
         IF(INITIAL_p)THEN
            INITIAL_p =.FALSE.
            ALLOCATE(tli_sum(na))
            ALLOCATE(tli_ave(na))
            ALLOCATE(vli_sum(na,ndim),pi_sum(na),cboroni_sum(na))
            ALLOCATE(vli_ave(na,ndim),pi_ave(na),cboroni_ave(na))
            tli_sum_num   =0.d0
            tli_sum(:)    =0.d0
            cboroni_sum(:)=0.d0
            pi_sum(:)     =0.d0
            vli_sum(:,:)  =0.d0         
         ENDIF
!      
!........Initialize the variable for time-averaging
!
         IF(time.ge.stavg.and.time.le.ftavg)THEN
!           time_avg_opt=1
            tli_sum_num=tli_sum_num+1.d0
            DO i=1,ncell_fluid_all
               tli_sum(i)=tli_sum(i)+tlall(i)-273.15d0
               vli_sum(i,:)=vli_sum(i,:)+vlp(i,:)
               pi_sum(i)=pi_sum(i)+pr(i)/1.0d6
               cboroni_sum(i)=cboroni_sum(i)+cboron(i)
               tli_ave(i)=tli_sum(i)/tli_sum_num
               pi_ave(i)=pi_sum(i)/tli_sum_num
              cboroni_ave(i)=cboroni_sum(i)/tli_sum_num
              vli_ave(i,:)=vli_sum(i,:)/tli_sum_num
            ENDDO
         ELSE
!           time_avg_opt=0   
!           tli_sum_num=1.d0
            DO i=1,ncell_fluid_all
               tli_sum(i)=tlall(i)-273.15d0
               vli_sum(i,:)=vlp(i,:)
               pi_sum(i)=pr(i)/1.0d6
               cboroni_sum(i)=cboron(i)
!              tli_ave(i)=tli_sum(i)/tli_sum_num
!              pi_ave(i)=pi_sum(i)/tli_sum_num
!              cboroni_ave(i)=cboroni_sum(i)/tli_sum_num
!              vli_ave(i,:)=vli_sum(i,:)/tli_sum_num
               tli_ave(i)=tli_sum(i)
               pi_ave(i)=pi_sum(i)
               cboroni_ave(i)=cboroni_sum(i)
               vli_ave(i,:)=vli_sum(i,:)
            ENDDO   
         ENDIF     
!      
!........Initialize the variables for space-averaging
!
         tl_vol_sum(:)=0.d0 
         vol_sum(:)   =0.d0
         tl_min(:)    =10000.d0
!
!........For the 2 radius
!
         DO level=1,Nradius 
!
!...........For the i-th cell
!         
            DO i=1,ncell_fluid_all
!            
!..............Check whether assigned radius is between maximum and minimum radii of current cell
!
               radiusmax=radiusmax_all(i)
               radiusmin=radiusmin_all(i)
               IF(radius(level).lt.radiusmin .or. radius(level).gt.radiusmax) CYCLE
!
!..............Check the height
!..............1.82-0.03 becaUSE of rectangular leg, if circular leg, it should be 1.8  
!
               height=xloc_tmp(i,ndim)-hzero 
               IF(height.lt.h1.or.height.gt.h2) CYCLE
!
!..............Sum of tl for the 2-d plane                 
!
               tl_vol_sum(level)=tl_vol_sum(level)+(tlall(i)-273.15d0)*volall(i)
               vol_sum(level)   =vol_sum(level)+volall(i)
               IF(tl_min(level).gt.(tlall(i)-273.15d0))tl_min(level)=tlall(i)-273.15d0            
            ENDDO
         ENDDO
!
!........Averaging the tl for the 2-d plane
!        
         DO i=1,nradius
            IF(vol_sum(i).gt.0.d0)THEN
               tl_vol_avg(i)=tl_vol_sum(i)/vol_sum(i)
            ELSE
               PRINT*,'vol_sum(i) is zero for',vol_sum(i),'level in SUBROUTINE rocom_output_cylindeer'           
            ENDIF           
         ENDDO    
         tl_space_avg(2)=tl_vol_avg(1)     
         tl_space_avg(3)=tl_vol_avg(2) 
         tl_space_min(2)=tl_min(1)    
         tl_space_min(3)=tl_min(2)   
!
      ENDIF
!   
      DEALLOCATE(vlp)
      DEALLOCATE(pr)
      DEALLOCATE(tlall)
      DEALLOCATE(volall)
      DEALLOCATE(cboron)  
      DEALLOCATE(radiusmax_all,radiusmin_all)
!    
      END SUBROUTINE rocom_output_cylinder_user
!      
!------------------------------------------------------------------------------------------------
!
      SUBROUTINE rocom_output_plane_user(time,stavg,ftavg)
!
!     Save lateral output for Origin
!
      USE VOL_DATA        , ONLY: cell                  
      USE Zmpi            , ONLY: ncell_fp
      USE Zzone           , ONLY: ncell_fluid,ncell_fluid_all  
      USE Zcore           , ONLY: myrank      
      USE Zparam          , ONLY: ndim
      USE Zcoord1         , ONLY: xloc_tmp
      USE Zcoord2         , ONLY: xloc_xfc_min,xloc_xfc_max
      USE Zcoord3         , ONLY: vol
      USE Zrocom_specific , ONLY: tl_space_avg,tl_space_min
      USE Zvector         , ONLY: vl_n
!
      IMPLICIT NONE 
!      
!.....Input
      REAL(8) time
!.....Local variables
      INTEGER :: i,na,wopt,loop,level     
!     INTEGER :: time_avg_opt
      INTEGER :: Nheight
      LOGICAL,SAVE::INITIAL_p =.TRUE.
      REAL(8),SAVE::tli_sum_num
      REAL(8) :: hzloc,lzloc         
      REAL(8) :: stavg,ftavg         
      REAL(8) :: alphag_sum,ia_sum,vg_sum,dsm_sum
      REAL(8) :: rdistance    
!.....Local arrays
      REAL(8) :: height(3),r(3*2)
      REAL(8) :: tl_vol_sum(3),vol_sum(3),tl_min(3),tl_vol_avg(3)
      REAL(8) :: hzloc_all(ncell_fluid_all),lzloc_all(ncell_fluid_all)
!.....Local allocatable arrays
      REAL(8),DIMENSION(:),ALLOCATABLE :: pr,tlall,cboron,volall
      REAL(8),DIMENSION(:),SAVE,ALLOCATABLE :: tli_sum,pi_sum,cboroni_sum
      REAL(8),DIMENSION(:),SAVE,ALLOCATABLE :: tli_ave,pi_ave,cboroni_ave
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: vlp
      REAL(8),DIMENSION(:,:),SAVE,ALLOCATABLE :: vli_sum,vli_ave
!            
      DATA Nheight,height/3,0.25,0.376,1.755/    !drum, corein,leg
      DATA r/0,0.437,0,0.437,0,1.0/
!      
      na=ncell_fluid_all
      IF(myrank.eq.0) THEN
         ALLOCATE(vlp(na,ndim))
         ALLOCATE(pr(na),tlall(na),cboron(na),volall(na))
      ELSE
         ALLOCATE(vlp(1,ndim))
         ALLOCATE(pr(1),tlall(1),cboron(1),volall(1))
      ENDIF
!            
      CALL gatherv_r_2d(vl_n,ncell_fp,vlp    ,ncell_fluid,na,0)
      CALL gatherv_r(cell%p     ,ncell_fluid,pr    ,na,0)
      CALL gatherv_r(cell%tl    ,ncell_fluid,tlall ,na,0)
      CALL gatherv_r(cell%cboron,ncell_fluid,cboron,na,0)
      CALL gatherv_r(vol        ,ncell_fluid,volall,na,0)
!
      IF(ndim.eq.2) THEN
         CALL gatherv_r(xloc_xfc_max(1,2),ncell_fluid,hzloc_all,na,0)
         CALL gatherv_r(xloc_xfc_min(1,2),ncell_fluid,lzloc_all,na,0)
      ELSE
         CALL gatherv_r(xloc_xfc_max(1,3),ncell_fluid,hzloc_all,na,0)
         CALL gatherv_r(xloc_xfc_min(1,3),ncell_fluid,lzloc_all,na,0)
      ENDIF
!        
      IF(myrank.eq.0) THEN
!      
         IF(ndim.ne.3) RETURN 
!
!........DO just one time      
!
         IF(INITIAL_p)THEN
            INITIAL_p=.FALSE.
            ALLOCATE(tli_sum(na))
            ALLOCATE(tli_ave(na))
            ALLOCATE(vli_sum(na,ndim),pi_sum(na),cboroni_sum(na))
            ALLOCATE(vli_ave(na,ndim),pi_ave(na),cboroni_ave(na))
            tli_sum_num   =0.d0
            tli_sum(:)    =0.d0
            cboroni_sum(:)=0.d0
            pi_sum(:)     =0.d0
            vli_sum(:,:)  =0.d0 
         ENDIF
!
!........Initialize the variable for time-averaging
!
         IF(time.ge.stavg .and. time.le.ftavg)THEN
!           time_avg_opt=1
            tli_sum_num=tli_sum_num+1.d0
            DO i=1,ncell_fluid_all
               tli_sum(i)=tli_sum(i)+tlall(i)-273.15d0
               vli_sum(i,:)=vli_sum(i,:)+vlp(i,:)
               pi_sum(i)=pi_sum(i)+pr(i)/1.0d6
               cboroni_sum(i)=cboroni_sum(i)++cboron(i)
               tli_ave(i)=tli_sum(i)/tli_sum_num
               pi_ave(i)=pi_sum(i)/tli_sum_num
               cboroni_ave(i)=cboroni_sum(i)/tli_sum_num
               vli_ave(i,:)=vli_sum(i,:)/tli_sum_num
            ENDDO
         ELSE
!           time_avg_opt=0   
!           tli_sum_num=1.d0
            DO i=1,ncell_fluid_all
               tli_sum(i)=tlall(i)-273.15d0
               vli_sum(i,:)=vlp(i,:)
               pi_sum(i)=pr(i)/1.0d6
               cboroni_sum(i)=cboron(i)
!              tli_ave(i)=tli_sum(i)/tli_sum_num
!              pi_ave(i)=pi_sum(i)/tli_sum_num
!              cboroni_ave(i)=cboroni_sum(i)/tli_sum_num
!              vli_ave(i,:)=vli_sum(i,:)/tli_sum_num
               tli_ave(i)=tli_sum(i)
               pi_ave(i)=pi_sum(i)
               cboroni_ave(i)=cboroni_sum(i)
               vli_ave(i,:)=vli_sum(i,:)
            ENDDO   
         ENDIF    
!
!........Initialize the variable for space-averaging
!
         tl_vol_sum(:)=0.d0
         vol_sum(:)   =0.d0
         tl_min(:)    =10000.d0  
!
!........For the three level
!
         DO level=1,Nheight 
            loop = 0
            alphag_sum=0.d0
            ia_sum    =0.d0
            vg_sum    =0.d0
            dsm_sum   =0.d0
!
!...........For the i-th cell
!      
            DO i = 1, ncell_fluid_all
!
!..............Check the height
!           
               hzloc=hzloc_all(i) 
               lzloc=lzloc_all(i) 
               IF(hzloc.le.height(level).or.lzloc.gt.height(level)) CYCLE
!
!..............Check the radius
!
               wopt=0
               rdistance=ABS(xloc_tmp(i,1))
               IF(ndim .eq. 3)rdistance=SQRT(xloc_tmp(i,1)**2+xloc_tmp(i,2)**2)
               IF(rdistance.ge.r((level-1)*2+1).and.rdistance.le.r((level-1)*2+2)) wopt=1
!
!..............Extract the DATA to be PRINTed IF wopt ==1
!
               IF(wopt.ne.1) CYCLE !when height and radius is in proper value
!
!..............Sum of tl for the 2-d plane                 
!
               tl_vol_sum(level)=tl_vol_sum(level)+(tlall(i)-273.15d0)*volall(i)
               vol_sum(level)=vol_sum(level)+volall(i)
               IF(tl_min(level) .gt. (tlall(i)-273.15d0))tl_min(level)=tlall(i)-273.15d0            
           ENDDO
         ENDDO
!
!.......Averaging the tl for the 2-d plane
!        
         DO i=1,nheight
            IF(vol_sum(i).gt.0.d0)THEN
               tl_vol_avg(i)=tl_vol_sum(i)/vol_sum(i)
            ELSE
               PRINT*,'vol_sum(i) is zero for',vol_sum(i),'level in SUBROUTINE rocom_output_plane'           
            ENDIF           
         ENDDO   
!
         tl_space_avg(1)=tl_vol_avg(1)     
         tl_space_min(1)=tl_min(1)    
!        
      ENDIF
!  
      DEALLOCATE(vlp)
      DEALLOCATE(pr,tlall)
      DEALLOCATE(cboron)  
      DEALLOCATE(volall)
!
      END SUBROUTINE rocom_output_plane_user
!     
!--------------------------------------------------------------------------------------------------------
!
      SUBROUTINE make_rectangle(x,y,node)
!
      USE Zparam          , ONLY: pi
      USE Zio_unit        , ONLY: unit_log
!
      IMPLICIT NONE 
!
!.....Input
      INTEGER :: node(4)
      REAL(8) :: x(4),y(4)
!.....Local variables
      INTEGER :: i,j
      REAL(8) :: xmid,ymid,theta(4)
!.....Local arrays 
      INTEGER :: indx(4),itmp(4)
      REAL(8) :: tmpx(4),tmpy(4)
!      
      xmid=(x(1)+x(2)+x(3)+x(4))/4.d0
      ymid=(y(1)+y(2)+y(3)+y(4))/4.d0
!
!.....Make 4 angles
!
      DO i=1,4
         x(i)=x(i)-xmid
         y(i)=y(i)-ymid
         IF(ABS(x(i)).gt.0.d0)THEN
            theta(i)=ATAN(y(i)/x(i))/pi*180.d0
            IF(theta(i).lt.0.d0)theta(i)=theta(i)+180.d0
         ELSE
            IF(y(i).gt.0.d0)THEN
               theta(i)=90.d0
            ELSEIF(y(i).lt.0.d0)THEN
               theta(i)=270.d0
            ELSE
               WRITE(*,*)'Error in make_rectangle:atan!!!'
               WRITE(unit_log,*)'Error in make_rectangle:atan!!!'
               STOP
            ENDIF 
         ENDIF
         IF(x(i).gt.0.d0 .and. y(i).ge.0.d0)theta(i)=theta(i)     
         IF(x(i).lt.0.d0 .and. y(i).ge.0.d0)theta(i)=theta(i)     
         IF(x(i).lt.0.d0 .and. y(i).lt.0.d0)theta(i)=theta(i)+180.d0
         IF(x(i).gt.0.d0 .and. y(i).lt.0.d0)theta(i)=theta(i)+180.d0
      ENDDO
!
!.....Sort 4 angles
!
      CALL sortx_r(theta,indx,4)
      DO i=1,4
         itmp(i)=node(i)
         tmpx(i)=x(i)
         tmpy(i)=y(i)
      ENDDO
      DO i=1,4
         j=indx(i)
         node(i)=itmp(j)
         x(i)=tmpx(j)+xmid
         y(i)=tmpy(j)+ymid
      ENDDO
!
      END  SUBROUTINE make_rectangle
