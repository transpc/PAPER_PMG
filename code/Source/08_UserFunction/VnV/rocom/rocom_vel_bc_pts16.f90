!
      SUBROUTINE rocom_vel_bc_pts16(time) 
!
!     This routine defines velocity boundary condition
!
      USE Zcore           , ONLY: myrank      
      USE Zb_condition    , ONLY: vin_liq
      USE Zbc_index       , ONLY: nvin
      USE Zboron          , ONLY: cboronb,cboronb_liq        
      USE Zio_unit        , ONLY: unit_log
!
      IMPLICIT NONE
!
      INTEGER i,j,err,err1
!
      REAL(8),SAVE, ALLOCATABLE:: t_boron(:),t_vel(:),cboron(:),velin(:,:)
      REAL(8) time
!
      LOGICAL, SAVE::initial
!
      DATA INITIAL /.TRUE./                
!      
      IF(initial)THEN
         initial=.FALSE.
         ALLOCATE(t_boron(1200),t_vel(601),cboron(1200),velin(4,601))
!         
         IF(myrank.eq.0)WRITE(*,*)'ROCOM_12 Exp. data were loaded'
         IF(myrank.eq.0)WRITE(unit_log,*)'IROCOM_12 Exp. data were loaded'
!
         OPEN(201,file='boron_in.dat',status='old',iostat=err)
         OPEN(202,file='vel_in.dat',status='old',iostat=err1)
         IF(myrank.eq.0)OPEN(203,file='rocom_inlet.dat')
         IF(err.ne.0.or.err1.ne.0)then
            IF(myrank.eq.0)WRITE(*,*)'Input of ROCOM Boron Problem is missing'
            IF(myrank.eq.0)WRITE(unit_log,*)'Input of ROCOM Boron Problem is missing'
            STOP
         ENDIF
!
         DO i=1,1200     ! Read boron concentration data
            READ(201,*) t_boron(i),cboron(i)
            cboron(i)=DMIN1(DMAX1(cboron(i),0.0d0),1.0d0)
         ENDDO
!
         DO i=1,601      ! Read inlet velocity data
            READ(202,*) t_vel(i),(velin(j,i),j=1,4)  
         ENDDO                  
!
      ENDIF
!
      DO i=1,1200-1
         IF(time.ge.t_boron(i).and.time.lt.t_boron(i+1))THEN
            cboronb(1)=cboron(i)+(time-t_boron(i))*(cboron(i+1)-cboron(i))/(t_boron(i+1)-t_boron(i))
            cboronb(1)=DMIN1(DMAX1(cboronb(1),0.0d0),1.0d0)            
            cboronb_liq(1)=cboronb(1)
            CYCLE
         ENDIF
      ENDDO  
!
      DO i=1,601-1
         IF(time.ge.t_vel(i).and.time.lt.t_vel(i+1))THEN
            DO j=1,nvin
               vin_liq(j)=-velin(j,i)-(time-t_vel(i))*(velin(j,i+1)-velin(j,i))/(t_vel(i+1)-t_vel(i))    !for vin_norm=1
            ENDDO
            CYCLE
         ENDIF
      ENDDO 
!
      IF(myrank.eq.0)WRITE(203,1001)time,cboronb_liq(1),(vin_liq(i),i=1,nvin)
1001     FORMAT(1x,6f15.5) 
!
100   continue  
!
      RETURN
      END SUBROUTINE rocom_vel_bc_pts16 
