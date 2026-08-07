!
      SUBROUTINE udfn_H2P1_outputTE(vAxis)
!
!     Output : FOV at TE for selected time
!
      USE Zconst1         , ONLY: vv_prob
      USE Ztimecon        , ONLY: time
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,vAxis
      INTEGER,SAVE :: TEcount
      LOGICAL,SAVE :: initialTE=.true.
      REAL(8),SAVE :: TEtime(10)
!
      IF(initialTE)THEN
         TEcount=0
         TEtime(:)=0.d0
         initialTE=.false.
!         
         IF(vv_prob.eq.'h2p1_0'.or.vv_prob.eq.'h2p1_0x'.or.vv_prob.eq.'VD_h2p1_0')THEN
            TEtime=(/200.d0, 400.d0, 600.d0, 800.d0, 1000.d0, 1200.d0, 1400.d0, &
                     1600.d0, 1800.d0, 1.d4/)
         ELSEIF(vv_prob.eq.'h2p1_1'.or.vv_prob.eq.'h2p1_1x')THEN
            TEtime=(/200.d0, 400.d0, 600.d0, 800.d0, 1000.d0, 1200.d0, 1400.d0, &
                     1600.d0, 1800.d0, 1.d4/)
         ELSEIF(vv_prob.eq.'h2p1_2'.or.vv_prob.eq.'h2p1_2x')THEN
            TEtime=(/200.d0, 400.d0, 600.d0, 800.d0, 1000.d0, 1200.d0, 1400.d0, &
                     1600.d0, 1800.d0, 1.d4/)
         ELSEIF(vv_prob.eq.'h2p1_3'.or.vv_prob.eq.'h2p1_3x')THEN
            TEtime=(/200.d0, 400.d0, 600.d0, 800.d0, 1000.d0, 1200.d0, 1400.d0, &
                     1600.d0, 1800.d0, 1.d4/)
         ELSEIF(vv_prob.eq.'h2p1_4'.or.vv_prob.eq.'h2p1_4x')THEN
            TEtime=(/200.d0, 400.d0, 600.d0, 800.d0, 1000.d0, 1200.d0, 1400.d0, &
                     1600.d0, 1800.d0, 1.d4/)
         ENDIF
      ENDIF
!
      DO i=1,10
         IF(time.ge.TEtime(i).and.TEcount.lt.i)THEN
            TEcount=i
            CALL udfn_H2P1_outputTE2(TEcount,vAxis)
         ENDIF
      ENDDO
!
      END SUBROUTINE udfn_H2P1_outputTE
!
!------------------------------------------------------------------------
!
      SUBROUTINE udfn_H2P1_outputTE2(TEcount,vAxis)
!
      USE Zcore           , ONLY: np,myrank
      USE Zcoord1         , ONLY: xloc
      USE Zcoord2         , ONLY: xfc_min,xfc_max
      USE Ztimecon        , ONLY: time
      USE Zturb           , ONLY: turb_keg
      USE Zvector         , ONLY: vg_n
      USE Zzone           , ONLY: ncell_fluid,ncell_fluid_all
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,j,jj,na
      INTEGER :: count,numtab
      INTEGER :: filenumb,TEcount,TEchar,vAxis
      INTEGER,SAVE :: vvx,vvy,vvz
      LOGICAL,SAVE :: initialmm=.TRUE.
      REAL(8) :: xmin,xmax,ymin,zmin,zmax
      REAL(8) :: Vzj,TKEj
!.....Local allocatable arrays
      INTEGER,DIMENSION(:),ALLOCATABLE :: checkY0,checkY1
      REAL(8),DIMENSION(:),ALLOCATABLE :: ytemp0,Vztemp0,TKEtemp0
      REAL(8),DIMENSION(:),ALLOCATABLE :: ytemp1,Vztemp1,TKEtemp1
      REAL(8),DIMENSION(:),ALLOCATABLE :: ytab0,Vztab0,TKEtab0
      REAL(8),DIMENSION(:),ALLOCATABLE :: ytab1,Vztab1,TKEtab1
      REAL(8),DIMENSION(:),SAVE,ALLOCATABLE :: xmin1,xmax1,zmin1,zmax1
!
      CHARACTER(20) filename
!
      na=ncell_fluid_all
!
!.....find cells at specified x,z range
      ALLOCATE(ytemp0(ncell_fluid),Vztemp0(ncell_fluid))
      ALLOCATE(TKEtemp0(ncell_fluid),checkY0(ncell_fluid))
      ytemp0(:)=0.d0
      Vztemp0(:)=0.d0
      TKEtemp0(:)=0.d0
      checkY0(:)=0
!
      IF(initialmm)THEN
         initialmm=.FALSE.
         IF(.not.ALLOCATED(xmin1)) ALLOCATE(xmin1(ncell_fluid))
         IF(.not.ALLOCATED(xmax1)) ALLOCATE(xmax1(ncell_fluid))
         IF(.not.ALLOCATED(zmin1)) ALLOCATE(zmin1(ncell_fluid))
         IF(.not.ALLOCATED(zmax1)) ALLOCATE(zmax1(ncell_fluid))
      IF(vAxis.eq.3)THEN         !height:z
         vvx=1
         vvy=2
         vvz=3
         DO i=1,ncell_fluid
            xmin1(i)=xfc_min(i,1)
            xmax1(i)=xfc_max(i,1)
            zmin1(i)=xfc_min(i,3)
            zmax1(i)=xfc_max(i,3)
         ENDDO
      ELSEIF(vAxis.eq.2)THEN    !height:y
         vvx=3
         vvy=1
         vvz=2
         DO i=1,ncell_fluid
            xmin1(i)=xfc_min(i,3)
            xmax1(i)=xfc_max(i,3)
            zmin1(i)=xfc_min(i,2)
            zmax1(i)=xfc_max(i,2)
         ENDDO
      ENDIF
      ENDIF
!      
      count=0
      DO i=1,ncell_fluid
         xmin=xmin1(i)
         xmax=xmax1(i)
         zmin=zmin1(i)
         zmax=zmax1(i)
!
         IF(xmin.le.0.d0.and.xmax.gt.0.d0)THEN
            IF(zmin.le.4.041d0.and.zmax.gt.4.041d0)THEN
               count=count+1
               ytemp0(i)=xloc(i,vvy)
               Vztemp0(i)=vg_n(i,vvz)
               TKEtemp0(i)=turb_keg(i)
               checkY0(i)=1
            ENDIF
         ENDIF
      ENDDO
!
!.....gathering selected ztemp0 and Hetemp0
      ALLOCATE(ytemp1(na),vztemp1(na),TKEtemp1(na),checkY1(na))
      IF(np.gt.1) CALL allreducei_i1(count)
      CALL allgatherv_r(ytemp0,ytemp1,ncell_fluid,na,0)
      CALL allgatherv_r(Vztemp0,Vztemp1,ncell_fluid,na,0)
      CALL allgatherv_r(TKEtemp0,TKEtemp1,ncell_fluid,na,0)
      CALL allgatherv_i(checkY0,checkY1,ncell_fluid,na,0)
!
!.....make z and He/Tg table : not arranged
      ALLOCATE(ytab0(count),Vztab0(count),TKEtab0(count))
      ytab0(:)=0.d0
      Vztab0(:)=0.d0
      TKEtab0(:)=0.d0
!
      numtab=1
      DO i=1,ncell_fluid_all
         IF(checkY1(i).gt.0)THEN
            ytab0(numtab)=ytemp1(i)
            Vztab0(numtab)=Vztemp1(i)
            TKEtab0(numtab)=TKEtemp1(i)
            numtab=numtab+1
         ENDIF
      ENDDO
!
!.....make z and He/Tg table : arrnaged (finish)
      ALLOCATE(ytab1(count),Vztab1(count),TKEtab1(count))
      DO i=1,count
         DO j=1,count
            IF(j.eq.1) ymin=ytab0(j)
            ymin=DMIN1(ymin,ytab0(j))
            IF(ymin.eq.ytab0(j))THEN
               jj=j
               Vzj=Vztab0(j)
               TKEj=TKEtab0(j)
            ENDIF
         ENDDO
         ytab0(jj)=10.d0
         ytab1(i)=ymin
         Vztab1(i)=Vzj
         TKEtab1(i)=TKEj
      ENDDO
!
!.....output
      IF(myrank.eq.0)THEN
         filenumb=570+TEcount
         TEchar=48+TEcount
         filename='z_VzTKE_TE'//ACHAR(TEchar)//'.dat'
!
         OPEN(filenumb,file=filename)
         WRITE(filenumb,4441) 'time = ', time
         DO i=1,count
            WRITE(filenumb,4440) ytab1(i),Vztab1(i),TKEtab1(i)
         ENDDO
         CLOSE(filenumb)
      ENDIF
!
    4440 FORMAT(18(f22.16,1x))
    4441 FORMAT(a,1x,f22.16)
!
      DEALLOCATE(ytemp0,Vztemp0,TKEtemp0,ytemp1,Vztemp1,TKEtemp1)
      DEALLOCATE(checkY0,checkY1)
      DEALLOCATE(ytab0,Vztab0,TKEtab0,ytab1,Vztab1,TKEtab1)
!
      END SUBROUTINE udfn_H2P1_outputTE2
