!
      SUBROUTINE udfn_H2P1_outputEF(vA)
!
!     Output : erosion front progression from temperature 110'C
!
      USE Zconst1         , ONLY: vv_prob
      USE Ztimecon        , ONLY: time
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,vA
      INTEGER,SAVE :: counttime,ntime
      LOGICAL,SAVE :: initialEF=.true.
!.....Local allocatable arrays
      REAL(8),SAVE,ALLOCATABLE :: rtime(:)
      IF(initialEF)THEN
         counttime=1
         initialEF=.false.
!
         IF(vv_prob.eq.'h2p1_0'.or.vv_prob.eq.'h2p1_0x'.or.vv_prob.eq.'VD_h2p1_0')THEN
            ntime=120
         ELSEIF(vv_prob.eq.'h2p1_1'.or.vv_prob.eq.'h2p1_1x')THEN
            ntime=150
         ELSEIF(vv_prob.eq.'h2p1_2'.or.vv_prob.eq.'h2p1_2x')THEN
            ntime=170
         ELSEIF(vv_prob.eq.'h2p1_3'.or.vv_prob.eq.'h2p1_3x')THEN
            ntime=260
         ELSEIF(vv_prob.eq.'h2p1_4'.or.vv_prob.eq.'h2p1_4x')THEN
            ntime=270
         ENDIF
!
         ALLOCATE(rtime(ntime))
         DO i=1,ntime
            rtime(i)=i*10
         ENDDO
      ENDIF
!
      DO i=1,ntime
         IF(time.ge.rtime(i).and.counttime.eq.i)THEN
            CALL udfn_H2P1_outputEF2(vA)
            counttime=counttime+1
         ENDIF
      ENDDO
!
      END SUBROUTINE udfn_H2P1_outputEF
!
!------------------------------------------------------------------------
!
      SUBROUTINE udfn_H2P1_outputEF2(vAxis)
!
      USE VOL_DATA        , ONLY: cell
      USE Zzone           , ONLY: ncell_fluid,ncell_fluid_all
      USE Zcore           , ONLY: np,myrank
      USE Ztimecon        , ONLY: time
      USE Zcoord1         , ONLY: xloc
      USE Zcoord2         , ONLY: xfc_min,xfc_max
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,j,jj,na
      INTEGER :: count,numtab,vAxis,efef
      INTEGER,SAVE :: vvx,vvy,vvz
      LOGICAL,SAVE :: initial=.TRUE.
      LOGICAL,SAVE :: initialmm=.TRUE.
      REAL(8) :: xmin,xmax,ymin,ymax,zmin
      REAL(8) :: Tgj,zX,Ztg
!.....Local allocatable arrays
      REAL(8),DIMENSION(:),ALLOCATABLE :: ztemp0,Tgtemp0
      REAL(8),DIMENSION(:),ALLOCATABLE :: ztemp1,Tgtemp1
      REAL(8),DIMENSION(:),ALLOCATABLE :: ztab0,Tgtab0
      REAL(8),DIMENSION(:),ALLOCATABLE :: ztab1,Tgtab1
      REAL(8),DIMENSION(:),SAVE,ALLOCATABLE :: xmin1,xmax1,ymin1,ymax1
!
      na=ncell_fluid_all
!
!.....find cells at specified x,y range
      ALLOCATE(ztemp0(ncell_fluid),Tgtemp0(ncell_fluid))
      ztemp0(:)=0.d0
      Tgtemp0(:)=0.d0
!
      IF(initialmm)THEN
         initialmm=.FALSE.
         IF(.not.ALLOCATED(xmin1)) ALLOCATE(xmin1(ncell_fluid))
         IF(.not.ALLOCATED(xmax1)) ALLOCATE(xmax1(ncell_fluid))
         IF(.not.ALLOCATED(ymin1)) ALLOCATE(ymin1(ncell_fluid))
         IF(.not.ALLOCATED(ymax1)) ALLOCATE(ymax1(ncell_fluid))
      IF(vAxis.eq.3)THEN         !height:z
         vvx=1
         vvy=2
         vvz=3
         DO i=1,ncell_fluid
            xmin1(i)=xfc_min(i,1)
            xmax1(i)=xfc_max(i,1)
            ymin1(i)=xfc_min(i,2)
            ymax1(i)=xfc_max(i,2)
         ENDDO
      ELSEIF(vAxis.eq.2)THEN    !height:y
         vvx=3
         vvy=1
         vvz=2
         DO i=1,ncell_fluid
            xmin1(i)=xfc_min(i,3)
            xmax1(i)=xfc_max(i,3)
            ymin1(i)=xfc_min(i,2)
            ymax1(i)=xfc_max(i,2)
         ENDDO
      ENDIF
      ENDIF
!
      count=0
      DO i=1,ncell_fluid
         xmin=xmin1(i)
         xmax=xmax1(i)
         ymin=ymin1(i)
         ymax=ymax1(i)
!
         IF(xmin.le.0.d0.and.xmax.gt.0.d0)THEN
            IF(ymin.le.0.d0.and.ymax.gt.0.d0)THEN
               count=count+1
               ztemp0(i)=xloc(i,vvz)
               Tgtemp0(i)=cell%tg(i)
            ENDIF
         ENDIF
      ENDDO
!
!.....gathering selected ztemp0 and Hetemp0/Tgtemp0
      ALLOCATE(ztemp1(na),Tgtemp1(na))
      IF(np.gt.1) CALL allreducei_i1(count)
      CALL allgatherv_r(ztemp0,ztemp1,ncell_fluid,na,0)
      CALL allgatherv_r(Tgtemp0,Tgtemp1,ncell_fluid,na,0)
!
!.....make z and He/Tg table : not arranged
      ALLOCATE(ztab0(count),Tgtab0(count))
      ztab0(:)=0.d0
      Tgtab0(:)=0.d0
!
      numtab=1
      DO i=1,ncell_fluid_all
         IF(ztemp1(i).gt.0.d0)THEN
            ztab0(numtab)=ztemp1(i)
            Tgtab0(numtab)=Tgtemp1(i)
            numtab=numtab+1
         ENDIF
      ENDDO
!
!.....make z and He/Tg table : arrnaged (finish)
      ALLOCATE(ztab1(count),Tgtab1(count))
      DO i=1,count
         DO j=1,count
            IF(j.eq.1) zmin=ztab0(j)
            zmin=DMIN1(zmin,ztab0(j))
            IF(zmin.eq.ztab0(j))THEN
               jj=j
               Tgj=Tgtab0(j)
            ENDIF
         ENDDO
         ztab0(jj)=10.d0
         ztab1(i)=zmin
         Tgtab1(i)=Tgj-273.16d0
      ENDDO
!
      IF(myrank.eq.0)THEN
         IF(initial)THEN
            OPEN(451,file='z_Tg_erosionfront.dat')
            WRITE(451,4441) 'time     height     Tg(oC)'
            initial=.false.
         ENDIF
!
         efef=1
         zX=0.d0
         zTg=0.d0
         DO i=1,count
            IF(ztab1(i).ge.5.8d0.and.ztab1(i).le.8.2d0)THEN
               IF(efef.eq.1)THEN
                  zX=ztab1(i)
                  zTg=Tgtab1(i)
                  efef=2
               ELSE
                  IF(Tgtab1(i).ge.110.d0.and.Tgtab1(i+1).lt.110.d0)THEN
                     zX=ztab1(i)
                     zTg=Tgtab1(i)
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         WRITE(451,4440) time, zX, zTg
      ENDIF

 4440 FORMAT(18(f22.16,1x))
 4441 FORMAT(a,1x,i12,1x,f22.16,1x,i2)
!
      END SUBROUTINE udfn_H2P1_outputEF2
