!
      SUBROUTINE udfn_CUBE_outputTs(yy)
!
!     Output : transient solid temperature
!
      USE SOLID_DATA      , ONLY: solid
      USE Zcoord1         , ONLY: xloc_c
      USE Zcoord2         , ONLY: xfc_c
      USE Zcore           , ONLY: np,myrank
      USE Ztimecon        , ONLY: time
      USE Znum_cell       , ONLY: i_neigh_c
      USE Zzone           , ONLY: ncell_cond,ncell_cond_all
!
      IMPLICIT NONE
!
      INTEGER i,j,jj,yy,na
      INTEGER count,numtab,ipos
      INTEGER,SAVE::vvx,vvy,vvz
      REAL(8) xmin,xmax,ymin,ymax,zmin, Tgj
!
      LOGICAL, SAVE::initialTs,initialmm
!
      REAL(8),ALLOCATABLE::ztemp0(:),Tgtemp0(:)
      REAL(8),ALLOCATABLE::ztemp1(:),Tgtemp1(:)
      REAL(8),ALLOCATABLE::ztab0(:),Tgtab0(:)
      REAL(8),ALLOCATABLE::ztab1(:),Tgtab1(:)
      REAL(8) solidT
      REAL(8) zz
      REAL(8),DIMENSION(:),SAVE,ALLOCATABLE :: xmin1,xmax1,ymin1,ymax1
!
      DATA initialTs,initialmm/.true.,.TRUE./
!
      na=ncell_cond_all
!
!.....find cells at specified x,y range
      ALLOCATE(ztemp0(ncell_cond),Tgtemp0(ncell_cond))
      ztemp0(:)=0.0d0
      Tgtemp0(:)=0.0d0
!
      IF(initialmm)THEN
         initialmm=.FALSE.
         IF(.not.ALLOCATED(xmin1)) ALLOCATE(xmin1(ncell_cond))
         IF(.not.ALLOCATED(xmax1)) ALLOCATE(xmax1(ncell_cond))
         IF(.not.ALLOCATED(ymin1)) ALLOCATE(ymin1(ncell_cond))
         IF(.not.ALLOCATED(ymax1)) ALLOCATE(ymax1(ncell_cond))
         vvx=1
         vvy=2
         vvz=3
         DO i=1,ncell_cond
            xmin=huge(0.d0)
            xmax=-xmin
            ymin=huge(0.d0)
            ymax=-ymin
            DO j=i_neigh_c(i),i_neigh_c(i+1)-1
               xmin=MIN(xmin,xfc_c(j,1))
               xmax=MAX(xmax,xfc_c(j,1))
               ymin=MIN(ymin,xfc_c(j,2))
               ymax=MAX(ymax,xfc_c(j,2))
            ENDDO
            xmin1(i)=xmin
            xmax1(i)=xmax
            ymin1(i)=ymin
            ymax1(i)=ymax
         ENDDO
      ENDIF
!
      count=0
!
      IF(yy.eq.1)THEN          !PZR
         DO i=1,ncell_cond
            xmin=xmin1(i)
            xmax=xmax1(i)
            ymin=ymin1(i)
            ymax=ymax1(i)
!            IF(xmin.le.-2.2d0.and.xmax.gt.-2.2d0)THEN
            IF(xmin.le.-2.22d0.and.xmax.gt.-2.22d0)THEN
               IF(ymin.le.1.05d0.and.ymax.gt.1.05d0)THEN
                  count=count+1
                  ztemp0(i)=xloc_c(i,vvz)
                  Tgtemp0(i)=solid%tsol(i)
               ENDIF
            ENDIF
         ENDDO
      ELSEIF(yy.eq.2)THEN      !SG01
         DO i=1,ncell_cond
            xmin=xmin1(i)
            xmax=xmax1(i)
            ymin=ymin1(i)
            ymax=ymax1(i)
            IF(xmin.le.0.0d0.and.xmax.gt.0.0d0)THEN
!               IF(ymin.le.-2.25d0.and.ymax.gt.-2.25d0)THEN
               IF(ymin.le.-2.3d0.and.ymax.gt.-2.3d0)THEN
                  count=count+1
                  ztemp0(i)=xloc_c(i,vvz)
                  Tgtemp0(i)=solid%tsol(i)
               ENDIF
            ENDIF
         ENDDO
      ELSEIF(yy.eq.3)THEN      !SG02
         DO i=1,ncell_cond
            xmin=xmin1(i)
            xmax=xmax1(i)
            ymin=ymin1(i)
            ymax=ymax1(i)
            IF(xmin.le.0.0d0.and.xmax.gt.0.0d0)THEN
!               IF(ymin.le.2.25d0.and.ymax.gt.2.25d0)THEN
               IF(ymin.le.2.3d0.and.ymax.gt.2.3d0)THEN
                  count=count+1
                  ztemp0(i)=xloc_c(i,vvz)
                  Tgtemp0(i)=solid%tsol(i)
               ENDIF
            ENDIF
         ENDDO
      ELSEIF(yy.eq.4)THEN      !RFP01
         DO i=1,ncell_cond
            xmin=xmin1(i)
            xmax=xmax1(i)
            ymin=ymin1(i)
            ymax=ymax1(i)
            IF(xmin.le.0.0d0.and.xmax.gt.0.0d0)THEN
               IF(ymin.le.-0.7d0.and.ymax.gt.-0.7d0)THEN
                  count=count+1
                  ztemp0(i)=xloc_c(i,vvz)
                  Tgtemp0(i)=solid%tsol(i)
               ENDIF
            ENDIF
         ENDDO
      ELSEIF(yy.eq.5)THEN      !RFP02
         DO i=1,ncell_cond
            xmin=xmin1(i)
            xmax=xmax1(i)
            ymin=ymin1(i)
            ymax=ymax1(i)
            IF(xmin.le.0.0d0.and.xmax.gt.0.0d0)THEN
               IF(ymin.le.0.7d0.and.ymax.gt.0.7d0)THEN
                  count=count+1
                  ztemp0(i)=xloc_c(i,vvz)
                  Tgtemp0(i)=solid%tsol(i)
               ENDIF
            ENDIF
         ENDDO
      ELSEIF(yy.eq.6)THEN      !SSW01
         DO i=1,ncell_cond
            xmin=xmin1(i)
            xmax=xmax1(i)
            ymin=ymin1(i)
            ymax=ymax1(i)
            IF(xmin.le.0.0d0.and.xmax.gt.0.0d0)THEN
!               IF(ymin.le.-2.25d0.and.ymax.gt.-2.25d0)THEN
               IF(ymin.le.-2.3d0.and.ymax.gt.-2.3d0)THEN
                  count=count+1
                  ztemp0(i)=xloc_c(i,vvz)
                  Tgtemp0(i)=solid%tsol(i)
               ENDIF
            ENDIF
         ENDDO
      ELSEIF(yy.eq.7)THEN      !SSW02
         DO i=1,ncell_cond
            xmin=xmin1(i)
            xmax=xmax1(i)
            ymin=ymin1(i)
            ymax=ymax1(i)
            IF(xmin.le.0.0d0.and.xmax.gt.0.0d0)THEN
!               IF(ymin.le.2.25d0.and.ymax.gt.2.25d0)THEN
               IF(ymin.le.2.3d0.and.ymax.gt.2.3d0)THEN
                  count=count+1
                  ztemp0(i)=xloc_c(i,vvz)
                  Tgtemp0(i)=solid%tsol(i)
               ENDIF
            ENDIF
         ENDDO
      ELSEIF(yy.eq.8)THEN      !PSW
         DO i=1,ncell_cond
            xmin=xmin1(i)
            xmax=xmax1(i)
            ymin=ymin1(i)
            ymax=ymax1(i)
!            IF(xmin.le.0.8d0.and.xmax.gt.0.8d0)THEN
            IF(xmin.le.0.75d0.and.xmax.gt.0.75d0)THEN   
               IF(ymin.le.0.0d0.and.ymax.gt.0.0d0)THEN
                  count=count+1
                  ztemp0(i)=xloc_c(i,vvz)
                  Tgtemp0(i)=solid%tsol(i)
               ENDIF
            ENDIF
         ENDDO
      ENDIF
!
!.....gathering selected ztemp0 and Tgtemp0
      ALLOCATE(ztemp1(na),Tgtemp1(na))
      ztemp1(:)=0.0d0
      Tgtemp1(:)=0.0d0
      IF(np.gt.1)THEN
         CALL allreducei_i1(count)
         CALL allgatherv_r(ztemp0,ztemp1,ncell_cond,na,1)
         CALL allgatherv_r(Tgtemp0,Tgtemp1,ncell_cond,na,1)
      ELSE
         ztemp1(:)=ztemp0(:)
         Tgtemp1(:)=Tgtemp0(:)
      ENDIF
!
!.....make z and He/Tg table : not arranged
      ALLOCATE(ztab0(count),Tgtab0(count))
      ztab0(:)=0.0d0
      Tgtab0(:)=0.0d0
!
      numtab=1
      DO i=1,ncell_cond_all
         IF(ztemp1(i).gt.0.0d0)THEN
            ztab0(numtab)=ztemp1(i)
            Tgtab0(numtab)=Tgtemp1(i)
            numtab=numtab+1
         ENDIF
      ENDDO
!
!.....make z and He/Tg table : arranged (finish)
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
         ztab0(jj)=20.0d0
         ztab1(i)=zmin
         Tgtab1(i)=Tgj
      ENDDO
!
!.....make output files
      IF(yy.eq.1)THEN
         zz=7.45d0
      ELSEIF(yy.eq.2)THEN
         zz=7.15d0
      ELSEIF(yy.eq.3)THEN
         zz=7.15d0
      ELSEIF(yy.eq.4)THEN
         zz=5.55d0
      ELSEIF(yy.eq.5)THEN
         zz=5.55d0
      ELSEIF(yy.eq.6)THEN
         zz=4.98d0
      ELSEIF(yy.eq.7)THEN
         zz=4.98d0
      ELSEIF(yy.eq.8)THEN
         zz=4.7d0
      ENDIF
!
      solidT=0.0d0
      ipos=1
      CALL INTERP1(count,zz,ztab1,Tgtab1,ipos,solidT)
      solidT=solidT-273.15d0
!
      IF(initialTs.and.myrank.eq.0)THEN
         IF(yy.eq.1)THEN
            OPEN(1721,file='z_Ts_transient_PZR.dat')
         ELSEIF(yy.eq.2)THEN
            OPEN(1722,file='z_Ts_transient_SG01.dat')
         ELSEIF(yy.eq.3)THEN
            OPEN(1723,file='z_Ts_transient_SG02.dat')
         ELSEIF(yy.eq.4)THEN
            OPEN(1724,file='z_Ts_transient_RFP01.dat')
         ELSEIF(yy.eq.5)THEN
            OPEN(1725,file='z_Ts_transient_RFP02.dat')
         ELSEIF(yy.eq.6)THEN
            OPEN(1726,file='z_Ts_transient_SSW01.dat')
         ELSEIF(yy.eq.7)THEN
            OPEN(1727,file='z_Ts_transient_SSW02.dat')
         ELSEIF(yy.eq.8)THEN
            OPEN(1728,file='z_Ts_transient_PSW.dat')
            initialTs=.false.  
         ENDIF
      ENDIF
!
      IF(yy.eq.1)THEN
         IF(myrank.eq.0) WRITE(1721,4440) time,solidT
      ELSEIF(yy.eq.2)THEN
         IF(myrank.eq.0) WRITE(1722,4440) time,solidT
      ELSEIF(yy.eq.3)THEN
         IF(myrank.eq.0) WRITE(1723,4440) time,solidT
      ELSEIF(yy.eq.4)THEN
         IF(myrank.eq.0) WRITE(1724,4440) time,solidT
      ELSEIF(yy.eq.5)THEN
         IF(myrank.eq.0) WRITE(1725,4440) time,solidT
      ELSEIF(yy.eq.6)THEN
         IF(myrank.eq.0) WRITE(1726,4440) time,solidT
      ELSEIF(yy.eq.7)THEN
         IF(myrank.eq.0) WRITE(1727,4440) time,solidT
      ELSEIF(yy.eq.8)THEN
         IF(myrank.eq.0) WRITE(1728,4440) time,solidT
      ENDIF
!
   4440 FORMAT(3(f22.16,1x))
!
      DEALLOCATE(ztemp0,Tgtemp0,ztemp1,Tgtemp1,ztab0,Tgtab0,ztab1,Tgtab1)
!
      RETURN
      END SUBROUTINE udfn_CUBE_outputTs
