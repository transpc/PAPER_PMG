!
      SUBROUTINE udfn_CUBE_outputTg(yy)
!
!     Output : transient gas temperature
!
      USE VOL_DATA        , ONLY: cell
      USE Zzone           , ONLY: ncell_fluid
      USE Zcore           , ONLY: np,myrank
      USE Ztimecon        , ONLY: time
      USE Zcoord1         , ONLY: xloc
      USE Zcoord2         , ONLY: xfc_min,xfc_max
!
      IMPLICIT NONE
!
      INTEGER,PARAMETER :: m=7
!.....Input
      INTEGER :: yy,vAxis
!.....Local variables
      INTEGER :: i,j
      INTEGER :: ip,i0,j0,loop
      INTEGER :: count,count_l
      INTEGER,SAVE :: vvx,vvy,vvz
      LOGICAL,SAVE :: initialmm=.TRUE.
      LOGICAL,SAVE :: initialHeTg=.TRUE.
!.....Local arrays
      INTEGER :: count_all(np),count_dsp(np)
      INTEGER :: count_all3(np),count_dsp3(np)
      INTEGER :: ia(np+1)
      INTEGER :: icell_index(ncell_fluid)
      INTEGER,SAVE :: INTERPnum(m,13)
      REAL(8) :: ztemp_l(ncell_fluid)
      REAL(8) :: gasTemp(7)
      REAL(8),SAVE :: zz(7)
!.....Local allocatable arrays
      INTEGER,DIMENSION(:),ALLOCATABLE :: indx
      REAL(8),DIMENSION(:),ALLOCATABLE :: ztab1,Tgtab1
      REAL(8),DIMENSION(:),SAVE,ALLOCATABLE :: xmin1,xmax1,ymin1,ymax1
      REAL(8),DIMENSION(:),ALLOCATABLE :: prnvar_all
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: prnvar_l
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: prnvar1
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: prnvar
!
!.....find cells at specified x,y range
!
      IF(initialmm)THEN
         INTERPnum(:,1)=(/0,2,3,4,5,6,7/)
         INTERPnum(:,2)=(/0,2,3,4,5,6,7/)
         INTERPnum(:,3)=(/0,2,3,4,5,6,7/)
         INTERPnum(:,4)=(/0,0,3,4,5,6,0/)
         INTERPnum(:,5)=(/0,2,3,4,5,6,7/)
         INTERPnum(:,6)=(/0,2,3,4,5,6,7/)
         INTERPnum(:,7)=(/0,0,3,4,5,6,0/)
         INTERPnum(:,8)=(/0,0,0,4,5,6,7/)
         INTERPnum(:,9)=(/1,0,3,4,5,6,7/)
         INTERPnum(:,10)=(/0,0,0,4,5,6,0/)
         INTERPnum(:,11)=(/0,0,0,4,5,6,7/)
         INTERPnum(:,12)=(/1,0,3,4,5,6,7/)
         INTERPnum(:,13)=(/0,0,0,4,5,6,0/)
!
!----------    07       06       05       04       03          02          01
         zz=(/4.98d0, 5.55d0, 7.15d0, 8.9592d0, 10.9592d0, 12.95936d0, 13.9592d0/)
!----------   1.580   2.150   3.750   5.5592     7.5592     9.55936    10.55920
!
         ALLOCATE(xmin1(ncell_fluid),xmax1(ncell_fluid))
         ALLOCATE(ymin1(ncell_fluid),ymax1(ncell_fluid))
         vvx=1
         vvy=2
         vvz=3
         DO i=1,ncell_fluid
            xmin1(i)=xfc_min(i,1)
            xmax1(i)=xfc_max(i,1)
            ymin1(i)=xfc_min(i,2)
            ymax1(i)=xfc_max(i,2)
         ENDDO
         initialmm=.FALSE.
      ENDIF
!
      count_l=0
      IF(yy.eq.1)THEN        ! A
         DO i=1,ncell_fluid
            IF(xmin1(i).le.0.d0.and.xmax1(i).gt.0.d0)THEN
               IF(ymin1(i).le.0.d0.and.ymax1(i).gt.0.d0)THEN
                  count_l=count_l+1
                  icell_index(count_l)=i
                  ztemp_l(count_l)=xloc(i,vvz)
               ENDIF
            ENDIF
         ENDDO
      ELSEIF(yy.eq.2)THEN    ! B
         DO i=1,ncell_fluid
            IF(xmin1(i).le.0.75d0.and.xmax1(i).gt.0.75d0)THEN
               IF(ymin1(i).le.0.0d0.and.ymax1(i).gt.0.0d0)THEN
                  count_l=count_l+1
                  icell_index(count_l)=i
                  ztemp_l(count_l)=xloc(i,vvz)
               ENDIF
            ENDIF
         ENDDO
      ELSEIF(yy.eq.3)THEN    ! C
         DO i=1,ncell_fluid
            IF(xmin1(i).le.1.5d0.and.xmax1(i).gt.1.5d0)THEN
               IF(ymin1(i).le.0.0d0.and.ymax1(i).gt.0.0d0)THEN
                  count_l=count_l+1
                  icell_index(count_l)=i
                  ztemp_l(count_l) =xloc(i,vvz)
               ENDIF
            ENDIF
         ENDDO
      ELSEIF(yy.eq.4)THEN    ! D
         DO i=1,ncell_fluid
            IF(xmin1(i).le.2.1d0.and.xmax1(i).gt.2.1d0)THEN
               IF(ymin1(i).le.0.0d0.and.ymax1(i).gt.0.0d0)THEN
                  count_l=count_l+1
                  icell_index(count_l)=i
                  ztemp_l(count_l) =xloc(i,vvz)
               ENDIF
            ENDIF
         ENDDO
      ELSEIF(yy.eq.5)THEN    ! E
         DO i=1,ncell_fluid
            IF(xmin1(i).le.-0.75d0.and.xmax1(i).gt.-0.75d0)THEN
               IF(ymin1(i).le.0.0d0.and.ymax1(i).gt.0.0d0)THEN
                  count_l=count_l+1
                  icell_index(count_l)=i
                  ztemp_l(count_l) =xloc(i,vvz)
               ENDIF
            ENDIF
         ENDDO
      ELSEIF(yy.eq.6)THEN    ! F
         DO i=1,ncell_fluid
            IF(xmin1(i).le.-1.5d0.and.xmax1(i).gt.-1.5d0)THEN
               IF(ymin1(i).le.0.0d0.and.ymax1(i).gt.0.0d0)THEN
                  count_l=count_l+1
                  icell_index(count_l)=i
                  ztemp_l(count_l) =xloc(i,vvz)
               ENDIF
            ENDIF
         ENDDO
      ELSEIF(yy.eq.7)THEN    ! G
         DO i=1,ncell_fluid
            IF(xmin1(i).le.-2.1d0.and.xmax1(i).gt.-2.1d0)THEN
               IF(ymin1(i).le.0.0d0.and.ymax1(i).gt.0.0d0)THEN
                  count_l=count_l+1
                  icell_index(count_l)=i
                  ztemp_l(count_l) =xloc(i,vvz)
               ENDIF
            ENDIF
         ENDDO
      ELSEIF(yy.eq.8)THEN    ! H
         DO i=1,ncell_fluid
            IF(xmin1(i).le.0.d0.and.xmax1(i).gt.0.d0)THEN
               IF(ymin1(i).le.-0.75d0.and.ymax1(i).gt.-0.75d0)THEN
                  count_l=count_l+1
                  icell_index(count_l)=i
                  ztemp_l(count_l) =xloc(i,vvz)
               ENDIF
            ENDIF
         ENDDO
      ELSEIF(yy.eq.9)THEN    ! I
         DO i=1,ncell_fluid
            IF(xmin1(i).le.0.d0.and.xmax1(i).gt.0.d0)THEN
               IF(ymin1(i).le.-1.5d0.and.ymax1(i).gt.-1.5d0)THEN
                  count_l=count_l+1
                  icell_index(count_l)=i
                  ztemp_l(count_l) =xloc(i,vvz)
               ENDIF
            ENDIF
         ENDDO
      ELSEIF(yy.eq.10)THEN    ! J
         DO i=1,ncell_fluid
            IF(xmin1(i).le.0.d0.and.xmax1(i).gt.0.d0)THEN
               IF(ymin1(i).le.-2.1d0.and.ymax1(i).gt.-2.1d0)THEN
                  count_l=count_l+1
                  icell_index(count_l)=i
                  ztemp_l(count_l) =xloc(i,vvz)
               ENDIF
            ENDIF
         ENDDO
      ELSEIF(yy.eq.11)THEN    ! K
         DO i=1,ncell_fluid
            IF(xmin1(i).le.0.d0.and.xmax1(i).gt.0.d0)THEN
               IF(ymin1(i).le.0.75d0.and.ymax1(i).gt.0.75d0)THEN
                  count_l=count_l+1
                  icell_index(count_l)=i
                  ztemp_l(count_l) =xloc(i,vvz)
               ENDIF
            ENDIF
         ENDDO
      ELSEIF(yy.eq.12)THEN    ! L
         DO i=1,ncell_fluid
            IF(xmin1(i).le.0.d0.and.xmax1(i).gt.0.d0)THEN
               IF(ymin1(i).le.1.5d0.and.ymax1(i).gt.1.5d0)THEN
                  count_l=count_l+1
                  icell_index(count_l)=i
                  ztemp_l(count_l) =xloc(i,vvz)
               ENDIF
            ENDIF
         ENDDO
      ELSEIF(yy.eq.13)THEN    ! M
         DO i=1,ncell_fluid
            IF(xmin1(i).le.0.d0.and.xmax1(i).gt.0.d0)THEN
               IF(ymin1(i).le.2.1d0.and.ymax1(i).gt.2.1d0)THEN
                  count_l=count_l+1
                  icell_index(count_l)=i
                  ztemp_l(count_l) =xloc(i,vvz)
               ENDIF
            ENDIF
         ENDDO         
      ENDIF
!
      ALLOCATE(prnvar_l(count_l,2))
      DO loop=1,count_l
         i=icell_index(loop)
         prnvar_l(loop,1)=ztemp_l(loop)
         prnvar_l(loop,2)=cell%tg(i)
      ENDDO
!
!.....Build count,disp info to call gather_vec
!
      CALL allgather_i(count_l,count_all)
      ia(1)=1
      DO ip=1,np
         ia(ip+1)=ia(ip)+count_all(ip)
         count_all3(ip)=count_all(ip)*3
      ENDDO
      count_dsp(1)=0
      count_dsp3(1)=0
      DO ip=2,np
         count_dsp(ip)=count_dsp(ip-1)+count_all(ip-1)
         count_dsp3(ip)=count_dsp3(ip-1)+count_all3(ip-1)
      ENDDO
      count=ia(np+1)-1
!
      IF(myrank.eq.0) THEN
         ALLOCATE(prnvar_all(count*3))
      ELSE
         ALLOCATE(prnvar_all(1))
      ENDIF
      CALL gather_vec_r(prnvar_l,count_l*2,prnvar_all,count*2,count_all3,count_dsp3)
      DEALLOCATE(prnvar_l)
      IF(myrank.eq.0) THEN
         ALLOCATE(prnvar(count,2))
         DO ip=1,np
            j0=count_dsp3(ip)
            i0=1
            DO i=ia(ip),ia(ip+1)-1
               prnvar(i,1)=prnvar_all(i0+j0)
               prnvar(i,2)=prnvar_all(i0+j0+  count_all(ip))
               i0=i0+1
            ENDDO
         ENDDO
      ENDIF
      DEALLOCATE(prnvar_all)
!
!.....Sort the extracted data
!
      IF(myrank.eq.0) THEN
         DO i=1,count
            IF(prnvar(i,1).gt.0.d0) CYCLE
            prnvar(i,1)=0.d0
            prnvar(i,2)=0.d0
         ENDDO
         ALLOCATE(indx(count),prnvar1(count,1))
         CALL sortex_r(prnvar(1,1),indx,count)
         DO i=1,count
            j=indx(i)
            prnvar1(i,1)=prnvar(j,2)
         ENDDO
         DO i=1,count
            prnvar(i,2)=prnvar1(i,1)
         ENDDO
         DEALLOCATE(indx,prnvar1)
!
         ALLOCATE(ztab1(count),Tgtab1(count))
         DO i=1,count
            ztab1(i )=prnvar(i,1)
            Tgtab1(i)=prnvar(i,2)
         ENDDO
         DEALLOCATE(prnvar)
!
!........Make output files
!
         gasTemp(:)=0.d0
!
         CALL INTERPV(count,m,zz,ztab1,Tgtab1,gasTemp,INTERPnum(1,yy))
         DO i=1,7
            IF(INTERPnum(m-i+1,yy).ne.0) gasTemp(i)=gasTemp(i)-273.15d0
         ENDDO
!
         IF(initialHeTg)THEN
            IF(yy.eq.1)THEN         !A
               OPEN(1701,file='z_Tg_transientA.dat')
            ELSEIF(yy.eq.2)THEN     !B
               OPEN(1702,file='z_Tg_transientB.dat')
            ELSEIF(yy.eq.3)THEN     !C
               OPEN(1703,file='z_Tg_transientC.dat')
            ELSEIF(yy.eq.4)THEN     !D
               OPEN(1704,file='z_Tg_transientD.dat')
            ELSEIF(yy.eq.5)THEN     !E
               OPEN(1705,file='z_Tg_transientE.dat')
            ELSEIF(yy.eq.6)THEN     !F
               OPEN(1706,file='z_Tg_transientF.dat')
            ELSEIF(yy.eq.7)THEN     !G
               OPEN(1707,file='z_Tg_transientG.dat')
            ELSEIF(yy.eq.8)THEN     !H
               OPEN(1708,file='z_Tg_transientH.dat')
            ELSEIF(yy.eq.9)THEN     !I
               OPEN(1709,file='z_Tg_transientI.dat')
            ELSEIF(yy.eq.10)THEN     !J
               OPEN(1710,file='z_Tg_transientJ.dat')
            ELSEIF(yy.eq.11)THEN     !K
               OPEN(1711,file='z_Tg_transientK.dat')
            ELSEIF(yy.eq.12)THEN     !L
               OPEN(1712,file='z_Tg_transientL.dat')
            ELSEIF(yy.eq.13)THEN     !M
               OPEN(1713,file='z_Tg_transientM.dat')
               initialHeTg=.false.  
            ENDIF
         ENDIF
!
         IF(yy.eq.1)THEN        !A
            WRITE(1701,4440) time,gasTemp
         ELSEIF(yy.eq.2)THEN    !B
           WRITE(1702,4440) time,gasTemp
         ELSEIF(yy.eq.3)THEN    !C
            WRITE(1703,4440) time,gasTemp
         ELSEIF(yy.eq.4)THEN    !D
            WRITE(1704,4440) time,gasTemp
         ELSEIF(yy.eq.5)THEN    !E
            WRITE(1705,4440) time,gasTemp
         ELSEIF(yy.eq.6)THEN    !F
            WRITE(1706,4440) time,gasTemp
         ELSEIF(yy.eq.7)THEN    !G
            WRITE(1707,4440) time,gasTemp
         ELSEIF(yy.eq.8)THEN    !H
            WRITE(1708,4440) time,gasTemp
         ELSEIF(yy.eq.9)THEN    !I
            WRITE(1709,4440) time,gasTemp
         ELSEIF(yy.eq.10)THEN    !J
            WRITE(1710,4440) time,gasTemp
         ELSEIF(yy.eq.11)THEN    !K
            WRITE(1711,4440) time,gasTemp
         ELSEIF(yy.eq.12)THEN    !L
            WRITE(1712,4440) time,gasTemp
         ELSEIF(yy.eq.13)THEN    !M
            WRITE(1713,4440) time,gasTemp
         ENDIF
         DEALLOCATE(ztab1,Tgtab1)
      ENDIF
!
4440  FORMAT(18(f22.16,1x))
4441  FORMAT(a,1x,i12,1x,f22.16)
!
      END SUBROUTINE udfn_CUBE_outputTg
!