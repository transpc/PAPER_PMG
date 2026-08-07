!
      SUBROUTINE udfn_H2P1_outputFOV(vAxis)
!
!     Output : FOV for selected time
!
      USE Zconst1         , ONLY: vv_prob
      USE Ztimecon        , ONLY: time
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,yy,avgOpt,vAxis
      INTEGER,SAVE :: FOVpos(10),FOVcount,FOVcount2
      LOGICAL,SAVE :: initialFOV=.true.
!.....Local arrays
      REAL(8),SAVE :: FOVtime(10),FOVs(10),FOVe(10)
!
      IF(initialFOV)THEN
         FOVcount=0
         FOVcount2=0
         FOVtime(:)=0.0d0
         FOVpos(:)=0
         initialFOV=.false.
!         
         IF(vv_prob.eq.'h2p1_0'.or.vv_prob.eq.'h2p1_0x'.or.vv_prob.eq.'VD_h2p1_0')THEN
            FOVtime=(/175.0d0, 321.0d0, 655.0d0, 1.0d4, 1.0d4, 1.0d4, 1.0d4, &
                     1.0d4, 1.0d4, 1.0d4/)
            FOVpos=(/2,2,3,0,0,0,0,0,0,0/)
         ELSEIF(vv_prob.eq.'h2p1_1'.or.vv_prob.eq.'h2p1_1x')THEN
            FOVtime=(/133.0d0, 423.0d0, 752.0d0, 947.0d0, 1078.0d0, 1739.0d0, &
                     1.0d4, 1.0d4, 1.0d4, 1.0d4/)
            FOVpos=(/1,2,3,3,3,1,0,0,0,0/)
         ELSEIF(vv_prob.eq.'h2p1_2'.or.vv_prob.eq.'h2p1_2x')THEN
            FOVtime=(/207.0d0, 509.0d0, 658.0d0, 954.0d0, 1114.0d0, 1767.0d0, &
                     2027.0d0, 1.0d4, 1.0d4, 1.0d4/)
            FOVpos=(/1,2,2,3,3,1,2,0,0,0/)
         ELSEIF(vv_prob.eq.'h2p1_3'.or.vv_prob.eq.'h2p1_3x')THEN
            FOVtime=(/157.0d0, 635.0d0, 865.0d0, 1130.0d0, 1520.0d0, 1750.0d0, &
                     2072.0d0, 1.0d4, 1.0d4, 1.0d4/)
            FOVpos=(/1,2,2,2,3,3,3,0,0,0/)
         ELSEIF(vv_prob.eq.'h2p1_4'.or.vv_prob.eq.'h2p1_4x')THEN
            FOVtime=(/317.0d0, 663.0d0, 983.0d0, 1269.0d0, 1602.0d0, 1858.0d0, &
                     2052.0d0, 2878.0d0, 3282.0d0, 1.0d4/)
            FOVpos=(/1,2,2,2,3,3,3,1,2,0/)
         ENDIF
         FOVs(:)=FOVtime(:)-102.4d0*0.5d0
         FOVe(:)=FOVtime(:)+102.4d0*0.5d0
      ENDIF
!
      DO i=1,10
         IF(time.ge.FOVtime(i).and.FOVcount.lt.i)THEN
            FOVcount=i
            yy=FOVpos(i)
            avgOpt=0
            CALL udfn_H2P1_outputFOV2(yy,FOVcount,avgOpt,vAxis)
            CALL udfn_H2P1_outputFOV3(FOVcount,avgOpt,vAxis)
         ENDIF
      ENDDO
!
      DO i=1,10
         IF(time.ge.FOVs(i).and.time.le.FOVe(i))THEN
            FOVcount2=i
            yy=FOVpos(i)
            avgOpt=1
            CALL udfn_H2P1_outputFOV2(yy,FOVcount2,avgOpt,vAxis)
            CALL udfn_H2P1_outputFOV3(FOVcount2,avgOpt,vAxis)
         ENDIF
      ENDDO
!
      END SUBROUTINE udfn_H2P1_outputFOV
!
!------------------------------------------------------------------------
!
      SUBROUTINE udfn_H2P1_outputFOV2(yy,FOVcount,avgOpt,vAxis)
!
!     FOV : horizontal
!
      USE Zzone           , ONLY: ncell_fluid
      USE Zcore           , ONLY: np,myrank
      USE Zcoord1         , ONLY: xloc
      USE Zcoord2         , ONLY: xfc_min,xfc_max
      USE Zturb           , ONLY: turb_keg
      USE Zvector         , ONLY: vg_n
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,j,k,yy
      INTEGER :: ip,i0,j0,loop
      INTEGER :: count_l,count
      INTEGER :: filenumb,FOVcount,FOVchar,avgOpt,vAxis
      INTEGER,SAVE :: avgCount(10),OUTcount
      INTEGER,SAVE :: vvx,vvy,vvz
      LOGICAL,SAVE :: initialavg=.TRUE.
      LOGICAL,SAVE :: initialmm=.TRUE.
      CHARACTER(20),SAVE :: filename,filename2
      REAL(8) :: xmin,xmax,zmin,zmax
!.....Local arrays
      INTEGER :: count_all(np),count_dsp(np)
      INTEGER :: count_all3(np),count_dsp3(np)
      INTEGER :: ia(np+1)
      INTEGER :: icell_index(ncell_fluid)
      REAL(8) :: ytemp10(ncell_fluid),Vztemp10(ncell_fluid)
!.....Local allocatable arrays
      INTEGER,DIMENSION(:),ALLOCATABLE :: indx
      REAL(8),DIMENSION(:),ALLOCATABLE :: ytab1,Vztab1,TKEtab1
      REAL(8),DIMENSION(:),SAVE,ALLOCATABLE :: Vztab1a,TKEtab1a
      REAL(8),DIMENSION(:),SAVE,ALLOCATABLE :: xmin1,xmax1,zmin1,zmax1
      REAL(8),DIMENSION(:),ALLOCATABLE :: prnvar1
      REAL(8),DIMENSION(:),ALLOCATABLE :: prnvar_all
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: prnvar_l
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: prnvar
      REAL(8),DIMENSION(:,:),SAVE,ALLOCATABLE :: Vztab1s,TKEtab1s
!
!.....find cells at specified x,z range
!
      IF(initialmm)THEN
         ALLOCATE(xmin1(ncell_fluid),xmax1(ncell_fluid))
         ALLOCATE(zmin1(ncell_fluid),zmax1(ncell_fluid))
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
         initialmm=.FALSE.
      ENDIF
!      
      count_l=0
      IF(yy.eq.1)THEN        !PosA
         DO i=1,ncell_fluid
            xmin=xmin1(i)
            xmax=xmax1(i)
            zmin=zmin1(i)
            zmax=zmax1(i)
            IF(xmin.le.0.0d0.and.xmax.gt.0.0d0)THEN
               IF(zmin.le.5.404d0.and.zmax.gt.5.404d0)THEN
                  count_l=count_l+1
                  icell_index(count_l)=i
                  ytemp10(count_l)=xloc(i,vvy)
                  Vztemp10(count_l)=vg_n(i,vvz)
               ENDIF
            ENDIF
         ENDDO
      ELSEIF(yy.eq.2)THEN    !PosB
         DO i=1,ncell_fluid
            xmin=xmin1(i)
            xmax=xmax1(i)
            zmin=zmin1(i)
            zmax=zmax1(i)
            IF(xmin.le.0.0d0.and.xmax.gt.0.0d0)THEN
               IF(zmin.le.6.150d0.and.zmax.gt.6.150d0)THEN
                  count_l=count_l+1
                  icell_index(count_l)=i
                  ytemp10(count_l)=xloc(i,vvy)
                  Vztemp10(count_l)=vg_n(i,vvz)
               ENDIF
            ENDIF
         ENDDO
      ELSEIF(yy.eq.3)THEN    !PosC
         DO i=1,ncell_fluid
            xmin=xmin1(i)
            xmax=xmax1(i)
            zmin=zmin1(i)
            zmax=zmax1(i)
            IF(xmin.le.0.0d0.and.xmax.gt.0.0d0)THEN
               IF(zmin.le.6.702d0.and.zmax.gt.6.702d0)THEN
                  count_l=count_l+1
                  icell_index(count_l)=i
                  ytemp10(count_l)=xloc(i,vvy)
                  Vztemp10(count_l)=vg_n(i,vvz)
               ENDIF
            ENDIF
         ENDDO
      ENDIF
!
      ALLOCATE(prnvar_l(count_l,3))
      DO loop=1,count_l
         i=icell_index(loop)
         prnvar_l(loop,1)=ytemp10(loop)
         prnvar_l(loop,2)=Vztemp10(loop)
         prnvar_l(loop,3)=turb_keg(i)
      ENDDO
!
      CALL allgather_i(count_l,count_all)
      ia(1)=1
      DO ip=1,np
         ia(ip+1)=ia(ip)+count_all(ip)
         count_all3(ip)=3*count_all(ip)
      ENDDO
      count_dsp(1)=0
      count_dsp3(1)=0
      DO ip=2,np
         count_dsp(ip)=count_dsp(ip-1)+count_all(ip-1)
         count_dsp3(ip)=count_dsp3(ip-1)+count_all3(ip-1)
      ENDDO
      count=ia(np+1)-1
!
      IF(myrank.eq.0)THEN
         ALLOCATE(prnvar_all(count*3))
      ELSE
         ALLOCATE(prnvar_all(1))
      ENDIF
      CALL gather_vec_r(prnvar_l,3*count_l,prnvar_all,count*3,count_all3,count_dsp3)
      DEALLOCATE(prnvar_l)
!
      IF(myrank.eq.0)THEN
         ALLOCATE(prnvar(count,3))
         DO ip=1,np
            j0=count_dsp3(ip)
            i0=1
            DO i=ia(ip),ia(ip+1)-1
               prnvar(i,1)=prnvar_all(i0+j0)
               prnvar(i,2)=prnvar_all(i0+j0+  count_all(ip))
               prnvar(i,3)=prnvar_all(i0+j0+2*count_all(ip))
               i0=i0+1
            ENDDO
         ENDDO
      ENDIF
      DEALLOCATE(prnvar_all)
!
      IF(myrank.eq.0)THEN
         ALLOCATE(indx(count))
         CALL sortex_r(prnvar(1,1),indx,count)
         ALLOCATE(prnvar1(count))
         DO k=2,3
            DO i=1,count
               j=indx(i)
               prnvar1(i)=prnvar(j,k)
            ENDDO
            DO i=1,count
               prnvar(i,k)=prnvar1(i)
            ENDDO
         ENDDO
         DEALLOCATE(prnvar1)
      ENDIF
!
!.....gathering selected ztemp0 and Hetemp0
      IF(myrank.eq.0)THEN
         ALLOCATE(ytab1(count),Vztab1(count),TKEtab1(count))
         DO i=1,count
            ytab1(i)=prnvar(i,1)
            Vztab1(i)=prnvar(i,2)
            TKEtab1(i)=prnvar(i,3)
         ENDDO
!
!.....output
         IF(initialavg)THEN
            ALLOCATE(Vztab1s(10,count),TKEtab1s(10,count),Vztab1a(count),TKEtab1a(count))
            avgCount(:)=0
            OUTcount=0
            Vztab1s(:,:)=0.0d0
            TKEtab1s(:,:)=0.0d0
            initialavg=.FALSE.
         ENDIF
!
         IF(avgOpt.eq.0)THEN          !center time
            filenumb=470+FOVcount
            FOVchar=48+FOVcount
            IF(yy.eq.1)THEN
               filename='z0'//ACHAR(FOVchar)//'_VzTKE_PosA.dat'
            ELSEIF(yy.eq.2)THEN
               filename='z0'//ACHAR(FOVchar)//'_VzTKE_PosB.dat'
            ELSEIF(yy.eq.3)THEN
               filename='z0'//ACHAR(FOVchar)//'_VzTKE_PosC.dat'
            ENDIF
!
            OPEN(filenumb,file=filename)
            DO i=1,count
               WRITE(filenumb,4440) ytab1(i),Vztab1(i),TKEtab1(i)
            ENDDO
            CLOSE(filenumb)
!
         ELSEIF(avgOpt.eq.1)THEN     !averaged time
            Vztab1a(:)=0.0d0
            TKEtab1a(:)=0.0d0
            avgCount(FOVcount)=avgCount(FOVcount)+1
            Vztab1s(FOVcount,:)=Vztab1s(FOVcount,:)+Vztab1(:)
            TKEtab1s(FOVcount,:)=TKEtab1s(FOVcount,:)+TKEtab1(:)
            Vztab1a(:)=Vztab1s(FOVcount,:)/avgCount(FOVcount)
            TKEtab1a(:)=TKEtab1s(FOVcount,:)/avgCount(FOVcount)
!
            filenumb=770+FOVcount
            FOVchar=48+FOVcount
            IF(yy.eq.1.and.OUTcount.lt.FOVcount)THEN
               filename2='z1'//ACHAR(FOVchar)//'_VzTKE_PosAa.dat'
               OUTcount=OUTcount+1
            ELSEIF(yy.eq.2.and.OUTcount.lt.FOVcount)THEN
               filename2='z1'//ACHAR(FOVchar)//'_VzTKE_PosBa.dat'
               OUTcount=OUTcount+1
            ELSEIF(yy.eq.3.and.OUTcount.lt.FOVcount)THEN
               filename2='z1'//ACHAR(FOVchar)//'_VzTKE_PosCa.dat'
               OUTcount=OUTcount+1
            ENDIF
!
            OPEN(filenumb,file=filename2)
            DO i=1,count
               WRITE(filenumb,4440) ytab1(i),Vztab1a(i),TKEtab1a(i)
            ENDDO
            CLOSE(filenumb)
         ENDIF
         DEALLOCATE(ytab1,Vztab1,TKEtab1)
         DEALLOCATE(indx)
         DEALLOCATE(prnvar)
      ENDIF
!
    4440 FORMAT(18(f22.16,1x))
!
      END SUBROUTINE udfn_H2P1_outputFOV2
!
!------------------------------------------------------------------------
!
      SUBROUTINE udfn_H2P1_outputFOV3(FOVcount,avgOpt,vAxis)
!
!     FOV : vertical
!
      USE Zzone           , ONLY: ncell_fluid
      USE Zcore           , ONLY: np,myrank
      USE Zcoord1         , ONLY: xloc
      USE Zcoord2         , ONLY: xfc_min,xfc_max
      USE Zvector         , ONLY: vg_n
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: FOVcount,avgOpt,vAxis
!.....Local variables
      INTEGER :: i,j
      INTEGER :: ip,i0,j0,loop
      INTEGER :: count_l,count
      INTEGER :: filenumb,FOVchar
      INTEGER,SAVE :: avgCount2(10),OUTcount2
      INTEGER,SAVE::vvx,vvy,vvz
      LOGICAL,SAVE :: initialavg2=.TRUE.
      LOGICAL,SAVE :: initialmm=.TRUE.
      REAL(8) :: xmin,xmax,ymin,ymax
!.....Local arrays
      INTEGER :: count_all(np),count_dsp(np)
      INTEGER :: count_all2(np),count_dsp2(np)
      INTEGER :: ia(np+1)
      REAL(8) :: ztemp10(ncell_fluid),Vytemp10(ncell_fluid)
!.....Local allocatable arrays
      INTEGER,DIMENSION(:),ALLOCATABLE :: indx
      REAL(8),DIMENSION(:),ALLOCATABLE :: ztab1,Vytab1
      REAL(8),DIMENSION(:),SAVE,ALLOCATABLE :: xmin1,xmax1,ymin1,ymax1
      REAL(8),DIMENSION(:),ALLOCATABLE :: prnvar1
      REAL(8),DIMENSION(:),ALLOCATABLE :: prnvar_all
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: prnvar_l
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: prnvar
      REAL(8),DIMENSION(:),SAVE,ALLOCATABLE :: Vytab1a
      REAL(8),DIMENSION(:,:),SAVE,ALLOCATABLE :: Vytab1s
!
      CHARACTER(20),SAVE :: filename11,filename12
!
!
!.....find cells at specified x,y range
!
!.....Get xfc scalar fron vector xn_nf first time or nbcon_change time
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
            ymin1(i)=xfc_min(i,1)
            ymax1(i)=xfc_max(i,1)
         ENDDO
      ENDIF
      ENDIF
!
      count_l=0
      DO i=1,ncell_fluid
         xmin=xmin1(i)
         xmax=xmax1(i)
         ymin=ymin1(i)
         ymax=ymax1(i)
!
         IF(xmin.le.0.0d0.and.xmax.gt.0.0d0)THEN
            IF(ymin.le.0.0d0.and.ymax.gt.0.0d0)THEN
               count_l=count_l+1
               ztemp10(count_l)=xloc(i,vvz)
               Vytemp10(count_l)=vg_n(i,vvz)
            ENDIF
         ENDIF
      ENDDO
!
      ALLOCATE(prnvar_l(count_l,3))
      DO loop=1,count_l
         prnvar_l(loop,1)=ztemp10(loop)
         prnvar_l(loop,2)=Vytemp10(loop)
      ENDDO
!
      CALL allgather_i(count_l,count_all)
      ia(1)=1
      DO ip=1,np
         ia(ip+1)=ia(ip)+count_all(ip)
         count_all2(ip)=2*count_all(ip)
      ENDDO
      count_dsp(1)=0
      count_dsp2(1)=0
      DO ip=2,np
         count_dsp(ip)=count_dsp(ip-1)+count_all(ip-1)
         count_dsp2(ip)=count_dsp2(ip-1)+count_all2(ip-1)
      ENDDO
      count=ia(np+1)-1
!
      IF(myrank.eq.0)THEN
         ALLOCATE(prnvar_all(count*2))
      ELSE
         ALLOCATE(prnvar_all(1))
      ENDIF
      CALL gather_vec_r(prnvar_l,2*count_l,prnvar_all,count*2,count_all2,count_dsp2)
      DEALLOCATE(prnvar_l)
!
      IF(myrank.eq.0)THEN
         ALLOCATE(prnvar(count,2))
         DO ip=1,np
            j0=count_dsp2(ip)
            i0=1
            DO i=ia(ip),ia(ip+1)-1
               prnvar(i,1)=prnvar_all(i0+j0)
               prnvar(i,2)=prnvar_all(i0+j0+  count_all(ip))
               i0=i0+1
            ENDDO
         ENDDO
         DO i=1,count
            IF(prnvar(i,1).gt.0.d0) CYCLE
            prnvar(i,1)=0.d0
            prnvar(i,2)=0.d0
         ENDDO
      ENDIF
      DEALLOCATE(prnvar_all)
!
      IF(myrank.eq.0)THEN
         ALLOCATE(indx(count))
         CALL sortex_r(prnvar(1,1),indx,count)
         ALLOCATE(prnvar1(count))
            DO i=1,count
               j=indx(i)
               prnvar1(i)=prnvar(j,2)
            ENDDO
            DO i=1,count
               prnvar(i,2)=prnvar1(i)
            ENDDO
         DEALLOCATE(prnvar1)
      ENDIF
!
!.....gathering selected ztemp0 and Hetemp0/Tgtemp0
      IF(myrank.eq.0)THEN
         ALLOCATE(ztab1(count),Vytab1(count))
         DO i=1,count
            ztab1(i)=prnvar(i,1)
            Vytab1(i)=prnvar(i,2)
         ENDDO
!
!.....output
         IF(initialavg2)THEN
            ALLOCATE(Vytab1s(10,count),Vytab1a(count))
            avgCount2(:)=0
            OUTcount2=0
            Vytab1s(:,:)=0.0d0
            initialavg2=.FALSE.
         ENDIF
!
         IF(avgOpt.eq.0)THEN          !center time
            filenumb=670+FOVcount
            FOVchar=48+FOVcount
            filename11='z2'//ACHAR(FOVchar)//'_Vy.dat'
!
            OPEN(filenumb,file=filename11)
            DO i=1,count
               WRITE(filenumb,4440) ztab1(i),Vytab1(i)
            ENDDO
            CLOSE(filenumb)
!
         ELSEIF(avgOpt.eq.1)THEN     !averaged time
            Vytab1a(:)=0.0d0
            avgCount2(FOVcount)=avgCount2(FOVcount)+1
            Vytab1s(FOVcount,:)=Vytab1s(FOVcount,:)+Vytab1(:)
            Vytab1a(:)=Vytab1s(FOVcount,:)/avgCount2(FOVcount)
!
            filenumb=870+FOVcount
            FOVchar=48+FOVcount
            IF(OUTcount2.lt.FOVcount)THEN
               filename12='z3'//ACHAR(FOVchar)//'_VyAvg.dat'
               OUTcount2=OUTcount2+1
            ENDIF
!
            OPEN(filenumb,file=filename12)
            DO i=1,count
               WRITE(filenumb,4440) ztab1(i),Vytab1a(i)
            ENDDO
            CLOSE(filenumb)
         ENDIF
         DEALLOCATE(ztab1,Vytab1)
         DEALLOCATE(indx)
         DEALLOCATE(prnvar)
      ENDIF
!
    4440 FORMAT(18(f22.16,1x))
!
      END SUBROUTINE udfn_H2P1_outputFOV3
