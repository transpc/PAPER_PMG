!
      SUBROUTINE udfn_H2P1_outputHeTg(yy,vAxis)
!
!     Output : transient He molar fraction & gas temperature
!
      USE VOL_DATA        , ONLY: cell
      USE Zzone           , ONLY: ncell_fluid
      USE Zcore           , ONLY: np,myrank
      USE Ztimecon        , ONLY: time,itim
      USE Zcoord1         , ONLY: xloc
      USE Zcoord2         , ONLY: xfc_min,xfc_max
      USE Zmodel          , ONLY: molefr
!
      IMPLICIT NONE
!
      INTEGER,PARAMETER :: m=15
!.....Input
      INTEGER :: yy,vAxis
!.....Local variables
      INTEGER :: i,j
      INTEGER :: ip,i0,j0,loop
      INTEGER :: count,count_l
      INTEGER,SAVE :: vvx,vvy,vvz
      LOGICAL,SAVE :: initialmm=.TRUE.
      LOGICAL,SAVE :: initialHeTg=.TRUE.
      REAL(8) :: xi_he,xi_steam,mol_he,mol_steam,ni_he,ni_steam
!.....Local arrays
      INTEGER :: count_all(np),count_dsp(np)
      INTEGER :: count_all3(np),count_dsp3(np)
      INTEGER :: ia(np+1)
      INTEGER :: icell_index(ncell_fluid)
      INTEGER,SAVE :: INTERPnum(m,3)
      REAL(8) :: ztemp_l(ncell_fluid)
      REAL(8) :: He_mole(15),gasTemp(15)
      REAL(8),SAVE :: zz(15)
!.....Local allocatable arrays
      INTEGER,DIMENSION(:),ALLOCATABLE :: indx
      REAL(8),DIMENSION(:),ALLOCATABLE :: ztab1,Hetab1,Tgtab1
      REAL(8),DIMENSION(:),SAVE,ALLOCATABLE :: xmin1,xmax1,ymin1,ymax1
      REAL(8),DIMENSION(:),ALLOCATABLE :: prnvar_all
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: prnvar_l
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: prnvar1
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: prnvar
!
      mol_he=4.003d0
      mol_steam=18.02d0
      DO i=1,ncell_fluid
         xi_he=cell%quala(i)
         xi_steam=1.d0-cell%quala(i)
         ni_he=xi_he/mol_he
         ni_steam=xi_steam/mol_steam
         molefr(i)=ni_he/(ni_he+ni_steam)
      ENDDO
!
!.....find cells at specified x,y range
!
      IF(initialmm)THEN
         INTERPnum(:,1)=(/1,2,3,0,0,6,7,8,9,10,11,12,13,14,15/)
         INTERPnum(:,2)=(/0,2,3,4,5,6,7,8,9,10,11,12,13, 0, 0/)
         INTERPnum(:,3)=(/0,2,3,4,5,6,7,8,9,10,11,12,13, 0, 0/)
!
!----------   a        b        c        d        f        g       gh        h
         zz=(/0.538d0, 1.076d0, 1.726d0, 2.376d0, 3.036d0, 4.100d0, 4.326d0,          &
              4.976d0, 5.301d0, 5.626d0, 6.000d0, 6.276d0, 6.926d0, 7.478d0, 8.030d0/)
!----------   i       j        m        n        r        s       t
!
         ALLOCATE(xmin1(ncell_fluid),xmax1(ncell_fluid))
         ALLOCATE(ymin1(ncell_fluid),ymax1(ncell_fluid))
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
         initialmm=.FALSE.
      ENDIF
!
      count_l=0
      IF(yy.eq.1)THEN        ! 20
         DO i=1,ncell_fluid
            IF(xmin1(i).le.0.d0.and.xmax1(i).gt.0.d0)THEN
               IF(ymin1(i).le.0.d0.and.ymax1(i).gt.0.d0)THEN
                  count_l=count_l+1
                  icell_index(count_l)=i
                  ztemp_l(count_l)=xloc(i,vvz)
               ENDIF
            ENDIF
         ENDDO
      ELSEIF(yy.eq.2)THEN    ! 14
         IF(vAxis.eq.3)THEN
            DO i=1,ncell_fluid
               IF(xmin1(i).le.0.d0.and.xmax1(i).gt.0.d0)THEN
                  IF(ymin1(i).le.1.430d0.and.ymax1(i).gt.1.430d0)THEN
                     count_l=count_l+1
                     icell_index(count_l)=i
                     ztemp_l(count_l)=xloc(i,vvz)
                  ENDIF
               ENDIF
            ENDDO
         ELSEIF(vAxis.eq.2)THEN
            DO i=1,ncell_fluid
               IF(xmin1(i).le.0.d0.and.xmax1(i).gt.0.d0)THEN
                  IF(ymin1(i).le.-1.430d0.and.ymax1(i).gt.-1.430d0)THEN
                     count_l=count_l+1
                     icell_index(count_l)=i
                     ztemp_l(count_l)=xloc(i,vvz)
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ELSEIF(yy.eq.3)THEN    ! 26
         IF(vAxis.eq.3)THEN
            DO i=1,ncell_fluid
               IF(xmin1(i).le.0.d0.and.xmax1(i).gt.0.d0)THEN
                  IF(ymin1(i).le.-1.430d0.and.ymax1(i).gt.-1.430d0)THEN
                     count_l=count_l+1
                     icell_index(count_l)=i
                     ztemp_l(count_l) =xloc(i,vvz)
                  ENDIF
               ENDIF
            ENDDO
         ELSEIF(vAxis.eq.2)THEN
            DO i=1,ncell_fluid
               IF(xmin1(i).le.0.d0.and.xmax1(i).gt.0.d0)THEN
                  IF(ymin1(i).le.1.430d0.and.ymax1(i).gt.1.430d0)THEN
                     count_l=count_l+1
                     icell_index(count_l)=i
                     ztemp_l(count_l)=xloc(i,vvz)
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDIF
!
      ALLOCATE(prnvar_l(count_l,3))
      DO loop=1,count_l
         i=icell_index(loop)
         prnvar_l(loop,1)=ztemp_l(loop)
         prnvar_l(loop,2)=molefr(i)
         prnvar_l(loop,3)=cell%tg(i)
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
      CALL gather_vec_r(prnvar_l,count_l*3,prnvar_all,count*3,count_all3,count_dsp3)
      DEALLOCATE(prnvar_l)
      IF(myrank.eq.0) THEN
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
!.....Sort the extracted data
!
      IF(myrank.eq.0) THEN
         DO i=1,count
            IF(prnvar(i,1).gt.0.d0) CYCLE
            prnvar(i,1)=0.d0
            prnvar(i,2)=0.d0
            prnvar(i,3)=0.d0
         ENDDO
         ALLOCATE(indx(count),prnvar1(count,2))
         CALL sortex_r(prnvar(1,1),indx,count)
         DO i=1,count
            j=indx(i)
            prnvar1(i,1)=prnvar(j,2)
            prnvar1(i,2)=prnvar(j,3)
         ENDDO
         DO i=1,count
            prnvar(i,2)=prnvar1(i,1)
            prnvar(i,3)=prnvar1(i,2)
         ENDDO
         DEALLOCATE(indx,prnvar1)
!
         ALLOCATE(ztab1(count),Hetab1(count),Tgtab1(count))
         DO i=1,count
            ztab1(i )=prnvar(i,1)
            Hetab1(i)=prnvar(i,2)
            Tgtab1(i)=prnvar(i,3)
         ENDDO
         DEALLOCATE(prnvar)
!
!........Make output files
!
         He_mole(:)=0.d0
         gasTemp(:)=0.d0
!
         CALL INTERPV(count,m,zz,ztab1,Hetab1,He_mole,INTERPnum(1,yy))
         CALL INTERPV(count,m,zz,ztab1,Tgtab1,gasTemp,INTERPnum(1,yy))
         DO i=1,15
            IF(INTERPnum(m-i+1,yy).ne.0) gasTemp(i)=gasTemp(i)-273.16d0
         ENDDO
!
         IF(initialHeTg)THEN
            IF(yy.eq.1)THEN         !20
               OPEN(442,file='z_He_initial.dat')
               OPEN(443,file='z_Tg_initial.dat')
               OPEN(444,file='z_He_transient20.dat')
               OPEN(447,file='z_Tg_transient20.dat')
!
               WRITE(442,4441) 'itim & time =',itim, time     !initial He
               DO i=1,count
                  WRITE(442,4440) ztab1(i),Hetab1(i)
               ENDDO
               WRITE(443,4441) 'itim & time =',itim, time     !initial Tg
               DO i=1,count
                  WRITE(443,4440) ztab1(i),Tgtab1(i)-273.16d0
               ENDDO
               WRITE(444,4441) 'Time     A20     B20     C20     D20     F20     G20     GH20     H20     I20     J20     M20     N20     R20     S20     T20'
               WRITE(447,4441) 'Time     A20     B20     C20     D20     F20     G20     GH20     H20     I20     J20     M20     N20     R20     S20     T20'
            ELSEIF(yy.eq.2)THEN    !14
               OPEN(445,file='z_He_transient14.dat')
               OPEN(448,file='z_Tg_transient14.dat')
               WRITE(445,4441) 'Time     A14     B14     C14     D14     F14     G14     GH14     H14     I14     J14     M14     N14     R14     S14     T14'
               WRITE(448,4441) 'Time     A14     B14     C14     D14     F14     G14     GH14     H14     I14     J14     M14     N14     R14     S14     T14'
            ELSEIF(yy.eq.3)THEN    !26
               OPEN(446,file='z_He_transient26.dat')
               OPEN(449,file='z_Tg_transient26.dat')
               WRITE(446,4441) 'Time     A26     B26     C26     D26     F26     G26     GH26     H26     I26     J26     M26     N26     R26     S26     T26'
               WRITE(449,4441) 'Time     A26     B26     C26     D26     F26     G26     GH26     H26     I26     J26     M26     N26     R26     S26     T26'
               initialHeTg=.false.  
            ENDIF
         ENDIF
!
         IF(yy.eq.1)THEN        !20
            WRITE(444,4440) time,He_mole
            WRITE(447,4440) time,gasTemp
         ELSEIF(yy.eq.2)THEN    !14
            WRITE(445,4440) time,He_mole
            WRITE(448,4440) time,gasTemp
         ELSEIF(yy.eq.3)THEN    !26
            WRITE(446,4440) time,He_mole
            WRITE(449,4440) time,gasTemp
         ENDIF
         DEALLOCATE(ztab1,Hetab1,Tgtab1)
      ENDIF
!
4440  FORMAT(18(f22.16,1x))
4441  FORMAT(a,1x,i12,1x,f22.16)
!
      END SUBROUTINE udfn_H2P1_outputHeTg
!
      SUBROUTINE interpv(n,m,x,xtbl,ytbl,y,INTERPnum)
!
!     This routine does linear interpolation on m selected points out of
!     a list of n points.
!     x,xtbl are ordered in ascending order
!     Ouput y is in reverse order
!     INTERPnum(i)=0 skips the interpolation for point i 
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n,m
      INTEGER :: INTERPnum(m)
      REAL(8) :: x(m)
      REAL(8) :: xtbl(n),ytbl(n)
!.....Output
      REAL(8) :: y(m)
!.....Local variables
!     INTEGER :: i
      INTEGER :: ix,ix1,ix2
      INTEGER :: ip,ip1,ip2
      REAL(8) :: fact
!.....Get the bounds
!
      ix1=0
101   CONTINUE
      ix1=ix1+1
      IF(ix1.gt.m) RETURN
      IF(INTERPnum(ix1).eq.0) GOTO 101
      IF(x(ix1).le.xtbl(1)) THEN
         y(m-ix1+1)=ytbl(1)
         GOTO 101
      ENDIF
      ix2=m+1
102   CONTINUE
      ix2=ix2-1
      IF(ix2.lt.1) RETURN
      IF(INTERPnum(ix2).eq.0) GOTO 102
      IF(x(ix2).ge.xtbl(n)) THEN
         y(m-ix2+1)=ytbl(n)
         GOTO 102
      ENDIF
!
!.....Sarch for first eligible element in list
!
      ix=ix1
      ip1=1
      ip2=n
103   CONTINUE
      ip=(ip1+ip2)/2
      IF(ip.eq.ip1) GOTO 200
      IF(x(ix).lt.xtbl(ip)) THEN
         ip2=ip
      ELSEIF(x(ix).gt.xtbl(ip)) THEN
         ip1=ip
      ELSE
         y(m-ix+1)=ytbl(ip)
         GOTO 300
      ENDIF
      GOTO 103
200   CONTINUE
!
!.....Interpolate first eligible element
!
      fact=(x(ix)-xtbl(ip))/(xtbl(ip+1)-xtbl(ip))
      y(m-ix+1)=ytbl(ip)+fact*(ytbl(ip+1)-ytbl(ip))
300   CONTINUE
!
!.....Perform matching of 2 sorted lists and interpolate
!
      ix=ix+1
      IF(ix.gt.ix2) GOTO 400
      IF(INTERPnum(ix).eq.0) GOTO 300
!     write(*,*) '0=>',ix,x(ix),ip,xtbl(ip+1)
      IF(x(ix).lt.xtbl(ip+1)) THEN
         fact=(x(ix)-xtbl(ip))/(xtbl(ip+1)-xtbl(ip))
         y(m-ix+1)=ytbl(ip)+fact*(ytbl(ip+1)-ytbl(ip))
      ELSEIF(x(ix).gt.xtbl(ip+1)) THEN
301      CONTINUE
         ip=ip+1
         IF(ip.ge.n) GOTO 302
         IF(x(ix).gt.xtbl(ip+1)) THEN
            GOTO 301
         ELSEIF(x(ix).lt.xtbl(ip+1)) THEN
            fact=(x(ix)-xtbl(ip))/(xtbl(ip+1)-xtbl(ip))
            y(m-ix+1)=ytbl(ip)+fact*(ytbl(ip+1)-ytbl(ip))
         ELSE
            y(m-ix+1)=ytbl(ip+1)
         ENDIF
302      CONTINUE
      ELSE 
         y(m-ix+1)=ytbl(ip+1)
      ENDIF
      GOTO 300
400   CONTINUE
!
      END SUBROUTINE interpv
