!
      SUBROUTINE calc_geo
!
      USE Zinterface
      USE Zzone       , ONLY: ncell_fluid
      USE Zcore       , ONLY: np
      USE Zparam      , ONLY: ndim,nin_max,nb_max,nb_sym,nb_mars
      USE Znum_cell   , ONLY: i_neigh,neigh
      USE Zbc_index   , ONLY: nbcon,ngrad
      USE Zcoord1     , ONLY: xloc,xloc_m
      USE Zcoord2     , ONLY: fac,fac1,xfc,cell_leng, &
                              xfc_min,xfc_max,xloc_xfc_min,xloc_xfc_max, &
                              xloc_xfc_radius_min,xloc_xfc_radius_max
      USE Zcoord3     , ONLY: sv,svp,vol,volr,volp,volpr,porosity,permeability
      USE Zcoord4     , ONLY: sa,saa,sap,dji,dji_x,dji_a,sad,dnj
      USE Znormal     , ONLY: xn
      USE Ztimecon    , ONLY: dxmin
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,j,k,ix
      REAL(8) :: dr,dr1,dr2
      REAL(8) :: xloc_k1,xloc_k2,xloc_k3
      REAL(8) :: theta,theta_max
      REAL(8) :: xmin,ymin,zmin
      REAL(8) :: xmax,ymax,zmax
      REAL(8) :: radiusn,radiusmin,radiusmax
!
      theta_max=0.d0
      dxmin=1.d10
!
      DO i=1,ncell_fluid
         volp(i)=vol(i)*porosity(i)
         volr(i) =1.d0/vol(i)
         volpr(i)=1.d0/volp(i)
      ENDDO 
!
      IF(ndim.eq.2) THEN
         DO i=1,ncell_fluid
            xmin=huge(0.d0)
            xmax=-xmin
            ymin=huge(0.d0)
            ymax=-ymin
            radiusn=ABS(xloc(i,1))
            radiusmax=radiusn
            radiusmin=radiusn
!DIR$ NOVECTOR
            DO j=i_neigh(i),i_neigh(i+1)-1
!
!..............Get min max xfc and min max xloc,xfc
!
               xmin=MIN(xmin,xfc(j,1))
               xmax=MAX(xmax,xfc(j,1))
               ymin=MIN(ymin,xfc(j,2))
               ymax=MAX(ymax,xfc(j,2))
               radiusn=SQRT(xfc(j,1)**2+xfc(j,2)**2)
               radiusmax=MAX(radiusmax,radiusn)
               radiusmin=MIN(radiusmin,radiusn)
!
!..............Primary geo data
!
               sa(j)=sv(j,1)**2+sv(j,2)**2
               sa(j)=sqrt(sa(j))
               xn(j,1)=sv(j,1)/sa(j)
               xn(j,2)=sv(j,2)/sa(j)
               IF(nbcon(j).eq.0)THEN
                  k=neigh(j)
                  xloc_k1=xloc(k,1)
                  xloc_k2=xloc(k,2)
               ELSE
                  xloc_k1=xfc(j,1)
                  xloc_k2=xfc(j,2)
               ENDIF
               dji_x(j,1)=xloc_k1-xloc(i,1)
               dji_x(j,2)=xloc_k2-xloc(i,2)
               dji(j)=dji_x(j,1)**2+dji_x(j,2)**2
               dr1=(xfc(j,1)-xloc(i,1))**2+(xfc(j,2)-xloc(i,2))**2
               dr2=(xfc(j,1)-xloc_k1  )**2+(xfc(j,2)-xloc_k2  )**2
!
               dji(j)=sqrt(dji(j))
               dr1=sqrt(dr1)
               dr2=sqrt(dr2)
               fac1(j)=dr2/(dr1+dr2)
               fac(j) =dr1/(dr1+dr2)
!              fac(j)=1.d0-fac1(j)
!..............Minimum length
               dxmin=min(dxmin,dji(j))
!
!..............Secondary variables
!
               saa(j)=sa(j)*permeability(j)
               dji_a(j)=dji_x(j,1)*xn(j,1)+dji_x(j,2)*xn(j,2)
               sap(j)=saa(j)/dji_a(j)
               sad(j)=saa(j)/dji(j)
!
!.............Modified cell center coordinates for non-orthogonal grid
!
              dr=(xfc(j,1)-xloc(i,1))*xn(j,1)+(xfc(j,2)-xloc(i,2))*xn(j,2)
              xloc_m(j,1)=xfc(j,1)-dr*xn(j,1)-xloc(i,1)
              xloc_m(j,2)=xfc(j,2)-dr*xn(j,2)-xloc(i,2)
!.............Why do this compute ??????
              dr2=(xfc(j,1)-xloc(i,1))**2+(xfc(j,2)-xloc(i,2))**2
              dr2=sqrt(dr2)
              dr2=abs(dr/dr2)
              dr2=max(dr2,0.d0)
              dr2=min(dr2,1.d0)
              theta=DACOSD(dr2)
              theta_max=max(theta_max,theta)
            ENDDO
            xfc_min(i,1)=xmin
            xfc_max(i,1)=xmax
            xfc_min(i,2)=ymin
            xfc_max(i,2)=ymax
            xloc_xfc_min(i,1)=MIN(xmin,xloc(i,1))
            xloc_xfc_max(i,1)=MAX(xmax,xloc(i,1))
            xloc_xfc_min(i,2)=MIN(ymin,xloc(i,2))
            xloc_xfc_max(i,2)=MAX(ymax,xloc(i,2))
            xloc_xfc_radius_max(i)=radiusmax
            xloc_xfc_radius_min(i)=radiusmin
         ENDDO
      ELSE
         DO i=1,ncell_fluid
            xmin=huge(0.d0)
            xmax=-xmin
            ymin=huge(0.d0)
            ymax=-ymin
            zmin=huge(0.d0)
            zmax=-zmin
            radiusn=SQRT(xloc(i,1)**2+xloc(i,2)**2)
            radiusmax=radiusn
            radiusmin=radiusn
!DIR$ NOVECTOR
            DO j=i_neigh(i),i_neigh(i+1)-1
!
!..............Get min max xfc and min max xloc,xfc
!
               xmin=MIN(xmin,xfc(j,1))
               xmax=MAX(xmax,xfc(j,1))
               ymin=MIN(ymin,xfc(j,2))
               ymax=MAX(ymax,xfc(j,2))
               zmin=MIN(zmin,xfc(j,3))
               zmax=MAX(zmax,xfc(j,3))
               radiusn=SQRT(xfc(j,1)**2+xfc(j,2)**2)
               radiusmax=MAX(radiusmax,radiusn)
               radiusmin=MIN(radiusmin,radiusn)
!
!..............Primary geo data
!
               sa(j)=sv(j,1)**2+sv(j,2)**2+sv(j,3)**2
               sa(j)=sqrt(sa(j))
               xn(j,1)=sv(j,1)/sa(j)
               xn(j,2)=sv(j,2)/sa(j)
               xn(j,3)=sv(j,3)/sa(j)
               IF(nbcon(j).eq.0)THEN
                  k=neigh(j)
                  xloc_k1=xloc(k,1)
                  xloc_k2=xloc(k,2)
                  xloc_k3=xloc(k,3)
               ELSE
                  xloc_k1=xfc(j,1)
                  xloc_k2=xfc(j,2)
                  xloc_k3=xfc(j,3)
               ENDIF
               dji_x(j,1)=xloc_k1-xloc(i,1)
               dji_x(j,2)=xloc_k2-xloc(i,2)
               dji_x(j,3)=xloc_k3-xloc(i,3)
               dji(j)=dji_x(j,1)**2+dji_x(j,2)**2+dji_x(j,3)**2
               dr1=(xfc(j,1)-xloc(i,1))**2+(xfc(j,2)-xloc(i,2))**2+(xfc(j,3)-xloc(i,3))**2
               dr2=(xfc(j,1)-xloc_k1  )**2+(xfc(j,2)-xloc_k2  )**2+(xfc(j,3)-xloc_k3  )**2
!
               dji(j)=sqrt(dji(j))
               dr1=sqrt(dr1)
               dr2=sqrt(dr2)
               fac1(j)=dr2/(dr1+dr2)
               fac(j) =dr1/(dr1+dr2)
!              fac(j)=1.d0-fac1(j)
!..............Minimum length
               dxmin=min(dxmin,dji(j))
!
!..............Secondary variables
!
               saa(j)=sa(j)*permeability(j)
               dji_a(j)=dji_x(j,1)*xn(j,1)+dji_x(j,2)*xn(j,2)+dji_x(j,3)*xn(j,3)
               sap(j)=saa(j)/dji_a(j)
               sad(j)=saa(j)/dji(j)
!
!..............Modified cell center coordinates for non-orthogonal grid
!
               dr=(xfc(j,1)-xloc(i,1))*xn(j,1)+(xfc(j,2)-xloc(i,2))*xn(j,2)+(xfc(j,3)-xloc(i,3))*xn(j,3)
               xloc_m(j,1)=xfc(j,1)-dr*xn(j,1)-xloc(i,1)
               xloc_m(j,2)=xfc(j,2)-dr*xn(j,2)-xloc(i,2)
               xloc_m(j,3)=xfc(j,3)-dr*xn(j,3)-xloc(i,3)
!..............Why do this compute ??????
               dr2=(xfc(j,1)-xloc(i,1))**2+(xfc(j,2)-xloc(i,2))**2+(xfc(j,3)-xloc(i,3))**2
               dr2=sqrt(dr2)
               dr2=abs(dr/dr2)
               dr2=max(dr2,0.d0)
               dr2=min(dr2,1.d0)
               theta=DACOSD(dr2)
               theta_max=max(theta_max,theta)
             ENDDO
            xfc_min(i,1)=xmin
            xfc_max(i,1)=xmax
            xfc_min(i,2)=ymin
            xfc_max(i,2)=ymax
            xfc_min(i,3)=zmin
            xfc_max(i,3)=zmax
            xloc_xfc_min(i,1)=MIN(xmin,xloc(i,1))
            xloc_xfc_max(i,1)=MAX(xmax,xloc(i,1))
            xloc_xfc_min(i,2)=MIN(ymin,xloc(i,2))
            xloc_xfc_max(i,2)=MAX(ymax,xloc(i,2))
            xloc_xfc_min(i,3)=MIN(zmin,xloc(i,3))
            xloc_xfc_max(i,3)=MAX(zmax,xloc(i,3))
            xloc_xfc_radius_max(i)=radiusmax
            xloc_xfc_radius_min(i)=radiusmin
         ENDDO
      ENDIF
!            
!.....Define index for pressure gradient calculation
!
      DO i=1,ncell_fluid
!DIR$ NOVECTOR
         DO j=i_neigh(i),i_neigh(i+1)-1
            IF(nbcon(j).gt.nin_max.and.nbcon(j).ne.nb_sym.and.nbcon(j).le.nb_max)THEN
!
!..............Pressure boundary
!
               ngrad(i)=2         
            ELSEIF(ngrad(i).ne.2.and.nbcon(j).ne.0.and.nbcon(j).ne.nb_sym.and.nbcon(j).lt.nb_mars)THEN
!
!..............Except pressure, cell, symmetric boundary
!  
               ngrad(i)=1
           ENDIF
         ENDDO 
      ENDDO
!
!.....Cell's physical length
!
      DO ix=1,ndim
         DO i=1,ncell_fluid
            xmin= 1.d10
            xmax=-1.d10
!DIR$ NOVECTOR
            DO j=i_neigh(i),i_neigh(i+1)-1
                svp(j,ix)=sv(j,ix)*permeability(j)
                dnj(j,ix)=xn(j,ix)-dji_x(j,ix)/dji(j)
                xmax=max(xmax,xfc(j,ix))
                xmin=min(xmin,xfc(j,ix))
            ENDDO 
            cell_leng(i,ix)=xmax-xmin
         ENDDO
      ENDDO
!
!.....Make fac and fac1 symmetric
!
!     DO i=1,ncell_fluid
!!DIR$ NOVECTOR
!        DO j=i_neigh(i),i_neigh(i+1)-1
!           IF(nbcon(j).eq.0)THEN
!              k=neigh(j)
!              IF(k.lt.i)THEN
!                 jk=nji(j)
!                 k0=i_neigh(k)-1
!                 fac(jk+k0)=fac1(j)
!                 fac1(jk+k0)=fac(j)
!              ENDIF
!           ENDIF
!        ENDDO
!     ENDDO
!
      IF(np.gt.1) THEN
         CALL communicate_1d(volp)
         CALL communicate_1d_csr(fac,i_neigh)
         CALL communicate_1d_csr(fac1,i_neigh)
         CALL communicate_2d_csr(xloc_m,i_neigh)
!........No need to communicate xn, if not activate here        
!        CALL communicate_2d_csr(xn,i_neigh)
!
         CALL allreducei_min_r1(dxmin)
      ENDIF
!
      END SUBROUTINE calc_geo
!
!*************************************************************************************
!*************************************************************************************
!
      SUBROUTINE calc_walln
!
      USE Zzone       , ONLY: ncell_fluid
      USE Zparam      , ONLY: ndim,nin_max,nb_max,nb_sym
      USE Znum_cell   , ONLY: i_neigh
      USE Zbc_index   , ONLY: nbcon,icell_type,iface_wall
      USE Zcoord2     , ONLY: fac,fac1
      USE Zcoord4     , ONLY: dji,dji_x
      USE Znormal     , ONLY: xn
      USE Zturb       , ONLY: walln,walln2,wallnr
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER i,j,j0,j1
      INTEGER itype
      REAL*8  a,awalln
!
      IF(ndim.eq.2) THEN
         DO i=1,ncell_fluid
            j0=iface_wall(i)
            itype=icell_type(i)
            IF(j0.ne.0) THEN
               j1=j0+i_neigh(i)-1
               a=dji_x(j1,1)*xn(j1,1)+dji_x(j1,2)*xn(j1,2)
               awalln   =fac (j1)*a
               walln2(i)=fac1(j1)*a
!DIR$ NOVECTOR
               DO j=i_neigh(i),i_neigh(i+1)-1
                  IF(j.ne.j1) THEN
                     IF((itype.eq.1.and. nbcon(j).lt.0)                               .or. &
                        (itype.eq.2.and.(nbcon(j).ge.1.and.nbcon(j).le.nin_max))      .or. &
                        (itype.eq.3.and.(nbcon(j).gt.nin_max.and.nbcon(j).le.nb_max)) .or. &
                        (itype.eq.4.and. nbcon(j).eq.nb_sym)) THEN
                        awalln=MIN(dji(j),awalln)
                     ENDIF
                  ENDIF
               ENDDO
               walln(i)=awalln
               wallnr(i)=1.d0/awalln
            ELSE
               walln(i) =0.d0
               walln2(i)=0.d0
            ENDIF
         ENDDO
      ELSE
         DO i=1,ncell_fluid
            itype=icell_type(i)
            j0=iface_wall(i)
            IF(j0.ne.0) THEN
               j1=j0+i_neigh(i)-1
               a=dji_x(j1,1)*xn(j1,1)+dji_x(j1,2)*xn(j1,2)+dji_x(j1,3)*xn(j1,3)
               awalln   =fac (j1)*a       
               walln2(i)=fac1(j1)*a
!DIR$ NOVECTOR
               DO j=i_neigh(i),i_neigh(i+1)-1
                  IF(j.ne.j1) THEN
                     IF((itype.eq.1.and. nbcon(j).lt.0)                                 .or. &
                        (itype.eq.2.and.(nbcon(j).ge.1.and.nbcon(j).le.nin_max))      .or. &
                        (itype.eq.3.and.(nbcon(j).gt.nin_max.and.nbcon(j).le.nb_max)) .or. &
                        (itype.eq.4.and. nbcon(j).eq.nb_sym)) THEN
                        awalln=MIN(dji(j),awalln)
                     ENDIF
                  ENDIF                
               ENDDO
               walln(i)=awalln
               wallnr(i)=1.d0/awalln
            ELSE
               walln(i) =0.d0
               walln2(i)=0.d0
            ENDIF
         ENDDO
      ENDIF
!
      END SUBROUTINE calc_walln
!
