!
      SUBROUTINE find_walls !(nbcon_all,xloc_all,xfc_all)
!
!     This routine find the group number of each wall according to the its nbcon/xn/xloc. Also find the closest wall and its distance 
!
      USE VOL_DATA            
      USE Zzone       , ONLY: ncell_fluid
      USE Zparam      , ONLY: ndim  
      USE Zbc_index   , ONLY: num_wallcells,wall_cell,num_wall_group
      USE Zcoord1     , ONLY: xloc      
      USE Zcoord2     , ONLY: xfc_wallcell
      USE Zndforce    , ONLY: d_bfc,c_bface,c_bface_indx
      USE Znormal     , ONLY: xn_wallcell
!
      IMPLICIT NONE
!
      INTEGER i,k,m,n
      INTEGER iblk
      PARAMETER(iblk=2000)
!            
      REAL(8) distance_m
      REAL(8) dr1,dr2,dr3
      REAL(8) dr10,dr20,dr30
      REAL(8) dr11,dr21,dr31
      REAL(8) dr12,dr22,dr32
      REAL(8) dr13,dr23,dr33
      REAL(8) wa1,wa2,wa3
      REAL(8) distance0,distance1,distance2,distance3
      REAL(8) distance_m0,distance_m1,distance_m2,distance_m3
      INTEGER m0,m1,m2,m3
      INTEGER n0
      REAL(8) distance_mi(ncell_fluid)
      INTEGER mi(ncell_fluid)
      INTEGER ns0,ne0
!
      !INTEGER, ALLOCATABLE::face_wall_group_all(:,:),cell_closewall_all(:,:)
!
!.....if num_wall_group=0, use the d_bfc which is the closest wall distance in all walls
!
      IF(num_wall_group.eq.0)THEN
!        d_bfc(:)=0.0d0
!!!!!!!!!!
         IF(ndim.eq.2) THEN
!
         DO i=1,ncell_fluid
            n=1
            m=n
            wa1=xfc_wallcell(n,1)
            wa2=xfc_wallcell(n,2)
            dr1=wa1-xloc(i,1)
            dr2=wa2-xloc(i,2)
            distance_mi(i)=dr1**2+dr2**2
            mi(i)=m
         ENDDO   
         DO ns0=1,num_wallcells,iblk
            ne0=min(ns0+iblk-1,num_wallcells)
         DO i=1,ncell_fluid-3,4
            distance_m0=distance_mi(i  )
            distance_m1=distance_mi(i+1)
            distance_m2=distance_mi(i+2)
            distance_m3=distance_mi(i+3)
            m0=mi(i  )
            m1=mi(i+1)
            m2=mi(i+2)
            m3=mi(i+3)
            DO n=ns0,ne0
               wa1=xfc_wallcell(n,1)
               wa2=xfc_wallcell(n,2)
               dr10=wa1-xloc(i  ,1)
               dr20=wa2-xloc(i  ,2)
               dr11=wa1-xloc(i+1,1)
               dr21=wa2-xloc(i+1,2)
               dr12=wa1-xloc(i+2,1)
               dr22=wa2-xloc(i+2,2)
               dr13=wa1-xloc(i+3,1)
               dr23=wa2-xloc(i+3,2)
               distance0=dr10**2+dr20**2
               distance1=dr11**2+dr21**2
               distance2=dr12**2+dr22**2
               distance3=dr13**2+dr23**2
               IF(distance0.le.distance_m0) THEN
                  distance_m0=distance0
                  m0=n
               ENDIF
               IF(distance1.le.distance_m1) THEN
                  distance_m1=distance1
                  m1=n
               ENDIF
               IF(distance2.le.distance_m2) THEN
                  distance_m2=distance2
                  m2=n
               ENDIF
               IF(distance3.le.distance_m3) THEN
                  distance_m3=distance3
                  m3=n
               ENDIF
            ENDDO
            distance_mi(i  )=distance_m0
            distance_mi(i+1)=distance_m1
            distance_mi(i+2)=distance_m2
            distance_mi(i+3)=distance_m3
            mi(i  )=m0
            mi(i+1)=m1
            mi(i+2)=m2
            mi(i+3)=m3
         ENDDO   
         IF(mod(ncell_fluid,4).eq.1) THEN
            i=ncell_fluid
            distance_m0=distance_mi(i  )
            m0=mi(i  )
            DO n=ns0,ne0
               wa1=xfc_wallcell(n,1)
               wa2=xfc_wallcell(n,2)
               dr10=wa1-xloc(i  ,1)
               dr20=wa2-xloc(i  ,2)
               distance0=dr10**2+dr20**2
               IF(distance0.le.distance_m0) THEN
                  distance_m0=distance0
                  m0=n
               ENDIF
            ENDDO
            distance_mi(i  )=distance_m0
            mi(i  )=m0
         ELSEIF(mod(ncell_fluid,4).eq.2) THEN
            i=ncell_fluid-1
            distance_m0=distance_mi(i  )
            distance_m1=distance_mi(i+1)
            m0=mi(i  )
            m1=mi(i+1)
            DO n=ns0,ne0
               wa1=xfc_wallcell(n,1)
               wa2=xfc_wallcell(n,2)
               dr10=wa1-xloc(i  ,1)
               dr20=wa2-xloc(i  ,2)
               dr11=wa1-xloc(i+1,1)
               dr21=wa2-xloc(i+1,2)
               distance0=dr10**2+dr20**2
               distance1=dr11**2+dr21**2
               IF(distance0.le.distance_m0) THEN
                  distance_m0=distance0
                  m0=n
               ENDIF
               IF(distance1.le.distance_m1) THEN
                  distance_m1=distance1
                  m1=n
               ENDIF
            ENDDO
            distance_mi(i  )=distance_m0
            distance_mi(i+1)=distance_m1
            mi(i  )=m0
            mi(i+1)=m1
         ELSEIF(mod(ncell_fluid,4).eq.3) THEN
            i=ncell_fluid-2
            distance_m0=distance_mi(i  )
            distance_m1=distance_mi(i+1)
            distance_m2=distance_mi(i+2)
            m0=mi(i  )
            m1=mi(i+1)
            m2=mi(i+2)
            DO n=ns0,ne0
               wa1=xfc_wallcell(n,1)
               wa2=xfc_wallcell(n,2)
               dr10=wa1-xloc(i  ,1)
               dr20=wa2-xloc(i  ,2)
               dr11=wa1-xloc(i+1,1)
               dr21=wa2-xloc(i+1,2)
               dr12=wa1-xloc(i+2,1)
               dr22=wa2-xloc(i+2,2)
               distance0=dr10**2+dr20**2
               distance1=dr11**2+dr21**2
               distance2=dr12**2+dr22**2
               IF(distance0.le.distance_m0) THEN
                  distance_m0=distance0
                  m0=n
               ENDIF
               IF(distance1.le.distance_m1) THEN
                  distance_m1=distance1
                  m1=n
               ENDIF
               IF(distance2.le.distance_m2) THEN
                  distance_m2=distance2
                  m2=n
               ENDIF
            ENDDO
            distance_mi(i  )=distance_m0
            distance_mi(i+1)=distance_m1
            distance_mi(i+2)=distance_m2
            mi(i  )=m0
            mi(i+1)=m1
            mi(i+2)=m2
         ENDIF
         ENDDO   
!
         DO i=1,ncell_fluid
            distance_m=distance_mi(i)
            m=mi(i)
            k=wall_cell(m)
            c_bface(i)=k                                          ! Cell number that have the closest wall, k is global cell number
            c_bface_indx(i)=m                                          ! Cell number that have the closest wall, k is global cell number
            dr1=xfc_wallcell(m,1)-xloc(i,1)
            dr2=xfc_wallcell(m,2)-xloc(i,2)
            d_bfc(i)=DABS(xn_wallcell(m,1)*dr1+xn_wallcell(m,2)*dr2)
         ENDDO   
!
         ELSE ! ndim
!
         DO i=1,ncell_fluid
            n=1
            m=n
            wa1=xfc_wallcell(n,1)
            wa2=xfc_wallcell(n,2)
            wa3=xfc_wallcell(n,3)
            dr1=wa1-xloc(i,1)
            dr2=wa2-xloc(i,2)
            dr3=wa3-xloc(i,3)
            distance_mi(i)=dr1**2+dr2**2+dr3**2
            mi(i)=m
         ENDDO   
         DO ns0=1,num_wallcells,iblk
            ne0=min(ns0+iblk-1,num_wallcells)
         DO i=1,ncell_fluid-3,4
            distance_m0=distance_mi(i  )
            distance_m1=distance_mi(i+1)
            distance_m2=distance_mi(i+2)
            distance_m3=distance_mi(i+3)
            m0=mi(i  )
            m1=mi(i+1)
            m2=mi(i+2)
            m3=mi(i+3)
            DO n=ns0,ne0
               wa1=xfc_wallcell(n,1)
               wa2=xfc_wallcell(n,2)
               wa3=xfc_wallcell(n,3)
               dr10=wa1-xloc(i  ,1)
               dr20=wa2-xloc(i  ,2)
               dr30=wa3-xloc(i  ,3)
               dr11=wa1-xloc(i+1,1)
               dr21=wa2-xloc(i+1,2)
               dr31=wa3-xloc(i+1,3)
               dr12=wa1-xloc(i+2,1)
               dr22=wa2-xloc(i+2,2)
               dr32=wa3-xloc(i+2,3)
               dr13=wa1-xloc(i+3,1)
               dr23=wa2-xloc(i+3,2)
               dr33=wa3-xloc(i+3,3)
               distance0=dr10**2+dr20**2+dr30**2
               distance1=dr11**2+dr21**2+dr31**2
               distance2=dr12**2+dr22**2+dr32**2
               distance3=dr13**2+dr23**2+dr33**2
               IF(distance0.le.distance_m0) THEN
                  distance_m0=distance0
                  m0=n
               ENDIF
               IF(distance1.le.distance_m1) THEN
                  distance_m1=distance1
                  m1=n
               ENDIF
               IF(distance2.le.distance_m2) THEN
                  distance_m2=distance2
                  m2=n
               ENDIF
               IF(distance3.le.distance_m3) THEN
                  distance_m3=distance3
                  m3=n
               ENDIF
            ENDDO
            distance_mi(i  )=distance_m0
            distance_mi(i+1)=distance_m1
            distance_mi(i+2)=distance_m2
            distance_mi(i+3)=distance_m3
            mi(i  )=m0
            mi(i+1)=m1
            mi(i+2)=m2
            mi(i+3)=m3
         ENDDO   
         IF(mod(ncell_fluid,4).eq.1) THEN
            i=ncell_fluid
            distance_m0=distance_mi(i  )
            m0=mi(i  )
            DO n=ns0,ne0
               wa1=xfc_wallcell(n,1)
               wa2=xfc_wallcell(n,2)
               wa3=xfc_wallcell(n,3)
               dr10=wa1-xloc(i  ,1)
               dr20=wa2-xloc(i  ,2)
               dr30=wa3-xloc(i  ,3)
               distance0=dr10**2+dr20**2+dr30**2
               IF(distance0.le.distance_m0) THEN
                  distance_m0=distance0
                  m0=n
               ENDIF
            ENDDO
            distance_mi(i  )=distance_m0
            mi(i  )=m0
         ELSEIF(mod(ncell_fluid,4).eq.2) THEN
            i=ncell_fluid-1
            distance_m0=distance_mi(i  )
            distance_m1=distance_mi(i+1)
            m0=mi(i  )
            m1=mi(i+1)
            DO n=ns0,ne0
               wa1=xfc_wallcell(n,1)
               wa2=xfc_wallcell(n,2)
               wa3=xfc_wallcell(n,3)
               dr10=wa1-xloc(i  ,1)
               dr20=wa2-xloc(i  ,2)
               dr30=wa3-xloc(i  ,3)
               dr11=wa1-xloc(i+1,1)
               dr21=wa2-xloc(i+1,2)
               dr31=wa3-xloc(i+1,3)
               distance0=dr10**2+dr20**2+dr30**2
               distance1=dr11**2+dr21**2+dr31**2
               IF(distance0.le.distance_m0) THEN
                  distance_m0=distance0
                  m0=n
               ENDIF
               IF(distance1.le.distance_m1) THEN
                  distance_m1=distance1
                  m1=n
               ENDIF
            ENDDO
            distance_mi(i  )=distance_m0
            distance_mi(i+1)=distance_m1
            mi(i  )=m0
            mi(i+1)=m1
         ELSEIF(mod(ncell_fluid,4).eq.3) THEN
            i=ncell_fluid-2
            distance_m0=distance_mi(i  )
            distance_m1=distance_mi(i+1)
            distance_m2=distance_mi(i+2)
            m0=mi(i  )
            m1=mi(i+1)
            m2=mi(i+2)
            DO n=ns0,ne0
               wa1=xfc_wallcell(n,1)
               wa2=xfc_wallcell(n,2)
               wa3=xfc_wallcell(n,3)
               dr10=wa1-xloc(i  ,1)
               dr20=wa2-xloc(i  ,2)
               dr30=wa3-xloc(i  ,3)
               dr11=wa1-xloc(i+1,1)
               dr21=wa2-xloc(i+1,2)
               dr31=wa3-xloc(i+1,3)
               dr12=wa1-xloc(i+2,1)
               dr22=wa2-xloc(i+2,2)
               dr32=wa3-xloc(i+2,3)
               distance0=dr10**2+dr20**2+dr30**2
               distance1=dr11**2+dr21**2+dr31**2
               distance2=dr12**2+dr22**2+dr32**2
               IF(distance0.le.distance_m0) THEN
                  distance_m0=distance0
                  m0=n
               ENDIF
               IF(distance1.le.distance_m1) THEN
                  distance_m1=distance1
                  m1=n
               ENDIF
               IF(distance2.le.distance_m2) THEN
                  distance_m2=distance2
                  m2=n
               ENDIF
            ENDDO
            distance_mi(i  )=distance_m0
            distance_mi(i+1)=distance_m1
            distance_mi(i+2)=distance_m2
            mi(i  )=m0
            mi(i+1)=m1
            mi(i+2)=m2
         ENDIF
         ENDDO   
!
         DO i=1,ncell_fluid
            distance_m=distance_mi(i)
            m=mi(i)
            k=wall_cell(m)
            c_bface(i)=k                                          ! Cell number that have the closest wall, k is global cell number
            c_bface_indx(i)=m                                          ! Cell number that have the closest wall, k is global cell number
            dr1=xfc_wallcell(m,1)-xloc(i,1)
            dr2=xfc_wallcell(m,2)-xloc(i,2)
            dr3=xfc_wallcell(m,3)-xloc(i,3)
            d_bfc(i)=DABS(xn_wallcell(m,1)*dr1+xn_wallcell(m,2)*dr2+xn_wallcell(m,3)*dr3)
         ENDDO   
         ENDIF
!
         RETURN
      ENDIF
!!
!!.....Define the group of wall according to the user-defined ranges of xn/xloc/nbcon
!!
!      face_wall_group_all(:,:)=0
!      DO n=1, num_wallcells
!         i=wall_cell(n)                                              ! i is global cell number 
!         DO j=1, num_neigh_tmp(i)
!            IF(nbcon_all(j,i).lt.0)THEN
!               SELECTCASE(select_wall_group)
!               CASE(1)                                               ! Divide the group of wall by surface vector
!!                  IF((xn_all(1,j,i).le.-DSQRT(0.5d0)+1.0d-4.and.xn_all(1,j,i).ge.-1.0d0-1.0d-4) .and. (xn_all(2,j,i).le.0.0d0+1.0d-4.and.xn_all(2,j,i).ge.-DSQRT(0.5d0)-1.0d-4)) face_wall_group_all(1,i)=j        !for group1
!!                  IF((xn_all(1,j,i).ge. DSQRT(0.5d0)-1.0d-4.and.xn_all(1,j,i).le. 1.0d0+1.0d-4) .and. (xn_all(2,j,i).le.0.0d0+1.0d-4.and.xn_all(2,j,i).ge.-DSQRT(0.5d0)-1.0d-4))  face_wall_group_all(2,i)=j     !for group2
!!                  IF(xn_all(2,j,i).eq.-1.0d0.and.xn_all(1,j,i).eq.0.0d0)  face_wall_group_all(3,i)=j     !for group3 
!                  IF((xn_all(1,j,i).le.-DSQRT(0.5d0)+1.0d-4.and.xn_all(1,j,i).ge.-1.0d0-1.0d-4) ) face_wall_group_all(1,i)=j        !for group1
!                  IF((xn_all(1,j,i).ge. DSQRT(0.5d0)-1.0d-4.and.xn_all(1,j,i).le. 1.0d0+1.0d-4) )  face_wall_group_all(2,i)=j     !for group2
!                  IF(xn_all(1,j,i).gt.0.0d0-1.0d-4.and.xn_all(1,j,i).lt.0.0d0+1.0d-4)  face_wall_group_all(3,i)=j     !for group3                   
!!
!               CASE(2)                                               ! Divide the group of wall by xyz coordinate,
!                  IF((xloc_all(1,i).ge.0.0d0.and.xloc_all(1,i).le.0.0075d0) .and. (xloc_all(2,i).ge.0.00155d0)) face_wall_group_all(1,i)=j
!                  IF((xloc_all(1,i).gt.0.0075d0) .and. (xloc_all(2,i).ge.0.0d0.and.xloc_all(2,i).le.0.00155d0)) face_wall_group_all(2,i)=j
!!
!               CASE(3)                                               ! Divide the group of wall by value of nbcon(j,i)
!                  IF(nbcon_all(j,i).eq.-1)face_wall_group_all(1,i)=j
!                  IF(nbcon_all(j,i).eq.-5)face_wall_group_all(2,i)=j
!               END SELECT
!            ENDIF
!         ENDDO        
!      ENDDO
!!
!!.....Calculate the multiple distances to the closest wall for each wall group if num_wall_group is not zero.
!!
!      dis_closewall(:,:)=0.0d0
!      cell_closewall_all(:,:)=0
!      DO m=1, num_wall_group
!         DO i=1,ncell_fluid
!            distance_m=huge(0.d0)
!            DO n=1,num_wallcells
!               ii=wall_cell(n)                                       ! ii is global cell number
!               IF(face_wall_group_all(m,ii).ne.0)THEN
!                  distance=0.0d0
!                  DO ix=1,ndim
!                     dr=xfc_all(ix,face_wall_group_all(m,ii),ii)-xloc(i,ix)
!                     distance=distance+dr*dr
!                     IF(distance.gt.distance_m) CYCLE
!                     IF(ix.eq.ndim)THEN
!                        distance_m=distance
!                        k=ii
                         n0=n
!                     ENDIF
!                  ENDDO
!               ENDIF
!            ENDDO
!!       
!            cell_closewall(m,i)=k                                   ! k is global cell number
!            cell_closewall_indx(m,i)=n0
!            j=face_wall_group_all(m,k)
!            DO ix=1,ndim
!               xbfcv(ix)=xfc_all(ix,j,k)-xloc(i,ix)
!            ENDDO
!            dis_closewall(m,i)=DABS(DOT_PRODUCT(xn_all(:,j,k),xbfcv(:)))
!         ENDDO    
!      ENDDO
!!      
!!.....Save sub-domain values from global values      
!!
!      DO ii=1,ncell_fluid
!         i=jperm(ii)    
!         DO n=1,num_wall_group
!            face_wall_group(n,ii)=face_wall_group_all(n,i)                 
!         ENDDO
!      ENDDO      
!!
!!.....Store d_bfc for turbulence model. Find the minimum value of dis_closewall
!!
!      IF(num_wall_group.gt.0)THEN
!         d_bfc(:)=1.0d4    
!         DO i=1,ncell_fluid   
!            DO m=1, num_wall_group
!               IF(d_bfc(i).gt.dis_closewall(m,i)) d_bfc(i)=dis_closewall(m,i)
!            ENDDO
!            IF(d_bfc(i).ge.1.0d4)THEN
!               WRITE(*,*) 'Minimum distance to the closest wall is large than 10000. Stop in find_walls subroutine.'
!               WRITE(97,*) 'Minimum distance to the closest wall is large than 10000. Stop in find_walls subroutine.'
!               STOP
!            ENDIF
!         ENDDO
!      ENDIF
!!   
!      DEALLOCATE(temp,tempall,xn_all,iface_wall_all)  
!      DEALLOCATE(face_wall_group_all,cell_closewall_all)           
!      
      RETURN
      END SUBROUTINE find_walls
