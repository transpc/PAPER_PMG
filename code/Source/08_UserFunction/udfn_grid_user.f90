!
      SUBROUTINE udfn_grid_user(ncell,xloc_tmp,npb_tmp)
!
!     Define user-defined geometry data for fluidic device
!   
      USE Zmpi       , ONLY: maxmt_cell
      USE Zcore      , ONLY: myrank
      USE Zconst1    , ONLY: vv_prob
      USE Zparam     , ONLY: nn,ndim,pi,nin_max
      USE Znum_cell  , ONLY: i_neigh_tmp,j_neigh_tmp,j_nbcon_tmp, &
                             sv_tmp1,xfc_tmp1
      USE Zio_unit   , ONLY: unit_log
      USE unitManager, ONLY: createUnit      

!      
      IMPLICIT NONE
!
!     input
      INTEGER ncell
      REAL(8) xloc_tmp(nn,ndim)
!     output
      INTEGER npb_tmp(nn)
!     local variables
      INTEGER i,j,k,m,IFd,j1,j2,k1
      INTEGER j0,k0,m0,k10
      REAL(8) sa
      REAL(8) r,rmax
!     local arrays
      INTEGER si_cell(4),si_nbcon(4),si_jth(4)
      REAL(8) si_z(4),si_r(4),si_theta(4),si_x(4),si_y(4),dr(4),dr_min(4),dx(4),dy(4),dz(4),x,y,z,hleg_z
      REAL(8) xn_tmp(maxmt_cell,ndim)
!     apr1400_lbloca      
      INTEGER si_opt,err,runit
      REAL(8) a_offset,si_z0,si_r0,si_ang(4)
!
      DATA a_offset/90/
!
      IF(myrank.eq.0) THEN
         IF(ndim.eq.2) THEN
            DO i=1,ncell
               j0=i_neigh_tmp(i)-1
               DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
                  sa=sv_tmp1(j,1)**2+sv_tmp1(j,2)**2
                  sa=1.d0/DSQRT(sa)
                  xn_tmp(j,1)=sv_tmp1(j,1)*sa
                  xn_tmp(j,2)=sv_tmp1(j,2)*sa
               ENDDO
            ENDDO
         ELSE
            DO i=1,ncell
               j0=i_neigh_tmp(i)-1
               DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
                  sa=sv_tmp1(j,1)**2+sv_tmp1(j,2)**2+sv_tmp1(j,3)**2
                  sa=1.d0/DSQRT(sa)
                  xn_tmp(j,1)=sv_tmp1(j,1)*sa
                  xn_tmp(j,2)=sv_tmp1(j,2)*sa
                  xn_tmp(j,3)=sv_tmp1(j,3)*sa
               ENDDO
            ENDDO
         ENDIF
!
!.....Fluidic_device
!     
         IF(vv_prob.eq.'fluidic_device') THEN
!
!........Change boundary condition to determine flow resistance
!  
            IF(.false.)THEN
               DO i=1,ncell
                  j0=i_neigh_tmp(i)-1
                  IF(xloc_tmp(i,ndim).gt.11.2d0)THEN  ! Inlet boundary for SIT top
                     DO j=1,6
                        IF(xn_tmp(j+j0,ndim).gt.0.5d0) THEN
                           j_nbcon_tmp(j+j0)=1
                        ENDIF  
                     ENDDO
                  ENDIF
!
                  IF(.false.)THEN ! wall boundary for stand pipe top
                     IF(xloc_tmp(i,ndim).gt.4.42d0.and.xloc_tmp(i,ndim).lt.4.58d0)THEN
                        IF(xloc_tmp(i,1).gt.1.04d0.and.xloc_tmp(i,1).lt.1.39d0.and.  &
                           xloc_tmp(i,2).gt.1.0412d0.and.xloc_tmp(i,2).lt.1.39d0)THEN
                           DO j=1,6
                              IF(xn_tmp(j+j0,ndim).gt.0.5d0)THEN
                                 j_nbcon_tmp(j+j0)=-1
                                 k=j_neigh_tmp(j+j0)
                                 k0=i_neigh_tmp(k)-1
                                 DO m=1,6
                                    IF(j_neigh_tmp(m+k0).eq.i) j_nbcon_tmp(m+k0)=-1
                                 ENDDO
                              ENDIF
                           ENDDO
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
!
!........Modify grid to meet the water volume at the top of FD inlet
!
            DO i=1,ncell
               j0=i_neigh_tmp(i)-1
               ifd=0
               IF(xloc_tmp(i,ndim).gt.0.25.and.xloc_tmp(i,ndim).lt.0.32)THEN
                  IF(xloc_tmp(i,1).gt.1.13.and.xloc_tmp(i,1).lt.1.3)THEN
                     IF(xloc_tmp(i,2).gt.0.444.and.xloc_tmp(i,2).lt.0.626) IFd=1
                     IF(xloc_tmp(i,2).gt.1.82.and.xloc_tmp(i,2).lt.1.99) IFd=1
                  ENDIF
                  IF(xloc_tmp(i,2).gt.1.13.and.xloc_tmp(i,2).lt.1.3)THEN
                     IF(xloc_tmp(i,1).gt.0.444.and.xloc_tmp(i,1).lt.0.626) IFd=1
                     IF(xloc_tmp(i,1).gt.1.82.and.xloc_tmp(i,1).lt.1.99) IFd=1
                  ENDIF
               ENDIF
!
               IF(ifd.gt.0)THEN
!
!...............First layer
!
                   DO j=1,6
                      IF(xn_tmp(j+j0,ndim).gt.0.5d0)THEN
                         k=j_neigh_tmp(j+j0)
                         k0=i_neigh_tmp(k)-1
                         DO j1=1,6
                            IF(dabs(xn_tmp(j1+k0,ndim)).lt.0.5)THEN
                               j_nbcon_tmp(j1+k0)=-1
                               m=j_neigh_tmp(j1+k0)
                               m0=i_neigh_tmp(m)-1
                               DO j2=1,6
                                  IF(j_neigh_tmp(j2+m0).eq.k) j_nbcon_tmp(j2+m0)=-1
                               ENDDO
                            ENDIF
                         ENDDO
                      ENDIF
                   ENDDO
!
!...............Second layer
!
                   DO j=1,6
                      k=j_neigh_tmp(j+j0)
                      k0=i_neigh_tmp(k)-1
                      IF(xn_tmp(j+k0,ndim).gt.0.5d0)THEN
                         k1=j_neigh_tmp(j+k0)
                         k10=i_neigh_tmp(k1)-1
                         DO j1=1,6
                            IF(dabs(xn_tmp(j1+k10,ndim)).lt.0.5)THEN
                               j_nbcon_tmp(j1+k10)=-1
                               m=j_neigh_tmp(j1+k10)
                               m0=i_neigh_tmp(m)-1
                               DO j2=1,6
                                  IF(j_neigh_tmp(j2+m0).eq.k1) j_nbcon_tmp(j2+m0)=-1
                               ENDDO
                            ENDIF
                         ENDDO
                      ENDIF
                   ENDDO
!
               ENDIF
!
            ENDDO
!
         ENDIF
!
!.....moving_wall, 2D_boiling, 3D_boiling
!     
         IF(vv_prob.eq.'moving_wall' .or. vv_prob.eq.'2D_boiling'.or.vv_prob.eq.'3D_boiling') THEN
            DO i=1,ncell
               DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
                  IF(j_nbcon_tmp(j).ne.0) j_neigh_tmp(j)=0
               ENDDO
            ENDDO
         ENDIF
!
!........OPR1000 Containment
!     
         IF(vv_prob.eq.'OPR1000_cont') THEN
!
!...........Change boundary condition to inlet for Hot Leg break region
!  
            sa=0.0d0
            DO i=1,ncell
               IF(xloc_tmp(i,1).gt.5.0d0 .and. xloc_tmp(i,1).lt.5.3d0)then
                  IF(xloc_tmp(i,2).gt.-0.94d0 .and. xloc_tmp(i,2).lt.0.94d0)then
                     IF(xloc_tmp(i,3).gt.0.55d0 .and. xloc_tmp(i,3).lt.2.45d0)then
                        j0=i_neigh_tmp(i)-1
                        DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
                           IF(j_nbcon_tmp(j).eq.-1)THEN
                              j_nbcon_tmp(j)=1
                              sa=sa+dsqrt(sv_tmp1(j,1)**2+sv_tmp1(j,2)**2+sv_tmp1(j,3)**2)               
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
            WRITE(* ,*)'Break area is ',sa,' m2'
            WRITE(unit_log,*)'Break area is ',sa,' m2'
         ENDIF
!
!....set 4 si pbouns at apr1400 lbloca
!
         IF(vv_prob.eq.'apr1400_lbloca'   .or.&
            vv_prob.eq.'opr1000_rv_lbloca' )THEN
            runit=createUnit('somaAddition.in')
            runit=813
            OPEN(runit,file='somaAddition.in',status='old',iostat=err)
            IF(err.ne.0.and.myrank.eq.0)THEN
               WRITE(*,"(11x,a)")'somaAddition.in is missing!'
               WRITE(unit_log,"(11x,a)")'somaAddition.in is missing!'
               PAUSE
               STOP
            ELSE
               READ(runit,*)hleg_z,si_z0,si_r0,si_opt
               READ(runit,*)a_offset,si_ang(1),si_ang(2),si_ang(3),si_ang(4)
               CLOSE(runit)
            ENDIF
            IF(si_opt.eq.0)THEN
               RETURN
            ELSEIF(si_opt.eq.1)THEN            
               WRITE(*,"(11x,a)")'4 SIs set to pbouns of 7,8,9,10...'
               si_z(1:4)=si_z0+hleg_z !from the center of hot leg
               si_r(1:4)=si_r0        !from center
               DO k=1,4
                  si_theta(k)=(a_offset+si_ang(k))/180.0d0*pi
               ENDDO                  
               DO k=1,4
                  si_x(k)=si_r(k)*DCOS(si_theta(k))
                  si_y(k)=si_r(k)*DSIN(si_theta(k))         
               ENDDO
            ELSEIF(si_opt.eq.2)THEN 
               si_z(1:4)=2.5d0+0.5d0/2.0d0
               si_x(1)=+1.0d0+0.5d0/2.0d0
               si_x(2)=-1.0d0-0.5d0/2.0d0
               si_x(3)=-1.0d0-0.5d0/2.0d0           
               si_x(4)=+1.0d0+0.5d0/2.0d0           
               si_y(1)=+2.5d0+0.5d0/2.0d0
               si_y(2)=+2.5d0+0.5d0/2.0d0
               si_y(3)=-2.5d0-0.5d0/2.0d0
               si_y(4)=-2.5d0-0.5d0/2.0d0
            ENDIF   
!
            si_nbcon(1)=7
            si_nbcon(2)=8
            si_nbcon(3)=9
            si_nbcon(4)=10
            dr_min(1:4)=100.0d0
            DO i=1,ncell
               x=xloc_tmp(i,1)
               y=xloc_tmp(i,2)
               z=xloc_tmp(i,3)
               DO k=1,4
                  dx(k)=(x-si_x(k))
                  dy(k)=(y-si_y(k))
                  dz(k)=(z-si_z(k))
                  dr(k)=dx(k)*dx(k)+dy(k)*dy(k)+dz(k)*dz(k)
                  dr(k)=dr(k)**0.5d0
                  IF(dr(k).lt.dr_min(k))THEN
                     dr_min(k)=dr(k)
                     si_cell(k)=i
                  ENDIF   
               ENDDO
            ENDDO 
            DO k=1,4
               i=si_cell(k)
               rmax=0.0d0
               DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
                  IF(j_nbcon_tmp(j).eq.-1)THEN
                     r=SQRT(xfc_tmp1(j,1)**2+xfc_tmp1(j,2)**2)
                     IF(r.gt.rmax)THEN
                        rmax=r
                        si_jth(k)=j
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO   
            DO k=1,4
               i=si_cell(k)
               j=si_jth(k)
               IF(j_nbcon_tmp(j).eq.-1)THEN
                  j_nbcon_tmp(j)=si_nbcon(k)
                  npb_tmp(i)=j_nbcon_tmp(j)-nin_max
               ENDIF
            ENDDO
         ENDIF 
!
      ENDIF  !myrank
!
!.....MASTER mapping, ASSEMBLY connectivity (1x1 assembly ONLY)
!
      !original mapping - ljr
      IF(vv_prob.eq.'opr1000_mc_rv'          .or. &
         vv_prob.eq.'opr1000_rv_lbloca'      .or. &         
         vv_prob.eq.'opr1000_rv'                   )then
         CALL GeoData_user_opr1000_1x1(ncell)
      ENDIF
!            
      !new mapping - yhy
      IF(vv_prob.eq.'OPR1000_fullvessel_1x1')then
         CALL GeoData_user_opr1000_1x1(ncell)
      ENDIF
!
!     APR1400 Assembly-wise      
      IF(vv_prob.eq.'apr1400_mc_rv' .or. &
         vv_prob.eq.'apr1400_rv')then
         CALL GeoData_user_apr1400_1x1(ncell)
      ENDIF
!
      IF(vv_prob.eq.'OPR1000_fullcore_modmesh02_rv'        .or. &
         vv_prob.eq.'OPR1000_fullcore_modmesh02_rv_vessel' .or. &
         vv_prob.eq.'OPR1000_single_assem'                       ) then
         CALL GeoData_user_opr1000_rod01_ser(ncell)
         ENDIF
!
!.....KSMR subchannel
!
      IF(vv_prob.eq.'KSMR') THEN
         CALL GeoData_user_KSMR(ncell) 
      ENDIF           
!
      RETURN
      END SUBROUTINE udfn_grid_user
