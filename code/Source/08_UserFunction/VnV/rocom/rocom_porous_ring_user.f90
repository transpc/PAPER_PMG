!     
      SUBROUTINE rocom_porous_ring_user
!
!     Set the porous into rocom mesh
!
      USE Zmpi         , ONLY: maxmt_fluid
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: myrank       
      USE Zparam       , ONLY: ndim
      USE Zcoord1      , ONLY: xloc
      USE Zcoord3      , ONLY: floss,porosity,vol,volpr
      USE Znum_cell    , ONLY: i_neigh
      USE Zvec_geo     , ONLY: xn_nf,sv_nf,svp_nf,xfc_nf,  &
                               djia_nf,sap_nf,sa_nf,saa_nf
      USE Zio_unit     , ONLY: unit_log
!
      IMPLICIT NONE 
!
!.....External functions
      REAL(8) cosine_func,sine_func
!.....Local variables
      INTEGER i,j,j0,i1
      INTEGER level,ropt,hopt
      REAL(8) heightn,heightmax,heightmin,radiusn
      REAL(8) rdiffer,rdiffermin,radiusnmin,radiusnmax
      REAL(8) a_drum,AA_drum,K_drum                                   !area of hole, area of orifice, foam loss factor
      REAL(8) f_pipe,n_pipe,a_core,a_pipe_side_1m                     !friction factor of 1 pipe per 1m, number of pipe, cross-area of core, area of 1m pipe side
!      
!.....Local arrays
      REAL(8) :: ring(9,4)
      REAL(8) :: radiusi(ncell_fluid)
      REAL(8) :: radius(maxmt_fluid)
      REAL(8) :: xfc(maxmt_fluid,ndim)
! 
      DATA a_drum,AA_drum,K_drum/1.256e-5,3.329d-5,9.961d0/                          !Dr.yoon_cherl-2.5 for ke model
      DATA f_pipe,n_pipe,a_core,a_pipe_side_1m/0.027d0, 193.d0, 0.4776d0, 0.02355d0/ !30mm pipe * #193, 0.437m !Dr. Jim.j.w.
!      
      IF(myrank.eq.0)then
         WRITE(*,*) '           Porosity for ROCOM in Rocom_porous_ring2.f90!'
         WRITE(unit_log,*)'           Porosity for ROCOM in Rocom_porous_ring2.f90!'
      ENDIF
!         
!.....Mixing drum
!
      ring(1,1)=0.092d0                                       !height_low,height_high,radius_small,radius_large,loss_coeff,porosity,permeability
      ring(2,1)=0.092d0+0.095d0+0.113d0                       !height_high
      ring(3,1)=0.290d0                                       !radius_small
      ring(4,1)=0.290d0                                       !radius_large
      ring(5,1)=a_drum/AA_drum*K_drum/2.0d0                   !foam loss constant  fx=a_hole/A_orIFice*K/2 *rho*u**2.*A_cell*cos_theta
      ring(6,1)=1.0d0                                         !porosity
      ring(7,1)=1.0d0                                         !x-permeability
      ring(8,1)=1.0d0                                         !y-permeability
      ring(9,1)=1.0d0                                         !z-permeability
!       
!.....Lower plate
!
      ring(1,2)=0.092d0+0.095d0+0.113d0                       !height_low,height_high,radius_small,radius_large,loss_coeff
      ring(2,2)=0.092d0+0.095d0+0.113d0+0.010d0               !height_low
      ring(3,2)=0.0d0                                         !radius_small
      ring(4,2)=0.338d0                                       !radius_large
      ring(5,2)=1.0d0                                         !foam loss costant
      ring(6,2)=(1.0d0-0.66d0)                                !porosity 
      ring(7,2)=1.0d0                                         !x-permeability
      ring(8,2)=1.0d0                                         !y-permeability
      ring(9,2)=1.0d0                                         !z-permeability
!
!.....Pipe
!
      ring(1,3)=0.092d0+0.095d0+0.113d0+0.010d0               !height_low,height_high,radius_small,radius_large,loss_coeff,porosity,permeability
      ring(2,3)=0.092d0+0.095d0+0.113d0+0.010d0+1.004d0       !height_high
      ring(3,3)=0.0d0                                         !radius_small
      ring(4,3)=0.338d0                                       !radius_large
      ring(5,3)=f_pipe/2.0d0*a_pipe_side_1m*n_pipe/a_core     !friction factor
      ring(6,3)=0.286d0                                       !0.1364/0.4776
      ring(7,3)=0.0d0                                         !x-permeability
      ring(8,3)=0.0d0                                         !y-permeability
      ring(9,3)=1.0d0                                         !z-permeability
!
!.....Upper plate
!
      ring(1,4)=0.092d0+0.095d0+0.113d0+0.010d0+1.004d0       !height_low,height_high,radius_small,radius_large,loss_coeff
      ring(2,4)=0.092d0+0.095d0+0.113d0+0.010d0+1.004d0 
      ring(3,4)=0.0d0
      ring(4,4)=0.338d0
      ring(5,4)=0.5d0                                         !upper support
      ring(6,4)=(1.0d0-0.66d0) 
      ring(7,4)=1.0d0                                         !x-permeability
      ring(8,4)=1.0d0                                         !y-permeability
      ring(9,4)=1.0d0                                         !z-permeability
!
!.....Get scalar variable xfc from vector xfc_nf
!
      CALL get_scalar_variable_p_ndim(xfc_nf,xfc)
!
!.....mixing drum: find best-estimated radius of the node, revise ring_radial_position from cell to face
!
      IF(ndim.eq.2) THEN
         DO i=1,ncell_fluid
            DO j=i_neigh(i),i_neigh(i+1)-1
               radius(j)=ABS(xfc(j,1)) 
            ENDDO
         ENDDO
      ELSE
         DO i=1,ncell_fluid
            DO j=i_neigh(i),i_neigh(i+1)-1
               radius(j)=SQRT(xfc(j,1)**2+xfc(j,2)**2)             
            ENDDO
         ENDDO
      ENDIF
!
      radiusnmin= 1000.d0
      radiusnmax=-1000.d0
      rdiffermin= 1000.d0
      DO i=1,ncell_fluid
         DO j=i_neigh(i),i_neigh(i+1)-1
            radiusn=radius(j)
            radiusnmin=MIN(radiusn,radiusnmin) 
            radiusnmax=MAX(radiusn,radiusnmax) 
            rdiffer=ring(3,1)-radiusn
            IF(rdiffer.ge.0.d0) rdiffermin=MIN(rdiffer,rdiffermin)
         ENDDO
      ENDDO
      ring(3,1)=ring(3,1)-rdiffermin
!      
!.....find the neighbor whose radius is the same as revised ring(3,1)
!
      IF(ndim.eq.2) THEN
         DO i=1,ncell_fluid
            floss(i,1)=0.d0
            floss(i,2)=0.d0
         ENDDO
      ELSE
         DO i=1,ncell_fluid
            floss(i,1)=0.d0
            floss(i,2)=0.d0
            floss(i,3)=0.d0
         ENDDO
      ENDIF
!
      level=4
!
      DO i=1,ncell_fluid
         j0=i_neigh(i)-1
         DO j=i_neigh(i),i_neigh(i+1)-1
            ropt=0
            hopt=0
            heightn=ABS(xfc(j,ndim))
            radiusn=radius(j)
            IF(heightn.ge.ring(1,1).and.heightn.le.ring(2,1)) hopt=1
            IF(ABS(radiusn-ring(3,1)).le.1.0d-3) ropt=1
!
            IF(hopt*ropt.eq.1)THEN
               CALL get_vector_disp(j-j0,i,i1)
               IF(i1.gt.0) THEN
                  floss(i,1)=a_drum/AA_drum*K_drum/2.0d0*sa_nf(i1)*cosine_func(xloc(i,1),xloc(i,2)) !*rho*u**2
                  floss(i,2)=a_drum/AA_drum*K_drum/2.0d0*sa_nf(i1)*sine_func  (xloc(i,1),xloc(i,2)) !*rho*u**2
                  floss(i,ndim)=0.0d0
                  porosity(i)=ring(6,1)
                  svp_nf(i1,1)=sv_nf(i1,1)*ring(7,level)
                  svp_nf(i1,2)=sv_nf(i1,2)*ring(8,level)
                  svp_nf(i1,3)=sv_nf(i1,3)*ring(ndim,level)
                  saa_nf(i1)=sa_nf(i1)*(ring(7,level)*ABS(xn_nf(i1,1))+ring(8,level)*DABS(xn_nf(i1,2))+ring(9,level)*DABS(xn_nf(i1,3)))
                  sap_nf(i1)=saa_nf(i1)/djia_nf(i1)        
               ENDIF
            ENDIF
         ENDDO
      ENDDO
!
!.....Pipe  
!
      level=3                      
      IF(ndim.eq.2) THEN
         DO i=1,ncell_fluid
            radiusi(i)=ABS(xloc(i,1)) 
         ENDDO
      ELSE
         DO i=1,ncell_fluid
            radiusi(i)=SQRT(xloc(i,1)**2+xloc(i,2)**2)
         ENDDO
      ENDIF
!
      DO i = 1, ncell_fluid
         heightn=ABS(xloc(i,ndim))
         radiusn=radiusi(i)
!         
!........Check the radius from ith cell to the center 
!
         ropt=0
         IF(ring(3,level).eq.ring(4,level))THEN                 !for a plate ring
            IF(ring(3,level).ge.radiusnmin.and.ring(3,level).le.radiusnmax) ropt=1
         ELSE !for the thick region
            IF(radiusn.ge.ring(3,level).and.radiusn.le.ring(4,level)) ropt=1
         ENDIF
!         
!........Check the height from ith cell to the bottom or y=0
!
         hopt=0
         IF(ring(1,level).eq.ring(2,level))THEN !for a disc 
            IF(ring(1,level).ge.heightmin.and.ring(1,level).le.heightmax) hopt=1
         ELSE  !for the thick region
            IF(heightn.ge.ring(1,level).and.heightn.le.ring(2,level)) hopt=1
         ENDIF
!
!........Set the geometrical information
!
         IF(ropt*hopt.eq.1)THEN
            floss(i,1)=0.0d0
            floss(i,2)=0.0d0
            IF(ndim.eq.3) floss(i,ndim)=f_pipe/2.0d0*a_pipe_side_1m*vol(i)*n_pipe/a_core
            porosity(i)=ring(6,level)
            j0=i_neigh(i)-1
            DO j=i_neigh(i),i_neigh(i+1)-1
               CALL get_vector_disp(j-j0,i,i1)
               IF(i1.gt.0) THEN
                  svp_nf(i1,1)=sv_nf(i1,1)*ring(7,level)
                  svp_nf(i1,2)=sv_nf(i1,2)*ring(8,level)
                  svp_nf(i1,3)=sv_nf(i1,3)*ring(ndim,level)
                  saa_nf(i1)=sa_nf(i1)*(ring(7,level)*ABS(xn_nf(i1,1))+ring(8,level)*dabs(xn_nf(i1,2))+ring(9,level)*dabs(xn_nf(i1,3)))
                  sap_nf(i1)=saa_nf(i1)/djia_nf(i1)        
               ENDIF
            ENDDO
         ENDIF
      ENDDO 
!
      IF(ndim.eq.2) THEN
         DO i=1,ncell_fluid
            floss(i,1)=floss(i,1)*volpr(i)
            floss(i,2)=floss(i,2)*volpr(i)
         ENDDO
      ELSE
         DO i=1,ncell_fluid
            floss(i,1)=floss(i,1)*volpr(i)
            floss(i,2)=floss(i,2)*volpr(i)
            floss(i,3)=floss(i,3)*volpr(i)
         ENDDO
      ENDIF
!         
      END SUBROUTINE rocom_porous_ring_user
!     
!----------------------------------------------------------------------
! 
      FUNCTION cosine_func(x,y)
      REAL*8 cosine_func,x,y,r
      r=SQRT(x**2+y**2)
      cosine_func=x/r
      END FUNCTION cosine_func
!      
      FUNCTION sine_func(x,y)
      REAL*8 sine_func,x,y,r
      r=SQRT(x**2+y**2)
      sine_func=y/r
      END FUNCTION sine_func
