!
      SUBROUTINE read_openfoam(n_face,n_bface,n_bc,n_zone,                       &
                               npb_tmp,vol_tmp,poro_tmp,nzone_tmp,nmaterial_tmp, &
                               celem,dznull)
!
      USE Zmpi         , ONLY: maxmt_nfluid,maxmt_cell
      USE Zparam       , ONLY: nn,nin_max,nb_max,ndim
      USE Zzone        , ONLY: ncell_fluid_all,ncell_cond_all
      USE Znum_cell    , ONLY: num_neigh_tmp,i_neigh_tmp,j_neigh_tmp,j_nbcon_tmp,     &
                               perm_tmp1,sv_tmp1,xfc_tmp1
      USE Zcoord1      , ONLY: xloc_tmp
      USE Znode        , ONLY: nd,node_face,xnode,neigh_face_tmp1
      USE Zio_unit     , ONLY: unit_grid,unit_log
!
      IMPLICIT NONE 
!
!     input
      INTEGER :: n_face,n_bface,n_bc,n_zone
!     output
      INTEGER :: npb_tmp(nn)
      REAL(8) :: vol_tmp(nn)
      INTEGER :: nzone_tmp(nn),nmaterial_tmp(nn),celem(nn)
      REAL(8) :: poro_tmp(nn),dznull
!     local variables
      INTEGER :: i,j,j0,j1,ix,k,m1,m2,m,n,i1,i2,i3,ii
      INTEGER :: ji,jk
      INTEGER :: n_owner,n_neighbour
      INTEGER :: nullbc_opt
      REAL(8) :: area,area_sum,vol
      REAL(8) :: sa_tmp
      REAL(8) :: vol_tmp0
      REAL(8) :: xfc_tmp10,xfc_tmp20,xfc_tmp30
      REAL(8) :: xc_tmp1,xc_tmp2,xc_tmp3
      REAL(8) :: xloc_tmp1,xloc_tmp2,xloc_tmp3
      REAL(8) :: sv_tmp10,sv_tmp20,sv_tmp30
      REAL(8) :: xn_tmp1,xn_tmp2,xn_tmp3
!     local arrays
      REAL(8) :: sv_t(3),xc(3),xc_tmp(3) 
      REAL(8) :: xnode1(3),xnode2(3),xnode3(3) 
      INTEGER :: nb_condition(n_bc),nb_number(n_bc),nb_start(n_bc)
!      INTEGER :: owner(n_face),neighbour(n_face)
      INTEGER :: nz_index(n_zone),nz_material(n_zone),nz_number(n_zone),nz_start(n_zone),nz_subdomain(n_zone)
!     REAL(8) :: xface(n_face,ndim)
      REAL(8) :: nz_porosity(n_zone),nz_perm(n_zone,ndim)          
!      REAL(8) :: sv_t_all(nn,ndim)
      
! sang test
      INTEGER(4),DIMENSION(:),ALLOCATABLE::owner,neighbour
      REAL(8),DIMENSION(:,:),ALLOCATABLE::sv_t_all
      
      ALLOCATE(sv_t_all(nn,ndim))
      ALLOCATE(owner(n_face),neighbour(n_face))
!
!.....Initialize local variables 
! remove abuse zeroing later
!     celem(:)=0
!
!     xloc_tmp(:,:)=0.d0
!
!
!     nzone_tmp(:)=0
!     nmaterial_tmp(:)=0
!     vol_tmp(:)=0.d0
!     poro_tmp(:)=0.d0
!
!     perm_tmp(:,:)=0.d0
!!!!!!!!!!!!!!
      n_owner=n_face
      n_neighbour=n_face-n_bface
!
!     xface(:,:)=0.0d0
!     owner(:)=0
!     neighbour(:)=0
!     nb_condition(:)=0
!     nb_number(:)=0
!     nb_start(:)=0
!     nz_index(:)=0
!     nz_material(:)=0
!     nz_number(:)=0
!     nz_start(:)=0
!     nz_subdomain(:)=0      
!     nz_porosity(:)=0
!     nz_perm(:,:)=0
      
      DO i=1,n_owner
         READ(unit_grid,*) owner(i)
         owner(i)=owner(i)+1
      ENDDO

      DO i=1,n_neighbour
         READ(unit_grid,*) neighbour(i)
         neighbour(i)=neighbour(i)+1
      ENDDO
!
      DO i=1,n_bc
         READ(unit_grid,*)nb_condition(i),nb_number(i),nb_start(i)
         nb_start(i)=nb_start(i)+1
      ENDDO
!
      DO i=1,n_zone
         READ(unit_grid,*)nz_index(i),nz_material(i),nz_number(i),nz_start(i),nz_subdomain(i),nz_porosity(i),(nz_perm(i,ix),ix=1,ndim)
         nz_start(i)=nz_start(i)+1
      ENDDO           
!
      CLOSE(unit_grid)      
!
      IF(n_zone.gt.0)THEN
         DO k=1,n_zone
            DO i=nz_start(k),nz_start(k)+nz_number(k)-1
               nmaterial_tmp(i)=nz_material(k)
            ENDDO
         ENDDO
      ELSE
         nmaterial_tmp(:)=0
      ENDIF
!
!     get ncell_fluid_all,ncell_cond_all to all tasks
!
      ncell_fluid_all=0
      ncell_cond_all=0
      DO i=1,nn
         IF(nmaterial_tmp(i).le.0) ncell_fluid_all=ncell_fluid_all+1
         IF(nmaterial_tmp(i).ne.0) ncell_cond_all=ncell_cond_all+1
      ENDDO
!
!     Get num_neigh_tmp only first
!
      num_neigh_tmp(:)=0
      IF(ndim.eq.2)THEN
         DO ii=1,n_owner
            nullbc_opt=0
            IF(ii.gt.n_neighbour)THEN
               DO m=1,n_bc
                  IF(nb_condition(m).eq.100) THEN
                     m1=nb_start(m)
                     m2=nb_start(m)+nb_number(m)-1
                     IF(ii.ge.m1.and.ii.le.m2) nullbc_opt=1
                  ENDIF
                ENDDO
            ENDIF
            IF(nullbc_opt.eq.1)CYCLE
            i=owner(ii)
            ji=num_neigh_tmp(i)
            ji=ji+1
            num_neigh_tmp(i)=ji
            IF(ii.le.n_neighbour)THEN
               k=neighbour(ii)
               jk=num_neigh_tmp(k)
               jk=jk+1
               num_neigh_tmp(k)=jk
            ENDIF
         ENDDO
      ELSE
         DO ii=1,n_owner
            i=owner(ii)
            ji=num_neigh_tmp(i)
            ji=ji+1
            num_neigh_tmp(i)=ji
            IF(ii.le.n_neighbour)THEN
               k=neighbour(ii)
               jk=num_neigh_tmp(k)
               jk=jk+1
               num_neigh_tmp(k)=jk
            ENDIF
         ENDDO
      ENDIF
!
!     We need to go csr format
!
!     write(*,*) 'ncell_fluid_all',ncell_fluid_all,ncell_cond_all,ncell,ns
      i_neigh_tmp(1)=1
      DO i=1,ncell_fluid_all
         i_neigh_tmp(i+1)=i_neigh_tmp(i)+num_neigh_tmp(i)
      ENDDO
      maxmt_nfluid=i_neigh_tmp(ncell_fluid_all+1)-1
!     write(*,*) 'maxmt_nfluid',maxmt_nfluid
!           
      maxmt_cell=maxmt_nfluid
      DO i=ncell_fluid_all+1,nn
         i_neigh_tmp(i+1)=i_neigh_tmp(i)+num_neigh_tmp(i)
      ENDDO
      maxmt_cell=i_neigh_tmp(nn+1)-1
!     write(*,*) 'maxmt_cell',maxmt_cell
!
      ALLOCATE(j_neigh_tmp(maxmt_cell),j_nbcon_tmp(maxmt_cell))
      ALLOCATE(perm_tmp1(maxmt_cell),sv_tmp1(maxmt_cell,ndim),xfc_tmp1(maxmt_cell,ndim))
      ALLOCATE(neigh_face_tmp1(maxmt_cell))
      num_neigh_tmp(:)=0
      j_neigh_tmp(:)=0
      j_nbcon_tmp(:)=0
      npb_tmp(:)=0
      sv_tmp1(:,:)=0.d0
      xfc_tmp1(:,:)=0.d0
      IF(ndim.eq.2)THEN
         DO ii=1,n_owner
            nullbc_opt=0
            IF(ii.gt.n_neighbour)THEN
               DO m=1,n_bc
                  IF(nb_condition(m).eq.100) THEN
                     m1=nb_start(m)
                     m2=nb_start(m)+nb_number(m)-1
                     IF(ii.ge.m1.and.ii.le.m2) nullbc_opt=1
                  ENDIF
               ENDDO
            ENDIF
            IF(nullbc_opt.eq.1)CYCLE
            i=owner(ii)
            ji=num_neigh_tmp(i)
            ji=ji+1
            num_neigh_tmp(i)=ji
!
! Calculate surface vectors (Ref. Comp. Methods for Fluid Dyn. pp.240)
!
            j0=i_neigh_tmp(i)-1
            neigh_face_tmp1(ji+j0)=ii
            sv_tmp10=sv_tmp1(ji+j0,1)
            sv_tmp20=sv_tmp1(ji+j0,2)
            DO m=3,nd(ii)
               i1=node_face(1,ii)
               i2=node_face(m-1,ii)
               i3=node_face(m,ii)
               xnode1(:)=xnode(i1,:)
               xnode2(:)=xnode(i2,:)
               xnode3(:)=xnode(i3,:)
!              CALL sv_tri_2d(xnode(i1,:),xnode(i2,:),xnode(i3,:),sv_t)
               CALL sv_tri_2d(xnode1,xnode2,xnode3,sv_t)
               sv_tmp10=sv_tmp10+sv_t(1)/dznull
               sv_tmp20=sv_tmp20+sv_t(2)/dznull
            ENDDO
            sv_tmp1(ji+j0,1)=sv_tmp10
            sv_tmp1(ji+j0,2)=sv_tmp20
!
            xfc_tmp10=xfc_tmp1(ji+j0,1)
            xfc_tmp20=xfc_tmp1(ji+j0,2)
            DO m=1,nd(ii)
               n=node_face(m,ii)
               xfc_tmp10=xfc_tmp10+xnode(n,1)
               xfc_tmp20=xfc_tmp20+xnode(n,2)
            ENDDO
            xfc_tmp1(ji+j0,1)=xfc_tmp10/DBLE(nd(ii))
            xfc_tmp1(ji+j0,2)=xfc_tmp20/DBLE(nd(ii))
!
            IF(ii.le.n_neighbour)THEN
               k=neighbour(ii)
               j_neigh_tmp(ji+j0)=k
               j_nbcon_tmp(ji+j0)=0
               j1=i_neigh_tmp(k)-1
               jk=num_neigh_tmp(k)
               jk=jk+1
               num_neigh_tmp(k)=jk
               j_neigh_tmp(jk+j1)=i
               j_nbcon_tmp(jk+j1)=0
               neigh_face_tmp1(jk+j1)=ii
               sv_tmp1(jk+j1,1)=-sv_tmp1(ji+j0,1)
               sv_tmp1(jk+j1,2)=-sv_tmp1(ji+j0,2)
               xfc_tmp1(jk+j1,1)=xfc_tmp1(ji+j0,1)
               xfc_tmp1(jk+j1,2)=xfc_tmp1(ji+j0,2)
            ELSE
               j_neigh_tmp(ji+j0)=0
               DO m=1,n_bc
                  m1=nb_start(m)
                  m2=nb_start(m)+nb_number(m)-1
                  IF(ii.ge.m1.and.ii.le.m2)then
                     j_nbcon_tmp(ji+j0)=nb_condition(m)
                     IF(nb_condition(m).gt.nin_max.and.nb_condition(m).le.nb_max)THEN 
                        npb_tmp(i)=nb_condition(m)-nin_max
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ELSE
         DO ii=1,n_owner
            i=owner(ii)
            ji=num_neigh_tmp(i)
            ji=ji+1
            num_neigh_tmp(i)=ji
!
! Calculate surface vectors (Ref. Comp. Methods for Fluid Dyn. pp.240)
!
            j0=i_neigh_tmp(i)-1
            neigh_face_tmp1(ji+j0)=ii
            area_sum=0.0d0
            sv_tmp10=sv_tmp1(ji+j0,1)
            sv_tmp20=sv_tmp1(ji+j0,2)
            sv_tmp30=sv_tmp1(ji+j0,3)
            xfc_tmp10=xfc_tmp1(ji+j0,1)
            xfc_tmp20=xfc_tmp1(ji+j0,2)
            xfc_tmp30=xfc_tmp1(ji+j0,3)
            DO m=3,nd(ii)
               i1=node_face(1,ii)
               i2=node_face(m-1,ii)
               i3=node_face(m,ii)
               xnode1(:)=xnode(i1,:)
               xnode2(:)=xnode(i2,:)
               xnode3(:)=xnode(i3,:)
!              CALL sv_tri(xnode(i1,:),xnode(i2,:),xnode(i3,:),sv_t,xc,area)
               CALL sv_tri(xnode1,xnode2,xnode3,sv_t,xc,area)
               sv_tmp10=sv_tmp10+sv_t(1)
               sv_tmp20=sv_tmp20+sv_t(2)
               sv_tmp30=sv_tmp30+sv_t(3)
               xfc_tmp10=xfc_tmp10+xc(1)*area
               xfc_tmp20=xfc_tmp20+xc(2)*area
               xfc_tmp30=xfc_tmp30+xc(3)*area
               area_sum=area_sum+area
            ENDDO
            sv_tmp1(ji+j0,1)=sv_tmp10
            sv_tmp1(ji+j0,2)=sv_tmp20
            sv_tmp1(ji+j0,3)=sv_tmp30
            xfc_tmp1(ji+j0,1)=xfc_tmp10/area_sum
            xfc_tmp1(ji+j0,2)=xfc_tmp20/area_sum
            xfc_tmp1(ji+j0,3)=xfc_tmp30/area_sum
!
            IF(ii.le.n_neighbour)THEN
               k=neighbour(ii)
               j_neigh_tmp(ji+j0)=k
               j_nbcon_tmp(ji+j0)=0
               j1=i_neigh_tmp(k)-1
               jk=num_neigh_tmp(k)
               jk=jk+1
               num_neigh_tmp(k)=jk
               j_neigh_tmp(jk+j1)=i
               j_nbcon_tmp(jk+j1)=0
               neigh_face_tmp1(jk+j1)=ii
               sv_tmp1(jk+j1,1)=-sv_tmp1(ji+j0,1)
               sv_tmp1(jk+j1,2)=-sv_tmp1(ji+j0,2)
               sv_tmp1(jk+j1,3)=-sv_tmp1(ji+j0,3)
               xfc_tmp1(jk+j1,1)=xfc_tmp1(ji+j0,1)
               xfc_tmp1(jk+j1,2)=xfc_tmp1(ji+j0,2)
               xfc_tmp1(jk+j1,3)=xfc_tmp1(ji+j0,3)
            ELSE
               j_neigh_tmp(ji+j0)=0
               DO m=1,n_bc
                  m1=nb_start(m)
                  m2=nb_start(m)+nb_number(m)-1
                  IF(ii.ge.m1.and.ii.le.m2)then
                     j_nbcon_tmp(ji+j0)=nb_condition(m)
                     IF(nb_condition(m).gt.nin_max.and.nb_condition(m).le.nb_max)THEN 
                        npb_tmp(i)=nb_condition(m)-nin_max
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDIF
!
      IF(ndim.eq.2)THEN

         DO i=1,nn
            vol_tmp0=0.d0
            xloc_tmp1=0.d0
            xloc_tmp2=0.d0
            DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
                 xloc_tmp1=xloc_tmp1+xfc_tmp1(j,1)
                 xloc_tmp2=xloc_tmp2+xfc_tmp1(j,2)
                 vol_tmp0=vol_tmp0+xfc_tmp1(j,1)*sv_tmp1(j,1) &
                                  +xfc_tmp1(j,2)*sv_tmp1(j,2)
            ENDDO
            vol_tmp(i)=vol_tmp0/DBLE(ndim)
            xloc_tmp(i,1)=xloc_tmp1/DBLE(num_neigh_tmp(i))
            xloc_tmp(i,2)=xloc_tmp2/DBLE(num_neigh_tmp(i))
         ENDDO
!
      ELSEIF(ndim.eq.3) THEN
!
         DO i=1,nn
            j0=i_neigh_tmp(i)-1
            vol_tmp0=0.d0
            xloc_tmp1=0.d0
            xloc_tmp2=0.d0
            xloc_tmp3=0.d0
            xc_tmp1=0.0d0
            xc_tmp2=0.0d0
            xc_tmp3=0.0d0
            DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
               xc_tmp1=xc_tmp1+xfc_tmp1(j,1)
               xc_tmp2=xc_tmp2+xfc_tmp1(j,2)
               xc_tmp3=xc_tmp3+xfc_tmp1(j,3)
            ENDDO
            xc_tmp(1)=xc_tmp1/DBLE(num_neigh_tmp(i))
            xc_tmp(2)=xc_tmp2/DBLE(num_neigh_tmp(i))
            xc_tmp(3)=xc_tmp3/DBLE(num_neigh_tmp(i))
            DO j=1,num_neigh_tmp(i)
               ii=neigh_face_tmp1(j+j0)
               DO m=3,nd(ii)
                  i1=node_face(1,ii)
                  i2=node_face(m-1,ii)
                  i3=node_face(m,ii)
               xnode1(:)=xnode(i1,:)
               xnode2(:)=xnode(i2,:)
               xnode3(:)=xnode(i3,:)
!                 CALL xc_vol_tetra(xc_tmp(:),xnode(i1,:),xnode(i2,:),xnode(i3,:),xc,vol)
                  CALL xc_vol_tetra(xc_tmp,xnode1,xnode2,xnode3,xc,vol)
                  xloc_tmp1=xloc_tmp1+xc(1)*vol
                  xloc_tmp2=xloc_tmp2+xc(2)*vol
                  xloc_tmp3=xloc_tmp3+xc(3)*vol
                  vol_tmp0=vol_tmp0+vol
               ENDDO
            ENDDO
            vol_tmp(i)=vol_tmp0
            xloc_tmp(i,1)=xloc_tmp1/vol_tmp(i)
            xloc_tmp(i,2)=xloc_tmp2/vol_tmp(i)
            xloc_tmp(i,3)=xloc_tmp3/vol_tmp(i)
         ENDDO
      ENDIF
!
      IF(ndim.eq.2)THEN
         DO i=1,nn
            sv_tmp10=0.0d0
            sv_tmp20=0.0d0
            DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
               sv_tmp10=sv_tmp10+sv_tmp1(j,1)
               sv_tmp20=sv_tmp20+sv_tmp1(j,2)
            ENDDO
            sv_t_all(i,1)=sv_tmp10
            sv_t_all(i,2)=sv_tmp20
         ENDDO  
      ELSEIF(ndim.eq.3) THEN
         DO i=1,nn
            sv_tmp10=0.0d0
            sv_tmp20=0.0d0
            sv_tmp30=0.0d0
            DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
               sv_tmp10=sv_tmp10+sv_tmp1(j,1)
               sv_tmp20=sv_tmp20+sv_tmp1(j,2)
               sv_tmp30=sv_tmp30+sv_tmp1(j,3)
            ENDDO
            sv_t_all(i,1)=sv_tmp10
            sv_t_all(i,2)=sv_tmp20
            sv_t_all(i,3)=sv_tmp30
         ENDDO  
      ENDIF
!
!........Check negative volumes and wrong surface factors
!
      DO i=1,nn
         IF(vol_tmp(i).le.1.0e-20)THEN
            WRITE(*,*) '***Negative volume occured', i,vol_tmp(i) 
            WRITE(unit_log,*) '***Negative volume occured', i,vol_tmp(i) 
            PAUSE
            STOP
         ENDIF
         DO ix=1,ndim         
            IF(ABS(sv_t_all(i,ix)).ge.1.0e-12)THEN
               WRITE(*,*) '***Error in SV occured',i, ix, sv_t_all(i,ix)
               WRITE(unit_log,*) '***Error in SV occured',i, ix, sv_t_all(i,ix)
               PAUSE
               STOP
            ENDIF
         ENDDO         
      ENDDO  
      WRITE(*,*) '          Volume and surface vector checked.'
      WRITE(unit_log,*) '          Volume and surface vector checked.'      
!
     IF(n_zone.gt.0)THEN
        IF(ndim.eq.2) THEN
           DO k=1,n_zone
              DO i=nz_start(k),nz_start(k)+nz_number(k)-1
                 DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
                    sa_tmp=SQRT(sv_tmp1(j,1)**2+sv_tmp1(j,2)**2)
                    xn_tmp1=sv_tmp1(j,1)/sa_tmp
                    xn_tmp2=sv_tmp1(j,2)/sa_tmp
!                   perm_tmp(j,i)=xn_tmp1*nz_perm(k,1)*xn_tmp1*nz_perm(k,1)+xn_tmp2*nz_perm(k,2)*xn_tmp2*nz_perm(k,2)
                    perm_tmp1(j)=SQRT((xn_tmp1*nz_perm(k,1))**2+(xn_tmp2*nz_perm(k,2))**2)
                    perm_tmp1(j)=MIN(1.0d0,perm_tmp1(j))
                 ENDDO
                 nzone_tmp(i)=nz_index(k) 
                 celem(i)=nz_subdomain(k)
                 poro_tmp(i)=nz_porosity(k)
              ENDDO
           ENDDO  
        ELSE
           DO k=1,n_zone
              DO i=nz_start(k),nz_start(k)+nz_number(k)-1
                 DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
                    sa_tmp=DSQRT(sv_tmp1(j,1)**2+sv_tmp1(j,2)**2+sv_tmp1(j,3)**2)
                    xn_tmp1=sv_tmp1(j,1)/sa_tmp
                    xn_tmp2=sv_tmp1(j,2)/sa_tmp
                    xn_tmp3=sv_tmp1(j,3)/sa_tmp
!                   perm_tmp(j,i)=xn_tmp1*nz_perm(k,1)*xn_tmp1*nz_perm(k,1)+xn_tmp2*nz_perm(k,2)*xn_tmp2*nz_perm(k,2)
!                   IF(ndim.eq.3)perm_tmp(j,i)=perm_tmp(j,i)+xn_tmp(j,ndim)*nz_perm(k,ndim)*xn_tmp(j,ndim)*nz_perm(k,ndim)
                    perm_tmp1(j)=SQRT((xn_tmp1*nz_perm(k,1))**2+(xn_tmp2*nz_perm(k,2))**2+(xn_tmp3*nz_perm(k,3))**2)
                    perm_tmp1(j)=MIN(1.0d0,perm_tmp1(j))
                 ENDDO
                 nzone_tmp(i)=nz_index(k) 
                 celem(i)=nz_subdomain(k)
                 poro_tmp(i)=nz_porosity(k)
              ENDDO
            ENDDO  
         ENDIF
      ELSE
         nzone_tmp(:)=1 
         celem(:)=0
         poro_tmp(:)=1.0d0
         perm_tmp1(:)=1.0d0         
      ENDIF
!!!!
!!!      IF(udfl_read_openfoam)CALL udfn_read_openfoam(poro_tmp,xloc_tmp,celem)
!!!!
      WRITE(*,*) '          Successful grid generation from openfoam grid files.'
      
! sang test
  
      DEALLOCATE(owner,neighbour)
      DEALLOCATE(sv_t_all)
!
      END SUBROUTINE read_openfoam
!
!----------------------------------------------------------------------
!
      SUBROUTINE sv_tri_2d(x1,x2,x3,sv)
!
      IMPLICIT NONE
!
      REAL(8) x1(3),x2(3),x3(3),sv(3)
      REAL(8) dx1,dx2,dy1,dy2,dz1,dz2
!
      dx1=x2(1)-x1(1)
      dx2=x3(1)-x1(1)
      dy1=x2(2)-x1(2)
      dy2=x3(2)-x1(2)
      dz1=x2(3)-x1(3)
      dz2=x3(3)-x1(3)
!
      sv(1)=(dy1*dz2-dz1*dy2)*0.5d0
      sv(2)=(dz1*dx2-dx1*dz2)*0.5d0
      sv(3)=(dx1*dy2-dy1*dx2)*0.5d0
!
      END SUBROUTINE sv_tri_2d
!
!----------------------------------------------------------------------
!
      SUBROUTINE sv_tri(x1,x2,x3,sv,xc,a)
!
      IMPLICIT NONE
!
      REAL(8) x1(3),x2(3),x3(3),sv(3),xc(3),a
      REAL(8) dx1,dx2,dy1,dy2,dz1,dz2
!
      dx1=x2(1)-x1(1)
      dx2=x3(1)-x1(1)
      dy1=x2(2)-x1(2)
      dy2=x3(2)-x1(2)
      dz1=x2(3)-x1(3)
      dz2=x3(3)-x1(3)
!
      sv(1)=(dy1*dz2-dz1*dy2)*0.5d0
      sv(2)=(dz1*dx2-dx1*dz2)*0.5d0
      sv(3)=(dx1*dy2-dy1*dx2)*0.5d0
!
      xc(1)=(x1(1)+x2(1)+x3(1))/3.0d0
      xc(2)=(x1(2)+x2(2)+x3(2))/3.0d0
      xc(3)=(x1(3)+x2(3)+x3(3))/3.0d0
      a=sv(1)**2+sv(2)**2+sv(3)**2
!
!      IF(a.eq.0.0) STOP '### false triagle !'
!
      a=SQRT(a)
!
      END SUBROUTINE sv_tri
!
!----------------------------------------------------------------------
!
      SUBROUTINE xc_vol_tetra(x0,x1,x2,x3,xc,vol)
!
      IMPLICIT NONE
!
!     input
      REAL(8) x0(3),x1(3),x2(3),x3(3)
!     output
      REAL(8) xc(3),vol
!     local variables
      REAL(8) a1,a2,a3
      REAL(8) b1,b2,b3
      REAL(8) c1,c2,c3
!
      xc(1)=(x0(1)+x1(1)+x2(1)+x3(1))/4.0d0
      xc(2)=(x0(2)+x1(2)+x2(2)+x3(2))/4.0d0
      xc(3)=(x0(3)+x1(3)+x2(3)+x3(3))/4.0d0
!
      a1=x1(1)-x0(1)
      b1=x2(1)-x0(1)
      c1=x3(1)-x0(1)
      a2=x1(2)-x0(2)
      b2=x2(2)-x0(2)
      c2=x3(2)-x0(2)
      a3=x1(3)-x0(3)
      b3=x2(3)-x0(3)
      c3=x3(3)-x0(3)
!
      vol=a1*(b2*c3-b3*c2)-b1*(a2*c3-a3*c2)+c1*(a2*b3-a3*b2)
      vol=DABS(vol)/6.0d0
!
      RETURN
      ENDSUBROUTINE xc_vol_tetra
!
!----------------------------------------------------------------------
!
