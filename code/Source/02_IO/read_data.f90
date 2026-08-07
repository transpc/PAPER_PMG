!
      SUBROUTINE read_data(npb_tmp,nzone_tmp,nmaterial_tmp,celem,vol_tmp,poro_tmp, &
                           n_face,nb_cell_node)
!
!.....This routine reads geometry data from 'somaGrid.in'
!
!
      USE Zmpi         , ONLY: maxmt_nfluid,maxmt_cell
      USE Zzone        , ONLY: ncell_fluid_all,ncell_cond_all
      USE Zparam       , ONLY: nn,ns,ndim,nin_max,nb_max,mesh_openfoam,mesh_binary
      USE Znum_cell    , ONLY: ncell,num_neigh_tmp,i_neigh_tmp,j_neigh_tmp,j_nbcon_tmp, &
                               perm_tmp1,sv_tmp1,xfc_tmp1
      USE Zcoord1      , ONLY: xloc_tmp
      USE Znode        , ONLY: nmax_vertex,n_node,neigh_face_tmp1,nd,node_face,xnode, &
                               i_cell_node_tmp,num_cell_node_tmp,cell_node_tmp,  &
                               swcn,rwcn_tmp,                                    &
                               nd_max
      USE Zgradoption  , ONLY: ifrink
      USE Znode        , ONLY: nd_max
      USE Zio_unit      , ONLY: unit_grid
!
      IMPLICIT NONE
!
!     input
!     output
      INTEGER :: n_face
      INTEGER :: npb_tmp(nn)
      INTEGER :: nzone_tmp(nn)
      INTEGER :: nmaterial_tmp(nn)
      INTEGER :: celem(nn)
      REAL(8) :: vol_tmp(nn),poro_tmp(nn)
!     local variables
      INTEGER :: i,ii,j,jj,ix,n1,j0
      INTEGER :: m,n
      INTEGER :: nb_cell_node
      INTEGER :: n_bface,n_bc,n_zone
      REAL(8) :: dznull
      CHARACTER*50 :: filename
!     local arrays
      INTEGER :: nc(2*ns)
!     in modules
!     
         IF(mesh_openfoam.eq.1)THEN
            READ(unit_grid,*)filename
            READ(unit_grid,*)n_node,nn,n_face,n_bface,n_bc,n_zone,nmax_vertex,ns,ndim,dznull
!
            ALLOCATE(xnode(n_node,3),nd(n_face),node_face(nmax_vertex,n_face)) 
            xnode(:,:)=0.0d0
            nd(:)=0
            node_face(:,:)=0
!
            DO i=1,n_node
               READ(unit_grid,*) (xnode(i,ix),ix=1,3) 
            ENDDO
!      
            DO i=1,n_face
               READ(unit_grid,*) nd(i),(node_face(j,i),j=1,nd(i))
               DO j=1,nd(i)
                  node_face(j,i)=node_face(j,i)+1
               ENDDO
            ENDDO
            CALL read_openfoam(n_face,n_bface,n_bc,n_zone,                       &
                               npb_tmp,vol_tmp,poro_tmp,nzone_tmp,nmaterial_tmp, &
                               celem,dznull)
!
            IF(ifrink.ge.1)THEN
               nd_max=0
               DO i=1,nn
                  j0=i_neigh_tmp(i)-1
!                 nc(:)=0
                  n1=0
                  DO j=1,num_neigh_tmp(i)
!                 DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
                     ii=neigh_face_tmp1(j+j0)
                     DO m=1,nd(ii)
                        n=node_face(m,ii)
                        DO jj=1,n1
                           IF(n.eq.nc(jj)) goto 110
                        ENDDO
                        n1=n1+1
                        nc(n1)=n
110                     CONTINUE
                     ENDDO
                  ENDDO
                  IF(n1.gt.nd_max) nd_max=n1
                  DO jj=1,n1
                     nc(jj)=0
                  ENDDO
               ENDDO
!
               ALLOCATE(i_cell_node_tmp(nn+1))
               ALLOCATE(cell_node_tmp(nd_max,nn))
               ALLOCATE(rwcn_tmp(nd_max,nn))
               ALLOCATE(num_cell_node_tmp(nn))
               ALLOCATE(swcn(n_node))
               cell_node_tmp(:,:)=0         
               rwcn_tmp(:,:)=0.0d0
               swcn(:)=0.0d0
               CALL frink_weight(neigh_face_tmp1,nd,node_face,n_face,xloc_tmp,  &
                                 num_cell_node_tmp,cell_node_tmp,rwcn_tmp)
               i_cell_node_tmp(1)=1
               DO i=1,ncell_fluid_all
                  i_cell_node_tmp(i+1)=i_cell_node_tmp(i)+num_cell_node_tmp(i)
               ENDDO
               nb_cell_node=i_cell_node_tmp(ncell_fluid_all+1)-1
            ENDIF                   
!
         ELSE
!
            DO i=1,ncell
               poro_tmp(i)=1.0d0
               IF(mesh_binary.eq.1)THEN
                  READ(unit_grid)ii,(xloc_tmp(i,ix),ix=1,ndim),vol_tmp(i),num_neigh_tmp(i),npb_tmp(i),nzone_tmp(i),nmaterial_tmp(i),celem(i),poro_tmp(i)
               ELSE   
                  READ(unit_grid,*)ii,(xloc_tmp(i,ix),ix=1,ndim),vol_tmp(i),num_neigh_tmp(i),npb_tmp(i),nzone_tmp(i),nmaterial_tmp(i),celem(i),poro_tmp(i)
               ENDIF   
            ENDDO
!
!     get ncell_fluid_all,ncell_cond_all
!
            ncell_fluid_all=0
            ncell_cond_all=0
            DO i=1,ncell
               IF(nmaterial_tmp(i).le.0) ncell_fluid_all=ncell_fluid_all+1
               IF(nmaterial_tmp(i).ne.0) ncell_cond_all=ncell_cond_all+1
            ENDDO
!
!     We need to go csr format
!
!        write(*,*) 'ncell_fluid_all',ncell_fluid_all,ncell_cond_all,ncell,ns
            i_neigh_tmp(1)=1
            DO i=1,ncell_fluid_all
               i_neigh_tmp(i+1)=i_neigh_tmp(i)+num_neigh_tmp(i)
            ENDDO
           maxmt_nfluid=i_neigh_tmp(ncell_fluid_all+1)-1
!          write(*,*) 'maxmt_nfluid',maxmt_nfluid
!
           maxmt_cell=maxmt_nfluid
           DO i=ncell_fluid_all+1,ncell
              i_neigh_tmp(i+1)=i_neigh_tmp(i)+num_neigh_tmp(i)
           ENDDO
           maxmt_cell=i_neigh_tmp(ncell+1)-1
!          write(*,*) 'maxmt_cell',maxmt_cell
!            
!........Read face data and put directly in csr format
!
            ALLOCATE(j_neigh_tmp(maxmt_cell),j_nbcon_tmp(maxmt_cell))
            ALLOCATE(perm_tmp1(maxmt_cell),sv_tmp1(maxmt_cell,ndim),xfc_tmp1(maxmt_cell,ndim))
            DO i=1,ncell
               npb_tmp(i)=0 
               DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
                  IF(mesh_binary.eq.1)THEN
                     READ(unit_grid) ii,jj,j_neigh_tmp(j),j_nbcon_tmp(j),(sv_tmp1(j,ix),ix=1,ndim),(xfc_tmp1(j,ix),ix=1,ndim),perm_tmp1(j)
                  ELSE   
                     READ(unit_grid,*) ii,jj,j_neigh_tmp(j),j_nbcon_tmp(j),(sv_tmp1(j,ix),ix=1,ndim),(xfc_tmp1(j,ix),ix=1,ndim),perm_tmp1(j)
                  ENDIF   
                  IF(j_nbcon_tmp(j).gt.nin_max.and.j_nbcon_tmp(j).le.nb_max) THEN
                     npb_tmp(i)=j_nbcon_tmp(j)-nin_max 
                  ENDIF   
               ENDDO
            ENDDO
            CLOSE(4)
         ENDIF
!
      END SUBROUTINE read_data
