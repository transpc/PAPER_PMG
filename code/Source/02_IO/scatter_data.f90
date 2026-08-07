!
      SUBROUTINE scatter_data(npb_tmp,nzone_tmp,sl_tmp,sgap_tmp,nmaterial_tmp,icore_tmp, &
                              vol_tmp,poro_tmp,hydraulicd_tmp,                           &
                              nd_tmp1,cell_node_tmp1,rwcn_tmp1,dxr_tmp1,                 &
                              n_face,nb_cell_node)
!
!.....This routine distributes global variables across mpi tasks
!
      USE Zmpi         , ONLY: ncell_fp,maxmt_nfluid,maxmt_cell,maxmt_fluid,maxmt_fp
      USE Zzone        , ONLY: nzone,ncell_fluid,ncell_fluid_all,icore,nmaterial
      USE Zcore        , ONLY: np,myrank
      USE Zbc_index    , ONLY: npb
      USE Zparam       , ONLY: nn,ns,ndim,mesh_openfoam
      USE Znum_cell    , ONLY: i_neigh_tmp,perm_tmp1,sv_tmp1,xfc_tmp1, &
                               i_neigh, &
                               index_sort
      USE Zconst2      , ONLY: hydraulicd,sl
      USE Zcoord1      , ONLY: xloc,xloc_tmp
      USE Zcoord2      , ONLY: xfc
      USE Zcoord3      , ONLY: porosity,sv,vol,permeability
      USE Znode        , ONLY: nd_max,nmax_vertex,neigh_face_tmp1,num_nd,num_cell_node,cell_node,node_face_cell, &
                               rwcn,swcn,dxr,                                                                    &
                               node_face,n_node,                                                                 &
                               i_cell_node_tmp,num_cell_node_tmp
      USE Zgradoption  , ONLY: ifrink
      USE Zporous      , ONLY: sgap
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n_face,nb_cell_node
      INTEGER,DIMENSION(nn) :: npb_tmp,nzone_tmp,icore_tmp,nmaterial_tmp
      INTEGER,DIMENSION(maxmt_nfluid) :: nd_tmp1
      INTEGER,DIMENSION(nb_cell_node) :: cell_node_tmp1
      REAL(8),DIMENSION(nn) :: vol_tmp,poro_tmp,hydraulicd_tmp
      REAL(8),DIMENSION(nn,ndim) :: sl_tmp,sgap_tmp
      REAL(8),DIMENSION(nb_cell_node) :: rwcn_tmp1
      REAL(8),DIMENSION(nb_cell_node,ndim) :: dxr_tmp1
!.....Local variables
      INTEGER :: i,j,k,ii,ix,j0,j1,m
      INTEGER :: nb_cell_node_l
      INTEGER :: izone=0
!.....Local arrays
      REAL(8),DIMENSION(ns) :: tmp
!.....Allocatable arrays
      INTEGER,DIMENSION(:),ALLOCATABLE :: neigh_face_l
!
      CALL scatterv_i(npb      ,npb_tmp      ,ncell_fluid,ncell_fluid_all,izone)
      CALL scatterv_i(nzone    ,nzone_tmp    ,ncell_fluid,ncell_fluid_all,izone)
      CALL scatterv_i(icore    ,icore_tmp    ,ncell_fluid,ncell_fluid_all,izone)
      CALL scatterv_i(nmaterial,nmaterial_tmp,ncell_fluid,ncell_fluid_all,izone)
!
      CALL scatterv_r(vol       ,vol_tmp       ,ncell_fluid,ncell_fluid_all,izone)
      CALL scatterv_r(porosity  ,poro_tmp      ,ncell_fluid,ncell_fluid_all,izone)
      CALL scatterv_r(hydraulicd,hydraulicd_tmp,ncell_fluid,ncell_fluid_all,izone)
!
      CALL scatterv_ndim_r(sl,sl_tmp,ncell_fluid,nn,izone)
      CALL scatterv_ndim_fp_r(xloc,xloc_tmp,ncell_fluid,ncell_fp,nn,izone)
      CALL scatterv_ndim_fp_r(sgap,sgap_tmp,ncell_fluid,ncell_fp,nn,izone)
!
      CALL scatterv_csr_r(permeability,maxmt_fluid,perm_tmp1,maxmt_cell,ncell_fluid_all,maxmt_nfluid,i_neigh_tmp,izone)
!
!.....Get sorted permeability
!
      DO i=1,ncell_fluid
         j0=i_neigh(i)-1 
         DO j=i_neigh(i),i_neigh(i+1)-1
             k=index_sort(j)
             tmp(j-j0)=permeability(k+j0)
         ENDDO
         DO j=i_neigh(i),i_neigh(i+1)-1
            permeability(j)=tmp(j-j0)
         ENDDO
      ENDDO
!
      DO ix=1,ndim
         CALL scatterv_csr_r(sv(1,ix),maxmt_fluid,sv_tmp1(1,ix),maxmt_cell,ncell_fluid_all,maxmt_nfluid,i_neigh_tmp,izone)
!
!........Get sorted sv
!
         DO i=1,ncell_fluid
            j0=i_neigh(i)-1 
            DO j=i_neigh(i),i_neigh(i+1)-1
               k=index_sort(j)
               tmp(j-j0)=sv(k+j0,ix)
            ENDDO
            DO j=i_neigh(i),i_neigh(i+1)-1
               sv(j,ix)=tmp(j-j0)
            ENDDO
         ENDDO
      ENDDO
!
      DO ix=1,ndim
         CALL scatterv_csr_r(xfc(1,ix),maxmt_fluid,xfc_tmp1(1,ix),maxmt_cell,ncell_fluid_all,maxmt_nfluid,i_neigh_tmp,izone)
!
!........Get sorted xfc
!
         DO i=1,ncell_fluid
            j0=i_neigh(i)-1 
            DO j=i_neigh(i),i_neigh(i+1)-1
               k=index_sort(j)
               tmp(j-j0)=xfc(k+j0,ix)
            ENDDO
            DO j=i_neigh(i),i_neigh(i+1)-1
               xfc(j,ix)=tmp(j-j0)
            ENDDO
         ENDDO
      ENDDO
!
      IF(mesh_openfoam.eq.1.and.ifrink.ge.1)THEN
         CALL scatterv_csr_fp_i(num_nd,maxmt_fluid,nd_tmp1,maxmt_nfluid,ncell_fluid_all,maxmt_nfluid,i_neigh_tmp,izone, &
                                maxmt_fp,ncell_fluid,i_neigh)
!
!.....Communicate num_nd for grad_frink
!
         IF(np.gt.1) CALL communicate_1d_csr_int(num_nd,i_neigh)  
!
!........We need to access node_face,swcn that is defined on rank 0 only
!
         IF(np.gt.1) THEN
            CALL broadcast_i1(n_face)
            CALL broadcast_i1(nmax_vertex)
         ENDIF
         IF(myrank.ne.0) THEN
            ALLOCATE(node_face(nmax_vertex,n_face)) 
            ALLOCATE(swcn(n_node)) 
         ENDIF
         IF(np.gt.1) THEN
            CALL broadcast_i(node_face,nmax_vertex*n_face)
            CALL broadcast_r(swcn,n_node)
         ENDIF
         ALLOCATE(neigh_face_l(maxmt_fluid))
         CALL scatterv_csr_i(neigh_face_l,maxmt_fluid,neigh_face_tmp1,maxmt_nfluid,ncell_fluid_all,maxmt_nfluid,i_neigh_tmp,izone)
         j0=0
         DO i=1,ncell_fluid
            j1=i_neigh(i)-1
            DO j=i_neigh(i),i_neigh(i+1)-1
               j0=j0+1
               ii=neigh_face_l(j0)
               DO m=1,num_nd(j)
                  node_face_cell(m,j)=node_face(m,ii)
               ENDDO
            ENDDO
         ENDDO
         CALL scatterv_i(num_cell_node,num_cell_node_tmp,ncell_fluid,ncell_fluid_all,izone)
!
!........Communicate node_face_cell for grad_frink
!
         IF(np.gt.1) THEN
            DO m=1,nmax_vertex
               CALL communicate_1d_csr_int(node_face_cell(m,:),i_neigh)
            ENDDO
         ENDIF
!
!........Get the nb_cell_node local to allocate buffer for scatter communicate
!
         nb_cell_node_l=0
         DO i=1,ncell_fluid
            nb_cell_node_l=nb_cell_node_l+num_cell_node(i)
         ENDDO
!
!........Num_cell_node_tmp different from num_neigh_tmp 
!........Use i_cell_node_tmp instead of i_neigh_tmp
!
         CALL scatterv_node_csr_fp_i(cell_node,nb_cell_node_l,cell_node_tmp1,nb_cell_node,ncell_fluid_all,nd_max,ncell_fp,nb_cell_node,i_cell_node_tmp, &
                                     ncell_fluid,num_cell_node)
         CALL scatterv_node_csr_fp_r(rwcn,nb_cell_node_l,rwcn_tmp1,nb_cell_node,ncell_fluid_all,nd_max,ncell_fp,nb_cell_node,i_cell_node_tmp, &
                                     ncell_fluid,num_cell_node)
         DO ix=1,ndim
            CALL scatterv_node_csr_fp_r(dxr(1,1,ix),nb_cell_node_l,dxr_tmp1(1,ix),nb_cell_node,ncell_fluid_all,nd_max,ncell_fp,nb_cell_node,i_cell_node_tmp, &
                                        ncell_fluid,num_cell_node)
         ENDDO
         DEALLOCATE(i_cell_node_tmp)
      ENDIF
!
      END SUBROUTINE scatter_data 
