      MODULE Znode
!      
      IMPLICIT NONE
      SAVE
!
      INTEGER :: nd_max,n_node,nmax_vertex
!
      INTEGER, ALLOCATABLE::nd(:)
      INTEGER, ALLOCATABLE::neigh_face_tmp1(:)
      INTEGER, ALLOCATABLE::i_cell_node_tmp(:),num_cell_node_tmp(:)
      INTEGER, ALLOCATABLE::num_cell_node(:),num_nd(:)
      INTEGER, ALLOCATABLE::node_face(:,:)
      INTEGER, ALLOCATABLE::cell_node_tmp(:,:)
      REAL(8), ALLOCATABLE::swcn(:)
!
      INTEGER, ALLOCATABLE::cell_node(:,:),node_face_cell(:,:)
      REAL(8), ALLOCATABLE::rwcn(:,:),xnode(:,:),dxr(:,:,:)
      REAL(8), ALLOCATABLE::rwcn_tmp(:,:)
!      
      END MODULE Znode
      
