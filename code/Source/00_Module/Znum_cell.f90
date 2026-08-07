      MODULE Znum_cell
!      
      IMPLICIT NONE
      SAVE
!
      INTEGER,PARAMETER :: nb_nf=11,nbc_nf=4
!
!.....Scalar variables    
!
      INTEGER :: ncell,nbfacemsh,nnbnd,nbncel,nnode,ireadnode
      INTEGER :: nf_number_nb,lens
!
!.....Arrays variables
!
      INTEGER,DIMENSION(2,-1:nb_nf) :: istart_nf,istart_nf_old
      INTEGER,DIMENSION(2,-1:nb_nf) :: istart_nb1
      INTEGER,DIMENSION(-1:nb_nf) :: nf_number_id,istart_nfs
      INTEGER,DIMENSION(0:nb_nf) :: istart_nbcon_nf,istart_nbcon_nf_old,istart_svp_nf,istart_svp_nf_old
      INTEGER,DIMENSION(2,-1:nbc_nf) :: istartc_nf,istartc_nb1
!
!.....1D allocatable variables    
!
      INTEGER,DIMENSION(:),ALLOCATABLE :: i_neigh,neigh,n_fluid
      INTEGER,DIMENSION(:),ALLOCATABLE :: ia_nb,icell_nb,iptr_nb_k,right_nb_k
      INTEGER,DIMENSION(:),ALLOCATABLE :: i_neigh_c,neigh_c
      INTEGER,DIMENSION(:),ALLOCATABLE :: iac_nb,icellc_nb,iptrc_nb_k,rightc_nb_k
      INTEGER,DIMENSION(:),ALLOCATABLE :: num_neigh_o,num_neigh_n                            ! mcc
      INTEGER,DIMENSION(:),ALLOCATABLE :: index_sort,indexr_sort
      INTEGER,DIMENSION(:),ALLOCATABLE :: index_sort_c
      INTEGER,DIMENSION(:),ALLOCATABLE :: num_neigh_tmp,i_neigh_tmp,j_neigh_tmp,j_nbcon_tmp
      INTEGER,DIMENSION(:),ALLOCATABLE :: i_neigh_c_tmp,j_neigh_c_tmp,j_nbcon_c_tmp
      INTEGER,DIMENSION(:),ALLOCATABLE :: neigh_face_tmp1                                    ! csr form
      INTEGER,DIMENSION(:),ALLOCATABLE :: index_sort_tmp1
      INTEGER,DIMENSION(:),ALLOCATABLE :: node_ring                                          ! ibundle=1
      REAL(8),DIMENSION(:),ALLOCATABLE :: perm_tmp1                                          ! csr form
!
!.....2D allocatable variables    
!
      INTEGER,DIMENSION(:),ALLOCATABLE :: index_nbcon_nf               ! nbcon change index to get old value
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: sv_tmp1,xfc_tmp1           ! csr form
!
      END MODULE Znum_cell
