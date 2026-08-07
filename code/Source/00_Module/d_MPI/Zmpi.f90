      MODULE Zmpi
      
      IMPLICIT NONE
      SAVE
!      
!.....Serial decomposition
!
      INTEGER,DIMENSION(:),ALLOCATABLE :: celem,celem_c
!     
!.....For Communication
!
      INTEGER :: metis
      INTEGER :: niut,alstatus
      INTEGER :: ncell_fp,maxmt,maxmt_pad,max_neigh
      INTEGER,DIMENSION(:),ALLOCATABLE :: ri,si,iut,rintf,sintf,iperm,jperm
      INTEGER,DIMENSION(:),ALLOCATABLE :: mapping_ext        
      INTEGER, PARAMETER :: mxnbne=60
!      
!.....For CSR array
!
      INTEGER :: maxmt_nfluid,maxmt_cell
      INTEGER :: maxmt_fluid,maxmt_fp
      INTEGER,DIMENSION(:),ALLOCATABLE::ia_a,ja_a,ju_a,iEND
      INTEGER,DIMENSION(:),ALLOCATABLE::ia_nrhs
      REAL(8),DIMENSION(:),ALLOCATABLE::au
!.....For vector amux0
      INTEGER,DIMENSION(:),ALLOCATABLE::jap,jaa,jaar
      INTEGER,DIMENSION(:,:),ALLOCATABLE::nbgroup,iaa,iap
      REAL(8),DIMENSION(:),ALLOCATABLE::ap
      INTEGER :: ngroup
!.....For symbolic factorization
      INTEGER :: levt,lev_typet,lev_typedt,nnp
      INTEGER :: levt_c,lev_typet_c
      INTEGER :: maxmt2,maxmt_r
      INTEGER :: maxmt2_c,maxmt_r_c
      INTEGER,DIMENSION(:),ALLOCATABLE :: ia_r,ja_r,ju_r
      INTEGER,DIMENSION(:),ALLOCATABLE :: ia_r_c,ja_r_c,ju_r_c
      INTEGER,DIMENSION(:),ALLOCATABLE :: perm_r,permi_r,index_r
      INTEGER,DIMENSION(:),ALLOCATABLE :: perm_r_c,permi_r_c,index_r_c
      REAL(8),DIMENSION(:),ALLOCATABLE :: au_r
      REAL(8),DIMENSION(:),ALLOCATABLE :: au_r_c
!.....For split rhs
      INTEGER maxmt_lu0,maxmt_lu1
      INTEGER,DIMENSION(:),ALLOCATABLE :: ja0,ja1,ia0,ia1
      REAL(8),DIMENSION(:),ALLOCATABLE :: diag_lu,alu0,alu1
!     REAL(8),DIMENSION(:),ALLOCATABLE :: alu
!
      INTEGER count
! 
      INTEGER :: maxmt_c,maxmt_pad_c,niut_c,ncell_ps,max_neigh_c
      INTEGER :: maxmt_nncond
      INTEGER :: maxmt_ncond,maxmt_ps
      INTEGER,DIMENSION(:),ALLOCATABLE::ia_a_c,ja_a_c,ju_a_c,iend_c
      REAL(8),DIMENSION(:),ALLOCATABLE::au_c
!.....For fsw 
      INTEGER :: nbj_fsw
      INTEGER,DIMENSION(:),ALLOCATABLE :: cell_fsw,ia_fsw,ja_fsw,celem_fsw
!.....For vector amux0 solid processing
      INTEGER,DIMENSION(:),ALLOCATABLE :: jap_c,jaa_c,jaar_c
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: nbgroup_c,iaa_c,iap_c
      REAL(8),DIMENSION(:),ALLOCATABLE :: ap_c
      INTEGER :: ngroup_c
!.....For split rhs
      INTEGER :: maxmt_lu0_c,maxmt_lu1_c
      INTEGER,DIMENSION(:),ALLOCATABLE :: ja0_c,ja1_c,ia0_c,ia1_c
      REAL(8),DIMENSION(:),ALLOCATABLE :: diag_lu_c,alu0_c,alu1_c
!     REAL(8),DIMENSION(:),ALLOCATABLE :: alu_c
!
      INTEGER,DIMENSION(:),ALLOCATABLE :: ri_c,si_c,iut_c,rintf_c,sintf_c,jperm_c
!.....For ALLGATHERV
      INTEGER,DIMENSION(:),ALLOCATABLE :: ncell_fluid1,ncell_fp1,ncell_fluid1_dsp
      INTEGER,DIMENSION(:),ALLOCATABLE :: ncell_fluid1_2d,ncell_fluid1_2d_dsp
      INTEGER,DIMENSION(:),ALLOCATABLE :: ncell_fluid1_c,ncell_fp1_c,ncell_fluid1_dsp_c
      INTEGER,DIMENSION(:),ALLOCATABLE :: ncell_csr_sz,ncell_csr_dsp
      INTEGER,DIMENSION(:),ALLOCATABLE :: ncell_ndim_sz,ncell_ndim_dsp
      INTEGER,DIMENSION(:),ALLOCATABLE :: ncell_node_csr_sz,ncell_node_csr_dsp
      INTEGER,DIMENSION(:),ALLOCATABLE :: ncell_csr_sz_c,ncell_csr_dsp_c
      INTEGER,DIMENSION(:),ALLOCATABLE :: ncell_ndim_sz_c,ncell_ndim_dsp_c
      INTEGER,DIMENSION(:),ALLOCATABLE :: ncell_csr_sz_nbcon0,ncell_csr_dsp_nbcon0
!
!.....Separate pre-processing
!
      INTEGER,DIMENSION(:),ALLOCATABLE :: numiut,jjperm,jjperm_c
      INTEGER,DIMENSION(:),ALLOCATABLE :: i3perm
      INTEGER,DIMENSION(:),ALLOCATABLE :: iperm1,jperm1,num_neigh1
      INTEGER,DIMENSION(:),ALLOCATABLE :: npb1
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: iutp
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: ssi,rri
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: sint,rint,iiperm
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: neigh1,nbcon1
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: nji1
      REAL(8),DIMENSION(:),ALLOCATABLE :: vol1
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: xloc1
      REAL(8),DIMENSION(:,:,:),ALLOCATABLE :: sv1,xfc1
!
      END MODULE Zmpi
      
