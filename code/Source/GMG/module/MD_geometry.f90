!=====================================================================!
! Geometry arrays ----------------------------------------------------|
!---------------------------------------------------------------------+
      module md_geometry
      implicit none
!---
      integer(4):: nelem,nnode,nelem_mg,nnodegl,nelemgl,nnode_mg
!      
	  real(8),dimension(:,:), allocatable :: coord
      
      integer,dimension(:),allocatable:: num_neigh,num_neigh_mg
      integer,dimension(:,:),allocatable:: e_neigh,neigh_mg
      integer(kind=4),dimension(:),allocatable :: imap
!---
      save
    end module

    
      module MD_MG_Global_C
      implicit none 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!       
      INTEGER(4) i_dir,nlv_glo  ,nlv_glomax  
      integer nnodeG,nnzG,nnodeC,nnzC
      INTEGER,DIMENSION(:),ALLOCATABLE :: imapG,imapGZ
      integer,dimension(:),allocatable :: iaG,jaG,juG
      REAL*8,DIMENSION(:),ALLOCATABLE :: eG,eG0,rG,rG0
      real*8,dimension(:),allocatable :: auG,auG0,aluG
      REAL*8,DIMENSION(:,:),ALLOCATABLE :: Ainv
! NEW
        INTEGER(4)nnziG
        INTEGER(4),DIMENSION(:),ALLOCATABLE ::iaiG,jaiG,jarG,iarG
        INTEGER(4),DIMENSION(:),ALLOCATABLE ::nnodeGC,nnzGC,nnziGC
        INTEGER(4),DIMENSION(:,:),ALLOCATABLE ::iaGC,jaGC,juGC
        INTEGER(4),DIMENSION(:,:),ALLOCATABLE ::iaiGC,jaiGC,iarGC,jarGC
        REAL(8),DIMENSION(:),ALLOCATABLE ::XintpG,XrestG,eGt,rGt
        REAL(8),DIMENSION(:,:),ALLOCATABLE ::eGC,rGC,XrestGC,XintpGC,auGC
        REAL(8),DIMENSION(:,:),ALLOCATABLE :: coordG
        INTEGER (4) inmaxG
        INTEGER (4),DIMENSION(:),ALLOCATABLE :: inmaxGC
        INTEGER (4) nnodecm,nnzicm,nnzcm
        
! FOR gather
        INTEGER (4) igather, nsengatA, nsengatR
        INTEGER (4), DIMENSION(:), ALLOCATABLE ::  irevgatA, irevgatR
        INTEGER (4), DIMENSION(:), ALLOCATABLE :: imapgatA,imapgatR, idispA, idispR
!---
      save
      end module