      module MD_MG_matrix
      implicit none 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!         CSR matrix for coarse level : P, R, Ac  
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! matrix P & R:
    
      integer nnzi1,nnzi2
      integer,dimension(:),allocatable :: iai1,jai1,iar1,jar1
      integer,dimension(:),allocatable :: iai2,jai2,iar2,jar2
      
      real*8,dimension(:),allocatable:: Xintp1,Xrest1
      real*8,dimension(:),allocatable:: Xintp2,Xrest2
      
! matrix Ac: 
    
      integer nnz1,nnz2
      
      integer,dimension(:),allocatable :: ia1,ja1,ju1
      integer,dimension(:),allocatable :: ia2,ja2,ju2
      
      real*8,dimension(:),allocatable :: au1,au2,alu_mg
      real*8,dimension(:),allocatable :: diagr1,diagr2 
      real*8,dimension(:),allocatable :: diagt
      real*8,dimension(:),allocatable :: diagrc      
      
! NEW
      REAL(8),DIMENSION(:),ALLOCATABLE :: r,rt,rc,rs,e,et,es
      REAL(8),DIMENSION(:),ALLOCATABLE :: auc,aus,Xrest,Xintp
      INTEGER(4),DIMENSION(:),ALLOCATABLE :: iac,jac,juc,ias,jas,jus,               &
                                             iai,jai,iar,jar
      INTEGER(4):: nnzc0,nnzi,nnzr,nnzs
! temp
      INTEGER(4), DIMENSION(:), ALLOCATABLE ::  iat, jat
      REAL(8), DIMENSION(:), ALLOCATABLE :: aut
!---
      save
      end module
	  
