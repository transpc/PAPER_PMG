      MODULE Zvec_index
! 
      IMPLICIT NONE
      SAVE
!
      INTEGER,DIMENSION(:),ALLOCATABLE :: left_nf,jneigh_nf,nbcon_nf
      INTEGER,DIMENSION(:),ALLOCATABLE :: right_non,kneigh_non      
      INTEGER,DIMENSION(:),ALLOCATABLE :: right_fsw
!
      REAL(8),DIMENSION(:),ALLOCATABLE :: delta_npb
!
      END MODULE Zvec_index
