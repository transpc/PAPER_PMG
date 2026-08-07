      MODULE Znormal
!      
      IMPLICIT NONE
      SAVE
!
      INTEGER :: num_wallcells_l
      INTEGER :: maxmt_fluid_nbcon0
      INTEGER,DIMENSION(:),ALLOCATABLE :: i_neigh_nbcon0
      INTEGER,DIMENSION(:),ALLOCATABLE :: nji,nji_c,wall_cell_l
      REAL(8),DIMENSION(:),ALLOCATABLE :: sa_wallcell_l,sa_walll
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: xn,xn_wallcell_l,xn_wallcell
!
      END MODULE Znormal
      
