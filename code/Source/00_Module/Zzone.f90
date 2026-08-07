      MODULE Zzone
!      
      IMPLICIT NONE
!
      INTEGER,PARAMETER :: num_max_zone=20
      INTEGER :: num_fzone,num_szone,ncell_fluid,num_pzone,ncell_porous,num_rzone,  &
                 ncell_rod,ncell_fluid_all,ncell_fluid_f,ncell_cond,ncell_cond_all, &
                 ncell_fluid_pad,ncell_fluid_padv,ncell_cond_pad,ncell_cond_padv
      INTEGER,DIMENSION(:),ALLOCATABLE :: nzone,nmaterial,nmaterial_c,nzone_c,icore
!
      END MODULE Zzone
      
