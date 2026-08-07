      MODULE Zbc_index
!      
      IMPLICIT NONE
      SAVE
!
!.....Scalar variables    
!
      INTEGER nvin,npin,i_horizontal_outlet
      INTEGER num_wallcells,num_wall_group,select_wall_group
!
!.....1D allocatable variables    
!
      INTEGER,DIMENSION(:),ALLOCATABLE :: nbcon_c,npb,ngrad,icell_type,iface_wall,wall_cell
      INTEGER,DIMENSION(:),ALLOCATABLE :: cellvin,cellpin,vin_norm,vin_mfr
      INTEGER,DIMENSION(:),ALLOCATABLE :: index_flux,index_property
      INTEGER,DIMENSION(:),ALLOCATABLE :: iface_wall0,iface_wall1
      INTEGER,DIMENSION(:),ALLOCATABLE :: icell_type_tmp,iface_wall_tmp
      INTEGER,DIMENSION(:),ALLOCATABLE :: nbcon,j_nbcon_old,nbcon_old
!
      LOGICAL l_horizontal_outlet_init 
!
      END MODULE Zbc_index
      
