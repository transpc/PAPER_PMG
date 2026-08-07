MODULE Zpdrop
!
   IMPLICIT NONE
   SAVE
!
   INTEGER :: npid
   INTEGER, ALLOCATABLE :: num_dp_region(:,:)
!
   REAL(8) :: time_pid_on,time_pid_off
   REAL(8),ALLOCATABLE :: err(:),Pcon(:),Dcon(:)
   REAL(8),ALLOCATABLE :: err_o(:),Icon(:),dp_set(:),dp_control(:)
   REAL(8),ALLOCATABLE :: dp_region(:)
   
!
END MODULE Zpdrop
