MODULE Zvalve
!
   IMPLICIT NONE
   SAVE
!
   INTEGER :: valve_model
   INTEGER :: num_valveloc,num_valve                  !number of valve faces (local domain, global domain)
   INTEGER, ALLOCATABLE :: fzone_valve(:,:)           !fluid zones surrounding a valve face
   INTEGER, ALLOCATABLE :: num_valveface(:)
   INTEGER, ALLOCATABLE :: n_face_valve(:,:)          !non-face index of tt-th valve (1:num_valve,1:num_valveface(tt))   
   INTEGER, ALLOCATABLE :: mapping_valve(:)
   INTEGER, ALLOCATABLE :: icell_valve(:,:)           !fluid zones surrounding a valve face   
   INTEGER, ALLOCATABLE :: ocell_valve(:,:)           !fluid zones surrounding a valve face   
   INTEGER, ALLOCATABLE :: valve_closed(:)            !a valve is closed (1:num_valve) 1=closed, 0=open
   REAL(8), ALLOCATABLE :: time_valve_closed(:,:)     !ith valve is closed during time_valve_close(i,1)~time_valve_close(i,2)

   REAL(8) :: sa_nf_o(5),sap_nf_o(5),saa_nf_o(5),sad_non_o(5)
!
!..nonk array for nf_nonk array
!      
   INTEGER, ALLOCATABLE :: nonk_valve(:,:)
   
!
END MODULE Zvalve
