MODULE Zrv_choke
!
   IMPLICIT NONE
   SAVE
!
   INTEGER :: ick,ick_dir  
   INTEGER :: ichokingflow
   INTEGER :: choke_model_opt,choke_update
   INTEGER :: choke_cell_num,choke_axis
   INTEGER :: num_throatface     !total number of throat faces
   INTEGER :: cellnum_throatface
   INTEGER :: icell_throatface(1000)
   INTEGER, ALLOCATABLE :: choke_cell(:)
   INTEGER, ALLOCATABLE :: n_face_throat(:),dir_face_throat(:)
!
!..Finding throat faces
!  
   INTEGER, ALLOCATABLE :: fzone_throat(:,:) !fluid zones of each throat (1)=(left zone), (2)=(right zone)
   INTEGER, ALLOCATABLE :: icell_throat(:)   !inner cell number of throat faces (left zone)
   INTEGER, ALLOCATABLE :: ocell_throat(:)   !outter cell number of throat faces (right zone)
   INTEGER :: env_press_option   
!
!..nonk array for nf_nonk array
!      
   INTEGER, ALLOCATABLE :: nonk_throat(:)   
!
!..upstream value
!
   REAL(8) cpgas,cppf0,cppg0,gamma,ploss,pvzero,             &
            pzero,rgas,rnc,sgas,sliq,svap,tzero,vsubf0,vsubg0,xnc,xzero      
   REAL(8) :: choke_throat_area,choke_throat_ratio  !m2
   REAL(8) :: choke_pout
   REAL(8) :: throat_area
!
!..relaxation
!   
   REAL(8) relax_choke
   
   REAL(8) vl_choke,vg_choke,vl_choke_o,vg_choke_o
   REAL(8),ALLOCATABLE:: ug_throatface(:),ul_throatface(:),ud_throatface(:)
   REAL(8),ALLOCATABLE:: fluxl_throatface(:),fluxg_throatface(:),fluxd_throatface(:)
!
!..average values
!   
   REAL(8) cell_leng_avg,ar_liq_avg,ar_gas_avg,p_avg,p_avg_out,pps_avg,vl_o_avg,vg_o_avg,     &
            rhog_avg,eg_avg,el_avg,quala_avg,quals_avg,tl_avg,alphag_avg,alphal_avg, &
            cvao_cell_avg,dcva_cell_avg,uao_cell_avg,ra_cell_avg,theta_avg
!
!..choke
!       
   LOGICAL choke    
!   
!...input  
!   
   CHARACTER(30) s_critical_flow      
   REAL(8) :: time_cflow_on
!   
!...apr1400_lbloca
!
   INTEGER, ALLOCATABLE :: idx_throatface(:),dir_throatface(:)  
  
!
END MODULE Zrv_choke
