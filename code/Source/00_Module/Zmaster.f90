!
!     modules needed for linking to MASTER
!
!.....MASTER4 
!
      MODULE MASTER4
         INTEGER i_flag,iok
         INTEGER NCB,NXY_TH,NZ_TH
         INTEGER NXYF,NPINX
         INTEGER rv_model_master
         CHARACTER(4),ALLOCATABLE :: CBNAM(:)
         REAL(8),ALLOCATABLE :: TFC(:,:),TFS(:,:),TCOO(:,:),DCOO(:,:),BCOO(:,:)
         REAL(8),ALLOCATABLE :: ZCB(:)
         REAL(8),ALLOCATABLE :: P3D_TH(:,:),PIN3D_TH(:,:,:,:)
         REAL(8),ALLOCATABLE :: VOL_TH(:,:)
         REAL(8) TTIME,DTTR,PWTH,PPM,PPCT_MASTER
         REAL(8) power_master
         CHARACTER(132) RSTFN_INP,RSTFN_OUT,RSTFN_REV
!
!........scram
!
         INTEGER mas_rx_trip         
         INTEGER ncb_mas
         CHARACTER*3,ALLOCATABLE:: cbnam_mas(:)
         REAL(8),ALLOCATABLE:: zcb_mas(:) 
         REAL(8) ppm_mas
         REAL(8) time_cri
!         
!........master control from initialization mode to operation mode.
!
         INTEGER dtemp_pass,dpower_pass,dmass_pass,dmaster_pass       
         INTEGER mas_dtemp_opt,mas_dpower_opt,mas_dmass_opt      
         INTEGER count_dtemp,count_dmass,count_dpower 
         REAL(8) mas_wait,mas_interval,mas_delay,mas_interval_delay
         REAL(8) mas_init_duration
         REAL(8) mas_dtemp,mas_dpower,mas_dmass  
         !REAL(8) cpc_vopt
         REAL(8),ALLOCATABLE:: dmass(:),dpres(:)       
!
!........new mapping
!
         INTEGER nchn,num_asm
         INTEGER,ALLOCATABLE::mst_to_asmi(:)
!
      END MODULE
