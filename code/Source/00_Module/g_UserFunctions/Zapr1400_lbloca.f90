      Module Zapr1400_lbloca
! 
      IMPLICIT NONE
      SAVE
!    
      INTEGER,PARAMETER :: m_max=5 !4 DVIs and 1 cold leg as a break
!
      INTEGER,DIMENSION(:),ALLOCATABLE :: num_loca
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: icell_loca,jth_loca
      REAL(8),DIMENSION(:),ALLOCATABLE :: p_loca,rhom_loca,q_loca,h_loca,cpg_loca,vg_loca,area_loca
      REAL(8),DIMENSION(:),ALLOCATABLE :: vol_loca_tot
      REAL(8),DIMENSION(:),ALLOCATABLE :: vol_l_tot,vol_g_tot
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: vol_loca,vin_liq_loca
!
      REAL(8) :: mflux_sit(m_max-1),mflux_sip(m_max-1),mflux_dvi(m_max-1),pres_break,pct,pres_ambient      
!
      REAL(8) :: sit_mass,sit_mass_spipe,sit_pre
      INTEGER :: sit_detect,sit_avail(m_max)
      REAL(8) :: hpsip_pre,hpsip_delay
      INTEGER :: hpsip_detect,hpsip_avail(m_max)
      REAL(8) :: mflux_dvi_int,mflux_sit_int,sit_time,hpsip_time
      INTEGER :: break_opt
      REAL(8) :: flux_break
      REAL(8) :: tl_si,tg_si,mflux_sit_phase1,mflux_sit_phase2
!.....Henry-Fauske Critical Flow
      REAL(8),DIMENSION(:),ALLOCATABLE :: ar_liq_loca,ar_gas_loca,pps_loca ,rhog_loca  ,eg_loca    ,el_loca    
      REAL(8),DIMENSION(:),ALLOCATABLE :: quala_loca ,quals_loca ,tl_loca  ,alphag_loca,alphal_loca,cvao_loca ,dcva_loca
      REAL(8),DIMENSION(:),ALLOCATABLE :: uao_loca   ,ra_loca    ,ul_o_loca,ug_o_loca     
      REAL(8),DIMENSION(:),ALLOCATABLE :: vly_o_loca ,vgy_o_loca 
!.....a hole to ambient
      INTEGER rpv_status
      INTEGER,DIMENSION(:),ALLOCATABLE :: num_ambient
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: icell_ambient,jth_ambient
!.....user_def_inp
      REAL(8):: tbreak_s,tbreak_f,tbreak_m,dt_min,pbnd_s,pbnd_f,pbnd_m 
      REAL(8):: time_restart
      REAL(8):: topenleg_s,topenleg_f      
      REAL(8),DIMENSION(:),ALLOCATABLE :: vin_liq_init
!.....time dependent SIT mass flow rate
      INTEGER :: nsit
      REAL(8),DIMENSION(:),ALLOCATABLE :: sit_mrate_time,sit_mrate
!.....boundary condition
      INTEGER,ALLOCATABLE :: nbcon_cell(:)
!
      END MODULE Zapr1400_lbloca     
