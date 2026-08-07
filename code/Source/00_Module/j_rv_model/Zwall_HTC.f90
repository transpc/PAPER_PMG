      MODULE Zwall_HTC
!
      IMPLICIT NONE
      SAVE
!        
      REAL(8) tw,tsat_t,dt_sat
      REAL(8) qflux_t,qflux_l,qflux_g
      REAL(8) chf,chf_mul
      REAL(8) rho,tl,cond,viscos,cp,beta,sigma
      REAL(8) h_liq,hfg,hfg_p,sat_hfp,qual_eq
      REAL(8) mflux_liqa,mflux_gasa,mflux_tota,vfg
      REAL(8) zqf,zqf_min,zqf_top
      REAL(8) HTC_t,HTC_tl,HTC_tst,HTC_tg,HTC_tgp,HTC_tgt,HTC_cond,HTC_d,pvblk
      INTEGER, ALLOCATABLE::mode(:),h_index(:),mode_1d(:)
      INTEGER, ALLOCATABLE::hmode_rv(:,:)
      REAL(8), ALLOCATABLE::twall_rv(:,:)
      REAL(8), ALLOCATABLE::dt_sat_rv(:),chf_rv(:) !pik-rv-debug      
      REAL(8), ALLOCATABLE::Mul_o(:),Mul1d_o(:)
      REAL(8), ALLOCATABLE::chfr(:)
!
      INTEGER inline_bundle,reflood,incnd
!      INTEGER, PARAMETER :: inline_bundle=0,reflood=1,rezone=1
      INTEGER, PARAMETER :: iter_tw=0
      INTEGER, PARAMETER :: f_direc=3               ! Main flow direction for using 1-D based correl. (3=z direction)              (!!!cyj 향후 자동화 방법 고려)
      INTEGER, PARAMETER :: c_direc1=1,c_direc2=2   ! Cross flow direction                            (1&2=x,y for cross flow)     (!!!cyj 향후 자동화 방법 고려)   
!
!      REAL(8), PARAMETER :: l_plate = 3.657d0         ! Length of rectangle to calculate HTC for natural convection by Chulchill-Chu
!      REAL(8), PARAMETER :: dia_rod = 0.009496d0, pit_dia=1.3262d0     ! Ratio of pitch and Diameter of fuel rod                    
!      REAL(8), PARAMETER :: h_bundle =3.996d0       ! 1.0d0 Axial length of bundle for Reflood calculation                          
!      REAL(8), PARAMETER :: base_bundle =0.0d0      ! height of the botton of bundle; base height   
      REAL(8)  hyd_core,l_plate,dia_rod,pit_dia,h_bundle,base_bundle   
      REAL(8)  dia_rod_cfd
!
!.....parameters for CHF calc.
!
      INTEGER, PARAMETER :: horiz_chf=0             ! Horizontal bundle option for CHF              (!!!cyj 향후 입력으로 처리)
      REAL(8), PARAMETER :: angle_horiz= 0.0d0      ! angle for horizontal bundle option            (!!!cyj 향후 계산 서브루틴 추가) 
      REAL(8), PARAMETER :: k_grid = 0.0d0          ! Loss coefficient for grid spacer              (!!!cyj 향후 입력으로 처리) 
      REAL(8), PARAMETER :: dis_grid = 0.0d0        ! Distance to grid spacer                       (!!!cyj 향후 입력으로 처리) 
      REAL(8), PARAMETER :: f_pp_axial = 1.0d0      ! Power peaking factor for axial direction      (!!!cyj 향후 입력으로 처리)
!
      LOGICAL reflod      
      REAL(8) rey_reflod      
      ! OPR1000 rod-scale 
      REAL(8),ALLOCATABLE:: qflux_l0(:),qflux_g0(:)
      real(8),ALLOCATABLE:: qf00(:),qf01(:)
      REAL(8), ALLOCATABLE::gamma_wall_rod(:)
      REAL(8),allocatable:: qf0(:),qf1(:)
      REAL(8),allocatable:: qg0(:),qg1(:)
      REAL(8),allocatable:: gw0(:),gw1(:)      
!
      END MODULE Zwall_HTC