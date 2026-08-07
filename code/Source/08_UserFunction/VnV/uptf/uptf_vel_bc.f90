!
      SUBROUTINE UPTF_vel_bc(time) 
!
!     This routine defines velocity boundary condition
!
      USE Vol_Data        , ONLY: cell
      USE Zzone           , ONLY: ncell_fluid
      USE Zcore           , ONLY: myrank      
      USE Znum_cell       , ONLY: i_neigh,indexr_sort
      USE Ztimecon        , ONLY: itim
      USE Zbc_index       , ONLY: nbcon
      USE Zbc_index       , ONLY: nvin
      USE Zconst2         , ONLY: iprn
      USE Zb_condition    , ONLY: vb_liq,vb_gas,vin_liq,vin_gas,p_fb, &
                                  tb_gas,tb_liq,tb_drp,eb_gas,eb_liq,eb_drp,rhob_gas,rhob_liq,rhob_drp,qualab
      USE Zncg            , ONLY: tao,cvao_nvin,uao_nvin,dcva_nvin,ra_nvin,qn_cell,n_ncg_sp
      USE Zuptf           , ONLY: cwl_limit,cwl_transition,cwl_time1,cwl_time2,      &
                                  cwl_tmp1,cwl_tmp2,cwl_tmp3,cwl_level1,cwl_level2,  &
                                  cwl_ratio,cwl_ratio0,cwl,                          &
                                  A_bottom,A_drain,v_drain
      USE Zvec_geo        , ONLY: sa_nf
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER i,j,j0,i1,err,err1
      INTEGER,SAVE:: lsj
      INTEGER,SAVE:: cell_drain      
      LOGICAL, SAVE :: initial=.TRUE.
      LOGICAL, SAVE :: cwl_initial=.TRUE.
      REAL(8) time,tmp1,tmp2
      REAL(8) ag_vb
!.....Local arrays
      REAL(8) qn_cell0(n_ncg_sp)
      REAL(8),DIMENSION(:),SAVE,ALLOCATABLE :: t_steam,t_vel,steam_in,steam_p,steam_t
      REAL(8),DIMENSION(:,:),SAVE,ALLOCATABLE :: velin,Pin,Tin
!
      IF(initial)THEN
         initial=.FALSE.      
         ALLOCATE(t_steam(22),t_vel(22),steam_in(22),steam_p(22),steam_t(22))
         ALLOCATE(velin(4,22),Pin(4,22),Tin(4,22))
!         
         OPEN(201,file='steam_in.dat',status='old',iostat=err)
         OPEN(202,file='vel_in.dat',status='old',iostat=err1)
         IF(myrank.eq.0)WRITE(*,*) 'UPTF data were loaded'
         IF(err.ne.0.or.err1.ne.0)then
            IF(myrank.eq.0)WRITE(*,*)'Input of UPTF 203 is missing'
            STOP
         ENDIF         
!
         IF(myrank.eq.0)OPEN(203,file='UPTF_inlet.dat')
         IF(myrank.eq.0)OPEN(204,file='UPTF_drain.dat')
!
         lsj=5
         DO i=1,lsj     ! Read boron concentration data
            READ(201,*) t_steam(i),steam_in(i) ,steam_p(i),steam_t(i)
         ENDDO
         DO i=1,lsj      ! Read inlet velocity data
            READ(202,*) t_vel(i),(velin(j,i),j=1,4),(Pin(j,i),j=1,4),(Tin(j,i),j=1,4)  
         ENDDO     
!
         cell_drain=0   
         DO i=1,ncell_fluid
            DO j=i_neigh(i),i_neigh(i+1)-1 
               IF(nbcon(j).eq.4) THEN
                  cell_drain=i
               ENDIF   
            ENDDO
         ENDDO          
!                        
      ENDIF     
!
!.....Interpolation
!
      vin_liq(:)=0.0d0
      vin_gas(:)=0.0d0
      p_fb(:)=0.0d0
      tb_liq(:)=0.0d0
      tb_gas(:)=0.0d0
      DO i=1,lsj-1
         IF(time.ge.t_steam(i).and.time.lt.t_steam(i+1))THEN
            j=5
            vin_gas(j)=-steam_in(i)-(time-t_steam(i))*(steam_in(i+1)-steam_in(i))/(t_steam(i+1)-t_steam(i))
            p_fb(j)=steam_p(i)+(time-t_steam(i))*(steam_p(i+1)-steam_p(i))/(t_steam(i+1)-t_steam(i))
            tb_liq(j)=steam_t(i)+(time-t_steam(i))*(steam_t(i+1)-steam_t(i))/(t_steam(i+1)-t_steam(i))
            tb_gas(j)=tb_liq(j)
            CYCLE
         ENDIF
      ENDDO        
!
      DO i=1,lsj-1
         IF(time.ge.t_vel(i).and.time.lt.t_vel(i+1))THEN
            DO j=1,nvin-1
               vin_liq(j)=-velin(j,i)-(time-t_vel(i))*(velin(j,i+1)-velin(j,i))/(t_vel(i+1)-t_vel(i))    !for vin_norm=1
               p_fb(j)=Pin(j,i)+(time-t_vel(i))*(Pin(j,i+1)-Pin(j,i))/(t_vel(i+1)-t_vel(i))
               tb_liq(j)=Tin(j,i)+(time-t_vel(i))*(Tin(j,i+1)-Tin(j,i))/(t_vel(i+1)-t_vel(i))
               tb_gas(j)=tb_liq(j)
!
!..............Drain control                       
!
               IF(cell_drain.ne.0.and.j.eq.4) THEN
                  p_fb(j)=cell%p(cell_drain)
                  tb_liq(j)=cell%tl(cell_drain)
                  tb_gas(j)=cell%tg(cell_drain)
               ENDIF               
!                  
            ENDDO
            CYCLE
         ENDIF
      ENDDO 
!
!.....Drain control        
!
      i=4   !4 means drainning velocity inlet index (only for the UPTF)
!      
!.....Initialize      
!
      vb_liq(i,1)=0.d0
      vb_liq(i,2)=0.d0
      vb_liq(i,3)=0.d0
!
      cwl_limit=0.4d0      !transient start [m]    for case1,2,3,6,8,9,10
      cwl_transition=0.4d0 !transient length [m]          
      cwl_time1=2.d0       !time1 to calculate water level incresing rate
      cwl_time2=3.d0       !time2 to calculate water level incresing rate    
!
!      cwl_limit=1.0d0      !transient start [m]    for case7
!      cwl_transition=0.5d0 !transient length [m]          
!      cwl_time1=10.d0       !time1 to calculate water level incresing rate
!      cwl_time2=11.d0       !time2 to calculate water level incresing rate    
!      
!      cwl_limit=0.2d0      !transient start [m]    for case4,5
!      cwl_transition=0.6d0 !transient length [m]          
!      cwl_time1=3.d0       !time1 to calculate water level incresing rate
!      cwl_time2=4.d0       !time2 to calculate water level incresing rate                
!
      IF(cwl_initial) THEN
         cwl_tmp1=1
         cwl_tmp2=1
         cwl_tmp3=1
         cwl_initial=.false.
         cwl=0.d0
      ENDIF   
      IF(time.gt.cwl_time1.and.cwl_tmp1.eq.1) THEN
         cwl_level1=cwl
         cwl_time1=time
         cwl_tmp1=0
      ELSEIF(time.gt.cwl_time2.and.cwl_tmp2.eq.1) THEN
         cwl_level2=cwl
         cwl_time2=time
         cwl_tmp2=0
         cwl_initial=.false.
         cwl_ratio=(cwl_level2-cwl_level1)/(cwl_time2-cwl_time1)
         cwl_ratio0=cwl_ratio         
      ENDIF
!
      IF(cwl.gt.cwl_limit) THEN
!........A_bottom, A_drain
         A_bottom=3.14d0*2.43d0*2.43d0 !radius=2.43m
         A_drain=0.32782065581820d0    !by pre-calculation
         ag_vb=0.d0
! 
!........Drain area   
         DO i=1,ncell_fluid
!..............neighbor have been sorted use indexr_sort to retrieve j=6
               j0=i_neigh(i)-1
               j=indexr_sort(6+j0)
            IF(nbcon(j+j0).eq.4) THEN  !j=6 in this problem only
               CALL get_vector_disp(j,i,i1)
               i1=ABS(i1)
               A_drain=sa_nf(i1)
               ag_vb=cell%alphag(i)
               EXIT
            ENDIF
         ENDDO
!
         IF(cwl.le.cwl_limit+cwl_transition) THEN
            cwl_ratio=cwl_ratio0*(cwl-cwl_limit)/cwl_transition
         ENDIF   
!
         v_drain=cwl_ratio*A_bottom/A_drain
         v_drain=v_drain*DMIN1(cwl_transition,cwl-cwl_limit)/cwl_transition/DMAX1(1.d-1,1.d0-ag_vb)  
         vb_liq(4,3)=-v_drain     !3=z-direction, 4=inlet number   
         vb_gas(4,3)=-v_drain     !3=z-direction, 4=inlet number   
      ENDIF
!      
      IF(MOD(itim,iprn).eq.0.and.myrank.eq.0) WRITE(204,1001) time,cwl,vb_liq(4,3)
!
      IF(myrank.eq.0) WRITE(203,1001)time,(vin_liq(i),i=1,nvin-1),vin_gas(5)
1001     FORMAT(1x,6f15.5) 
!
!.....Calculate time dependant rho & e for inlet conditions
!
      DO i=1,nvin
         IF(cell_drain.eq.0) CYCLE
!      
!........Liquid, steam/gas      
!
         qn_cell0(:)=qn_cell(i,:)
         CALL convert_temp2erg(p_fb(i),tb_liq(i),tb_gas(i),qualab(i),eb_liq(i),eb_gas(i),rhob_liq(i),rhob_gas(i),tmp1,tmp2, &
                               tao,cvao_nvin(i),uao_nvin(i),dcva_nvin(i),ra_nvin(i),qn_cell0)         
!
!........Droplet property is same as liquid
!
         tb_drp(i)  =tb_liq(i)
         rhob_drp(i)=rhob_liq(i)
         eb_drp(i)  =eb_liq(i)         
      ENDDO 
!
      END SUBROUTINE UPTF_vel_bc     
