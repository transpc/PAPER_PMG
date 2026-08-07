!
      SUBROUTINE user_def_output(nout)
!
!     This routine controls user/problem specific output variables & format
!
      USE VOL_DATA        , ONLY: cell
      USE Wall_DATA       , ONLY: face      
      USE Zmpi            , ONLY: ncell_fp
      USE Zzone           , ONLY: ncell_fluid,ncell_fluid_all,ncell_cond_all
      USE Zcore           , ONLY: np,myrank,cupid_mars  
      USE Zparam          , ONLY: nn,ndim,nb_max,nin_max
      USE Zb_condition    , ONLY: pbnd,alphab_liq,rhob_liq,alphab_gas,rhob_gas,alphab_drp,rhob_drp
      USE Zbc_index       , ONLY: index_flux,index_property,npb
      USE Zconst1         , ONLY: vv_prob,restart
      USE Zconst2         , ONLY: dt,iprn,grav
      USE Zcoord1         , ONLY: xloc,xloc_tmp
      USE Zcoord2         , ONLY: prnvar
      USE Zcoord3         , ONLY: vol,volp
      USE Zmass_conv      , ONLY: tot_mass
      USE Zporous         , ONLY: chn_type_tmp
      USE Zpress          , ONLY: p
      USE Zrocom_specific , ONLY: tl_space_min,tl_space_avg,in_mdot,out_mdot
      USE Ztimecon        , ONLY: time,t_end,toutstep,t_end_ctrl,itim
      USE Zmodel          , ONLY: qconden,wVertical
      USE Zvector         , ONLY: vl_n,vg_n,vl_o,vg_o
      USE Zwall_HTC       , ONLY: twall_rv
      USE Zapr1400_lbloca , ONLY: tbreak_s,pbnd_f,pres_break,pct,mflux_sit,mflux_sip,mflux_dvi,m_max,mflux_sit_int,mflux_dvi_int 
      USE Zsbloca         , ONLY: break_flow,si_flow,break_flow_int,si_flow_int 
      USE Zrv_mpi         , ONLY: ncell_fuel_rod_p
      USE Zrv_ncell       , ONLY: ncell_fluid_core,ncell_fuel_rod,ncell_fuel_rod_all
      USE Zrv_ncell       , ONLY: chn_nx,chn_ny,chn_nz
      USE Zrv_hts_2d      , ONLY: power_2d,t_fuel,nr_2d,ri_2d
      USE Zrv_choke       , ONLY: vl_choke,vg_choke,vl_o_avg,choke,choke_update
      USE Zrv_gap_cond    , ONLY: h_gap,width_gap      
      USE Zvec_major      , ONLY: flux_l_nf,flux_g_nf,flux_d_nf
      USE Zio_unit        , ONLY: unit_log
      USE Zmars           , ONLY: mass_nf_mcc
      USE Ztplot          , ONLY: tplot_cell_loc
      USE Zvector         , ONLY: ul_o
!      
      IMPLICIT NONE
!
      INCLUDE '../10_LinkToMARS/c3com.h'
!
!.....Input
      INTEGER :: nout
!.....Local variables
      INTEGER :: i,j,k,ii1,ii2
      INTEGER :: prnopt,fdim,sdim     
      INTEGER :: na,na_c,na_r
      INTEGER :: ii,num                                                  !wh14x28
      INTEGER :: ichok  
      INTEGER,SAVE :: nloop
      INTEGER,SAVE :: nplot           
      INTEGER,SAVE :: n       
      LOGICAL,SAVE :: initial=.true.
! 
      REAL(8) :: choking_onoff
      REAL(8) :: x1
      REAL(8) :: level1,level2,level
      REAL(8) :: vel_an, temp_an,temp_ave,tw                                      !VFS7
      REAL(8) :: vel                                                              !pnl7x7
      REAL(8) :: height
      REAL(8) :: pct_tmp,pres_break_tmp,v1_tmp,v2_tmp
      REAL(8) :: temp
      REAL(8),SAVE :: time1
      REAL(8),SAVE :: print_time,print_interval,CALL_time,stavg,ftavg
!.....Local arrays
      CHARACTER*20 :: varname(25)
      INTEGER,SAVE :: cellindex(100,2)
      INTEGER,SAVE :: pnl_height(25),ce_index(6,8),wh_index(7,29)                 !pnl7x7,ce15x15
      INTEGER,SAVE :: cellNum(20)
      REAL(8) :: vel_ce(10),vel_wh(7)                                             !ce15x15,wh14x28
      REAL(8) :: mflux_sit_tmp(m_max-1),mflux_sip_tmp(m_max-1),mflux_dvi_tmp(m_max-1)
      REAL(8) :: eq_av,eq_corner,eq_side,eq_center,mf_av,mf_corner,mf_side,mf_center          !eq=Exit quality
      REAL(8) :: volsum(4),qualsum(4),mfluxg_sum(4),mfluxl_sum(4),mflux_g, mflux_l            !GE3x3
      REAL(8),SAVE :: vy(500,52)
      REAL(8),SAVE :: temp_krane(100),v_krane(100),u_krane(100)                   !for VFTX
      REAL(8),SAVE :: zzkuhn(13)
      REAL(8),SAVE :: data_gap(19)
!     
!     CUBE ST2-CT-01
      INTEGER,SAVE:: ttnum
      REAL(8),SAVE:: ttrec(9000)
      LOGICAL, SAVE::initialtt
      DATA initialtt /.true./
!
!.....Local allocatable arrays
      !OPR1000 fullvessel
!     REAL(8) :: pave(num_fzone),tave(num_fzone)
!     REAL(8) :: tave_assem(15,15,26),pave_assem(15,15,26)
      INTEGER,DIMENSION(:),ALLOCATABLE :: tempi1
      REAL(8),DIMENSION(:),ALLOCATABLE :: temp1,temp2,temp3,temp4,temp5
      REAL(8),DIMENSION(:),ALLOCATABLE :: temp6,temp7,temp8,temp9
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: temp0
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: t_fuel_all                            !gap conductance test
     !PNNL
      INTEGER cx,cy,cz
      REAL(8) ul(ncell_fluid)
!
      DATA varname/'time','p_dvi','pct','break_flow','break_flow_int','si_flow','si_flow_int',                        &
                                         'mflux_sit1','mflux_sit2','mflux_sit3','mflux_sit4',                         &
                                         'mflux_sip1','mflux_sip2','mflux_sip3','mflux_sip4',                         &
                                         'mflux_dvi1','mflux_dvi2','mflux_dvi3','mflux_dvi4',                         &
                                         'power','speed_choke','speed_pboun','mflux_sit_int','mflux_dvi_int','ichok'/
!        
      DATA data_gap/12.2,23.6,37.0,47.7,55.7,61.6,65.9,66.2,67.3,69.1, &
                    66.7,63.0,58.6,53.8,48.4,41.5,32.5,20.9,13.6/
!
      na=ncell_fluid_all
      na_c=ncell_cond_all
      na_r=ncell_fuel_rod_all
!
!.....Default: write the restart file and screen output
!
      CALL save_output(nout)
      CALL prn_m_err(tot_mass)
!
  11  FORMAT(i10,3e15.7)
  20  FORMAT(i10,3e15.7)
  300 FORMAT(3(e14.7,1x))
  400 FORMAT(10(e14.7,1x))  
!  
!-----VnV output-----------------------------------------------------------------------------  
!
!
!.........LSJ
!
!.....Gap conductance model Test: calculate only for single core
!
      IF(vv_prob.eq.'gap_conductance') THEN
!      
         IF(myrank.eq.0) THEN
            ALLOCATE(t_fuel_all(na_r,nr_2d)) 
         ELSE
            ALLOCATE(t_fuel_all(1,nr_2d)) 
         ENDIF
         CALL gatherv_r_2d(t_fuel,ncell_fuel_rod_p,t_fuel_all,ncell_fuel_rod,na_r,3)
!
         IF(time.gt.t_end) THEN
!
            IF(myrank.eq.0) THEN
               ALLOCATE(temp0(na_r,2)) 
            ELSE
               ALLOCATE(temp0(1,2)) 
            ENDIF
            CALL gatherv_r(h_gap    ,ncell_fuel_rod,temp0(1,1),na_r,3)
            CALL gatherv_r(width_gap,ncell_fuel_rod,temp0(1,2),na_r,3)
!         
            IF(myrank.eq.0) THEN
               OPEN(433, file='VFTX_ref-centerT.dat')
               OPEN(434, file='VFTX_ref-htcGap.dat')
               OPEN(435, file='VFTX_ref-widthGap.dat')
               OPEN(436, file='VFTX_ref-radialNode2.dat')         
               OPEN(437, file='VFTX_ref-radialNode9.dat')   
!
               j=1
               DO i=1,19
                  temp=(t_fuel_all(j,1)+t_fuel_all(j+1,1)+t_fuel_all(j+2,1)) /3.d0    !centerline temperature 
                  WRITE(433,334) data_gap(i),temp
                  j=j+3
               ENDDO
!
               j=1
               DO i=1,19
                  temp=(temp0(j,1)+temp0(j+1,1)+temp0(j+2,1)) /3.d0
                  WRITE(434,334) data_gap(i),temp
                  j=j+3
               ENDDO
!
               j=1
               DO i=1,19
                  temp=(temp0(j,2)+temp0(j+1,2)+temp0(j+2,2)) /3.d0
                  WRITE(435,334) data_gap(i),temp
                  j=j+3
               ENDDO
!            
               DO i=1,nr_2d
                  temp=(t_fuel_all( 4,i)+t_fuel_all( 5,i)+t_fuel_all( 6,i))/3.d0      !node2 with 3divs
                  WRITE(436,334) ri_2d(i),temp
                  temp=(t_fuel_all(25,i)+t_fuel_all(26,i)+t_fuel_all(27,i))/3.d0      !node9 with 3divs
                  WRITE(437,334) ri_2d(i),temp
               ENDDO  
!
               CLOSE(433)            
               CLOSE(434)             
               CLOSE(435)             
               CLOSE(436)
               CLOSE(437)
            ENDIF   
            DEALLOCATE(temp0)           
         ENDIF   
! 
!........For comparison to np4 result   
!
         IF(myrank.eq.0) THEN
            IF(np.eq.1) THEN
               OPEN(431, file='VFTX_ref-compare_np1.dat') 
            ELSE
               OPEN(431, file='VFTX_ref-compare_np4.dat')
            ENDIF    
            DO i=1,na_r
               WRITE(431,311) i, t_fuel_all(i,1) 
            ENDDO 
            CLOSE(431)
         ENDIF   
         DEALLOCATE(t_fuel_all)  
311      FORMAT(i3,2x,e14.7)          
!
      ELSEIF(vv_prob.eq.'choking_edward') THEN
         IF(initial) THEN
            IF(myrank.eq.0) THEN
               OPEN(333, file='VFTX_ref.dat')
               WRITE(333,"(a)") 'time    alphag_11    p_11    vl_29     vg_29   vl_choke  vg_choke alphag_30   choking_onoff vl_30  vg_30'
            ENDIF
            print_time=0.d0
            initial=.false.
         ENDIF
! 
!...........Compare
!
         IF(myrank.eq.0) THEN
            IF(np.eq.1) THEN
               OPEN(232,file='VFTX_ref-compare_np1.dat')
            ELSE
               OPEN(232,file='VFTX_ref-compare_np4.dat')
            ENDIF
         ENDIF   
         
         IF(time.ge.print_time)THEN
            print_time=print_time+0.001d0
            IF(myrank.eq.0) THEN
               ALLOCATE(temp0(na,4))
            ELSE
               ALLOCATE(temp0(1,4))
            ENDIF
            CALL gatherv_r(cell%alphag,ncell_fluid,temp0(1,1),na,0)
            CALL gatherv_r(cell%p     ,ncell_fluid,temp0(1,2),na,0)
            CALL gatherv_r(vl_o(1,1)  ,ncell_fluid,temp0(1,3),na,0)
            CALL gatherv_r(vg_o(1,1)  ,ncell_fluid,temp0(1,4),na,0)
!
!...........Compare
!
            IF(myrank.eq.0) THEN
               DO i=1,ncell_fluid_all
                  WRITE(232,311) i,temp0(i,1)   
               ENDDO         
               CLOSE(232)
!         
               IF(choke) THEN
                  choking_onoff=1.d0
               ELSE
                  choking_onoff=0.d0
               ENDIF
               WRITE(333,334) time,temp0(11,1),temp0(11,2),temp0(29,3),temp0(29,4),vl_choke,vg_choke,temp0(30,1),choking_onoff,temp0(30,3),temp0(30,4)
            ENDIF         
            DEALLOCATE(temp0)
         ENDIF
         
      ELSEIF(vv_prob.eq.'mcp_ftype'.or.vv_prob.eq.'mcp_ftype2') THEN
         IF(initial) THEN
            IF(myrank.eq.0) THEN
               OPEN(333, file='VFT_ref.dat')
            ENDIF
            print_time=0.d0
            initial=.false.
         ENDIF
         IF(myrank.eq.0) THEN
            ALLOCATE(temp0(na,4))
         ELSE
            ALLOCATE(temp0(1,4))
         ENDIF         
         CALL gatherv_r(cell%p     ,ncell_fluid,temp0(1,2),na,0)
         CALL gatherv_r(ul_o       ,ncell_fluid,temp0(1,3),na,0)         
         IF(myrank.eq.0) WRITE(333,990) time,temp0(1,2),temp0(36,2),temp0(1,3),temp0(36,3)
990      FORMAT(4x,100(e15.7,1x))            
         DEALLOCATE(temp0)  
!
      ELSEIF(vv_prob.eq.'blowtest2') THEN
         IF(initial) THEN
            IF(myrank.eq.0) THEN
               OPEN(333, file='VFTX_ref.dat')
            ENDIF
            print_time=0.d0
            initial=.false.
         ENDIF
!
         IF(time.ge.print_time)THEN
            print_time=print_time+0.0005d0
            IF(myrank.eq.0) THEN
               ALLOCATE(temp0(na,4))
            ELSE
               ALLOCATE(temp0(1,4))
            ENDIF
            CALL gatherv_r(cell%alphag,ncell_fluid,temp0(1,1),na,0)
            CALL gatherv_r(cell%p     ,ncell_fluid,temp0(1,2),na,0)
            CALL gatherv_r(vl_o(1,1)  ,ncell_fluid,temp0(1,3),na,0)
            CALL gatherv_r(vg_o(1,1)  ,ncell_fluid,temp0(1,4),na,0)
!
            IF(myrank.eq.0) THEN
                WRITE(333,444) time,temp0(95,1),temp0(100,1),temp0(271,1),temp0(295,1), &
                                  temp0(95,2),temp0(100,2),temp0(271,2),temp0(295,2), &
                                  temp0(95,3),temp0(100,3),temp0(271,3),temp0(295,3), &
                                  temp0(95,4),temp0(100,4),temp0(271,4),temp0(295,4), &
                                  vl_choke,vg_choke,choke_update
            ENDIF         
            DEALLOCATE(temp0)
444          FORMAT(19(e22.15,1x),(i3,3x))
         ENDIF           
!        
      ELSEIF(vv_prob.eq.'nat_conv_krane') THEN                       
!
         IF(time.gt.t_end) THEN
            IF(myrank.eq.0) THEN
               ALLOCATE(temp1(na),temp0(na,ndim))
            ELSE
               ALLOCATE(temp1(1),temp0(1,ndim))
            ENDIF         
            CALL gatherv_r(cell%tg  ,ncell_fluid,temp1,na,0)
            CALL gatherv_r_2d(vg_n,ncell_fp,temp0,ncell_fluid,na,0)
!
            IF(myrank.eq.0) THEN
               OPEN(335, file='list.txt')
               DO i=1,100
                  READ(335,*) cellindex(i,1),cellindex(i,2)
               ENDDO
               CLOSE(335)           
               OPEN(333, file='VFTX_ref.dat')
               DO i=1,100
                  ii1=cellindex(i,1)
                  ii2=cellindex(i,2)
                  temp_krane(i)=temp1(ii1)
                  v_krane(i)=temp0(ii1,2)
                  u_krane(i)=temp0(ii2,1)   
!
!             non-dimensionalization
!
                  temp_krane(i)=(temp_krane(i)-283.d0)/20.d0
                  v_krane(i)=v_krane(i)/dsqrt(9.81d0*3.41*0.001d0*20.d0*0.04493d0)
                  u_krane(i)=u_krane(i)/dsqrt(9.81d0*3.41*0.001d0*20.d0*0.04493d0)
                  WRITE(333,100) 0.01d0*(i-1), temp_krane(i),v_krane(i),u_krane(i)
               ENDDO
               CLOSE(333)
            ENDIF
            DEALLOCATE(temp1,temp0)
         ENDIF
!                   
      ELSEIF(vv_prob.eq.'PAFS-POOL') THEN
         IF(MOD(itim,iprn).eq.0)THEN
            IF(initial) THEN
               IF(myrank.eq.0)THEN 
!                 OPEN(30, file='plot.out')
                  OPEN(333, file='VD7_ref.dat')
               ENDIF   
               initial=.false.
            ENDIF   
!
            IF(myrank.eq.0)THEN 
               ALLOCATE(temp0(na,4))
            ELSE
               ALLOCATE(temp0(1,4))
            ENDIF
            CALL gatherv_r(cell%alphal,ncell_fluid,temp0(1,1),na,0)
            CALL gatherv_r(cell%rhol  ,ncell_fluid,temp0(1,2),na,0)
            CALL gatherv_r(vol        ,ncell_fluid,temp0(1,3),na,0)
            CALL gatherv_r(cell%tl    ,ncell_fluid,temp0(1,4),na,0)
!
            IF(myrank.eq.0) THEN
               level1=0.d0
               level2=0.d0
               DO i=1,na
                  level1 = level1 + temp0(i,1)*temp0(i,2)*temp0(i,3)
                  level2 = level2 + temp0(i,2)*temp0(i,3)
               ENDDO
               level = level1/level2*11.d0+0.0515d0 
               WRITE(333,27) time,level,temp0(274,4),temp0(176,4) !collapsed water level,274=(0.225,2.8),176=(0.225,1.41)
            ENDIF
            DEALLOCATE(temp0)  
         ENDIF         
   27    FORMAT(4x,4e15.7)    
!
!.........YHY
!
      ELSEIF(vv_prob.eq.'cavitation') THEN
         IF(initial) THEN
            IF(myrank.eq.0) THEN
!              OPEN(30, file='plot.out')
!              OPEN(50, file='Tec2d.dat')
               OPEN(333, file='VFT4_ref.dat')
            ENDIF
            initial=.false.
         ENDIF
         IF(time.gt.0.001d0*nplot) THEN
            nplot=nplot+1
         ENDIF
         IF(myrank.eq.0) THEN
            ALLOCATE(temp0(na,3))
         ELSE
            ALLOCATE(temp0(1,3))
         ENDIF
         CALL gatherv_r(cell%alphag,ncell_fluid,temp0(1,1),na,0)
         CALL gatherv_r(cell%tg    ,ncell_fluid,temp0(1,2),na,0)
         CALL gatherv_r(cell%tl    ,ncell_fluid,temp0(1,3),na,0)
!
         IF(myrank.eq.0) THEN
            WRITE(333,334) time,temp0(1582,1),temp0(1582,2),temp0(1582,3)   ! foamGrid : 1582 / somaGrid : 2000
         ENDIF
         DEALLOCATE(temp0)
!
      ELSEIF (vv_prob.eq.'plume') THEN
         IF(initial)THEN
            IF(myrank.eq.0)THEN
!              OPEN(30, file='plot.out')
!              OPEN(50, file='Tec2d.dat')
               OPEN(333, file='VFT5_ref.dat')
            ENDIF
            initial=.false.
         ENDIF
         IF(myrank.eq.0) THEN
            ALLOCATE(temp0(na,3))
         ELSE
            ALLOCATE(temp0(1,3))
         ENDIF
         CALL gatherv_r(cell%alphag,ncell_fluid,temp0(1,1),na,0)
         CALL gatherv_r(cell%tg    ,ncell_fluid,temp0(1,2),na,0)
         CALL gatherv_r(cell%tl    ,ncell_fluid,temp0(1,3),na,0)
!
         IF(myrank.eq.0) THEN
            WRITE(333,334) time, temp0(3141,1), temp0(3141,2), temp0(3141,3)  ! foamGrid : 3141 / somaGrid : 1620
         ENDIF
         DEALLOCATE(temp0)
         IF (time.gt.0.1*nplot) THEN
            nplot=nplot+1
         ENDIF
!
      ELSEIF(vv_prob.eq.'flashing') THEN
         IF(initial) THEN
            IF(myrank.eq.0) THEN
!               OPEN(111, file='sa.dat')
               OPEN(333, file='VFT3_ref.dat')
            ENDIF
            initial=.false.
         ENDIF
!
         IF(myrank.eq.0)THEN 
            ALLOCATE(temp0(na,4))
         ELSE
            ALLOCATE(temp0(1,4))
         ENDIF
         CALL gatherv_r(cell%alphag,ncell_fluid,temp0(1,1),na,0)
         CALL gatherv_r(p          ,ncell_fluid,temp0(1,2),na,0)
         CALL gatherv_r(cell%tg    ,ncell_fluid,temp0(1,3),na,0)
         CALL gatherv_r(cell%tl    ,ncell_fluid,temp0(1,4),na,0)
         IF(myrank.eq.0)  THEN
            WRITE(333,334) time,temp0(243,1),temp0(243,3),temp0(243,4)  ! foamGrid : 243 / somaGrid : 140
         ENDIF
         DEALLOCATE(temp0)
!
      ELSEIF (vv_prob.eq.'manometric') THEN
         IF (initial.and.myrank.eq.0) THEN
            OPEN(333, file='VFT7_ref.dat')
!           OPEN(30, file='plot.out')
            initial=.false.
         ENDIF
         IF(myrank.eq.0)THEN 
            ALLOCATE(temp0(na,ndim))
         ELSE
            ALLOCATE(temp0(1,ndim))
         ENDIF
         CALL gatherv_r_2d(vl_n,ncell_fp,temp0,ncell_fluid,na,0)
         IF(myrank.eq.0) THEN
            WRITE(333,334) time, temp0(519,1), temp0(519,2)  ! foamGrid : 281 or 519 / somaGrid : 677
         ENDIF
         DEALLOCATE(temp0)
         IF(time.gt.0.5*nplot) THEN
            nplot=nplot+1
         ENDIF
!
      ELSEIF(vv_prob.eq.'2D_loca') THEN
         IF (initial) THEN
            IF (myrank.eq.0) THEN
!              OPEN(30, file='plot.out')
               OPEN(20, file='VFT10_ref.dat')
            ENDIF
            time1=0
            initial=.false.
         ENDIF
         CALL sbloca_out_yhy_user
         time1 = time1 + dt
         IF (time1 .gt. 0.2d0) THEN
            time1 = 0.
            CALL sbloca_out_user
         ENDIF
!
      ELSEIF(vv_prob.eq.'fluidic_device') THEN
         CALL fd_post_user
!
      ELSEIF (vv_prob.eq.'T_blowdown') THEN
         CALL T_blowdown_user
!        
      ELSEIF(vv_prob.eq.'sgp_separator')THEN
         CALL udfn_sg_post
!
      ELSEIF (vv_prob.eq.'multi_ncg') THEN
         CALL multi_ncg_out
!
      ELSEIF (vv_prob.eq.'multi_ncg_diff') THEN
         CALL multi_ncg_diff_out
!
!.........LJR
!
      ELSEIF(vv_prob.eq.'block_porous') THEN 
         IF (initial) THEN
            IF(myrank.eq.0) THEN
               OPEN(111, file='VFS9_ref.dat')
            ENDIF
            initial=.false.
         ENDIF
         IF(myrank.eq.0) THEN
            ALLOCATE(temp1(na))
         ELSE
            ALLOCATE(temp1(1))
         ENDIF
         CALL gatherv_r(vl_n(1,2),ncell_fluid,temp1,na,0)
!
         IF(myrank.eq.0) THEN
            WRITE(111,400) time,temp1(58),temp1(98)
         ENDIF
!
      ELSEIF(vv_prob.eq.'annul_porous') THEN
         IF (initial) THEN
            IF(myrank.eq.0) THEN
               OPEN(111, file='VFS10_ref.dat')
            ENDIF
            initial=.false.
         ENDIF
         IF(myrank.eq.0) THEN
            ALLOCATE(temp1(na))
         ELSE
            ALLOCATE(temp1(1))
         ENDIF
         CALL gatherv_r(vl_n(1,2),ncell_fluid,temp1,na,0)
         IF(myrank.eq.0) THEN
         ! WRITE(111,400) time,temp1(445),temp1(134) !somaGrid.in
            WRITE(111,400) time,temp1(383),temp1(134) !foamGrid.in
         ENDIF
         DEALLOCATE(temp1)
!
      ELSEIF(vv_prob.eq.'stern') THEN
         IF(initial) THEN
            IF(myrank.eq.0)THEN
!              OPEN(30, file='plot.out')
               OPEN(333, file='VD9_ref.dat')
            ENDIF
            initial=.false.
         ENDIF
         IF(MOD(itim,iprn).eq.0)THEN
            IF(myrank.eq.0) THEN
               ALLOCATE(temp1(na),temp0(na,ndim))
            ELSE
               ALLOCATE(temp1(1),temp0(1,ndim))
            ENDIF
            CALL gatherv_r(cell%tl  ,ncell_fluid,temp1,na,0)
            CALL gatherv_r_2d(vl_n,ncell_fp,temp0,ncell_fluid,na,0)
            IF(myrank.eq.0) THEN
               WRITE(333,400) time,temp1(1000),temp0(1000,1),temp0(1000,2)
            ENDIF
            DEALLOCATE(temp1,temp0)
         ENDIF
  35     FORMAT(4(e20.10,1x))
!
      ELSEIF(vv_prob.eq.'kuhn_111')THEN
         IF(initial) THEN
            zzkuhn=(/2.418d0, 2.315d0, 2.181d0, 2.047d0, 1.897d0, 1.709d0, &
                    1.531d0, 1.313d0, 1.097d0, 0.837d0, 0.569d0, 0.237d0, 0.0d0/)
            initial=.false.
         ENDIF
!
         IF(time.gt.t_end) THEN
            ALLOCATE(temp0(12,3))
            temp0(:,1)=0.d0
            temp0(:,2)=0.d0
            temp0(:,3)=0.d0
            DO i=1,ncell_fluid
               IF(wVertical(i).eq.1)THEN
               DO j=1,12
                  IF(xloc(i,3).lt.zzkuhn(j).and.xloc(i,3).gt.zzkuhn(j+1))THEN
                     temp0(j,1)=-qconden(i)
                     temp0(j,2)=2.418d0-xloc(i,3)
                     temp0(j,3)=face%twall_partition(i)
                  ENDIF
               ENDDO
               ENDIF
            ENDDO
!
            IF(np.gt.1)THEN
               CALL allreducei_max_r(temp0(1,1),12)
               CALL allreducei_max_r(temp0(1,2),12)
               CALL allreducei_max_r(temp0(1,3),12)
            ENDIF
!
            IF(myrank.eq.0)THEN
               OPEN(440, file='VFT14_ref.dat')
               OPEN(441, file='VFT14_Twall.dat')
               WRITE(440,777) 'time =', time
               WRITE(441,777) 'time =', time
               DO i=1,12
                  WRITE(440,334) temp0(i,2),temp0(i,1)
                  WRITE(441,334) temp0(i,2),temp0(i,3)
               ENDDO
               CLOSE(440)
               CLOSE(441)
            ENDIF
            DEALLOCATE(temp0)
         ENDIF
777   FORMAT(a,1x,f5.2)  
!        
!.........PIK
!
!
      ELSEIF(vv_prob.eq.'halden650_5') then
!         
         CALL halden650_lbloca_output  
!         
      ELSEIF(vv_prob.eq.'Horizontal_flow')THEN
         IF (initial) THEN
            IF(myrank.eq.0) THEN
               OPEN(333, file='VFS15_ref.dat')            
               WRITE(112,*) 'time    vx(1)'
            ENDIF
            time1=0
            initial=.false.
         ENDIF
         time1 = time1 + dt
         IF (time1 .gt. 0.1d0) THEN
            time1 = 0.
            IF(myrank.eq.0) THEN
               WRITE(333,334)time,vl_n(1,1)
            ENDIF
         ENDIF     
!
      ELSEIF(vv_prob.eq.'icarus2002') THEN
!         
          CALL icarus2002_reflood_output
!
      ELSEIF(vv_prob.eq.'rocom' .or. vv_prob.eq.'rocom_mc') THEN
         IF(initial) THEN
            IF(myrank.eq.0)THEN
               IF(vv_prob.eq.'rocom')THEN
                  OPEN(333, file='VD5_ref.dat')
               ELSEIF(vv_prob.eq.'rocom_mc')THEN
                  OPEN(333, file='MCC4_ROCOM_ref.dat')
               ENDIF       
               WRITE(333,*)'Time CoreMin  DcinMin  DcoutMin  CoreAve   DcinAve  DcoutAve'
            ENDIF        
            stavg=173.0d0
            ftavg=183.0d0
            stavg=0.0d0
            ftavg=1.0d0
            print_time=0.0d0
            CALL_time=toutstep*0.001d0
            initial=.false.
         ENDIF
         prnopt=0
         IF(time.ge.CALL_time)THEN
            CALL_time=CALL_time+toutstep*0.01d0
            prnopt=0
            IF(time.ge.print_time.and.time.lt.(ftavg+(ftavg-stavg)*0.1d0))THEN
               print_time=print_time+(ftavg-stavg)*0.1d0
               prnopt=1
            ENDIF
            CALL rocom_output_cylinder_user(time,stavg,ftavg) 
            CALL rocom_output_plane_user(time,stavg,ftavg) 
            IF(myrank.eq.0) THEN
               WRITE(333,1000)time-t_end_ctrl(1),(tl_space_min(i),i=1,3),(tl_space_avg(i),i=1,3)       
               WRITE(unit_log,1001)time,(pbnd(i),i=1,4),(in_mdot(i),i=1,4),(out_mdot(i),i=1,4)
            ENDIF
   1001     FORMAT(1x,'time&pbnd&in&out=',5e13.6,8e13.3)            
   1000     FORMAT(1x,1e20.10,50e14.4)             
         ENDIF
!
      ELSEIF(vv_prob.eq.'vawl'        .or. &
             vv_prob.eq.'bankoff'     .or. &
             vv_prob.eq.'subo'        .or. &
             vv_prob.eq.'turbulent'   .or. &
             vv_prob.eq.'turbulent3d' .or. &
             vv_prob.eq.'boron_trans')THEN
        IF(initial) THEN
            print_time=toutstep
            CALL_time=toutstep*0.1
            initial=.false.
        ENDIF
        IF(time.ge.call_time)THEN
           call_time=call_time+toutstep*0.1
           prnopt=0
           IF(time.ge.print_time)THEN
               print_time=print_time+toutstep
               prnopt=1
           ENDIF
           CALL print_radial_oned_user(time)
           CALL print_axial_oned_user(time)
           IF(vv_prob.eq.'boron_trans')CALL print_cboron_user(time)
        ENDIF
!
      ELSEIF(vv_prob.eq.'siphon')THEN
          CALL siphon_print_ref       
!                     
!.........CYJ
!
      ELSEIF (vv_prob.eq.'moving_wall') THEN
         IF(initial)THEN
            IF(myrank.eq.0)THEN
               OPEN(112, file='VFS2_ref.dat')            
               WRITE(112,101) 'x t=10.0 t=20.0 t=30.0 t=40.0 t=50.0 t=60.0 t=70.0 t=80.0 t=90.0 t=100.'
            ENDIF
            n=0
            time1=0
            initial=.false.
         ENDIF
  100   FORMAT(100(e14.7,1x))
  101   FORMAT(a150)        
!               
!........Additional output for V&V report
!              
         time1 = time1 + dt
         IF (time1 .ge. 10.d0) THEN
            time1 = 0.d0
            n=n+1
            height=0.02d0                                ! x-location       
            fdim=2                                       ! dimension number for finding (y direction for this problem)
            sdim=1                                       ! dimension number for sorting (x direction for this problem)
!
            CALL Get_extracted_data(fdim,sdim,height,nloop,n,vl_o(1,fdim))
!
            IF(myrank.eq.0)THEN 
               DO i=1,nloop
                  vy(i,n)=prnvar(i,3)
               ENDDO
               IF(time.ge.100.d0)THEN
                  DO i=1,nloop
                     WRITE(112,100)  prnvar(i,1),(vy(i,j),j=1,n)
                  ENDDO
               ENDIF
            ENDIF
         ENDIF   
  5010   FORMAT(1x,4e16.8)   
!
      ELSEIF (vv_prob.eq.'mass_diffusion') THEN
         IF(initial)THEN
            IF(myrank.eq.0)THEN
               OPEN(112, file='VFS14_ref.dat')            
               WRITE(112,*) '    x       VF_NC.'
            ENDIF   
            initial=.false.
         ENDIF   
!               
!........Additional output for V&V report
!              
        IF (time .gt. t_end) THEN
            height=0.15d0                                ! x-location       
            fdim=2                                       ! dimension number for finding (y direction for this problem)
            sdim=1                                       ! dimension number for sorting (x direction for this problem)
!
            CALL Get_extracted_data(fdim,sdim,height,nloop,1,cell%quala)
!
            IF(myrank.eq.0)THEN 
               DO i=1,nloop
                  WRITE(112,100) prnvar(i,1),prnvar(i,3)
               ENDDO
               DEALLOCATE(prnvar)
            ENDIF
         ENDIF   
!
      ELSEIF(vv_prob.eq.'2D_boiling') THEN
         IF(initial) THEN
            IF (myrank.eq.0) THEN
!            OPEN(111, file='sa.dat')
               OPEN(333, file='VFT1_ref.dat')
               OPEN(112, file='VFT1_vv.dat')
            ENDIF
            initial=.false.
         ENDIF
!                  
         IF(myrank.eq.0)THEN 
            ALLOCATE(temp0(na,3))
         ELSE
            ALLOCATE(temp0(1,3))
         ENDIF
         CALL gatherv_r(cell%alphag,ncell_fluid,temp0(1,1),na,0)
         CALL gatherv_r(cell%tg    ,ncell_fluid,temp0(1,2),na,0)
         CALL gatherv_r(cell%tl    ,ncell_fluid,temp0(1,3),na,0)
!         IF(time.gt.t_end.and.myrank.eq.0) THEN
!            nnx=5
!            nny=50
!            DO i=1,nny
!               avg_a=0.d0
!               DO j=1,nnx
!                  avg_a=avg_a+temp0(nnx*(i-1)+j,1)
!               ENDDO
!               avg_a=avg_a/nnx
!               k=nnx*(i-1)+3
!               WRITE(111,300) xloc(k,2),avg_a
!            ENDDO
!         ENDIF
         IF(myrank.eq.0) THEN
            WRITE(333,334) time, temp0(98,1),temp0(98,2),temp0(98,3)    ! x=0.1m, y=1.02m for cell#=98
         ENDIF
         DEALLOCATE(temp0)        
!               
!........Additional output for V&V report
!              
         IF(time.gt.t_end)THEN
            height=0.1d0                                ! x-location        
            fdim=1
            sdim=2
!
            CALL Get_extracted_data(fdim,sdim,height,nloop,1,cell%alphag)
!
            IF(myrank.eq.0)THEN 
               DO i=1,nloop
                  WRITE(112,5010) time,prnvar(i,1),prnvar(i,2),prnvar(i,3)
               ENDDO
               DEALLOCATE(prnvar)
            ENDIF
         ENDIF          
!
      ELSEIF(vv_prob.eq.'3D_boiling') THEN
!         
      ELSEIF(vv_prob.eq.'2D_laminar'.and.time.gt.1.0*nplot)THEN
         IF(initial) THEN
            IF(myrank.eq.0)THEN
               OPEN(111, file='VFS3_ref.dat')
               OPEN(112, file='VFS3_vv.dat')
            ENDIF
            initial=.false.
         ENDIF
         nplot=nplot+1   
         fdim=2
         sdim=1
!               
         IF(myrank.eq.0) THEN
            ALLOCATE(temp0(na,ndim))
         ELSE
            ALLOCATE(temp0(1,ndim))
         ENDIF
         CALL gatherv_r_2d(vl_n,ncell_fp,temp0,ncell_fluid,ncell_fluid_all,0)
         IF(myrank.eq.0) THEN
            WRITE(111,400) time, temp0(401,1), temp0(401,2)
         ENDIF
         DEALLOCATE(temp0)
!               
!........Additional output for V&V report
!              
         IF(time.gt.t_end)THEN
            height=0.9d0                        
!
            CALL Get_extracted_data(fdim,sdim,height,nloop,1,vl_n(1,fdim))
!
            IF(myrank.eq.0) THEN
               DO i=1,nloop         
                  WRITE(112,5009) time,prnvar(i,1),prnvar(i,2),prnvar(i,3)
               ENDDO
               DEALLOCATE(prnvar)
            ENDIF !myrank              
         ENDIF                                        
  5009   FORMAT(1x,4e16.8)   
!        
      ELSEIF(vv_prob.eq.'PSBT_sngl')THEN
         IF (time.gt.0.01*nplot.and.time.le.5.0) THEN
            nplot=nplot+1
            CALL psbt_post
!         ELSEIF (time.gt.0.5*(nplot-32).and.time.gt.4.0) THEN
!            nplot=nplot+1
!            CALL psbt_post                     
         ENDIF                   
!        IF(time.ge.t_end_ctrl(nctrl))THEN
!            CALL heat_partition_sens
!        ENDIF    
!
      ELSEIF(vv_prob.eq.'h2p1_0'.or.vv_prob.eq.'h2p1_0x'.or. &
             vv_prob.eq.'h2p1_1'.or.vv_prob.eq.'h2p1_1x'.or. &
             vv_prob.eq.'h2p1_2'.or.vv_prob.eq.'h2p1_2x'.or. &
             vv_prob.eq.'h2p1_3'.or.vv_prob.eq.'h2p1_3x'.or. &
             vv_prob.eq.'h2p1_4'.or.vv_prob.eq.'h2p1_4x'.or. &
             vv_prob.eq.'VD_h2p1_0')THEN
         IF(grav(2).gt.grav(3))THEN          !height:z
            CALL udfn_H2P1_outputHeTg(1,3)   ! 20
            CALL udfn_H2P1_outputHeTg(2,3)   ! 14
            CALL udfn_H2P1_outputHeTg(3,3)   ! 26
            CALL udfn_H2P1_outputFOV(3)
            CALL udfn_H2P1_outputTE(3)
            CALL udfn_H2P1_outputEF(3)
         ELSEIF(grav(2).lt.grav(3))THEN      !height:y
            CALL udfn_H2P1_outputHeTg(1,2)   ! 20
            CALL udfn_H2P1_outputHeTg(2,2)   ! 14
            CALL udfn_H2P1_outputHeTg(3,2)   ! 26
            CALL udfn_H2P1_outputFOV(2)
            CALL udfn_H2P1_outputTE(2)
            CALL udfn_H2P1_outputEF(2)
         ENDIF
!        
      ELSEIF(vv_prob.eq.'ST2-CT-01'.or.vv_prob.eq.'ST2-CT-02'.or. &
             vv_prob.eq.'ST2-CT-03')THEN
         IF(initialtt)THEN
            initialtt=.false.
            ttnum=1
            DO i=1,9000
               ttrec(i)=(i-1)*1.0d0
            ENDDO
         ENDIF
!
         DO i=1,9000
            IF(time.ge.ttrec(i).and.ttnum.eq.i)THEN
               DO j=1,13
                  CALL udfn_CUBE_outputTg(j)
               ENDDO
               DO j=1,8
                  CALL udfn_CUBE_outputTs(j)
               ENDDO
               ttnum=ttnum+1
            ENDIF
         ENDDO
!        
!.........CHK
!
      ELSEIF(vv_prob.eq.'2D_conduction')THEN
         IF(MOD(itim,iprn).eq.0)THEN
            IF(myrank.eq.0) THEN
               OPEN(111, file='VFS7_ref.dat')
               OPEN(112, file='VFS7_analy.dat')
            ENDIF
            IF(initial) THEN
               initial = .false.
            ENDIF
            IF(myrank.eq.0)THEN 
               ALLOCATE(temp0(na,4))
            ELSE
               ALLOCATE(temp0(1,4))
            ENDIF
!...........temp3 solid is useless ????
!           ALLOCATE(temp3(na_c)
            CALL gatherv_r(vl_n(1,2) ,ncell_fluid,temp0(1,1),na,0)
            CALL gatherv_r(cell%tl   ,ncell_fluid,temp0(1,2),na,0)
            CALL gatherv_r(cell%rhol ,ncell_fluid,temp0(1,3),na,0)
            CALL gatherv_r(cell%cpl  ,ncell_fluid,temp0(1,4),na,0)
!           CALL gatherv_r(solid%tsol,ncell_cond,temp3,na_c,1)
            IF(myrank.eq.0) THEN
               DO i=1,ncell_fluid_all
                  IF(xloc_tmp(i,2).ge.0.95d0.and.xloc_tmp(i,2).lt.0.96d0) THEN
                     x1=xloc_tmp(i,1)-0.15d0
                     WRITE(111,400) x1,temp0(i,1),temp0(i,2)
                     vel_an=0.15d0*(1.0d0-(x1/0.05d0)**2)
                     temp_ave=300.d0+1.0d6*0.95d0/(1.0d-2*temp0(i,3)*temp0(i,4)) !Tin=300,y=0.95, 2*q=1.0d6, 
                     tw=327.0251202d0  !temp3(1048)
                     temp_an=tw-(2.57353d-1*(5.0d0-6.0d0*(x1/0.05d0)**2+(x1/0.05d0)**4))*(tw-temp_ave) !35/136=2.57353d-1, 0.5*W=0.05d0
                     WRITE(112,400) x1,vel_an,temp_an
                  ENDIF
               ENDDO                   
               CLOSE(111)
               CLOSE(112)
            ENDIF                  
            DEALLOCATE(temp0)
         ENDIF
!         
       ELSEIF(vv_prob.eq.'separation')THEN
         IF(MOD(itim,iprn).eq.0)THEN
            IF(initial) THEN
               IF(myrank.eq.0) THEN
                  OPEN(111, file='VFT6_ref.dat')
                  OPEN(112, file='VFT6_vv.dat')
                  IF(nn.le.2500)THEN
                     WRITE(112,5012) 
                  ELSE
                     WRITE(112,5013)  
                  ENDIF   
 5012             FORMAT(7x,'height         time=0.0        time=1.0        time=2.0        time=3.0        time=4.0        time=5.0        time=6.0        time=7.0        time=8.0        time=9.0        time=10.0')  
 5013             FORMAT(7x,'height         time=0.0        time=1.0        time=2.0        time=3.0        time=4.0        time=5.0')                                     
               ENDIF
               nplot=0 
               CALL_time=0.0d0
               initial=.false.
            ENDIF
            IF(myrank.eq.0)THEN
               ALLOCATE(temp1(na))
            ELSE
               ALLOCATE(temp1(1))
            ENDIF
            CALL gatherv_r(cell%alphag,ncell_fluid,temp1,na,0)
            IF(myrank.eq.0)THEN
               WRITE(111,400) time,temp1(250)
            ENDIF
            DEALLOCATE(temp1)
               
!               
!...........Additional output for V&V report
!              
            fdim=1
            sdim=2
            IF(time.ge.CALL_time)THEN                           
               CALL_time=CALL_time+toutstep
               nplot=nplot+1
               IF(nn.le.500)THEN
                  height=0.05d0 
               ELSE
                  height=0.5d0 
               ENDIF
!
               CALL Get_extracted_data(fdim,sdim,height,nloop,nplot,cell%alphag)
!
               IF(myrank.eq.0)THEN
                  DO i=1,nloop
                     vy(i,nplot)=prnvar(i,3)
                  ENDDO
               ENDIF
!
!              IF((nn.gt.500.and.nplot.eq.6).or.(nn.le.500.and.nplot.eq.12))THEN 
               IF((nn.gt.2500.and.nplot.eq.6).or.(nn.gt.500.and.nn.le.2500.and.nplot.eq.11).or.(nn.le.500.and.nplot.eq.11))THEN 
                  IF(myrank.eq.0)THEN
                     DO i=1,nloop         
                        WRITE(112,5011) prnvar(i,2),(vy(i,k),k=1,nplot)
                     !WRITE(112,5011) vy(i,2),vy(i,3),vy(i,4),vy(i,5),vy(i,6),vy(i,7),vy(i,8)
                     ENDDO
                  ENDIF
               ENDIF
            ENDIF
         ENDIF 
  5011   FORMAT(1x,13e16.8)   
!           
      ELSEIF(vv_prob.eq.'dam_break')THEN
         IF(MOD(itim,iprn).eq.0)THEN
            IF(initial) THEN
               IF(myrank.eq.0)THEN
                  OPEN(111, file='VFT8_ref.dat')
               ENDIF
               initial=.false.
            ENDIF
            IF(myrank.eq.0)THEN
               ALLOCATE(temp1(na))
            ELSE
               ALLOCATE(temp1(1))
            ENDIF
            CALL gatherv_r(cell%alphag,ncell_fluid,temp1,na,0)
            IF(myrank.eq.0)THEN
               WRITE(111,400) time,temp1( 6),temp1( 8),temp1(10),temp1(12), &
                                   temp1(14),temp1(16),temp1(18),temp1(20)
            ENDIF
            DEALLOCATE(temp1)
            CALL dambreaking_front_user(time)
         ENDIF  
!           
      ELSEIF(vv_prob.eq.'3D_air/water')THEN
         IF(MOD(itim,iprn).eq.0)THEN
            IF(initial) THEN
               IF(myrank.eq.0)THEN
                  OPEN(111, file='VFT9_tetra_ref.dat')
               ENDIF
               initial=.false.
            ENDIF
            IF(myrank.eq.0)THEN
               ALLOCATE(temp0(na,4))
            ELSE
               ALLOCATE(temp0(1,4))
            ENDIF
            CALL gatherv_r(cell%alphag,ncell_fluid,temp0(1,1),na,0)
            CALL gatherv_r(cell%quala ,ncell_fluid,temp0(1,2),na,0)
            CALL gatherv_r(vg_n(1,3)  ,ncell_fluid,temp0(1,3),na,0)
            CALL gatherv_r(vl_n(1,3)  ,ncell_fluid,temp0(1,4),na,0)
            IF(myrank.eq.0) THEN
               WRITE(111,400) time,temp0(2090,1),temp0(2090,2),temp0(2090,3),temp0(2090,4)
            ENDIF
            DEALLOCATE(temp0)
         ENDIF            
!
!.....KJT
!
      ELSEIF(vv_prob.eq.'rbla'.or.vv_prob.eq.'hmta_wf'.or.vv_prob.eq.'copain_porous')THEN
         IF(initial)THEN
            IF(vv_prob.eq.'rbla')THEN
               cellNum=(/32501, 34601, 36701, 38801, 40901, 43001, 45101, 47201, 49301, 51401, &
                        11501, 13601, 15701, 17801, 19901, 22001, 24101, 26201, 28301, 30401/)
            ELSEIF(vv_prob.eq.'hmta_wf')THEN
               cellNum=(/13001, 13841, 14681, 15521, 16361, 17201, 18041, 18881, 19721, 20561, &
                        4601, 5441, 6281, 7121, 7961, 8801, 9641, 10481, 11321, 12161/)
            ELSE
               cellNum=(/3367, 3379, 3391, 3403, 3415, 3427, 3439, 3451, 3463, 3475, &
                        3487, 3499, 3511, 3523, 3535, 3547, 3559, 3571, 3583, 3595/)
            ENDIF
            initial=.false.
         ENDIF
!
         IF(time.gt.t_end) THEN
            IF(myrank.eq.0)THEN
               ALLOCATE(temp0(na,2))
            ELSE
               ALLOCATE(temp0(1,2))
            ENDIF
            CALL gatherv_r(xloc(1,3),ncell_fluid,temp0(1,1),na,0)
            CALL gatherv_r(qconden  ,ncell_fluid,temp0(1,2),na,0)
            IF(myrank.eq.0)THEN
               temp0(:,1)=2.5d0-temp0(:,1)
               temp0(:,2)=-temp0(:,2)
!
               OPEN(440, file='VD16_ref.dat')
               DO i=1,20
                  ii=cellNum(i)
                  WRITE(440,334) temp0(ii,1), temp0(ii,2)
               ENDDO
               CLOSE(440)
            ENDIF
            DEALLOCATE(temp0)
         ENDIF
!        
!.....run-SUB
!
      ELSEIF(vv_prob.eq.'cnen') then       
!
         IF(time.gt.t_end) THEN
            IF(myrank.eq.0) THEN
               OPEN(333, file='cnen_ref.dat')
               OPEN(334, file='cnen_result.dat')
            ENDIF  
            IF(myrank.eq.0) THEN
               ALLOCATE(temp1(na))
            ELSE
               ALLOCATE(temp1(1))
            ENDIF
            CALL gatherv_r(vl_n(1,ndim),ncell_fluid,temp1,na,0)
!
            IF(myrank.eq.0) THEN
               WRITE(333,*) 'vl_n(corner),vl_n(side),vl_n(center)'
               WRITE(333,100) temp1(50),temp1(100),temp1(400)
               WRITE(334,100)  4.2103, temp1(50)
               WRITE(334,100)  4.8369, temp1(100)
               WRITE(334,100)  5.5846, temp1(400)
               CLOSE(333)
               CLOSE(334)
            ENDIF  
            DEALLOCATE(temp1)
         ENDIF  
!        
      ELSEIF(vv_prob.eq.'pnl_70pro') then
!
         IF(initial) THEN
            IF(myrank.eq.0) THEN
               OPEN(334, file='index_channel1.dat')
               DO i=1,25
                  READ(334,*) pnl_height(i)
               ENDDO
               CLOSE(334)
            ENDIF    
            initial=.false.
         ENDIF    
!        
         IF(myrank.eq.0) THEN
            ALLOCATE(temp1(na))
         ELSE
            ALLOCATE(temp1(1))
         ENDIF
         CALL gatherv_r(vl_n(1,3),ncell_fluid,temp1,na,0)
!
         IF(time.gt.t_end) THEN
            IF(myrank.eq.0) THEN
               ALLOCATE(temp2(na))
            ELSE
               ALLOCATE(temp2(1))
            ENDIF
            CALL gatherv_r(xloc(1,3),ncell_fluid,temp2,na,0)
            IF(myrank.eq.0) THEN
               OPEN(333, file='pnl70pro_result.dat')
               DO i=1,25
                  k=pnl_height(i)
!                 vel=SQRT(temp1(k)**2)
                  vel=ABS(temp1(k))
                  WRITE(333,100) temp2(k),vel/1.73736d0
               ENDDO   
               CLOSE(333)
            ENDIF
            DEALLOCATE(temp2)         
         ENDIF
!
         IF(myrank.eq.0) THEN
            IF(np.eq.1) THEN
               OPEN(431, file='VFTX_ref-compare_np1.dat')
            ELSE
               OPEN(431, file='VFTX_ref-compare_np4.dat')
            ENDIF
            DO i=1,25
               k=pnl_height(i) 
               WRITE(431,311) i,temp1(k)
            ENDDO               
            CLOSE(431)
         ENDIF
         DEALLOCATE(temp1)         
!             
      ELSEIF(vv_prob.eq.'ce15x15') then
!
         IF(time.gt.t_end) THEN
            IF(myrank.eq.0) THEN
               ALLOCATE(temp0(na,ndim),temp1(na))
            ELSE
               ALLOCATE(temp0(1,ndim),temp1(1))
            ENDIF
            CALL gatherv_r_2d(vl_n,ncell_fp,temp0,ncell_fluid,na,0)
            CALL gatherv_r(xloc(1,1),ncell_fluid,temp1,na,0)
!
            IF(myrank.eq.0) THEN         
               OPEN(334, file='index_ce.dat')
               DO i=1,8
                  READ(334,*) (ce_index(j,i),j=1,6)
               ENDDO
               CLOSE(334)
               OPEN(333, file='ce_velocity.dat')    
               OPEN(335, file='ce_velocity_result.dat')                
               DO i=1,8
                  DO j=1,6
                     k=ce_index(j,i)
                     vel_ce(j)=dsqrt(temp0(k,1)**2+temp0(k,2)**2+temp0(k,3)**2)
                  ENDDO
                  WRITE(333,100) temp1(k),(vel_ce(j),j=1,6)
                  WRITE(335,100) temp1(k),vel_ce(1),vel_ce(3),vel_ce(5)  !center (if side 2,4,6)
               ENDDO              
               CLOSE(333)
               CLOSE(335)
            ENDIF   
            DEALLOCATE(temp0,temp1)
         ENDIF   
!
      ELSEIF(vv_prob.eq.'rpi2x2') then       
!
         IF(time.gt.t_end) THEN
            IF(myrank.eq.0) THEN
               ALLOCATE(temp1(na))
            ELSE
               ALLOCATE(temp1(1))
            ENDIF
            CALL gatherv_r(cell%alphag,ncell_fluid,temp1,na,0)
!
            IF(myrank.eq.0) THEN
               OPEN(333, file='rpi_ref.dat')
               OPEN(334, file='rpi_result.dat')
               WRITE(333,*) 'ag_n(corner),ag_n(side),ag_n(center)'
               WRITE(333,100) temp1(45),temp1(90),temp1(225)
               WRITE(334,100) 0.143478261, temp1(45)
               WRITE(334,100) 0.208695652, temp1(90)
               WRITE(334,100) 0.263768116, temp1(225)            
               CLOSE(333)
               CLOSE(334)
            ENDIF
            DEALLOCATE(temp1)     
         ENDIF
!           
      ELSEIF(vv_prob.eq.'VdR0') then       
          
          IF(MOD(itim,1000).eq.0)then
          DO i=1,29
              ii1=30-i
              ii2=29+i
              write(3090,*)i,cell%alphag(ii1),cell%alphag(ii2)
          ENDDO    
          CLOSE(3090)
          ENDIF
          
      ELSEIF(vv_prob.eq.'VdR') then       
         IF(initial) THEN
            initial=.false.
            CALL_time=0.d0
         ENDIF
         
         IF(time.gt.CALL_time)THEN
            CALL_time=CALL_time+toutstep
            nplot=nplot+1
            
            ALLOCATE(temp1(na))
            IF(np.gt.1) THEN
               CALL gatherv_r(cell%alphag,ncell_fluid,temp1,na,0)
            ELSE
               temp1(:)=cell%alphag(:)
            ENDIF
            IF(myrank.eq.0) THEN
               OPEN(333, file='VD27_VdR_result.dat')
               ALLOCATE(temp4(na/2),temp5(na/2),temp6(na/2))
               DO i=1,na
                  cx=chn_nx(i) 
                  cz=chn_nz(i)
                  IF(cx.eq.1)then
                     temp4(cz)=xloc_tmp(i,3)
                     temp5(cz)=temp1(i)
                  ELSEIF(cx.eq.2)then   
                     temp6(cz)=temp1(i)
                  ENDIF   
               ENDDO
               DO i=1,na/2
                  WRITE(333,*)temp4(i),temp5(i),temp6(i)
               ENDDO   
               CLOSE(333)
               DEALLOCATE(temp4,temp5,temp6)
            ENDIF      
            DEALLOCATE(temp1)
         ENDIF
          
!          
      ELSEIF(vv_prob.eq.'wh14x28') then
!
         num=29
         IF(initial) THEN
            IF(myrank.eq.0) THEN
               OPEN(334, file='index_wh.dat')
               DO i=1,num
                  READ(334,*) (wh_index(j,i),j=1,7)
               ENDDO
               CLOSE(334)
            ENDIF   
            initial=.false.
         ENDIF   
!         
         IF(myrank.eq.0) THEN
            ALLOCATE(temp0(na,ndim))
         ELSE
            ALLOCATE(temp0(1,ndim))
         ENDIF
         CALL gatherv_r_2d(vl_n,ncell_fp,temp0,ncell_fluid,na,0)
!
         IF(time.gt.t_end) THEN
            IF(myrank.eq.0) THEN
               ALLOCATE(temp1(na))
            ELSE
               ALLOCATE(temp1(1))
            ENDIF
            CALL gatherv_r(xloc(1,1),ncell_fluid,temp1,na,0)
!
            IF(myrank.eq.0) THEN
               OPEN(333, file='wh_velocity.dat')    
               OPEN(335, file='wh_level7_result.dat')  
               DO i=1,num
                  DO j=1,7
                     k=wh_index(j,i)
                     vel_wh(j)=sqrt(temp0(k,1)**2+temp0(k,2)**2+temp0(k,3)**2)
                  ENDDO
                  WRITE(333,100) temp1(k),(vel_wh(j),j=1,7)
               ENDDO              
               DO i=1,num
                  j=7
                  k=wh_index(j,i)
                  vel_wh(j)=sqrt(temp0(k,1)**2+temp0(k,2)**2+temp0(k,3)**2)
                  WRITE(335,100) temp1(k),vel_wh(j)
               ENDDO    
               CLOSE(333)
               CLOSE(335)
            ENDIF
            DEALLOCATE(temp1)                
         ENDIF
!        
!        compare np1 vs np4         
         IF(myrank.eq.0) THEN
            IF(np.eq.1) THEN
               OPEN(431, file='VFTX_ref-compare_np1.dat')
            ELSE
               OPEN(431, file='VFTX_ref-compare_np4.dat')
            ENDIF   
            DO i=1,num
               j=7
               k=wh_index(j,i)
               WRITE(431,311) i, temp0(k,3)
            ENDDO                  
            CLOSE(431)
         ENDIF         
         DEALLOCATE(temp0)                
!
      ELSEIF(vv_prob.eq.'PNNL2x6') THEN
         IF(initial) THEN
            CALL_time=0.0d0
            initial=.false.
         ENDIF
         
         IF(time.ge.CALL_time)THEN                           
            ALLOCATE(temp1(na),temp2(na),temp3(7),temp4(7),temp5(7),temp6(7),temp7(7))
            CALL_time=CALL_time+toutstep
            CALL gatherv_r(cell%tl,ncell_fluid,temp1,na,0)
            DO i=1,ncell_fluid
               ul(i)=dsqrt(vl_n(i,1)*vl_n(i,1)+vl_n(i,2)*vl_n(i,2)+vl_n(i,3)*vl_n(i,3))
            ENDDO   
            CALL gatherv_r(ul,ncell_fluid,temp2,na,0)
            IF(myrank.eq.0)THEN
               OPEN(333, file='VD26_PNNL2x6_Temp.dat')
               OPEN(334, file='VD26_PNNL2x6_Velo.dat')
               DO i=1,na
                  cx=chn_nx(i) 
                  cy=chn_ny(i)
                  IF(xloc_tmp(i,3).gt.0.75 .and. xloc_tmp(i,3).lt.0.83)then
                     IF(cy.eq.2)then
                        temp3(cx)=xloc_tmp(i,1)
                        temp4(cx)=temp1(i)
                        temp5(cx)=temp2(i)
                     ENDIF
                  ELSEIF(xloc_tmp(i,3).gt.1.30 .and. xloc_tmp(i,3).lt.1.38)then
                     IF(cy.eq.2)then
                        temp6(cx)=temp1(i)
                        temp7(cx)=temp2(i)
                     ENDIF
                  ENDIF   
               ENDDO
               DO i=1,7
                  WRITE(333,*)temp3(i),temp4(i),temp6(i)
                  WRITE(334,*)temp3(i),temp5(i),temp7(i)
               ENDDO
               CLOSE(333)
               CLOSE(334)
            ENDIF   
            DEALLOCATE(temp1,temp2,temp3,temp4)
         ENDIF   
!
      ELSEIF(vv_prob.eq.'Tapucu') then       
         IF(initial) THEN
            initial=.false.
            CALL_time=0.d0
         ENDIF
         
         IF(time.gt.CALL_time)THEN
            CALL_time=CALL_time+toutstep
            nplot=nplot+1
            
            ALLOCATE(temp1(na),temp2(na),temp3(na))
            IF(np.gt.1) THEN
               CALL gatherv_r(cell%alphag,ncell_fluid,temp1,na,0)
               CALL gatherv_r(cell%tg    ,ncell_fluid,temp2,na,0)
               CALL gatherv_r(cell%tl    ,ncell_fluid,temp3,na,0)
            ELSE
               temp1(:)=cell%alphag(:)
               temp2(:)=cell%tg(:)
               temp3(:)=cell%tl(:)
            ENDIF
            
            
            IF(myrank.eq.0) THEN
               OPEN(333, file='VD26_Tapucu_result.dat')
               ALLOCATE(temp4(na/2),temp5(na/2),temp6(na/2))
               DO i=1,na
                  cx=chn_nx(i) 
                  cz=chn_nz(i)
                  IF(cx.eq.1)then
                     temp4(cz)=xloc_tmp(i,3)
                     temp5(cz)=temp1(i)
                  ELSEIF(cx.eq.2)then   
                     temp6(cz)=temp1(i)
                  ENDIF   
               ENDDO
               DO i=1,na/2
                  WRITE(333,*)temp4(i),temp5(i),temp6(i)
               ENDDO   
               CLOSE(333)
               DEALLOCATE(temp4,temp5,temp6)
            ENDIF      
            DEALLOCATE(temp1,temp2,temp3)
         ENDIF
!
!.....run-2nd
      ELSEIF (vv_prob.eq.'xn_2nd_conv') THEN
         CALL xn_2nd_conv_out
!
      ELSEIF (vv_prob.eq.'mom_2nd_conv') THEN
         CALL mom_2nd_conv_out
!
      ELSEIF (vv_prob.eq.'mult_ncg_2nd_conv') THEN
         CALL mult_ncg_2nd_conv_out
!         
      ELSEIF (vv_prob.eq.'eng_2nd_conv') THEN
         CALL eng_2nd_conv_out
!
!.....run-RV
!
      ELSEIF (vv_prob.eq.'fs_31203' .or.       &
              vv_prob.eq.'fs_31302' .or.       &            
              vv_prob.eq.'fs_31701' .or.       &            
              vv_prob.eq.'fs_31805'     ) THEN
         CALL fs_1D_out
!
      ELSEIF (vv_prob.eq.'fs_31203_3D' .or.      &
              vv_prob.eq.'fs_31302_3D' .or.      &            
              vv_prob.eq.'fs_31701_3D' .or.      &            
              vv_prob.eq.'fs_31805_3D'    ) THEN
         CALL fs_1D_out
!         
      ELSEIF (vv_prob.eq.'rbht1196_1d'.or.vv_prob.eq.'rbht1196_3d') THEN
         CALL rbht_1D_out
!
      ELSEIF(vv_prob.eq.'UPTF-RV') THEN           
         CALL uptf_3D_out         
!
      ELSEIF(vv_prob.eq.'apr1400_lbloca'      .or.   &
             vv_prob.eq.'opr1000_rv_lbloca'   .or.   &
             vv_prob.eq.'opr1000_mc_rv_lbloca' )THEN
!      
         IF(initial) THEN
            initial=.false.
            DO i=1,4
               mflux_sit(i)=0.d0
               mflux_sip(i)=0.d0
               mflux_dvi(i)=0.d0
            ENDDO
            IF(myrank.eq.0)THEN
               OPEN(333, file='VD15_apr1400_lbloca_ref.dat')
               OPEN(334, file='VD15_apr1400_lbloca_vv.dat')
            ENDIF   
            print_time=time-0.001d0
            IF(myrank.eq.0) WRITE(333,5000)(varname(i),i=1,25)
         ENDIF
         CALL apr1400_lbloca_out_user_sub            
         IF(time.ge.print_time)THEN
            IF(vv_prob.eq.'opr1000_mc_rv_lbloca')THEN
               pres_break=0.0d0
               IF(tplot_cell_loc(1).gt.0)THEN
                 pres_break=cell%p(tplot_cell_loc(1))
               ELSE
                 pres_break=cell%p(1)
               ENDIF
               IF(np.gt.1) CALL allreducei_max_r1(pres_break)
            ENDIF
            CALL apr1400_lbloca_out         
            CALL  apr1400_lbloca_out_user            
            pct=0.0d0
            IF(ncell_fluid_core.gt.0)pct=MAXVAL(twall_rv(:,1))
            CALL reduce_max_r1(pct,pct_tmp)
            pres_break_tmp=pres_break
            CALL allreduce_r(mflux_sit,mflux_sit_tmp,4)
            CALL allreduce_r(mflux_sip,mflux_sip_tmp,4)
            CALL allreduce_r(mflux_dvi,mflux_dvi_tmp,4)
            CALL reduce_r1(vl_choke,v1_tmp)
            CALL reduce_r1(vl_o_avg,v2_tmp)
!
! mod. by lsj
            IF(choke) THEN
               ichok=1
            ELSE
               ichok=0
            ENDIF             
            IF(np.gt.1) CALL allreducei_i1(ichok)
!            
            IF(myrank.eq.0)THEN
               WRITE(333,"(24e20.10,1i3)")time,pres_break_tmp,pct_tmp,break_flow,break_flow_int,si_flow,si_flow_int,       &
                                               mflux_sit_tmp(1:m_max-1),mflux_sip_tmp(1:m_max-1),mflux_dvi_tmp(1:m_max-1), &
                                               power_2d,v1_tmp,v2_tmp,mflux_sit_int,mflux_dvi_int,ichok
            ENDIF
            IF(time.gt.tbreak_s.and.pres_break_tmp.gt.pbnd_f*2.0d0)THEN
               print_time=print_time+0.05d0 
            ELSE
               print_time=print_Time+0.05d0
            ENDIF
         ENDIF     
 5000    FORMAT(1x,25(A20))            
!             
!.....run-MCC
!     
      ELSEIF (vv_prob.eq.'manosteam_mcc') THEN
!
         IF (initial.and.myrank.eq.0) THEN
            OPEN(333, file='MCC3_manometer_ref.dat')
            WRITE(333,*) 'time  x_velocity'
            initial=.false.
            print_time=0.0d0
         ENDIF
         IF(time.ge.print_time)THEN
            print_time=print_time+0.01d0
            IF(myrank.eq.0) WRITE(333,334) time,vl_n(212,1) !,c3vl(1,1),c3vl(1,2),c3vg(1,1),c3vg(1,2)
         ENDIF
!             
     ELSEIF (vv_prob.eq.'check_couple_mcc') THEN
         IF (initial.and.myrank.eq.0) THEN
            OPEN(333, file='MCC1_check_couple_ref.dat')
            initial=.false.
            print_time=0.0d0
            IF(myrank.eq.0) WRITE(333,"(a)")'Time VL1 VL2 VL3 VL4 VL5 VL6 VL7 VL8 VL9 VL10 P1 P2 P3 P4 P5 P6 P7 P8 P9 P10 VG1 VG2 VG3 VG4 VG5 VG6 VG7 VG8 VG9 VG10 ag1 ag2 ag3 ag4 ag5 ag6 ag7 ag8 ag9 ag10 c3pa c3vl c3vg'  
         ENDIF
         IF(time.ge.print_time)THEN
            print_time=print_time+toutstep*0.1
            IF(myrank.eq.0) WRITE(333,333) time,vl_n(1,3),vl_n(2,3),vl_n(3,3),vl_n(4,3),vl_n(5,3),vl_n(6,3),vl_n(7,3),vl_n(8,3),vl_n(9,3),vl_n(10,3) &
                                               ,cell%p(1),cell%p(2),cell%p(3),cell%p(4),cell%p(5)                                                    &
                                               ,cell%p(6),cell%p(7),cell%p(8),cell%p(9),cell%p(10)                                                   &
                                               ,vg_n(1,3),vg_n(2,3),vg_n(3,3),vg_n(4,3),vg_n(5,3),vg_n(6,3),vg_n(7,3),vg_n(8,3),vg_n(9,3),vg_n(10,3) &
                                               ,cell%alphag(1),cell%alphag(2),cell%alphag(3),cell%alphag(4),cell%alphag(5) &            
                                               ,cell%alphag(6),cell%alphag(7),cell%alphag(8),cell%alphag(9),cell%alphag(10) &
                                               ,c3pa(1,1),c3vl(1,1),c3vg(1,1)
         ENDIF
!
      ELSEIF (vv_prob.eq.'mass_check_mcc'.or.&
              vv_prob.eq.'single_channel_mer') THEN
         IF (initial.and.myrank.eq.0) THEN
            IF(vv_prob.eq.'mass_check_mcc')THEN
               OPEN(333, file='MCC2_mass_check_ref.dat')
            ELSE
               OPEN(333, file='NSSC_channel_ref.dat')
            ENDIF 
            initial=.false.
            print_time=time
            IF(myrank.eq.0) WRITE(333,"(a)")'Time          MF3          MF4          VL_Bot1      V_Bot2       VL_Bot3      VL_Top1      VL_Top2      VL_Top3'  
         ENDIF
!
         IF(time.ge.print_time)THEN
            print_time=print_time+toutstep*0.1
            IF(myrank.eq.0) WRITE(333,333) time,mass_nf_mcc(1)*-1.0d0,mass_nf_mcc(2), &
                                           vl_n(nn-1,2),vl_n(nn,2)
                                            !vl_n(1,2),vl_n(2,2),vl_n(3,2),&
                                            !vl_n(148,2),vl_n(149,2),vl_n(150,2)
         ENDIF 
!
      ELSEIF (vv_prob.eq.'single_channel'.or.&
              vv_prob.eq.'single_channel_dif'.or.&
              vv_prob.eq.'single_channel_dif_3D'.or.&
              vv_prob.eq.'multi_channel_2x1'.or.&
              vv_prob.eq.'multi_channel_2x0') THEN
         IF (initial.and.myrank.eq.0) THEN
            OPEN(333, file='NSSC_channel_ref.dat')
            initial=.false.
            print_time=time
            IF(myrank.eq.0)THEN
               IF(vv_prob.eq.'single_channel')THEN
                  WRITE(333,"(a)")'Time          vl1          vl3          vl5          vl7          vl9          p1           p3           p       5    p7           p9' 
               ELSE
                  WRITE(333,"(a)")'Time          VL_Bot1      V_Bot2       VL_Bot3      VL_Top1      VL_Top2      VL_Top3      p1           p2           p3           p4           p5' 
               ENDIF            
            ENDIF
         ENDIF
         IF(myrank.eq.0)THEN
            IF(time.ge.print_time)THEN
               print_time=print_time+toutstep*0.1
               IF(vv_prob.eq.'single_channel')THEN
                  WRITE(333,333) time,vl_n(5,3),vl_n(3,3),vl_n(1,3),vl_n(9,3),vl_n(7,3),&
                                      cell%p(5),cell%p(3),cell%p(1),cell%p(9),cell%p(7)
               ELSEIF(vv_prob.eq.'single_channel_dif')THEN
                  WRITE(333,333) time,vl_n(11,2),vl_n(10,2),vl_n(3,2),vl_n(7,2),vl_n(2,2),vl_n(15,2),&
                                      cell%p(10),cell%p(1),cell%p(14),cell%p(4),cell%p(2),&
                                      vl_n(17,2),vl_n(9,2)
               ELSEIF(vv_prob.eq.'multi_channel_2x1')THEN
                  WRITE(333,333) time,vl_n(2,2),vl_n(13,2),vl_n(15,2),vl_n(11,2),vl_n(9,2),vl_n(3,2),&
                                      cell%p(13),cell%p(17),cell%p(1),cell%p(14),cell%p(9),&
                                      vl_n(10,2),vl_n(5,2),vl_n(16,2)
               ELSEIF(vv_prob.eq.'multi_channel_2x0')THEN
                  WRITE(333,333) time,vl_n(7,2),vl_n(14,2),vl_n(13,2),vl_n(3,2),vl_n(4,2),vl_n(10,2),&
                                      cell%p(14),cell%p(16),cell%p(17),cell%p(6),cell%p(4),&
                                      vl_n(15,2),vl_n(5,2),vl_n(6,2),vl_n(3,1),vl_n(10,1)
               ENDIF   
            ENDIF 
         ENDIF
!
      ELSEIF (vv_prob.eq.'GE3x3') THEN 
         IF (initial) THEN
            IF(myrank.eq.0) THEN
               OPEN(333, file='VD28_GE3x3.dat')
               WRITE(333,"(a)")'     Time    exit_quality  qual_corner   qual_side   qual_center  exit_mflux  mflux_corner  mflux_side  mflux_center ' 
            ENDIF
            initial=.false.
            print_time=time
         ENDIF
!
         IF(time.ge.print_time)THEN
            print_time=print_time+toutstep*0.1
!            
            IF(myrank.eq.0)THEN
               ALLOCATE(tempi1(na))
               ALLOCATE(temp1(na),temp2(na),temp4(na),temp5(na),temp6(na),temp7(na),temp8(na),temp9(na))
            ELSE
               ALLOCATE(tempi1(na))
               ALLOCATE(temp1(1),temp2(1),temp3(1),temp4(1),temp5(1),temp6(1),temp7(1),temp8(1),temp9(1))
            ENDIF
            CALL gatherv_r(cell%quals,ncell_fluid,temp1,na,0)
            CALL gatherv_r(volp,ncell_fluid,temp2,na,0)   
            CALL gatherv_i(npb,ncell_fluid,tempi1,na,0)   
  
            CALL gatherv_r(cell%rhog,ncell_fluid,temp4,na,0)   
            CALL gatherv_r(cell%rhol,ncell_fluid,temp5,na,0)   
            CALL gatherv_r(cell%alphag,ncell_fluid,temp6,na,0) 
            CALL gatherv_r(cell%alphal,ncell_fluid,temp7,na,0)   
            CALL gatherv_r(vg_n(1,3),ncell_fluid,temp8,na,0)       
            CALL gatherv_r(vl_n(1,3),ncell_fluid,temp9,na,0)                           
!            
            IF(myrank.eq.0) THEN
               volsum=0.0d0
               qualsum=0.0d0
               mfluxg_sum=0.0d0
               mfluxl_sum=0.0d0
               DO i=1,ncell_fluid_all
                  IF(tempi1(i).eq.1)THEN                        !pboun cell only            
                     mflux_g=DABS(temp8(i)*temp4(i)*temp6(i))
                     mflux_l=DABS(temp9(i)*temp5(i)*temp7(i))   
                     volsum(4)=volsum(4)+ temp2(i)
                     qualsum(4)=qualsum(4)+ temp1(i)*temp2(i)
                     mfluxg_sum(4)=mfluxg_sum(4)+ mflux_g*temp2(i) 
                     mfluxl_sum(4)=mfluxl_sum(4)+ mflux_l*temp2(i)
                     IF(chn_type_tmp(i).eq.1)THEN       !center
                        volsum(1)=volsum(1)+ temp2(i)
                        qualsum(1)=qualsum(1)+ temp1(i)*temp2(i)    
                        mfluxg_sum(1)=mfluxg_sum(1)+ mflux_g*temp2(i) 
                        mfluxl_sum(1)=mfluxl_sum(1)+ mflux_l*temp2(i)
                     ELSEIF(chn_type_tmp(i).eq.2)THEN   !side
                        volsum(2)=volsum(2)+ temp2(i)
                        qualsum(2)=qualsum(2)+ temp1(i)*temp2(i)  
                        mfluxg_sum(2)=mfluxg_sum(2)+ mflux_g*temp2(i) 
                        mfluxl_sum(2)=mfluxl_sum(2)+ mflux_l*temp2(i)
                     ELSEIF(chn_type_tmp(i).eq.3)THEN   !corner
                        volsum(3)=volsum(3)+ temp2(i)
                        qualsum(3)=qualsum(3)+ temp1(i)*temp2(i)  
                        mfluxg_sum(3)=mfluxg_sum(3)+ mflux_g*temp2(i) 
                        mfluxl_sum(3)=mfluxl_sum(3)+ mflux_l*temp2(i)
                     ELSE
                        WRITE(*,*) 'subchannel type was not designated!'
                        STOP
                     ENDIF  
                  ENDIF
               ENDDO
               eq_av=    qualsum(4)/volsum(4)
               eq_center=qualsum(1)/volsum(1)
               eq_side=  qualsum(2)/volsum(2)
               eq_corner=qualsum(3)/volsum(3)
               mf_av=    (mfluxg_sum(4)+mfluxl_sum(4))/volsum(4)
               mf_center=(mfluxg_sum(1)+mfluxl_sum(1))/volsum(1)
               mf_side=  (mfluxg_sum(2)+mfluxl_sum(2))/volsum(2)
               mf_corner=(mfluxg_sum(3)+mfluxl_sum(3))/volsum(3)         
               WRITE(333,333) time,eq_av,eq_corner,eq_side,eq_center,mf_av,mf_corner,mf_side,mf_center
            ENDIF
            DEALLOCATE(temp1,temp2,temp4,temp5,temp6,temp7,temp8,temp9)     
            DEALLOCATE(tempi1)
         ENDIF           
!         
      ELSEIF(vv_prob.eq.'opr1000_rv')THEN
         IF(cupid_mars.eq.0)THEN
            IF (initial) THEN
               initial=.false.
               CALL Get_face_boundary
               IF(myrank.eq.0) THEN
                  OPEN(333, file='opr1000_ref.dat')
                  WRITE(333,*)'time   inlet1   inlet2   inlet3   inlet4   pboun1   pboun2'
               ENDIF
               print_time=0.0d0
               print_interval=0.1d0
            ENDIF
            IF(time.gt.print_time)THEN
               print_time=print_time+print_interval
               ALLOCATE(temp1(nb_max),temp2(nb_max))
               temp1(:)=0.0d0
               temp2(:)=0.0d0
               DO i=1,nb_max
                  ii1=index_flux(i)
                  ii=index_property(i)
                  IF(ii1.gt.0.and.ii.gt.0)THEN
                     IF(i.le.nin_max)THEN
                        temp1(i)= alphab_liq(i)*rhob_liq(i)*flux_l_nf(ii1) &
                                 +alphab_gas(i)*rhob_gas(i)*flux_g_nf(ii1) &
                                 +alphab_drp(i)*rhob_drp(i)*flux_d_nf(ii1)
                        temp2(i)=cell%p(ii)            
                     ELSE
                        temp1(i)= cell%alphal(ii)*cell%rhol(ii)*flux_l_nf(ii1) &
                                 +cell%alphag(ii)*cell%rhog(ii)*flux_g_nf(ii1) &
                                 +cell%alphad(ii)*cell%rhod(ii)*flux_d_nf(ii1)
                        temp2(i)=cell%p(ii) 
                     ENDIF               
                  ENDIF            
               ENDDO                      
               IF(np.gt.1) THEN
                  CALL allreducei_r(temp1,nb_max)
                  CALL allreducei_r(temp2,nb_max)
               ENDIF
               IF(myrank.eq.0) THEN
                  WRITE(333,334)time,(temp1(i),i=1,6),(temp2(i),i=1,6)
               ENDIF
               DEALLOCATE(temp1,temp2)
            ENDIF
         ENDIF
!
      ELSEIF(vv_prob.eq.'OPR1000_fullcore_modmesh02_rv_vessel') then
         IF(cupid_mars.eq.0)THEN
            IF(restart.ne.0)then
               IF(MOD(itim,10).eq.0)then
                  CALL udfn_subchannel_CHF
               ENDIF
            ENDIF
         ENDIF
         IF(cupid_mars.eq.1)THEN
            IF(restart.ne.0)then
               IF(MOD(itim,10).eq.0)then
                  CALL udfn_subchannel_CHF
               ENDIF
            ENDIF
         ENDIF         
!
!.........DSJ
!          
!--------------------------------------------------------------------------------------------
!
      ELSE
         !IF(itim.eq.1.and.myrank.eq.0) WRITE(*,*)'          No user defined outout assigned'     
!      
      ENDIF
! 
 334  FORMAT(100(e22.15,1x)) 
 333  FORMAT(100(e12.5,1x)) 
!
      END SUBROUTINE user_def_output
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
      SUBROUTINE prn_m_err(tot_mass)
!
!     This routine  writes time control variables in 'his.dat' file and on the screen
!
      USE Zcore           ,ONLY:np,myrank
      USE Zio_unit        ,ONLY:unit_screen,unit_log
      USE Zconst2         ,ONLY:iprn,dt
      USE Ztimecon        ,ONLY:dp_max,time,itim
!
      IMPLICIT NONE
!      
!.....Local variables
      REAL(8) :: tot_mass
!      REAL(8) :: yp_min,yp_max
!
!!.....Local arrays
!      REAL(8) :: yplus_min(np),yplus_max(np)
!      REAL(8) :: yplus_min1,yplus_max1
!
      IF(np.gt.1) CALL allreducei_r1(tot_mass)
!  
!      ip=myrank+1
!      IF(MOD(itim,iprn).eq.0)THEN
!         yplus_min(:)=0.0d0
!         yplus_max(:)=0.0d0
!         yplus_min(ip)=1.0d5
!         DO i=1,ncell_fluid
!            IF(icell_type(i).eq.1)THEN
!               yplus_min(ip)=MIN(yplus_min(ip),yplus(i))
!               yplus_max(ip)=MAX(yplus_max(ip),yplus(i))
!            ENDIF 
!         ENDDO 
!         IF(np.gt.1) THEN
!            CALL allreducei_r(yplus_min,np)
!            CALL allreducei_r(yplus_max,np)            
!         ENDIF 
!         yp_min=MINVAL(yplus_min(1:np))
!         yp_max=MAXVAL(yplus_max(1:np))
!!!!!!!!!!!!!!!!!!!!!!!!!
!         yplus_max1=0.0d0
!         yplus_min1=1.0d5
!         DO i=1,ncell_fluid
!            IF(icell_type(i).eq.1)THEN
!               yplus_min1=MIN(yplus_min1,yplus(i))
!               yplus_max1=MAX(yplus_max1,yplus(i))
!            ENDIF 
!         ENDDO 
!         IF(np.gt.1) CALL allreducei_min_r1(yplus_min1)
!         IF(np.gt.1) CALL allreducei_max_r1(yplus_max1)
!         if(yp_min.ne.yplus_min1) stop 10
!         if(yp_max.ne.yplus_max1) stop 11
!      ENDIF
!      
      IF(MOD(itim,iprn).eq.0.and.myrank.eq.0)THEN
         !WRITE(unit_log,20)itim,time,dp_max,tot_mass,dt,yp_min,yp_max
         WRITE(unit_log,20)   itim,time,dp_max,tot_mass,dt
         WRITE(unit_screen,11)itim,time,dp_max,tot_mass,dt
      ENDIF
!      
   11 FORMAT('+',5x,i10,2(1e12.4),2(1e12.4),2(f6.1)) !OPEN(6,CARRIAGECONTROL='FORTRAN')-->+,1,2
   20 FORMAT(i10,4(1pe15.7),1x,2(1pe12.4))
!
      END SUBROUTINE prn_m_err
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
      SUBROUTINE prn_iter
!
      USE Zcore,       ONLY:np,myrank
      USE Zio_unit,    ONLY:unit_screen,unit_log
      USE Zconst2,     ONLY:dt
      USE Zmass_conv,  ONLY:tot_mass
      USE Ztimecon,    ONLY:dp_max,time,itim,iter_p
!
      IMPLICIT NONE
!
      IF(np.gt.1) CALL allreducei_r1(tot_mass)
!
      IF(myrank.eq.0) THEN
         WRITE(unit_log,11)itim,'-',iter_p,time,dp_max,tot_mass,dt
         WRITE(unit_screen,11)itim,'-',iter_p,time,dp_max,tot_mass,dt
      ENDIF
   11 FORMAT('+',5x,i10,a1,i3,2(1pe13.5),2(1pe13.5))
!
      END SUBROUTINE prn_iter 
!-------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
      SUBROUTINE Get_face_boundary
      USE Zparam       , ONLY: nin_max
      USE Zbc_index    , ONLY: index_flux,index_property,npb
!
      USE Znum_cell    , ONLY: istart_nf,istart_nbcon_nf
      USE Zvec_index   , ONLY: left_nf,nbcon_nf,right_non
!
      IMPLICIT NONE
!      
      INTEGER :: i,k
      INTEGER :: ii,kk
      INTEGER :: nf_number,istart,isize,istart2,i1,i2
!
!.....Define inlet and outlet face
!
      nf_number=2
      istart=istart_nf(1,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      isize =istart_nf(2,nf_number)
      DO i=1,isize
         i1=istart+i
         ii=left_nf(i1)
         i2=istart2+i
         k=nbcon_nf(i2)
         index_flux(k)=i1
         index_property(k)=ii
      ENDDO
!
      nf_number=0
      istart=istart_nf(1,nf_number)
      isize   =istart_nf(2,nf_number)
      DO i=1,isize  
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i) 
         IF((npb(ii).eq.0.and.npb(kk).gt.0))THEN
           k=npb(ii)+npb(kk)+nin_max
           index_flux(k)=i1
           index_property(k)=ii
         ELSEIF((npb(ii).gt.0.and.npb(kk).eq.0))THEN
           k=npb(ii)+npb(kk)+nin_max
           index_flux(k)=i1
           index_property(k)=kk
        ENDIF
      ENDDO  
!      
      END SUBROUTINE Get_face_boundary
       
      SUBROUTINE Get_extracted_data(fdim,sdim,height,nloop,nplot,v_temp1)
!      
      USE Zcore           , ONLY: np,myrank
      USE Zzone           , ONLY: ncell_fluid,ncell_fluid_all
      USE Zcoord1         , ONLY: xloc
      USE Zcoord2         , ONLY: xloc_xfc_min,xloc_xfc_max,prnvar
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: fdim,sdim,nplot
      REAL(8) :: height
      REAL(8) :: v_temp1(ncell_fluid)
!.....Output
      INTEGER :: nloop
!.....Local variable
      INTEGER :: i,j,k
      INTEGER :: ip,i0,j0
      INTEGER :: na
      INTEGER :: loop
      INTEGER,SAVE :: nloop_l
      REAL(8) :: hzloc,lzloc
!.....Local arrays
      INTEGER :: nloop_all3(np),nloop_dsp3(np)
!.....Local allocatablearrays
      INTEGER,DIMENSION(:),SAVE,ALLOCATABLE :: indx
      INTEGER,DIMENSION(:),SAVE,ALLOCATABLE :: icell_index
      INTEGER,DIMENSION(:),SAVE,ALLOCATABLE :: ia,nloop_all,nloop_dsp
      REAL(8),DIMENSION(:),ALLOCATABLE :: prnvar0
      REAL(8),DIMENSION(:),ALLOCATABLE :: prnvar2_all
      REAL(8),DIMENSION(:),ALLOCATABLE :: prnvar3_l,prnvar3
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: prnvar2
!
      na=ncell_fluid_all
!
      IF(nplot.eq.1) THEN
         IF(.not.ALLOCATED(icell_index)) ALLOCATE(icell_index(ncell_fluid))
         ip=myrank+1
         nloop_l=0                 
         DO i=1,ncell_fluid
             hzloc=xloc_xfc_max(i,fdim)
             lzloc=xloc_xfc_min(i,fdim)
             IF(lzloc.le.height.and.hzloc.gt.height)THEN
                nloop_l=nloop_l+1
                icell_index(nloop_l)=i
             ENDIF              
         ENDDO     
         ALLOCATE(prnvar2(nloop_l,3))
         DO loop=1,nloop_l
            i=icell_index(loop)
            prnvar2(loop,1)=xloc(i,1)
            prnvar2(loop,2)=xloc(i,2)
            prnvar2(loop,3)=v_temp1(i)         
         ENDDO     
!
         IF(.not.ALLOCATED(ia)) ALLOCATE(ia(np+1))
         IF(.not.ALLOCATED(nloop_all)) ALLOCATE(nloop_all(np))
         IF(.not.ALLOCATED(nloop_dsp)) ALLOCATE(nloop_dsp(np))
         CALL allgather_i(nloop_l,nloop_all)
         ia(1)=1
         do ip=1,np
            ia(ip+1)=ia(ip)+nloop_all(ip)
            nloop_all3(ip)=3*nloop_all(ip)
         enddo
         nloop_dsp(1)=0
         nloop_dsp3(1)=0
         DO ip=2,np
            nloop_dsp(ip)=nloop_dsp(ip-1)+nloop_all(ip-1)
            nloop_dsp3(ip)=nloop_dsp3(ip-1)+nloop_all3(ip-1)
         ENDDO         
         nloop=ia(np+1)-1
!
         IF(myrank.eq.0)THEN 
            ALLOCATE(prnvar2_all(nloop*3))
         ELSE
            ALLOCATE(prnvar2_all(1))
         ENDIF
         CALL gather_vec_r(prnvar2,3*nloop_l,prnvar2_all,nloop*3,nloop_all3,nloop_dsp3)
!
         IF(myrank.eq.0)THEN 
            IF(ALLOCATED(prnvar)) DEALLOCATE(prnvar)
            ALLOCATE(prnvar(nloop,3))
            DO ip=1,np
               j0=nloop_dsp3(ip)
               i0=1
               DO i=ia(ip),ia(ip+1)-1
                  prnvar(i,1)=prnvar2_all(i0+j0)
                  prnvar(i,2)=prnvar2_all(i0+j0+  nloop_all(ip))
                  prnvar(i,3)=prnvar2_all(i0+j0+2*nloop_all(ip))
                  i0=i0+1
               ENDDO
            ENDDO
!
!...........Sort the extracted DATA      
!
            IF(ALLOCATED(indx)) DEALLOCATE(indx)
            ALLOCATE(indx(nloop))
            CALL sortx_r(prnvar(1,sdim),indx,nloop)
            ALLOCATE(prnvar0(nloop))
            DO k=1,3
               IF(k.eq.sdim) CYCLE
                  DO i=1,nloop
                     j=indx(i)
                     prnvar0(i)=prnvar(j,k)
                  ENDDO
                  DO i=1,nloop
                     prnvar(i,k)=prnvar0(i)
                  ENDDO
            ENDDO
            DEALLOCATE(prnvar0)
         ENDIF
      ELSE
         ALLOCATE(prnvar3_l(nloop_l))
         DO loop=1,nloop_l
            i=icell_index(loop)
            prnvar3_l(loop)=v_temp1(i)         
         ENDDO     
         IF(myrank.eq.0)THEN 
            ALLOCATE(prnvar3(nloop))
         ELSE
            ALLOCATE(prnvar3(1))
         ENDIF
         CALL gather_vec_r(prnvar3_l,nloop_l,prnvar3,nloop,nloop_all,nloop_dsp)
         IF(myrank.eq.0)THEN 
            DO ip=1,np
               j0=nloop_dsp(ip)
               i0=1
               DO i=ia(ip),ia(ip+1)-1
                  prnvar(i,3)=prnvar3(i0+j0)
                  i0=i0+1
               ENDDO
            ENDDO
            ALLOCATE(prnvar0(nloop))
            DO i=1,nloop
               j=indx(i)
               prnvar0(i)=prnvar(j,3)
            ENDDO
            DO i=1,nloop
               prnvar(i,3)=prnvar0(i)
            ENDDO
            DEALLOCATE(prnvar0)
         ENDIF
         DEALLOCATE(prnvar3_l,prnvar3)
      ENDIF
!
      END SUBROUTINE Get_extracted_data
