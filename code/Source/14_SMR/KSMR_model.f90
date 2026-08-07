!   
   subroutine ReadZone()
      use Zzone,        only: ncell_fluid, ncell_fluid_all
      use Zmpi,         only: jperm
      use KSMR
      implicit none
      INTEGER :: fUnit = 1293, size, ix, i, k
      LOGICAL, ALLOCATABLE :: tmp_all(:)
      INTEGER :: zone_ix(n_zone)
      CHARACTER(LEN=20) :: filename(n_zone)

      filename(1:n_zone) = ['core_in.dat', 'core.dat','core_out.dat', 'cea.dat', 'crdm.dat', 'upper_cavity.dat', 'mcp.dat', 'pzr.dat', 'discharge.dat', 'sg.dat', 'dc.dat', 'lp.dat']
      zone_ix(1:n_zone) = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]

      if(ALLOCATED(zone_comp)) deallocate(zone_comp)
      ALLOCATE(zone_comp(ncell_fluid)); zone_comp = 0
      ALLOCATE(tmp_all(ncell_fluid_all),zone_comp_all(ncell_fluid_all))

      do k = 1, n_zone
         tmp_all = .false.
         OPEN(fUnit, file=trim(filename(k)))
         READ(fUnit, *) size
         do i = 1, size
            READ(fUnit, *) ix
            tmp_all(ix) = .true.
         end do
         CLOSE(fUnit)
   
         do i = 1, ncell_fluid
            if(tmp_all(jperm(i)) == .true.) then
               zone_comp(i) = zone_ix(k)
            end if
         end do
      end do
      
      do k = 1, n_zone
         OPEN(fUnit, file=trim(filename(k)))
         READ(fUnit, *) size
         do i = 1, size
            READ(fUnit, *) ix
            zone_comp_all(ix) = k
         end do
         CLOSE(fUnit)
      end do      

      zone_comp_area = 0.d0
      zone_comp_height = 0.d0

      OPEN(fUnit, file='zone_area.dat')
      do i = 1, n_zone
         READ(fUnit, *), zone_comp_area(i)
      end do
      CLOSE(fUnit)
      
      
      OPEN(fUnit, file='zone_height.dat')
      do i = 1, n_zone
         READ(fUnit, *), zone_comp_height(i)
      end do
      CLOSE(fUnit)
   end subroutine ReadZone
!============================================================================================
   subroutine AssignCoreHeat()
      use Zzone,        only: ncell_fluid
      use Zcoord3,      only: volp
      use Zqvol,        only: qvol_liq
      use Zcore,        only: np
      use Ztimecon,     only: time
      use KSMR

      implicit none
      real(8) :: core_vol !, core_vol_tot
      INTEGER :: i

      core_vol = 0.d0
      do i = 1, ncell_fluid
         if(zone_comp(i) == 2) then   !core
            core_vol = core_vol+volp(i)
         end if
      end do
      IF(np.gt.1) call allreducei_r1(core_vol)
!
!      decay heat: ANS73      
      IF(time.lt.t_trip) THEN
         power_ans73=core_power
      ELSE
         power_ans73=1.2*core_power*(58116.01+9769.69*(time-t_trip))**(-1./4.0108) 
      ENDIF           
!       
      do i = 1, ncell_fluid
         if(zone_comp(i) == 2) then
            qvol_liq(i) = power_ans73 / core_vol
            qvol_liq(i) = qvol_liq(i) * min(max(0.d0,time/time_power_ramping),1.d0)
         end if
      end do
   end subroutine AssignCoreHeat

!============================================================================================   
   subroutine AssignSGHeat()
      use Zzone,        only: ncell_fluid
      use Zcoord3,      only: volp
      use Zqvol,        only: qvol_liq
      use Zcore,        only: np
      use Ztimecon,     only: time
      use KSMR

      implicit none
      real(8) :: sg_vol !, sg_vol_tot
      INTEGER :: i
!
      sg_vol = 0.d0
      do i = 1, ncell_fluid
         if(zone_comp(i) == 10) then !sg
            sg_vol = sg_vol+volp(i)
         end if
      end do
      if(np.gt.1) call allreducei_r1(sg_vol)
!
!      decay heat: ANS73      
      IF(time.lt.t_trip) THEN
         power_ans73=core_power
      ELSE
         power_ans73=1.2*core_power*(58116.01+9769.69*(time-t_trip))**(-1./4.0108) 
      ENDIF      
!          
      do i = 1, ncell_fluid
         if(zone_comp(i) == 10) then !sg
            qvol_liq(i) = -power_ans73 / sg_vol
            qvol_liq(i) = qvol_liq(i) * min(max(0.d0,time/time_power_ramping),1.d0)
         end if
      end do
   end subroutine AssignSGHeat

!============================================================================================
   subroutine InitializeTemperature()
      use Zzone,            only: ncell_fluid
      use Vol_DATA,         only: cell
      USE STM_TBL_cupid,    ONLY: st_tbl,             & 
                                 nt,np,ns,ns2,ndxstd
      use KSMR
      implicit none
      LOGICAL erx   
      REAL betafs,betags,cpfs,cpgs,entfs,entgs,hsubfs,hsubgs,         &
           kapafs,kapags,psats,s(36),tsat,usubgs,usubfs,vsubgs,vsubfs
      EQUIVALENCE(s( 1),tsat),    &
                 (s(10),psats),   &
                 (s(11),vsubfs),  &
                 (s(12),vsubgs),  &
                 (s(13),usubfs),  &
                 (s(14),usubgs),  &
                 (s(15),hsubfs),  &
                 (s(16),hsubgs),  &
                 (s(17),betafs),  &
                 (s(18),betags),  &
                 (s(19),kapafs),  &
                 (s(20),kapags),  &
                 (s(21),cpfs),    &
                 (s(22),cpgs),    &
                 (s(25),entfs),   &
                 (s(26),entgs)
      real(8) :: temperature_cold = T_in
      real(8) :: temperature_hot  = T_out
      INTEGER :: i, it

      do i = 1, ncell_fluid
         if(zone_comp(i) >= 4 .and. zone_comp(i) <= 8) then
            cell%tl(i) = temperature_hot
            cell%tl_o(i) = temperature_hot
         else if(zone_comp(i) >= 9 .and. zone_comp(i) <= 12) then
            cell%tl(i) = temperature_cold
            cell%tl_o(i) = temperature_cold
         else if(zone_comp(i) >= 1 .and. zone_comp(i) <= 3) then
            cell%tl(i) = temperature_cold
            cell%tl_o(i) = temperature_cold
         end if
      end do

      DO i=1,ncell_fluid
         s(:)=0.d0
         s(1)=cell%tl(i)
         s(2)=cell%p(i)
         CALL sth2x3_cupid(s,it,erx,                          &
                           st_tbl(ndxstd),                    &
                           st_tbl(ndxstd+nt),                 &
                           st_tbl(ndxstd+nt+np+13*ns+13*ns2))
         
         cell%rhol(i)=1.0/vsubfs
         cell%el(i)=usubfs
      ENDDO
   
   end subroutine InitializeTemperature

!============================================================================================
   subroutine ForceCoreInletTemperatrue()
      use Zzone,            only: ncell_fluid
      use Vol_DATA,         only: cell
      USE STM_TBL_cupid,    ONLY: st_tbl,             & 
                                 nt,np,ns,ns2,ndxstd
      use KSMR
      implicit none
      LOGICAL erx   
      REAL betafs,betags,cpfs,cpgs,entfs,entgs,hsubfs,hsubgs,         &
           kapafs,kapags,psats,s(36),tsat,usubfs,usubgs,vsubfs,vsubgs
      EQUIVALENCE(s( 1),tsat),    &
                 (s(10),psats),   &
                 (s(11),vsubfs),  &
                 (s(12),vsubgs),  &
                 (s(13),usubfs),  &
                 (s(14),usubgs),  &
                 (s(15),hsubfs),  &
                 (s(16),hsubgs),  &
                 (s(17),betafs),  &
                 (s(18),betags),  &
                 (s(19),kapafs),  &
                 (s(20),kapags),  &
                 (s(21),cpfs),    &
                 (s(22),cpgs),    &
                 (s(25),entfs),   &
                 (s(26),entgs)
      INTEGER :: i, it
      real(8) :: temperature_cold = T_in

      do i = 1, ncell_fluid
         if(zone_comp(i) == 1) then           !core-in
            cell%tl(i) = temperature_cold
            cell%tl_o(i) = temperature_cold
            s(:)=0.d0
            s(1)=cell%tl(i)
            s(2)=cell%p(i)
            CALL sth2x3_cupid(s,it,erx,st_tbl(ndxstd),st_tbl(ndxstd+nt),st_tbl(ndxstd+nt+np+13*ns+13*ns2))
            cell%rhol(i)=1.0/vsubfs
            cell%el(i)=usubfs
         end if
      end do
   end subroutine ForceCoreInletTemperatrue

!============================================================================================
   subroutine PrintRcsFlowRate()
         use Zzone,         only: ncell_fluid
         use Vol_DATA,      only: cell
         use Zvec_geo,      only: sa_nf
         use Znum_cell,     only: istart_nf
         use Zvec_index,    only: left_nf, right_non
         use Zvec_major,    only: flux_l_nf
         use Zcore,         only: myrank, np
         use Ztimecon,      only: time
         use KSMR
         
         implicit none
         
!.....External function
         INTEGER :: get_global_cell         
         INTEGER :: i,ii,kk,ii1,kk1,i1
         INTEGER :: size, nf_number, istart1
!         
         REAL(8) :: area
         REAL(8) :: flow_rate1, flow_rate2, flow_rate3, flow_rate4, flow_rate5
         LOGICAL, save :: initial = .true.
         INTEGER :: funit = 123123

         if(initial) then
            initial = .false.
            if(myrank == 0) OPEN(funit, file='RCS_flow.dat')
         end if
!         
         flow_rate1 = 0.d0
         flow_rate2 = 0.d0
         flow_rate3 = 0.d0
         flow_rate4 = 0.d0 
         flow_rate5 = 0.d0 
!
         nf_number = 0
         istart1 = istart_nf(1,nf_number)
         size = istart_nf(2,nf_number)
         do i = 1, size
            i1 = istart1+i    !...nf_tot
            area = sa_nf(i1)
            ii = left_nf(i1)  
            kk = right_non(i) 

            ii1=get_global_cell(ii)
            kk1=get_global_cell(kk)  
!
! core_in (kg/s)
            if(kk <= ncell_fluid) then
               if(zone_comp_all(ii1) == 12 .or. zone_comp_all(ii1) == 1) then
                  if((zone_comp_all(ii1) == 12 .and. zone_comp_all(kk1) == 1) .or. (zone_comp_all(ii1) == 1 .and. zone_comp_all(kk1) == 12)) then
                     if(zone_comp_all(ii1) == 12) then
                        flow_rate1 = flow_rate1 + cell%rhol(ii)*flux_l_nf(i1)
                     else
                        flow_rate1 = flow_rate1 - cell%rhol(kk)*flux_l_nf(i1)
                     end if
                  end if
               end if
            else
               if(zone_comp_all(ii1) == 12) then
                  if((zone_comp_all(ii1) == 12 .and. zone_comp_all(kk1) == 1) .or. (zone_comp_all(ii1) == 1 .and. zone_comp_all(kk1) == 12)) then
                     if(zone_comp_all(ii1) == 12) then
                        flow_rate1 = flow_rate1 + cell%rhol(ii)*flux_l_nf(i1)
                     else
                        flow_rate1 = flow_rate1 - cell%rhol(kk)*flux_l_nf(i1)
                     end if
                  end if
               end if
            end if  
!
! sg (kg/s)
            if(kk <= ncell_fluid) then
               if(zone_comp_all(ii1) == 10 .or. zone_comp_all(ii1) == 11) then
                  if((zone_comp_all(ii1) == 10 .and. zone_comp_all(kk1) == 11) .or. (zone_comp_all(ii1) == 11 .and. zone_comp_all(kk1) == 10)) then
                     if(zone_comp_all(ii1) == 10) then
                        flow_rate2 = flow_rate2 + cell%rhol(ii)*flux_l_nf(i1)
                     else
                        flow_rate2 = flow_rate2 - cell%rhol(kk)*flux_l_nf(i1)
                     end if
                  end if
               end if
            else
               if(zone_comp_all(ii1) == 10) then
                  if((zone_comp_all(ii1) == 10 .and. zone_comp_all(kk1) == 11) .or. (zone_comp_all(ii1) == 11 .and. zone_comp_all(kk1) == 10)) then
                     if(zone_comp_all(ii1) == 10) then
                        flow_rate2 = flow_rate2 + cell%rhol(ii)*flux_l_nf(i1)
                     else
                        flow_rate2 = flow_rate2 - cell%rhol(kk)*flux_l_nf(i1)
                     end if
                  end if
               end if
            end if     
! mcp_out
            if(kk <= ncell_fluid) then
               if(zone_comp_all(ii1) == 8 .or. zone_comp_all(ii1) == 7) then
                  if((zone_comp_all(ii1) == 8 .and. zone_comp_all(kk1) == 7) .or. (zone_comp_all(ii1) == 7 .and. zone_comp_all(kk1) == 8)) then
                     if(zone_comp_all(ii1) == 8) then
                        flow_rate3 = flow_rate3 + cell%rhol(ii)*flux_l_nf(i1)
                     else
                        flow_rate3 = flow_rate3 - cell%rhol(kk)*flux_l_nf(i1)
                     end if
                  end if
               end if
            else
               if(zone_comp_all(ii1) == 8) then
                  if((zone_comp_all(ii1) == 8 .and. zone_comp_all(kk1) == 7) .or. (zone_comp_all(ii1) == 7 .and. zone_comp_all(kk1) == 8)) then
                     if(zone_comp_all(ii1) == 8) then
                        flow_rate3 = flow_rate3 + cell%rhol(ii)*flux_l_nf(i1)
                     else
                        flow_rate3 = flow_rate3 - cell%rhol(kk)*flux_l_nf(i1)
                     end if
                  end if
               end if
            end if
!
         end do         
!
         IF(np.gt.1) THEN
             CALL allreducei_r1(flow_rate1)
             CALL allreducei_r1(flow_rate2)
             CALL allreducei_r1(flow_rate3)
             CALL allreducei_r1(flow_rate4)
             CALL allreducei_r1(flow_rate5)
         ENDIF    
!         
         if(myrank == 0) write(funit,105) time, ABS(flow_rate1), ABS(flow_rate2), ABS(flow_rate3), ABS(flow_rate4), ABS(flow_rate5)
105      format(50(f15.7,1x))
         
    end subroutine PrintRcsFlowRate

!============================================================================================
      subroutine PrintPressure()
         use Zzone,         only:ncell_fluid
         use Vol_DATA,      only: cell
         use Znum_cell,     only: istart_nf
         use Zvec_index,    only: left_nf, right_non
         use Zcore,         only: myrank, np
         use Ztimecon,      only: time
         USE Zpdrop,        ONLY: npid,num_dp_region,dp_region
         use KSMR
         
         implicit none
         
!.....External function
         INTEGER :: get_global_cell         
         INTEGER :: i,ii,kk,ii1,kk1,i1,dd
         INTEGER :: size, nf_number, istart1
!         
         INTEGER :: region(4)
         REAL(8) :: p1(npid),p2(npid),icnt1(npid),icnt2(npid)
!         
         LOGICAL, save :: initial = .true.
         INTEGER :: funit = 123125
!         
         if(initial) then
            initial = .false.
            if(myrank == 0) OPEN(funit, file='del_pressure.dat')
            ALLOCATE(dp_region(npid))
         end if         
!
         p1(:)=0.d0
         p2(:)=0.d0
         icnt1(:)=0.d0
         icnt2(:)=0.d0
!
         nf_number = 0
         istart1 = istart_nf(1,nf_number)
         size = istart_nf(2,nf_number)
         do i = 1, size
            i1 = istart1+i    !...nf_tot
            ii = left_nf(i1)  
            kk = right_non(i) 

            ii1=get_global_cell(ii)
            kk1=get_global_cell(kk)  
            
            DO dd=1,npid
               region(:)=num_dp_region(dd,:) !from input file(rv_parameters.in)
               if(kk <= ncell_fluid) then
                  if(zone_comp_all(ii1) == region(1) .or. zone_comp_all(ii1) == region(2)) then
                     if((zone_comp_all(ii1) == region(1) .and. zone_comp_all(kk1) == region(2)) .or. (zone_comp_all(ii1) == region(2) .and. zone_comp_all(kk1) == region(1))) then
                        if(zone_comp_all(ii1) == region(2)) then
                           p1(dd) = p1(dd) + cell%p(ii)
                           icnt1(dd)=icnt1(dd)+1.d0
                        elseif(zone_comp_all(kk1) == region(2)) then
                           p1(dd) = p1(dd) + cell%p(kk)
                           icnt1(dd)=icnt1(dd)+1.d0
                        end if
                     end if
                  end if
               else
                  if(zone_comp_all(ii1) == region(2)) then
                     if((zone_comp_all(ii1) == region(1) .and. zone_comp_all(kk1) == region(2)) .or. (zone_comp_all(ii1) == region(2) .and. zone_comp_all(kk1) == region(1))) then
                        if(zone_comp_all(ii1) == region(2)) then
                           p1(dd) = p1(dd) + cell%p(ii)
                           icnt1(dd)=icnt1(dd)+1.d0
                        elseif(zone_comp_all(kk1) == region(2)) then
                           p1(dd) = p1(dd) + cell%p(kk)
                           icnt1(dd)=icnt1(dd)+1.d0
                        end if
                     end if
                  end if
               end if  
               
               if(kk <= ncell_fluid) then
                  if(zone_comp_all(ii1) == region(3) .or. zone_comp_all(ii1) == region(4)) then
                     if((zone_comp_all(ii1) == region(3) .and. zone_comp_all(kk1) == region(4)) .or. (zone_comp_all(ii1) == region(4) .and. zone_comp_all(kk1) == region(3))) then
                        if(zone_comp_all(ii1) == region(3)) then
                           p2(dd) = p2(dd) + cell%p(ii)
                           icnt2(dd)=icnt2(dd)+1.d0
                        elseif(zone_comp_all(kk1) == region(3)) then
                           p2(dd) = p2(dd) + cell%p(kk)
                           icnt2(dd)=icnt2(dd)+1.d0
                        end if
                     end if
                  end if
               else
                  if(zone_comp_all(ii1) == region(3)) then
                     if((zone_comp_all(ii1) == region(3) .and. zone_comp_all(kk1) == region(4)) .or. (zone_comp_all(ii1) == region(4) .and. zone_comp_all(kk1) == region(3))) then
                        if(zone_comp_all(ii1) == region(3)) then
                           p2(dd) = p2(dd) + cell%p(ii)
                           icnt2(dd)=icnt2(dd)+1.d0
                        elseif(zone_comp_all(kk1) == region(3)) then
                           p2(dd) = p2(dd) + cell%p(kk)
                           icnt2(dd)=icnt2(dd)+1.d0
                        end if
                     end if
                  end if
               end if                
            ENDDO
         end do
!
         IF(np.eq.1) THEN
            p1(:)=p1(:)/icnt1(:)
            p2(:)=p2(:)/icnt2(:)
            dp_region(:)=p1(:)-p2(:)
         ELSE
            CALL barrier_mpi
            CALL allreducei_r(p1,npid)
            CALL allreducei_r(p2,npid)
            CALL allreducei_r(icnt1,npid)
            CALL allreducei_r(icnt2,npid)
            p1(:)=p1(:)/icnt1(:)
            p2(:)=p2(:)/icnt2(:)
            dp_region(:)=p1(:)-p2(:)
         ENDIF    
!
         if(myrank == 0) write(funit,105) time,(dp_region(dd),dd=1,npid)
105      format(50(f15.7,1x))
         
      end subroutine PrintPressure

!============================================================================================      
      subroutine PrintTemperatureCoreOut()
         use Zzone,         only:ncell_fluid
         use Vol_DATA,      only: cell
         use Zcore,         only: myrank, np
         use Ztimecon,      only: time
         use KSMR
         
         implicit none
         INTEGER :: i,cnt,cnt_in !, cnt_tot
         REAL(8) :: T,T_inlet !, T_tot
         REAL(8), SAVE :: pre_value
         REAL(8) :: current_value,current_value_in
         LOGICAL, save :: initial = .true.
         INTEGER :: funit = 123124

         if(initial) then
            initial = .false.
            if(myrank == 0) OPEN(funit, file='T_core_out.dat')
            pre_value = 540.d0
         end if

         T      = 0.d0
         T_inlet=0.d0
         cnt    = 0
         cnt_in =0
         do i = 1, ncell_fluid
            if(zone_comp(i) == 3) then  !core-out
               T = T + cell%tl(i)
               cnt = cnt+1
            elseif(zone_comp(i) == 1) then  !core-in
                T_inlet=T_inlet + cell%tl(i)
                cnt_in = cnt_in+1
            end if
         end do

         !if(np > 1) then
         !   call allreduce_r(T, T_tot, 1)
         !   call allreduce_i(cnt, cnt_tot, 1)
         !else
         !   T_tot = T
         !   cnt_tot = cnt
         !end if
         IF(np.gt.1) THEN
             CALL allreducei_r1(T)
             CALL allreducei_i1(cnt)
             CALL allreducei_r1(T_inlet)
             CALL allreducei_i1(cnt_in)
         ENDIF
 
!         current_value = T_tot/dble(cnt_tot)
         current_value = T/dble(cnt)
         current_value_in=T_inlet/dble(cnt_in)
!        Print out: time, T_core out, (T_core-out - T_core-in)
         if(myrank == 0) write(funit, '(F15.7, x, F15.7, x, F15.7)') time, current_value, current_value-current_value_in
         pre_value = current_value
!         
    end subroutine PrintTemperatureCoreOut

!============================================================================================
      subroutine PressureDrop()
         use Zzone            ,ONLY: ncell_fluid
         use Vol_DATA         ,ONLY: cell
         use Zvector          ,ONLY: ul_o
         USE STM_TBL_cupid    ,ONLY: st_tbl,nt,np,ns,ns2,ndxstd
         USE Zmodel           ,ONLY: resist
         USE Zvector          ,ONLY: ul_o
         USE Zconst1          ,ONLY: vv_prob
         USE Zconst2          ,ONLY: dt
         USE Zpdrop           ,ONLY: npid,time_pid_on,time_pid_off,err_o,Icon,dp_set,dp_control,err,Pcon,Dcon,dp_region
         use Ztimecon         ,only: time    
         USE Zporous          ,ONLY: l_subchannel,chn_type   !PSH
         use KSMR
!         
         implicit none
         INTEGER :: i, zone_ix, it
!         INTEGER, PARAMETER :: npid=3 !control 3 regions
         REAL(8), PARAMETER :: rcs_flow = 4000.d0 !3719.d0 !...Target RCS mass flow rate
         REAL(8), PARAMETER :: time_steady = 0.d0 !...Target RCS mass flow rate
         REAL(8) :: vel, dpdz
         REAL(8) :: rho_hot, rho_cold, rho_ave
!pid control         
         REAL(8) :: Tu,Ti,Td,Kp,Ku
!         
         INTEGER, PARAMETER :: fUnit = 1293
         LOGICAL, save :: initial = .true., initial_dp=.true., initial_control=.true.
         LOGICAL erx   
         REAL betafs,betags,cpfs,cpgs,entfs,entgs,hsubfs,hsubgs,         &
            kapafs,kapags,psats,s(36),tsat,usubfs,usubgs,vsubfs,vsubgs
         EQUIVALENCE(s( 1),tsat),    &
                  (s(10),psats),   &
                  (s(11),vsubfs),  &
                  (s(12),vsubgs),  &
                  (s(13),usubfs),  &
                  (s(14),usubgs),  &
                  (s(15),hsubfs),  &
                  (s(16),hsubgs),  &
                  (s(17),betafs),  &
                  (s(18),betags),  &
                  (s(19),kapafs),  &
                  (s(20),kapags),  &
                  (s(21),cpfs),    &
                  (s(22),cpgs),    &
                  (s(25),entfs),   &
                  (s(26),entgs)

         if(initial) then
!            initial = .false.
             
            if(ALLOCATED(resist)) deallocate(resist)
            allocate(resist(ncell_fluid)); resist = 0.d0;
!
            s(:) = 0.d0  
            s(1) = T_in
            s(2) = 15.d6
            CALL sth2x3_cupid(s,it,erx,st_tbl(ndxstd),st_tbl(ndxstd+nt),st_tbl(ndxstd+nt+np+13*ns+13*ns2))
            rho_cold = 1.0/vsubfs
            s(1) = T_out
            s(2) = 15.d6
            CALL sth2x3_cupid(s,it,erx,st_tbl(ndxstd),st_tbl(ndxstd+nt),st_tbl(ndxstd+nt+np+13*ns+13*ns2))
            rho_hot = 1.0/vsubfs
            rho_ave = 0.5d0*(rho_cold + rho_hot)
!
!           pid control            
            IF(initial_control) THEN
               ALLOCATE(dp_control(npid),dp_set(npid))
               dp_control(:)=0.d0
               initial_control=.false.                
            ENDIF   
            IF(time.gt.time_pid_on.and.time.lt.time_pid_off) THEN
               IF(initial_dp) THEN
                  ALLOCATE(err(npid),Pcon(npid),Dcon(npid))
                  ALLOCATE(err_o(npid),Icon(npid))
                  Icon(:)=0.d0
!                  
                  err_o(:)=dp_region(:)-dp_set(:)
                  !err_o(1)=sg_dp-dp_set(1)
                  !err_o(2)=core_dp-dp_set(2)
                  !err_o(3)=dc_dp-dp_set(3)
                  initial_dp=.false.
               ENDIF   
                
               Ku=0.1d0
               Tu=0.1d0
                    
               Kp=0.6*Ku
               Ti=0.5*Tu
               Td=0.125*Tu

               err(:)=dp_region(:)-dp_set(:)
               !err(1)=sg_dp-dp_set(1)
               !err(2)=core_dp-dp_set(2)
               !err(3)=dc_dp-dp_set(3)
               Pcon(:)=Kp*err(:)
               Icon(:)=Icon(:)+Kp*err(:)*dt/Ti
               Dcon(:)=Kp*Td*(err(:)-err_o(:))/dt
               dp_control(:)=-(Pcon(:)+Icon(:)+Dcon(:))
               err_o(:)=err(:)
               write(*,5003) err,Pcon,Icon,Dcon,dp_control,time
5003           FORMAT(100(e15.5,1x))                    
            ENDIF
!
            do i = 1, ncell_fluid
               IF(l_subchannel)THEN     
                  IF(chn_type(i).ne.0)CYCLE
               ENDIF
               zone_ix = zone_comp(i)
               select case (zone_ix)
!old-dsj                   
               case(1, 2, 3)  !...Core active / in / out
                 dp_set(2)=69.75d3  
                 dpdz = (dp_set(2)+dp_control(2)) / (zone_comp_height(1) + zone_comp_height(2) + zone_comp_height(3)) !69.75d3 / (zone_comp_height(1) + zone_comp_height(2) + zone_comp_height(3))                   
                 vel = rcs_flow / rho_ave / zone_comp_area(2)
!KSMR
               !case(1, 10, 11)  !...Core active / in / out   !KSMR when porosity in core is implemented
               !  dpdz = 0.d0
               !  vel = rcs_flow / rho_ave / zone_comp_area(1)
!                 
               case(4, 5, 6) !...CEA / CRDM / UPPER_CAVITY 
                 dpdz = 29.2d3 / (zone_comp_height(4) + zone_comp_height(5) + zone_comp_height(6))
                 vel = rcs_flow / rho_hot / zone_comp_area(4)
               case(7) !...MCP
                 dpdz =0.d0 !77.1d3 / (zone_comp_height(7))
                 vel = rcs_flow / rho_cold / zone_comp_area(10) !vel = rcs_flow / rho_ave / zone_comp_area(7)
               case(8, 9) !...PZR, Discharger
                 dpdz = 0.d0 !77.1d3 / (zone_comp_height(7) + zone_comp_height(8) + zone_comp_height(9))
                 vel = rcs_flow / rho_cold / zone_comp_area(10) !vel = rcs_flow / rho_ave / zone_comp_area(8)                 
               case(10) !...SG
                 dp_set(1)=107.8d3 !0.d0 !107.8d3
                 dpdz =(dp_set(1)+dp_control(1)) / (zone_comp_height(10))
!                
                 IF(vv_prob.eq.'KSMR-SG-porous') dpdz=0.d0
!                 
                 vel = rcs_flow / rho_cold / zone_comp_area(10)
               case(11) !...DC
                 dpdz = 13.2d3 / (zone_comp_height(11))   
                 vel = rcs_flow / rho_cold / zone_comp_area(11)
               case(12) !...LP
                 dpdz = 20.0d3 / (zone_comp_height(12))
                 vel = rcs_flow / rho_cold / zone_comp_area(12)
               case default
                  cycle
               end select
   
               resist(i) = dpdz/vel/vel
               resist(i) = resist(i) !* 0.1d0               
            end do
         end if
!
         do i = 1, ncell_fluid
            IF(l_subchannel)THEN
               IF(chn_type(i).ne.0)CYCLE
            ENDIF
            cell%vfwl(i) = resist(i)*ul_o(i)
         end do
!         
      end subroutine PressureDrop
