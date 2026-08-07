   subroutine ReadZone_NuScale()
      use Zzone, only: ncell_fluid, ncell_fluid_all
      use Zmpi, only: jperm
      use NuScale
         
      implicit none
!      
      INTEGER :: fUnit = 1293, size, ix, i, k 
      LOGICAL, ALLOCATABLE :: tmp_all(:)
      INTEGER :: zone_ix(n_zone)
      CHARACTER(LEN=20) :: filename(n_zone)

      filename(1:n_zone) = ['core_in.dat','core_active.dat','core_out.dat', 'core_bypass.dat','riser.dat', 'up.dat', 'pre_sg.dat', 'sg.dat', 'downcomer.dat', 'lp.dat']
      zone_ix(1:n_zone) = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

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
!!
!      nf_number = 0
!      istart1 = istart_nf(1,nf_number)
!      size = istart_nf(2,nf_number)
!      do i = 1, size
!         i1 = istart1+i  
!         area = sa_nf(i1)
!         ii = left_nf(i1)  
!         kk = right_non(i1)
!!         
!         ig=get_global_cell(ii)
!         kg=get_global_cell(kk)
!         
!         IF((zone_comp_all(ig).eq.1.and.zone_comp_all(kg).eq.3).or.(zone_comp_all(kg).eq.1.and.zone_comp_all(ig).eq.3)) THEN
!             zone_comp_area(1)=zone_comp_area(1)+area !core_active
!         ENDIF
!!         
!         IF((zone_comp_all(ig).eq.2.and.zone_comp_all(kg).eq.3).or.(zone_comp_all(kg).eq.2.and.zone_comp_all(ig).eq.3)) THEN
!             zone_comp_area(2)=zone_comp_area(2)+area !core_bypass
!         ENDIF         
!!         
!         IF((zone_comp_all(ig).eq.5.and.zone_comp_all(kg).eq.6).or.(zone_comp_all(kg).eq.5.and.zone_comp_all(ig).eq.6)) THEN
!             zone_comp_area(5)=zone_comp_area(5)+area !pre_sg
!         ENDIF             
!!         
!         IF((zone_comp_all(ig).eq.6.and.zone_comp_all(kg).eq.7).or.(zone_comp_all(kg).eq.6.and.zone_comp_all(ig).eq.7)) THEN
!             zone_comp_area(6)=zone_comp_area(6)+area !sg
!         ENDIF  
!!         
!         IF((zone_comp_all(ig).eq.7.and.zone_comp_all(kg).eq.8).or.(zone_comp_all(kg).eq.7.and.zone_comp_all(ig).eq.8)) THEN
!             zone_comp_area(7)=zone_comp_area(7)+area !downcomer
!         ENDIF           
!         
!      ENDDO   
!      zone_comp_area(4)=zone_comp_area(7)+zone_comp_area(1)+zone_comp_area(2)                 !up=downcomer+core_active+core_bypass
!      zone_comp_area(8)=zone_comp_area(4)                                                 !lp=up
!      zone_comp_area(3)=zone_comp_area(4)-zone_comp_area(5)                                 !riser
!
!      do i = 1, n_zone
!         print*,i,zone_comp_area(i)
!      end do  
      
      zone_comp_area(1)=1.710325d0  !core_in
      zone_comp_area(2)=1.710325d0  !core_active
      zone_comp_area(3)=1.710325d0  !core_out
      zone_comp_area(4)=1.028660d0  !core_bypass
      zone_comp_area(5)=1.626232d0  !riser
      zone_comp_area(6)=5.865608d0  !up
      zone_comp_area(7)=4.239448d0  !pre_sg
      zone_comp_area(8)=4.239448d0  !sg
      zone_comp_area(9)=3.126694d0  !downcomer
      zone_comp_area(10)=5.865680d0 !lp
      
      zone_comp_height(1)=0.2d0     !core_in
      zone_comp_height(2)=1.6d0     !core_active
      zone_comp_height(3)=0.2d0     !core_out
      zone_comp_height(4)=2.d0      !core_bypass  
      zone_comp_height(5)=11.3d0    !riser
      zone_comp_height(6)=1.66d0    !up 
      zone_comp_height(7)=0.5d0     !pre_sg
      zone_comp_height(8)=6.6d0     !sg
      zone_comp_height(9)=6.d0      !downcomer 
      zone_comp_height(10)=0.5d0    !lp
      
      !OPEN(fUnit, file='zone_area.dat')
      !do i = 1, n_zone
      !   READ(fUnit, *), zone_comp_area(i)
      !end do
      !CLOSE(fUnit)
!      
      !OPEN(fUnit, file='zone_height.dat')
      !do i = 1, n_zone
      !   READ(fUnit, *), zone_comp_height(i)
      !end do
      !CLOSE(fUnit)
   end subroutine ReadZone_NuScale

   subroutine AssignCoreHeat_NuScale()
      use Zzone, only: ncell_fluid
      use Zcoord3, only: volp
      use Zqvol, only: qvol_liq
      use Zcore, only: np
      use Ztimecon, only:time
      use NuScale

      implicit none
      real(8) :: core_vol !, core_vol_tot
      INTEGER :: i

      core_vol = 0.d0

      do i = 1, ncell_fluid
         if(zone_comp(i) == 2) then   !core_active
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
      IF(time.gt.time_power_off) qvol_liq(:)=0.d0      
   end subroutine AssignCoreHeat_NuScale

   
   subroutine AssignSGHeat_NuScale()
      use Zzone, only: ncell_fluid
      use Zcoord3, only: volp
      use Zqvol, only: qvol_liq
      use Zcore, only: np
      use Ztimecon, only: time
      use NuScale

      implicit none
      real(8) :: sg_power
      real(8) :: sg_vol !, sg_vol_tot
      INTEGER :: i

      sg_power = -core_power
      sg_vol = 0.d0

      do i = 1, ncell_fluid
         if(zone_comp(i) == 8) then !sg
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
         if(zone_comp(i) == 8) then !sg
            qvol_liq(i) = -power_ans73 / sg_vol
            qvol_liq(i) = qvol_liq(i) * min(max(0.d0,time/time_power_ramping),1.d0)
         end if
      end do
      IF(time.gt.time_power_off) qvol_liq(:)=0.d0      
   end subroutine AssignSGHeat_NuScale

   subroutine InitializeTemperature_NuScale()
      use Zzone, only: ncell_fluid
      use Vol_DATA, only: cell
      USE STM_TBL_cupid  , ONLY: st_tbl,             & 
                                 nt,np,ns,ns2,ndxstd
      use NuScale
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
         if(zone_comp(i) >= 5 .and. zone_comp(i) <= 7) then     
            cell%tl(i) = temperature_hot
            cell%tl_o(i) = temperature_hot
         else if(zone_comp(i) >= 8 .and. zone_comp(i) <= 10) then
            cell%tl(i) = temperature_cold
            cell%tl_o(i) = temperature_cold
         else if(zone_comp(i) >= 1 .and. zone_comp(i) <= 4) then
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
   
   end subroutine InitializeTemperature_NuScale


   subroutine ForceCoreInletTemperatrue_NuScale()
      use Zzone, only: ncell_fluid
      use Vol_DATA, only: cell
      USE STM_TBL_cupid  , ONLY: st_tbl,             & 
                                 nt,np,ns,ns2,ndxstd
      use NuScale
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
      end subroutine ForceCoreInletTemperatrue_NuScale

      subroutine PrintRcsFlowRate_NuScale()
         use Zzone, only:ncell_fluid,ncell_fluid_all
         use Vol_DATA, only: cell
         use Zvec_geo, only: sa_nf
         use Znum_cell, only: istart_nf
         use Zvec_index, only: left_nf, right_non
         use Zvec_major, only: flux_l_nf
         use Zcore, only: myrank, np
         use Zcoord1, only: xloc
         use Ztimecon, only: time
         use NuScale
         
         implicit none
         
!.....External function
         INTEGER :: get_global_cell       
         INTEGER :: i,ii,kk,ii1,kk1,i1
         INTEGER :: size, nf_number, istart1
         REAL(8) :: area, flow_rate1, flow_rate2, flow_rate3, flow_rate4, flow_rate5, rho_tmp(ncell_fluid_all)
         LOGICAL, save :: initial = .true.
         INTEGER :: funit = 123123

         if(initial) then
            initial = .false.
            if(myrank == 0) OPEN(funit, file='RCS_flow.dat')
         end if
!         
         rho_tmp=0.d0
         IF(np.gt.1) THEN
            CALL gatherv_r(cell%rhol,ncell_fluid,rho_tmp,ncell_fluid_all,0) !0=fluid 
            CALL broadcast_r(rho_tmp,ncell_fluid_all)
         ELSE
            rho_tmp=cell%rhol
         ENDIF           
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
! core_in            
            if(kk <= ncell_fluid) then
               if(zone_comp_all(ii1) == 1 .or. zone_comp_all(ii1) == 2) then
                  if((zone_comp_all(ii1) == 1 .and. zone_comp_all(kk1) == 2) .or. (zone_comp_all(ii1) == 2 .and. zone_comp_all(kk1) == 1)) then
                     if(zone_comp_all(ii1) == 1) then
                        flow_rate1 = flow_rate1 + cell%rhol(ii)*flux_l_nf(i1)
                     else
                        flow_rate1 = flow_rate1 - cell%rhol(kk)*flux_l_nf(i1)
                     end if
                  end if
               end if
            else
               if(zone_comp_all(ii1) == 1) then
                  if((zone_comp_all(ii1) == 1 .and. zone_comp_all(kk1) == 2) .or. (zone_comp_all(ii1) == 2 .and. zone_comp_all(kk1) == 1)) then
                     if(zone_comp_all(ii1) == 1) then
                        flow_rate1 = flow_rate1 + cell%rhol(ii)*flux_l_nf(i1)
                     else
                        flow_rate1 = flow_rate1 - cell%rhol(kk)*flux_l_nf(i1)
                     end if
                  end if
               end if
            end if  
!
! sg
            if(kk <= ncell_fluid) then
               if(zone_comp_all(ii1) == 8 .or. zone_comp_all(ii1) == 9) then
                  if((zone_comp_all(ii1) == 8 .and. zone_comp_all(kk1) == 9) .or. (zone_comp_all(ii1) == 9 .and. zone_comp_all(kk1) == 8)) then
                     if(zone_comp_all(ii1) == 8) then
                        flow_rate2 = flow_rate2 + cell%rhol(ii)*flux_l_nf(i1)
                     else
                        flow_rate2 = flow_rate2 - cell%rhol(kk)*flux_l_nf(i1)
                     end if
                  end if
               end if
            else
               if(zone_comp_all(ii1) == 8) then
                  if((zone_comp_all(ii1) == 8 .and. zone_comp_all(kk1) == 9) .or. (zone_comp_all(ii1) == 9 .and. zone_comp_all(kk1) == 8)) then
                     if(zone_comp_all(ii1) == 8) then
                        flow_rate2 = flow_rate2 + cell%rhol(ii)*flux_l_nf(i1)
                     else
                        flow_rate2 = flow_rate2 - cell%rhol(kk)*flux_l_nf(i1)
                     end if
                  end if
               end if
            end if              
!
! pre-sg
            if(kk <= ncell_fluid) then
               if(zone_comp_all(ii1) == 7 .or. zone_comp_all(ii1) == 8) then
                  if((zone_comp_all(ii1) == 7 .and. zone_comp_all(kk1) == 8) .or. (zone_comp_all(ii1) == 8 .and. zone_comp_all(kk1) == 7)) then
                     if(zone_comp_all(ii1) == 7) then
                        flow_rate3 = flow_rate3 + cell%rhol(ii)*flux_l_nf(i1)
                     else
                        flow_rate3 = flow_rate3 - cell%rhol(kk)*flux_l_nf(i1)
                     end if
                  end if
               end if
            else
               if(zone_comp_all(ii1) == 7) then
                  if((zone_comp_all(ii1) == 7 .and. zone_comp_all(kk1) == 8) .or. (zone_comp_all(ii1) == 8 .and. zone_comp_all(kk1) == 7)) then
                     if(zone_comp_all(ii1) == 7) then
                        flow_rate3 = flow_rate3 + cell%rhol(ii)*flux_l_nf(i1)
                     else
                        flow_rate3 = flow_rate3 - cell%rhol(kk)*flux_l_nf(i1)
                     end if
                  end if
               end if
            end if  
!
! downcomer
            if(kk <= ncell_fluid) then
               if(zone_comp_all(ii1) == 9 .or. zone_comp_all(ii1) == 10) then
                  if((zone_comp_all(ii1) == 9 .and. zone_comp_all(kk1) == 10) .or. (zone_comp_all(ii1) == 10 .and. zone_comp_all(kk1) == 9)) then
                     if(zone_comp_all(ii1) == 9) then
                        flow_rate4 = flow_rate4 + cell%rhol(ii)*flux_l_nf(i1)
                     else
                        flow_rate4 = flow_rate4 - cell%rhol(kk)*flux_l_nf(i1)
                     end if
                  end if
               end if
            else
               if(zone_comp_all(ii1) == 9) then
                  if((zone_comp_all(ii1) == 9 .and. zone_comp_all(kk1) == 10) .or. (zone_comp_all(ii1) == 10 .and. zone_comp_all(kk1) == 9)) then
                     if(zone_comp_all(ii1) == 9) then
                        flow_rate4 = flow_rate4 + cell%rhol(ii)*flux_l_nf(i1)
                     else
                        flow_rate4 = flow_rate4 - cell%rhol(kk)*flux_l_nf(i1)
                     end if
                  end if
               end if
            end if        
!dsj-SG
            if(kk <= ncell_fluid) then
               if(zone_comp(ii) == 8 .or. zone_comp(ii) == 9) then
                  if((xloc(ii,3)-6.0d0)*(xloc(kk,3)-6.0d0) < 0.d0) then
                     if(xloc(ii,3) < 6.0d0) then
                        flow_rate5 = flow_rate5 + cell%rhol(kk)*flux_l_nf(i1)
                     else
                        flow_rate5 = flow_rate5 - cell%rhol(ii)*flux_l_nf(i1)
                     end if
                  end if
               end if
            else
               if(zone_comp(ii) == 8) then
                  if((xloc(ii,3)-6.0d0)*(xloc(kk,3)-6.0d0) < 0.d0) then
                     if(xloc(ii,3) < 6.0d0) then
                        flow_rate5 = flow_rate5 + cell%rhol(kk)*flux_l_nf(i1)
                     else
                        flow_rate5 = flow_rate5 - cell%rhol(ii)*flux_l_nf(i1)
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
         if(myrank == 0) write(funit, '(F15.7, x, F15.7, x, F15.7, x, F15.7, x, F15.7, x, F15.7)') time, ABS(flow_rate1), ABS(flow_rate2), ABS(flow_rate3), ABS(flow_rate4), ABS(flow_rate5)
!         
         
         
      end subroutine PrintRcsFlowRate_NuScale

      subroutine PrintTemperatureCoreOut_NuScale()
         use Zzone, only:ncell_fluid
         use Vol_DATA, only: cell
         use Zcore, only: myrank, np
         use Ztimecon, only: time
         use NuScale
         
         implicit none
         INTEGER :: i, cnt !, cnt_tot
         REAL(8) :: T !, T_tot
         REAL(8), SAVE :: pre_value
         REAL(8) :: current_value
         LOGICAL, save :: initial = .true.
         INTEGER :: funit = 123124

         if(initial) then
            initial = .false.
            if(myrank == 0) OPEN(funit, file='T_core_out.dat')
            pre_value = 540.d0
         end if

         T = 0.d0
         cnt = 0
         do i = 1, ncell_fluid
            if(zone_comp(i) == 3) then  !core-out
               T = T + cell%tl(i)
               cnt = cnt+1
            end if
         end do

         IF(np.gt.1) THEN
             CALL allreducei_r1(T)
             CALL allreducei_i1(cnt)
         ENDIF
 
         current_value = T/dble(cnt)
         if(myrank == 0) write(funit, '(F15.7, x, F15.7)') time, current_value
         pre_value = current_value
!         
      end subroutine PrintTemperatureCoreOut_NuScale

      subroutine PressureDrop_NuScale()
         use Zzone            ,ONLY: ncell_fluid
         use Vol_DATA         ,ONLY: cell
         use Zvector          ,ONLY: ul_o
         USE STM_TBL_cupid    ,ONLY: st_tbl,nt,np,ns,ns2,ndxstd
         USE Zmodel           ,ONLY: resist     
         USE Zvector          ,ONLY: ul_o
         use NuScale
!         
         implicit none
         INTEGER :: i, zone_ix, it
         REAL(8), PARAMETER :: rcs_flow = 3000.d0 !600.d0  !...Target RCS mass flow rate
         REAL(8), PARAMETER :: time_steady = 0.d0 !...Target RCS mass flow rate
         REAL(8) :: vel, dpdz
         REAL(8) :: rho_hot, rho_cold, rho_ave
         INTEGER, PARAMETER :: fUnit = 1293
         LOGICAL, save :: initial = .true.
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
            do i = 1, ncell_fluid
               zone_ix = zone_comp(i)
               select case (zone_ix)
!old-dsj                   
               case(1, 2, 3)  !...Core active / in / out
                 dpdz = 69.75d3 / (zone_comp_height(1) + zone_comp_height(2) + zone_comp_height(3)) 
                 vel = rcs_flow / rho_ave / zone_comp_area(2)
!KSMR
               !case(1, 10, 11)  !...Core active / in / out   
               !  dpdz = 0.d0
               !  vel = rcs_flow / rho_ave / zone_comp_area(1)
!                 
               case(4) !...core_bypass
                 dpdz = 29.2d10 / (zone_comp_height(4))
                 vel = rcs_flow / rho_hot / zone_comp_area(4)
               case(5) !...riser
                 dpdz =29.2d3 / (zone_comp_height(5))
                 vel = rcs_flow / rho_cold / zone_comp_area(5) 
               case(6) !...up
                 dpdz = 20.d3 !20.d3 
                 vel = rcs_flow / rho_cold / zone_comp_area(6)
               case(7) !...pre_sg
                 dpdz = 44.d3
                 vel = rcs_flow / rho_cold / zone_comp_area(7)
               case(8) !...sg
                 dpdz = 100.1d3 / (zone_comp_height(8))
                 vel = rcs_flow / rho_cold / zone_comp_area(8)
               case(9) !...downcomer
                 dpdz = 20.0d3 / (zone_comp_height(9))
                 vel = rcs_flow / rho_cold / zone_comp_area(9)
               case(10) !...LP
                 dpdz = 20.0d3 / (zone_comp_height(10))
                 vel = rcs_flow / rho_cold / zone_comp_area(10)
               case default
                  cycle
               end select
   
               resist(i) = dpdz/vel/vel
               resist(i) = resist(i) !* 0.1d0 !test
            end do
         end if
!
         do i = 1, ncell_fluid
            cell%vfwl(i) = resist(i)*ul_o(i)
         end do
!         
      end subroutine PressureDrop_NuScale
