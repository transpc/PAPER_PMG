!
      SUBROUTINE icarus2002_reflood_output 
!
      USE VOL_DATA        , ONLY: cell
      USE Zmpi            , ONLY: jperm
      USE Zzone           , ONLY: ncell_fluid,ncell_fluid_all
      USE Zcore           , ONLY: myrank
      USE Zcoord1         , ONLY: xloc_tmp
      USE Ztimecon        , ONLY: time
      USE Zwall_HTC       , ONLY: twall_rv
      USE Zrv_ncell       , ONLY: ncell_fluid_core,ncell_fuel_rod,ncell_fuel_rod_all,ncell_fluid_core_all,cupid_cell_channel
      USE Zwall_HTC       , ONLY: twall_rv        
!    
      USE Zrv_ncell       , ONLY: ncell_fuel_rod,cupid_cell_hts2d,nrod_fuel_rod
      USE Zrv_hts_2d      , ONLY: twall_fuel,n_ch_frap_input  
      !USE Zcupid_fraptran , ONLY: fuel_radius_F2C,CFuel,CFRAP_null      
      !USE Zrv_mpi         , ONLY: n_rods
      USE Zconst1         , ONLY: vv_prob
!      
      IMPLICIT NONE
!
      INTEGER :: i,j,ii,No_rod
      INTEGER,SAVE :: icell_bot,icell_mid,icell_top,icell_p4,icell_p3,icell_p2,icell_p1
      INTEGER,DIMENSION(:),ALLOCATABLE::itemp1,itemp2,itemp3,cupid_rv_jperm,cupid_rod_jperm
      LOGICAL,SAVE :: initial=.true.
      REAL(8),SAVE :: print_time,print_interval
      REAL(8),DIMENSION(:),ALLOCATABLE :: fuel_radius
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: temp0,temp1,temp2
      REAL(8) :: twall_tmp1,twall_tmp2,fuel_radius1,fuel_radius2,fuel_radius3 
      
!
!.....file open and print head
!      
      IF(initial)THEN
         IF(vv_prob.eq.'icarus2002')THEN   
            icell_bot=56 
            icell_mid=53 
            icell_top=50 
            icell_p4=80 
            icell_p3=58
            icell_p2=58 
            icell_p1=48 
         ENDIF
         print_time=0.0d0
         print_interval=1.d0         
         initial=.false.         
         IF(myrank.eq.0)THEN
            OPEN(333, file='VD30_icarus2002_ref.dat')
            OPEN(335, file='icarus_twall1.dat')              
         ENDIF   
      ENDIF  
!      
!.....print value of variables
!      
      IF(time.ge.print_time)THEN
      !IF(mod(itim,iprn).eq.0)THEN         
!         
         print_time=print_time+print_interval
!         
         IF(ncell_fluid_core_all.gt.0) THEN
            ALLOCATE(cupid_rv_jperm(ncell_fluid_core))
            DO i=1,ncell_fluid_core
               ii=cupid_cell_channel(i)
               cupid_rv_jperm(i)=jperm(ii)
            ENDDO   
         ELSE
            ALLOCATE(cupid_rv_jperm(1))
            cupid_rv_jperm(:)=0
         ENDIF 
!         
         IF(ncell_fuel_rod_all.gt.0) THEN
            ALLOCATE(cupid_rod_jperm(ncell_fuel_rod))
            DO i=1,ncell_fuel_rod
               ii=cupid_cell_hts2d(i)
               cupid_rod_jperm(i)=jperm(ii)
            ENDDO   
         ELSE
            ALLOCATE(cupid_rod_jperm(1))
            cupid_rod_jperm(:)=0
         ENDIF 
!         
         IF(myrank.eq.0) THEN
            ALLOCATE(temp0(ncell_fluid_all,3))
            ALLOCATE(itemp1(ncell_fluid_core_all))
            ALLOCATE(temp1(ncell_fluid_core_all,2))
            ALLOCATE(itemp2(ncell_fuel_rod_all))            
            ALLOCATE(itemp3(ncell_fuel_rod_all)) 
            ALLOCATE(temp2(ncell_fuel_rod_all,2))
         ELSE
            ALLOCATE(temp0(1,3))
            ALLOCATE(itemp1(1))
            ALLOCATE(temp1(1,2)) 
            ALLOCATE(itemp2(1))            
            ALLOCATE(itemp3(1))             
            ALLOCATE(temp2(1,2))             
         ENDIF
!
         ALLOCATE(fuel_radius(ncell_fuel_rod))
         fuel_radius(:)=0.0d0

         !IF(CFuel.eq.1.and.CFRAP_Null.eq.0) THEN         
         !   DO j=1,n_rods
         !      IF(n_ch_frap(j).eq.1)THEN
         !         nz_ff=nz0_2d*nz_fine
         !         idx_s=(j-1)*nz_ff+1
         !         idx_f=j*nz_ff               
         !         fuel_radius(idx_s:idx_f)=fuel_radius_F2C(idx_s:idx_f)
         !      ENDIF   
         !   ENDDO   
         !ENDIF         
!         
         CALL gatherv_i(cupid_rv_jperm,ncell_fluid_core,itemp1,ncell_fluid_core_all,2) !1~ncell_fluid_core=itemp1(1~ncell_fluid)
         CALL gatherv_r(twall_rv(1,1),ncell_fluid_core,temp1(1,1) ,ncell_fluid_core_all,2)
!         
         CALL gatherv_i(cupid_rod_jperm,ncell_fuel_rod,itemp2,ncell_fuel_rod_all,3) !1~ncell_fluid_core=itemp1(1~ncell_fluid)
         CALL gatherv_i(nrod_fuel_rod,ncell_fuel_rod,itemp3,ncell_fuel_rod_all,3) !1~ncell_fluid_core=itemp1(1~ncell_fluid)
         CALL gatherv_r(fuel_radius,ncell_fuel_rod,temp2(1,1),ncell_fuel_rod_all,3)
!
         CALL gatherv_r(cell%tg   ,ncell_fluid,temp0(1,1),ncell_fluid_all,0)
         CALL gatherv_r(twall_fuel,ncell_fluid,temp0(1,2),ncell_fluid_all,0)
         CALL gatherv_r(cell%p    ,ncell_fluid,temp0(1,3),ncell_fluid_all,0)
         
         IF(myrank.eq.0)THEN
            twall_tmp1=temp0(icell_bot,2)
            twall_tmp2=temp0(icell_top,2)  
            !IF(CFuel.eq.1.and.CFRAP_Null.eq.0) THEN
            !   fuel_radius1=temp2(itemp2(icell_bot),1)
            !   fuel_radius2=temp2(itemp2(icell_mid),1)
            !   fuel_radius3=temp2(itemp2(icell_top),1)  
            !ELSE
               fuel_radius1=0.d0
               fuel_radius2=0.d0
               fuel_radius3=0.d0
            !ENDIF             
         ENDIF
                
         IF(myrank.eq.0)THEN
            WRITE(333,155) time,temp0(icell_bot,1)-273.15,temp0(icell_mid,1)-273.15,temp0(icell_top,1)-273.15,&
                           twall_tmp1-273.15,twall_tmp2-273.15 ,&
                           fuel_radius1,fuel_radius2,fuel_radius3,&
!                          cell%p(icell_p4)-cell%p(icell_p3),cell%p(icell_p2)-cell%p(icell_p1)
                           temp0(icell_p4,3)-temp0(icell_p3,3),temp0(icell_p2,3)-temp0(icell_p1,3)
            DO i=1,ncell_fluid_all
               j=itemp2(i)
               No_rod=itemp3(j)
               IF(n_ch_frap_input(No_rod).eq.1)THEN
                  WRITE(335,155) time,dble(i), xloc_tmp(i,3),xloc_tmp(i,2),temp0(i,1)-273.15,temp0(i,2)-273.15
               ENDIF   
            ENDDO 
         ENDIF   
         DEALLOCATE(temp0)         
         DEALLOCATE(temp1,itemp1,cupid_rv_jperm)
         DEALLOCATE(temp2,itemp2,cupid_rod_jperm)   
         DEALLOCATE(itemp3)            
         DEALLOCATE(fuel_radius)
!         
      ENDIF !time.ge.print_time                
!
155   FORMAT(100(e14.7,1x)) 
156   FORMAT(i3,3x,100(e14.7,1x))   
!      
      RETURN
      ENDSUBROUTINE icarus2002_reflood_output
         

   
