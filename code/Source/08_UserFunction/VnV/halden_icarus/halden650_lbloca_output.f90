!
      SUBROUTINE halden650_lbloca_output !pik-halden
!
      USE VOL_DATA        , ONLY: cell
      USE SOLID_DATA      , ONLY: solid      
      USE Zmpi            , ONLY: jperm
      USE Zzone           , ONLY: ncell_fluid,ncell_fluid_all,ncell_cond_all,ncell_cond
      USE Zcore           , ONLY: myrank
      USE Ztimecon        , ONLY: time
      USE Zvector         , ONLY: vl_o
      USE Zwall_HTC       , ONLY: twall_rv
      USE Zrv_ncell       , ONLY: ncell_fluid_core,ncell_fluid_core_all,cupid_cell_channel
!
      USE Zqvol           , ONLY: qvol_ice_solid
      USE Zrad_comp      , ONLY: qrad_rod,qrad_sol,nproperty_out_tmp,ncell_out
      USE Zrad_comp      , ONLY: heat_rad,t1,t2,nset,nsize,cell_out_tmp,frac_out_tmp
!      
      USE viewData_common , ONLY: cupid_rv_jperm
!      
      IMPLICIT NONE
!
      CHARACTER*20:: tsol_name(50),hrad_name(50),lastname
      INTEGER :: i,j,ii
      INTEGER,DIMENSION(:),ALLOCATABLE::itemp1
      LOGICAL,SAVE :: initial=.true.
      REAL(8),SAVE :: print_time,print_interval
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: temp0,temp0_c,temp1
      REAL(8),ALLOCATABLE::tsol_out(:,:),qrad_out(:,:)
!
!.....file open and print head
!      
      IF(initial)THEN
         initial=.FALSE.
         print_time=0.0d0
         print_interval=1.d0
         IF(myrank.eq.0) THEN
            OPEN(331, file='t1_t2_heat_rad.dat')
            OPEN(332, file='2D_rod_twall.dat')
            OPEN(333, file='3D_cond_tsol.dat')
            OPEN(334, file='3D_cond_qsol.dat')
            OPEN(335, file='VD29_halden650_5_ref.dat')
            OPEN(336, file='qrod_qsol.dat')            
            WRITE(332,"(a15,200i8)") 'time',(i,i=1,ncell_fluid_core_all),(i,i=1,ncell_fluid_core_all)
            WRITE(333,"(a15,200i8)") 'time',(i,i=1,ncell_cond_all)
            WRITE(334,"(a15,200i12)") 'time',(i,i=1,ncell_cond_all)
            DO i=1,ncell_out
               WRITE(lastname,*)i !lastname size should 20
               lastname=ADJUSTL(lastname)
               tsol_name(i)='tsol'//TRIM(lastname)//'(C)'
               hrad_name(i)='qrad'//TRIM(lastname)//'(W)'
            ENDDO
            WRITE(335,"(1a15,100a12,5a12)") 'time(s) ',(tsol_name(i),i=1,ncell_out),(hrad_name(i),i=1,ncell_out)
            WRITE(336,"(a15,200i8)") 'time',(i,i=1,ncell_fluid_all),(i,i=1,ncell_cond_all)
         ENDIF            
      ENDIF   
!  
!.....print value of variables
!      
      IF(time.ge.print_time)THEN
!         
         print_time=print_time+print_interval
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
         IF(myrank.eq.0) THEN
            ALLOCATE(temp0(ncell_fluid_all,5))
            ALLOCATE(temp0_c(ncell_cond_all,3))
            ALLOCATE(itemp1(ncell_fluid_core_all))
            ALLOCATE(temp1(ncell_fluid_core_all,2))
         ELSE
            ALLOCATE(temp0(1,5))
            ALLOCATE(temp0_c(1,3))
            ALLOCATE(itemp1(1))
            ALLOCATE(temp1(1,2))            
         ENDIF
!         
         CALL gatherv_i(cupid_rv_jperm,ncell_fluid_core,itemp1,ncell_fluid_core_all,2) !1~ncell_fluid_core=itemp1(1~ncell_fluid)
         CALL gatherv_r(twall_rv(1,1),ncell_fluid_core,temp1(1,1) ,ncell_fluid_core_all,2) 
         CALL gatherv_r(twall_rv(1,2),ncell_fluid_core,temp1(1,2) ,ncell_fluid_core_all,2) 
!
         CALL gatherv_r(cell%alphag  ,ncell_fluid,temp0(1,1),ncell_fluid_all,0)
         CALL gatherv_r(vl_o(1,1)    ,ncell_fluid,temp0(1,2),ncell_fluid_all,0)
         CALL gatherv_r(qrad_rod     ,ncell_fluid,temp0(1,3),ncell_fluid_all,0)
!         
         CALL gatherv_r(solid%tsol    ,ncell_cond,temp0_c(1,1),ncell_cond_all,1)         
         CALL gatherv_r(qrad_sol      ,ncell_cond,temp0_c(1,2),ncell_cond_all,1)         
         CALL gatherv_r(qvol_ice_solid,ncell_cond,temp0_c(1,3),ncell_cond_all,1)         
!
         IF(myrank.eq.0) THEN
            WRITE(332,"(1e17.10,200f8.1)") time,(temp1(i,1),i=1,ncell_fluid_core_all),(temp1(i,2),i=1,ncell_fluid_core_all)
            WRITE(333,"(1e17.10,200f8.1)") time,(temp0_c(i,1),i=1,ncell_cond_all)
            WRITE(334,"(1e17.10,200e12.5)")time,(temp0_c(i,2),i=1,ncell_cond_all)
            WRITE(336,"(1e17.10,200f8.1)") time,(temp0(i,3),i=1,ncell_fluid_all),(temp0_c(i,2),i=1,ncell_cond_all)
            ALLOCATE(tsol_out(ncell_out,2))
            ALLOCATE(qrad_out(ncell_out,2))
            tsol_out=0.0d0
            qrad_out=0.0d0
            DO i=1,ncell_out
               IF(nproperty_out_tmp(i).eq.1)THEN
                  tsol_out(i,1)=temp1(itemp1(cell_out_tmp(i,1)),1)-273.15d0 !2d rod
                  tsol_out(i,2)=temp1(itemp1(cell_out_tmp(i,2)),1)-273.15d0 !2d rod
                  qrad_out(i,1)=temp0(cell_out_tmp(i,1),3)
                  qrad_out(i,2)=temp0(cell_out_tmp(i,2),3)
               ELSE
                  tsol_out(i,1)=temp0_c(cell_out_tmp(i,1),1)-273.15d0       !solid 
                  tsol_out(i,2)=temp0_c(cell_out_tmp(i,2),1)-273.15d0       !solid 
                  qrad_out(i,1)=temp0_c(cell_out_tmp(i,1),2)
                  qrad_out(i,2)=temp0_c(cell_out_tmp(i,2),2)                  
               ENDIF   
               qrad_out(2,1)=qrad_out(2,1)*frac_out_tmp(2,1)+qrad_out(2,2)*frac_out_tmp(2,2)
               tsol_out(2,1)=tsol_out(2,1)*frac_out_tmp(2,1)+tsol_out(2,2)*frac_out_tmp(2,2)
            ENDDO   
!            
            WRITE(335,"(1e17.10,100f12.3,100e12.4)") time,(tsol_out(i,1),i=1,ncell_out),(qrad_out(i,1),i=1,ncell_out)
            DO j=1,nset
                  WRITE(331,"(1e17.10,200f8.1)")time,(t1(i,j)-273.15d0,t2(i,j)-273.15d0,heat_rad(i,j),i=1,nsize(1,j))
            ENDDO
            DEALLOCATE(tsol_out,qrad_out)
         ENDIF   
!         
         DEALLOCATE(temp0)
         DEALLOCATE(temp0_c)            
         DEALLOCATE(temp1,itemp1,cupid_rv_jperm)
!         
      ENDIF !time.ge.print_time                
!      
      RETURN
      ENDSUBROUTINE halden650_lbloca_output
         

   
