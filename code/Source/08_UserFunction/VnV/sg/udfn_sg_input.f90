!
      SUBROUTINE udfn_sg_input 
!
      USE Zmpi         , ONLY: jperm
      USE Zzone        , ONLY: ncell_fluid,ncell_fluid_all
      USE Zcore        , ONLY: np
      USE Zsg          , ONLY: n_group,n_1d,sd_cell,mult_1d_group,mult_1d_cell,mult_3d_cell1,mult_3d_cell2,    &
                               ar_tube,h_tube,vol_1d,ht_area,f1_mult,f2_mult,tube_length,tin_1d,pin_1d,        &
                               do_tube,th_tube,dp_primary,pr_flow_rated,max_1d,nr_tube,q_rated,mult_cell,      &
                               ng_tube,time_sg_heat_tune,pitch,ih,iavb,idc,izp,ihp,igr,j1d,idc,z_econ,dz_fw,   &
                               relax_hb,tune_hb,iboil_model
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER num_mult
      INTEGER i,ie,g,m
!.....Local arrays
      INTEGER,ALLOCATABLE::igr_tmp(:),j1d_tmp(:),idc_tmp(:),m_cell1_tmp(:),m_cell2_tmp(:)
!
      ALLOCATE(igr_tmp(ncell_fluid_all),j1d_tmp(ncell_fluid_all),idc_tmp(ncell_fluid_all))
      igr_tmp(:)=0
      j1d_tmp(:)=0
      idc_tmp(:)=0
!
      OPEN(223,file='utube_geo.in')
      OPEN(224,file='sg_param.in')
!
      READ(223,*) n_group,max_1d
!
      CALL udfn_sg_allocate
!
      ALLOCATE(n_1d(n_group),sd_cell(n_group,max_1d))
      ALLOCATE(ar_tube(n_group),h_tube(n_group,max_1d),vol_1d(n_group,max_1d))
      ALLOCATE(mult_cell(n_group,max_1d),tube_length(n_group),ng_tube(n_group))
      tube_length(:)=0.0d0
!
      DO g=1,n_group
         READ(223,*) n_1d(g),(sd_cell(g,m),m=1,n_1d(g))
         READ(223,*) ar_tube(g),ng_tube(g)
         READ(223,*) (h_tube(g,m),m=1,n_1d(g))
         READ(223,*) (vol_1d(g,m),m=1,n_1d(g))
         READ(223,*) (ht_area(g,m),m=1,n_1d(g))
         READ(223,*) (mult_cell(g,m),m=1,n_1d(g))
         READ(223,*) (ih(g,m),m=1,n_1d(g))
         READ(223,*) (iavb(g,m),m=1,n_1d(g))
         READ(223,*) (izp(g,m),m=1,n_1d(g))
         READ(223,*) (ihp(g,m),m=1,n_1d(g))
         DO m=1,n_1d(g)
            tube_length(g)=tube_length(g)+h_tube(g,m)
         ENDDO
      ENDDO
!
      READ(223,*) num_mult
!
      ALLOCATE(mult_1d_group(num_mult),mult_1d_cell(num_mult),mult_3d_cell1(num_mult),mult_3d_cell2(num_mult))
      ALLOCATE(f1_mult(num_mult),f2_mult(num_mult),m_cell1_tmp(num_mult),m_cell2_tmp(num_mult))
!
      DO i=1,num_mult
         READ(223,*) mult_1d_group(i),mult_1d_cell(i),m_cell1_tmp(i),m_cell2_tmp(i)
         READ(223,*) f1_mult(i),f2_mult(i)
      ENDDO
!
      DO i=1,ncell_fluid_all
         READ(223,*) idc_tmp(i)
      ENDDO
!
      READ(223,*) z_econ,dz_fw          ! Height of economizer and feedwater inlet window
!
!.....Sg design input parameters
!
      READ(224,*) pin_1d                ! Inlet pressure
      READ(224,*) dp_primary            ! Pressure drop of the utube
      READ(224,*) tin_1d                ! Inlet temperature
      READ(224,*) pr_flow_rated         ! Primary coolant rated flow rate
      READ(224,*) q_rated               ! Rated thermal power of SG
      READ(224,*) nr_tube               ! Number of radial nodes for heat conduction calculation of utube
      READ(224,*) time_sg_heat_tune     ! Time for the initialization of heat transfer rate
      READ(224,*) relax_hb,tune_hb      ! Relaxation and tunning factor for the boiling heat transfer coeffient
      READ(224,*) do_tube,th_tube,pitch ! U-tube outer diameter, thickness, pitch
      READ(224,*) iboil_model           ! Option for boiling heat transfer model
!
!.....Initialize sg primary side variables
!
      CALL udfn_sg_primary_ini
!
      DO g=1,n_group
         DO m=1,n_1d(g)
            i=sd_cell(g,m)
            igr_tmp(i)=g
            j1d_tmp(i)=m
         ENDDO
      ENDDO
!
      IF(np.eq.1)THEN
         idc(:)=idc_tmp(:)
         igr(:)=igr_tmp(:)
         j1d(:)=j1d_tmp(:)
         mult_3d_cell1(:)=m_cell1_tmp(:)
         mult_3d_cell2(:)=m_cell2_tmp(:)
      ELSE
         DO i=1,ncell_fluid
            ie=jperm(i)
            idc(i)=idc_tmp(ie)
            igr(i)=igr_tmp(ie)
            j1d(i)=j1d_tmp(ie)
            DO m=1,num_mult
               IF(ie.eq.m_cell1_tmp(m)) mult_3d_cell1(m)=i
               IF(ie.eq.m_cell2_tmp(m)) mult_3d_cell2(m)=i
            ENDDO
         ENDDO
      ENDIF
!
      DEALLOCATE(igr_tmp,j1d_tmp,idc_tmp,m_cell1_tmp,m_cell2_tmp)
!
      END SUBROUTINE udfn_sg_input
