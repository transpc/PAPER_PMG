!
      SUBROUTINE rv_calc_geo_2d
!
      USE Zzone       , ONLY: ncell_fluid
      USE Zcore       , ONLY: myrank
      USE Zmpi        , ONLY: jperm,celem
      USE Zconst1     , ONLY: cplmaster
      USE Zrv_hts_2d  , ONLY: nrod_2d,nz0_2d,nr_2d,ht_geo_2d,nmat_2d,ri_2d,t_fuel,nqvol,           &
                               qvol_time,qvol_norm_2d,qvol_tot,height_fuel,qvol_2f0,qvol_2f,     &
                               rweight_power,l_ht_str_2d_qcell,qcell_to_qvol_mul
      USE Zrv_mpi     , ONLY: jperm_fuel_rod,ncell_fuel_rod_p,celem_fuel_rod,celem_fluid_core,     &
                               jperm_fluid_core
      USE Zrv_ncell   , ONLY: ncell_fuel_rod,cupid_cell_hts2d,ncell_fluid_core,n_channel_fluid,    &
                               nz_fluid,cupid_cell_channel,channel_cell_hts2d,num_ch,num_fuel_rod, &
                               nz_fine,nz_fuel_rod,nrod_fuel_rod,ncell_fuel_rod_all,  &
                               ncell_fluid_core_all
      USE Zporous     , ONLY: l_subchannel
      USE Zio_unit    , ONLY: unit_log
      USE Zrv_hts_2d     , ONLY: ri_2d_opt,ri_2d_input,qvol_norm_2d,qvol_tot_channel,n_ch_frap,n_ch_frap_input !icarus2002
!
      IMPLICIT NONE
!
      INTEGER i,j,k,m,n,ii,io
      INTEGER ch_opt,rod_opt, nmat_opt,temp_opt,qvolz_opt,qvolr_opt,qvolt_opt
      INTEGER read_stat,idummy,n_rweight_power,nrod_2d_frap
!
      INTEGER,ALLOCATABLE :: cupid_cell_hts2d_tmp(:),nf_input(:,:),n_ch(:),n_channel_fluid_tmp(:)
      INTEGER,ALLOCATABLE :: nz_fluid_tmp(:),nrod_fuel_rod_tmp(:),nz0_fuel_rod_tmp(:),nz_fuel_rod_tmp(:)
      INTEGER,ALLOCATABLE :: neigh_fuel_rod_tmp(:,:),cupid_cell_channel_tmp(:),channel_cell_hts2d_tmp(:)
!
      REAL(8) num_rod_tot
      REAL(8),ALLOCATABLE :: t_2d_input(:,:,:),t_2d_tmp(:,:),qvol_norm_z_input(:,:),qvol_norm_r_input(:,:),qvol_2d_tmp(:,:)
      REAL(8),ALLOCATABLE :: qcell_to_qvol_mul_tmp(:,:)

!
      OPEN(53,file='ht_str_2d_qcell.in',status='old',iostat=io)   
      IF(io.eq.0)THEN
         l_ht_str_2d_qcell=.true.
      ELSE
         l_ht_str_2d_qcell=.false.
      ENDIF  
      CLOSE(53)
!      
      OPEN(52,file='ht_str_2d.in',status='old',iostat=io)
      IF(io.gt.0) RETURN
!
!.....Number of fuel channel and axial cells
!
      READ(52,*) num_ch,nz0_2d,ch_opt,rod_opt, nmat_opt,temp_opt,qvolz_opt,qvolr_opt,qvolt_opt
      IF(nz0_2d.eq.0)RETURN
      IF(myrank.eq.0)WRITE(*,"(11x,a)")'Reading ht_str_2d.in...'
      IF(myrank.eq.0)WRITE(unit_log,"(11x,a)")'Reading ht_str_2d.in...'      
!
!.....Number of fuel rods for each channel
!
      ALLOCATE(num_fuel_rod(num_ch))
      IF(ch_opt.eq.1)THEN
         DO i=1,num_ch
            READ(52,*) num_fuel_rod(i)
         ENDDO
      ELSE
         READ(52,*) num_fuel_rod(1)
         DO i=1,num_ch
            num_fuel_rod(i)=num_fuel_rod(1)
         ENDDO          
      ENDIF     
!
!.....CUPID cell index for 2D hts (num_ch x nz0_2d)
!
      ALLOCATE(nf_input(num_ch,nz0_2d))
      DO i=1,num_ch
         READ(52,*) (nf_input(i,j),j=1,nz0_2d)
      ENDDO
!
!.....Number of radial cells
!
      READ(52,*) nrod_2d,nrod_2d_frap
!
!.....Coolant channel index of each fuel rod
!
      ALLOCATE(n_ch(nrod_2d))
      ALLOCATE(n_ch_frap_input(nrod_2d))
      IF(rod_opt.eq.1)THEN
         DO i=1,nrod_2d
            READ(52,*) n_ch(i),n_ch_frap_input(i)
         ENDDO
      ELSE
         DO i=1,nrod_2d
            n_ch(i)=i
            n_ch_frap_input(i)=nrod_2d_frap
         ENDDO         
      ENDIF   
      n_ch_frap(:)=0   
!
!.....Number of fine nodes in z-direction, Number of radial cells, Heat structure geometry
!
      READ(52,*) nz_fine,nr_2d,ht_geo_2d
!
!.....Height of fuel rod
!
      READ(52,*) height_fuel
!
!.....Material property and coordinates of radial cells
!
      !!!ALLOCATE(nmat_2d(nr_2d),ri_2d(nr_2d))
      !!!DO i=1,nr_2d
      !!!   READ(52,*) nmat_2d(i),ri_2d(i)
      !!!ENDDO
      
      ALLOCATE(nmat_2d(nrod_2d,nr_2d),ri_2d_input(nrod_2d,nr_2d))
      ALLOCATE(ri_2d(nr_2d))
      ri_2d_opt=nmat_opt
      IF(nmat_opt.eq.1)THEN
         DO i=1,nrod_2d
            DO j=1,nr_2d
               READ(52,*) nmat_2d(i,j),ri_2d_input(i,j)
            ENDDO      
         ENDDO
!
!........For a special case of 1 rod per 1 core         
         ri_2d(:)=ri_2d_input(1,:) 
         DO i=1,nrod_2d
            IF(myrank.eq.i-1)THEN !icarus2002
               ri_2d(:)=ri_2d_input(i,:)
            ENDIF   
         ENDDO   
      ELSE
         DO i=1,1
            DO j=1,nr_2d
               READ(52,*) nmat_2d(i,j),ri_2d_input(i,j)
            ENDDO      
         ENDDO    
         DO i=1,nrod_2d
            nmat_2d(i,:)=nmat_2d(1,:)
         ENDDO 
         ri_2d(:)=ri_2d_input(1,:)         
      ENDIF
!
!.....Initial temperature of fuel rods
!
      ALLOCATE(t_2d_input(nrod_2d,nz0_2d,nr_2d))
      ! temp_opt
      ! - 0: read ONLY 1 line and copy to whole region
      ! - 1: Specify whole region with input data
      IF(temp_opt.eq.1)THEN
         DO i=1,nrod_2d
            DO j=1,nz0_2d
               READ(52,*) (t_2d_input(i,j,k),k=1,nr_2d)
            ENDDO
         ENDDO
      ELSE
         DO i=1,1
            DO j=1,nz0_2d
               READ(52,*) (t_2d_input(i,j,k),k=1,nr_2d)
            ENDDO
         ENDDO
         DO i=1,nrod_2d
             t_2d_input(i,:,:)=t_2d_input(1,:,:)
         ENDDO                
      ENDIF        
!
!.....Normalized axial power shape
!
      ALLOCATE(qvol_norm_z_input(nrod_2d,nz0_2d))
      ! qvolz_opt
      ! - 0: read ONLY 1 line and copy to whole region
      ! - 1: Specify whole region with input data
      IF(qvolz_opt.eq.1)THEN
         DO i=1,nrod_2d
            READ(52,*) (qvol_norm_z_input(i,j),j=1,nz0_2d)
         ENDDO
      ELSE
         READ(52,*) (qvol_norm_z_input(1,j),j=1,nz0_2d)
         DO i=1,nrod_2d
            qvol_norm_z_input(i,:)=qvol_norm_z_input(1,:)
         ENDDO               
      ENDIF       
!
!.....Normalized radial power shape
!
      ALLOCATE(qvol_norm_r_input(nrod_2d,nr_2d))
      
      IF(qvolr_opt.eq.1)THEN
         DO i=1,nrod_2d
            READ(52,*) (qvol_norm_r_input(i,j),j=1,nr_2d)
         ENDDO
      ELSE   
         READ(52,*) (qvol_norm_r_input(1,j),j=1,nr_2d)
         DO i=1,nrod_2d
            qvol_norm_r_input(i,:)=qvol_norm_r_input(1,:)
         ENDDO            
      ENDIF   
!
!.....Number of time intervals and initial power
!
!     
!.....Transient time and normalized power
!
      ALLOCATE(qvol_tot_channel(nrod_2d))
      IF(qvolt_opt.eq.1)THEN
         READ(52,*) nqvol,(qvol_tot_channel(i),i=1,nrod_2d)
          ALLOCATE(qvol_time(nqvol),qvol_norm_2d(nrod_2d,nqvol))
         qvol_tot=0.0d0
         DO i=1,nrod_2d
            qvol_tot=qvol_tot+qvol_tot_channel(i)        
         ENDDO
         DO j=1,nqvol
            READ(52,*) qvol_time(j),(qvol_norm_2d(i,j),i=1,nrod_2d)
         ENDDO
      ELSE
         DO i=1,1
            READ(52,*) nqvol,qvol_tot
            ALLOCATE(qvol_time(nqvol),qvol_norm_2d(nrod_2d,nqvol))
            DO j=1,nqvol
               READ(52,*) qvol_time(j),qvol_norm_2d(1,j)
            ENDDO
         ENDDO  
         DO i=1,nrod_2d
            qvol_norm_2d(i,:)=qvol_norm_2d(1,:)
         ENDDO    
         qvol_tot_channel(:)=0.0d0        
      ENDIF
!
!.....Radial power distribution
!
     ALLOCATE(rweight_power(nrod_2d))
     READ(52,*,iostat=read_stat)n_rweight_power
     IF(read_stat.eq.0)THEN
        IF(myrank.eq.0)WRITE(*,"(11x,a,1i5)")'Reading rweight_power:n_rweight_power=',n_rweight_power
        IF(n_rweight_power.ne.nrod_2d)THEN
            WRITE(*,"(a,1i5,a)")'n_wreight should be equal to',nrod_2d,'.'
            PAUSE
            STOP          
        ENDIF
        DO i=1,n_rweight_power
           READ(52,*)idummy,rweight_power(i) 
        ENDDO  
     ELSE
        IF(myrank.eq.0)WRITE(*,"(11x,a)")'Setting rweight_power as 1.0.'
        rweight_power(:)=1.0d0
     ENDIF   
!
!.....Total number of fluid cells connected to 2D hts
!
      ncell_fluid_core_all=num_ch*nz0_2d
!
!.....Total number of 2D hts cells
!
      ncell_fuel_rod_all=nrod_2d*nz0_2d*nz_fine
!
      n=ncell_fuel_rod_all
      ALLOCATE(neigh_fuel_rod_tmp(2,n),cupid_cell_hts2d_tmp(n))
      ALLOCATE(nz_fuel_rod_tmp(n),nz0_fuel_rod_tmp(n),nrod_fuel_rod_tmp(n))
      ALLOCATE(t_2d_tmp(n,nr_2d),channel_cell_hts2d_tmp(n),qvol_2d_tmp(n,nr_2d))
      ALLOCATE(qcell_to_qvol_mul_tmp(n,nr_2d))
!
      neigh_fuel_rod_tmp(:,:)=0
      cupid_cell_hts2d_tmp(:)=0
      channel_cell_hts2d_tmp(:)=0
      nz0_fuel_rod_tmp(:)=0
      nz_fuel_rod_tmp(:)=0
      t_2d_tmp(:,:)=0.0d0
      qvol_2d_tmp(:,:)=0.0d0
!
      n=ncell_fluid_core_all
      ALLOCATE(celem_fluid_core(n),n_channel_fluid_tmp(n),nz_fluid_tmp(n),cupid_cell_channel_tmp(n))
      celem_fluid_core(:)=0
      n_channel_fluid_tmp(:)=0
      nz_fluid_tmp(:)=0
      cupid_cell_channel_tmp(:)=0
!
!.....Coloring of 2D hts fluid cells for domain decompisition(celem_fluid_core)
!.....Set the channel index(n_channel_fluid_tmp), axial node number(nz_fluid_tmp)
!     and CUPID cell index(cupid_cell_channel_tmp) of a 2D hts fluid cell
!
      m=0
      DO i=1,num_ch
         DO j=1,nz0_2d
            m=m+1
            ii=nf_input(i,j)
            celem_fluid_core(m)=celem(ii)
            n_channel_fluid_tmp(m)=i
            nz_fluid_tmp(m)=j
            cupid_cell_channel_tmp(m)=ii
         ENDDO
      ENDDO
!
!.....Define axial connectivity for 2D heat structure for the fuel rod
!
      ALLOCATE(celem_fuel_rod(ncell_fuel_rod_all))
      celem_fuel_rod(:)=0
!
!.....Set connectivity of 2D hts cells (neigh_fuel_rod_tmp)
!.....Coloring of 2D hts cells for domain decompisition (celem_fuel_rod)
!.....Set the cell index:
!     - 2D hts fluid cell index(channel_cell_hts2d_tmp)
!     - Axial coarse node number(nz0_fuel_rod_tmp), 
!     - Axial fine node number(nz_fuel_rod_tmp)
!     - Coolant channel index(nrod_fuel_rod_tmp),
!     - CUPID cell index of a 2D hts cell(cupid_cell_hts2d_tmp)
!.....Initialize 2D hts temperature(t_2d_tmp)
!
      m=0
      DO i=1,nrod_2d
         DO j=1,nz0_2d
            ii=nf_input(n_ch(i),j)
            DO k=1,nz_fine
!
               m=m+1
               IF(j.eq.1.and.k.eq.1)THEN
                  neigh_fuel_rod_tmp(1,m)=0
                  neigh_fuel_rod_tmp(2,m)=m+1
               ELSEIF(j.eq.nz0_2d.and.k.eq.nz_fine)THEN
                  neigh_fuel_rod_tmp(1,m)=m-1
                  neigh_fuel_rod_tmp(2,m)=0
               ELSE
                  neigh_fuel_rod_tmp(1,m)=m-1
                  neigh_fuel_rod_tmp(2,m)=m+1
               ENDIF
!
               celem_fuel_rod(m)=celem(ii)
               cupid_cell_hts2d_tmp(m)=ii
               channel_cell_hts2d_tmp(m)=(n_ch(i)-1)*nz0_2d+j
               nz0_fuel_rod_tmp(m)=j
               nz_fuel_rod_tmp(m)=nz_fine*(j-1)+k
               nrod_fuel_rod_tmp(m)=i
!
               t_2d_tmp(m,:)=t_2d_input(i,j,:)
!
            ENDDO
         ENDDO
      ENDDO
!
!.....Set heat generation rate for the fuel rods
!
      num_rod_tot=0.0d0
      DO i=1,num_ch
         num_rod_tot=num_rod_tot+num_fuel_rod(i)
      ENDDO
!
      DO i=1,ncell_fuel_rod_all
         j=nrod_fuel_rod_tmp(i)
         k=nz0_fuel_rod_tmp(i)
         DO m=1,nr_2d
            IF(cplmaster.gt.0)then
               qvol_norm_z_input(j,k)=1.0d0
               qvol_2d_tmp(i,m)=1.d0*(qvol_norm_z_input(j,k)/nz_fine)*qvol_norm_r_input(j,m)*rweight_power(j)
            ELSE
               IF(l_ht_str_2d_qcell)THEN
                  qcell_to_qvol_mul_tmp(i,m)=(1.0d0/num_fuel_rod(j))*qvol_norm_r_input(j,m) !*(1.0d0/nz_fine)*rweight_power(j)
               ELSE 
                  IF(qvolt_opt.eq.0) THEN 
                     IF(l_subchannel)then
                        qvol_2d_tmp(i,m)=qvol_tot/dble(nrod_2d)*(qvol_norm_z_input(j,k)/nz_fine)*qvol_norm_r_input(j,m)*rweight_power(j)
                     ELSE
                        qvol_2d_tmp(i,m)=qvol_tot/num_rod_tot*(qvol_norm_z_input(j,k)/nz_fine)*qvol_norm_r_input(j,m)*rweight_power(j)
                     ENDIF
                  ELSE
                     qvol_2d_tmp(i,m)=qvol_tot_channel(j)/num_fuel_rod(j)*(qvol_norm_z_input(j,k)/nz_fine)*qvol_norm_r_input(j,m)*rweight_power(j)
                  ENDIF
               ENDIF   
            ENDIF
         ENDDO
      ENDDO
!
!.....Domain decomposition for 2D hts cells
!
      CALL rv_subdomain_info_fuel_rod(ncell_fuel_rod_all,neigh_fuel_rod_tmp)
!
!.....Domain decomposition for 2D hts fluid cells
!
      CALL rv_subdomain_info_fluid_core(ncell_fluid_core_all)
!
!.....Allocate global variables 
!
      CALL rv_allocate_hts_2d
!
      CALL rv_allocate_var_2d
!
!.....Calculate geometric input
!
      CALL rv_geo_hts_2d
!
!.....Assign local values from global values for parallel calculation (2D hts fluid cells)
!.....For serial calculation, local values=global values
!
      n=ncell_fluid_core
      ALLOCATE(n_channel_fluid(n),nz_fluid(n),cupid_cell_channel(n))
      n_channel_fluid(:)=0
      nz_fluid(:)=0
      cupid_cell_channel(:)=0
!
      DO ii=1,ncell_fluid_core
         i=jperm_fluid_core(ii)
         n_channel_fluid(ii)=n_channel_fluid_tmp(i)
         nz_fluid(ii)=nz_fluid_tmp(i)
      ENDDO
!
!.....Assign local values from global values for parallel calculation (2D hts cells)
!.....For serial calculation, local values=global values
!
      n=ncell_fuel_rod
      m=ncell_fuel_rod_p
      ALLOCATE(cupid_cell_hts2d(n),channel_cell_hts2d(n),nz_fuel_rod(n),nrod_fuel_rod(n))
      ALLOCATE(qvol_2f0(n,nr_2d),qvol_2f(n,nr_2d),t_fuel(m,nr_2d))
      ALLOCATE(qcell_to_qvol_mul(n,nr_2d))
      cupid_cell_hts2d(:)=0
      channel_cell_hts2d(:)=0
      nz_fuel_rod(:)=0
      nrod_fuel_rod(:)=0
      t_fuel(:,:)=0.0d0
      qvol_2f0(:,:)=0.0d0
      qvol_2f(:,:)=0.0d0 
!
      DO ii=1,ncell_fuel_rod
         i=jperm_fuel_rod(ii)
         t_fuel(ii,:)=t_2d_tmp(i,:)
         qvol_2f0(ii,:)=qvol_2d_tmp(i,:)
         nz_fuel_rod(ii)=nz_fuel_rod_tmp(i)
         nrod_fuel_rod(ii)=nrod_fuel_rod_tmp(i)
!         
         qcell_to_qvol_mul(ii,:)=qcell_to_qvol_mul_tmp(i,:)         
      ENDDO
      qvol_2f(:,:)=qvol_2f0(:,:) 
!
!.....Find CUPID cell index connected to each 2D hts fluid cells
!
      DO i=1,ncell_fluid_core
         DO ii=1,ncell_fluid
            IF(cupid_cell_channel_tmp(jperm_fluid_core(i)).eq.jperm(ii)) cupid_cell_channel(i)=ii
         ENDDO         
         IF(cupid_cell_channel(i).eq.0)THEN
            STOP '### Fluid cell is not matched for a 2d rod 1 ###'
         ENDIF   
      ENDDO
!
      DO i=1,ncell_fuel_rod
!
!.....Find CUPID cell index connected to each 2D hts cells
!
         DO ii=1,ncell_fluid
            IF(cupid_cell_hts2d_tmp(jperm_fuel_rod(i)).eq.jperm(ii)) cupid_cell_hts2d(i)=ii
         ENDDO
         IF(cupid_cell_hts2d(i).eq.0) STOP '### Fluid cell is not matched for a 2d rod 2 ###'
!
!.....Find 2D hts fluid cell index connected to each 2D hts cells
!
         DO ii=1,ncell_fluid_core
            IF(channel_cell_hts2d_tmp(jperm_fuel_rod(i)).eq.jperm_fluid_core(ii)) channel_cell_hts2d(i)=ii
         ENDDO
         IF(channel_cell_hts2d(i).eq.0) STOP '### Fluid cell is not matched for a 2d rod 3 ###'
!
      ENDDO
!
      DEALLOCATE(neigh_fuel_rod_tmp,t_2d_input,nf_input,qvol_norm_z_input,qvol_norm_r_input)
      DEALLOCATE(cupid_cell_hts2d_tmp,n_ch,n_channel_fluid_tmp,nz_fluid_tmp,cupid_cell_channel_tmp)
      DEALLOCATE(nz_fuel_rod_tmp,qvol_2d_tmp,nrod_fuel_rod_tmp,nz0_fuel_rod_tmp)
      DEALLOCATE(qcell_to_qvol_mul_tmp)
!
      CLOSE(52)
!
      RETURN
      END SUBROUTINE rv_calc_geo_2d
