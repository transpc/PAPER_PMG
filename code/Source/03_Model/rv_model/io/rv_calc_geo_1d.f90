!
      SUBROUTINE rv_calc_geo_1d
!
      USE Zcore       , ONLY: myrank
      USE Zmpi        , ONLY: jperm,celem
      USE Zrv_mpi     , ONLY: jperm_hts_1d,celem_hts_1d
      USE Zrv_hts_1d  , ONLY: nr_1d,ng_hts_1d,imp_cond_1d,ht_geo_1d,nmat_1d,bcl_1d,bcr_1d,  &
                               ig_hts_1d,num_rod_1d,ri_1d,t_hts_1d,cupid_cell_1d,ncell_hts_1d,       &  
                               ncell_hts_1d_all,ig_hts_1d,slw_1d,srw_1d,ncell_hts_1d
      USE Zrv_hts_2d  , ONLY: nz0_2d
      USE Zrv_ncell   , ONLY: num_ch_1d
      USE Zzone       , ONLY: ncell_fluid
      USE Zio_unit    , ONLY: unit_log
!
      IMPLICIT NONE
!
      INTEGER i,j,k,ii,io,n,cupid_cell,g
!
      INTEGER,ALLOCATABLE::cupid_cell_1d_tmp(:),ig_hts_1d_tmp(:),inz0_2d_tmp(:)
      INTEGER,ALLOCATABLE::nf_input(:,:)
      INTEGER,ALLOCATABLE::num_ch_1dg(:),n_ch_1dg(:,:)
!
      REAL(8) num_rod_1dg_tot
      REAL(8),ALLOCATABLE::t_1d_tmp(:,:),t_1dg(:,:,:),num_rod_1d_tmp(:)  
      REAL(8),ALLOCATABLE::num_rod_1dg(:,:)
!
      OPEN(51,file='ht_str_1d.in',status='old',iostat=io)
      IF(io.gt.0) RETURN
!
!.....Number of coolant channels for 1d rods
!
      READ(51,*) num_ch_1d
      IF(num_ch_1d.eq.0) RETURN
      IF(myrank.eq.0)WRITE(*,"(11x,a)")'Reading ht_str_1d.in...'
      IF(myrank.eq.0)WRITE(unit_log,"(11x,a)")'Reading ht_str_1d.in...'      
!
!.....Cupid cell index connected to each 1d hts cell
!
      ALLOCATE(nf_input(num_ch_1d,nz0_2d))
      DO i=1,num_ch_1d
         READ(51,*) (nf_input(i,j),j=1,nz0_2d)
      ENDDO
!
!.....Total number of 1d hts cells
!
      READ(51,*) ncell_hts_1d_all
!
      n=ncell_hts_1d_all
      ALLOCATE(cupid_cell_1d_tmp(n),ig_hts_1d_tmp(n),celem_hts_1d(n))
      cupid_cell_1d_tmp(:)=0
      ig_hts_1d_tmp(:)=0
      celem_hts_1d(:)=0

!
!.....Number of radial cells, Number of rod groups, Numerical method for the 1D conduction
!
      READ(51,*) nr_1d,ng_hts_1d,imp_cond_1d

      n=ng_hts_1d
      ALLOCATE(ht_geo_1d(n),nmat_1d(n),bcl_1d(n),bcr_1d(n),ri_1d(n,nr_1d))
      ALLOCATE(slw_1d(n),srw_1d(n))
!        
      ht_geo_1d(:)=0
      nmat_1d(:)=0
      bcl_1d(:)=0
      bcr_1d(:)=0
      ri_1d(:,:)=0.0d0
      slw_1d(:)=0.0d0
      srw_1d(:)=0.0d0       
!
!.....Geometry, Material property, left surface b.c., right surface b.c., Number of rods in a rod group
!.....Radial coordinates of 1D cells
!
      ALLOCATE(num_ch_1dg(ng_hts_1d),n_ch_1dg(ng_hts_1d,num_ch_1d),num_rod_1dg(ng_hts_1d,num_ch_1d))
      num_ch_1dg(:)=0
      n_ch_1dg(:,:)=0
      num_rod_1dg(:,:)=0.0d0
!      
      DO i=1,ng_hts_1d
         READ(51,*) ht_geo_1d(i),nmat_1d(i),bcl_1d(i),bcr_1d(i),num_rod_1dg_tot,num_ch_1dg(i)
         READ(51,*) (ri_1d(i,j),j=1,nr_1d)
         IF(ht_geo_1d(i).eq.1) READ(51,*)slw_1d(i)
         srw_1d(i)=slw_1d(i)
         DO j=1,num_ch_1dg(i)
            READ(51,*)n_ch_1dg(i,j),num_rod_1dg(i,j)
         ENDDO
      ENDDO
!
      ALLOCATE(num_rod_1d_tmp(ncell_hts_1d_all),inz0_2d_tmp(ncell_hts_1d_all))      
      num_rod_1d_tmp(:)=0
      inz0_2d_tmp(:)=0
      k=0
!
!....Connectivities of 1d heat structure
!      
      DO g=1,ng_hts_1d
         DO j=1,num_ch_1d
            IF(n_ch_1dg(g,j).eq.0)CYCLE
            DO i=1,nz0_2d
               k=k+1
               ig_hts_1d_tmp(k)=g
               inz0_2d_tmp(k)=i
               num_rod_1d_tmp(k)=num_rod_1dg(g,j)
               cupid_cell=nf_input(n_ch_1dg(g,j),i)
               cupid_cell_1d_tmp(k)=cupid_cell
               celem_hts_1d(k)=celem(cupid_cell)
            ENDDO
         ENDDO      
      ENDDO 
!        
      DEALLOCATE(num_ch_1dg,n_ch_1dg,num_rod_1dg)
!      
      IF(k.ne.ncell_hts_1d_all) STOP '### Incorrect 1D heat structure input..'
!
      ALLOCATE(t_1d_tmp(ncell_hts_1d_all,nr_1d))
      t_1d_tmp(:,:)=0.0d0
!
!.....Initial temperature od 1d hts
!
      ALLOCATE(t_1dg(ng_hts_1d,nz0_2d,nr_1d))
!      
      DO g=1,ng_hts_1d
         DO j=1,nz0_2d
            READ(51,*) (t_1dg(g,j,i),i=1,nr_1d)
         ENDDO   
      ENDDO
!      
      DO i=1,ncell_hts_1d_all
         g=ig_hts_1d_tmp(i)
         j=inz0_2d_tmp(i)
         t_1d_tmp(i,:)=t_1dg(g,j,:)
      ENDDO  
!      
      DEALLOCATE(t_1dg,inz0_2d_tmp)
!
!.....Domain decomposition of 1D hts cells
!
      CALL rv_subdomain_info_hts_1d(ncell_hts_1d_all)
!      
      ALLOCATE(num_rod_1d(ncell_hts_1d))
      num_rod_1d(:)=0      
!
      CALL rv_allocate_hts_1d
!
!.....Assign local values from global values for parallel calculation
!.....For serial calculation, local values=global values
!
      DO ii=1,ncell_hts_1d
         i=jperm_hts_1d(ii)
         t_hts_1d(ii,:)=t_1d_tmp(i,:)
         ig_hts_1d(ii)=ig_hts_1d_tmp(i)
         num_rod_1d(ii)=num_rod_1d_tmp(i)
      ENDDO
!
!.....Calculate geometric input
!
      CALL rv_geo_hts_1d
!
!.....Find CUPID cell index connected to each 1D hts cells
!
      DO i=1,ncell_hts_1d
         DO ii=1,ncell_fluid
            IF(cupid_cell_1d_tmp(jperm_hts_1d(i)).eq.jperm(ii)) cupid_cell_1d(i)=ii
         ENDDO         
         IF(cupid_cell_1d(i).eq.0) STOP '### Fluid cell is not matched for a 1d rod ###'
      ENDDO
!
      DEALLOCATE(cupid_cell_1d_tmp,ig_hts_1d_tmp,nf_input,t_1d_tmp)
      DEALLOCATE(num_rod_1d_tmp)
!
      CLOSE(51)
!
      RETURN
      END SUBROUTINE rv_calc_geo_1d
