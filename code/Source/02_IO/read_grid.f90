      SUBROUTINE read_grid
!
!.....This routine reads geometry data from 'somaGrid.in'
!
      USE Zinterface
      USE Zmpi         , ONLY: ncell_fp,maxmt,maxmt_pad,                   &
                               ncell_ps,maxmt_nfluid,maxmt_cell,max_neigh, &
                               au,ia_a,ja_a,ju_a,iend,                     &
                               ap,iap,jap,jaa,jaar,iaa,nbgroup,ngroup,     &
                               celem,metis
      USE Zzone        , ONLY: num_fzone,num_pzone,num_szone,                  &
                               ncell_fluid_all,ncell_fluid,ncell_fluid_pad,ncell_fluid_padv, &
                               ncell_cond_all,ncell_cond 
      USE Zcore        , ONLY: cupid_mars,np,myrank
      USE Zparam       , ONLY: nn,nin_max,nb_max,ndim,mesh_openfoam,nb_mars
      USE Zbc_index    , ONLY: nbcon,cellvin,cellpin,iface_wall1
      USE Zconst1      , ONLY: cplmaster,iturb,Nwlf,iheatpart,rv_htmodel_forCFD
      USE Zcoord1      , ONLY: xloc,xloc_tmp
      USE Zcoord3      , ONLY: porosity
      USE Zcoord4      , ONLY: sa
      USE Znum_cell    , ONLY: ncell,neigh,num_neigh_tmp, &
                               i_neigh_tmp,j_neigh_tmp,   &
                               i_neigh, & 
                               perm_tmp1,sv_tmp1,xfc_tmp1
      USE Znormal      , ONLY: sa_walll
      USE Zgradoption  , ONLY: ifrink
      USE Zmars        , ONLY: n_marsbc
      USE Znode        , ONLY: nd_max,nmax_vertex,neigh_face_tmp1,nd,node_face, &
                               i_cell_node_tmp,num_cell_node_tmp,cell_node_tmp, &
                               xnode,rwcn_tmp
      USE Zuserdefined , ONLY: udfl_grid_user,udfl_porous_user,MG_solver
      USE Zporous      , ONLY: l_subchannel
      USE Zio_unit     , ONLY: unit_log
      USE Zmpi         , ONLY: jperm
      USE Zcoord3      , ONLY: vol
!
      USE MD_MPI       , ONLY: myrankt
      USE omp_lib
      USE MD_OpenMP
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,ii,ix,j,j0,jj,k,kk,celem_max
      INTEGER :: n_face
      INTEGER :: nb_cell_node
      INTEGER :: wUnit = 1001, colorin
      LOGICAL :: check
      REAL(8) :: permi
!.....Local arrays
      INTEGER :: itmp(4)
      INTEGER, ALLOCATABLE::npb_tmp(:)
!.....Global arrays
      INTEGER, ALLOCATABLE::nzone_tmp(:),nmaterial_tmp(:)
      INTEGER, ALLOCATABLE::icore_tmp(:)
!     INTEGER, ALLOCATABLE::ncell_zone(:)
      REAL(8), ALLOCATABLE::vol_tmp(:),poro_tmp(:)
      REAL(8), ALLOCATABLE::hydraulicd_tmp(:),sl_tmp(:,:),sgap_tmp(:,:)    !LSJ 161122 porous      
!     mesh_openfoam.eq.1.and.ifrink.ge.1
      INTEGER, ALLOCATABLE::nd_tmp1(:),cell_node_tmp1(:)
      REAL(8), ALLOCATABLE::rwcn_tmp1(:),dxr_tmp1(:,:)
!
! PMG
      IF (MG_solver) THEN
          CALL read_input_mg
      END IF
      
     i= OMP_get_max_threads()
     IF(myrank.EQ.0) THEN
     write(999,*)'maximum threads',i
     ENDIF

! set OMP here
     CALL omp_set_num_threads(nthre)
     
! = = = = = = = = = = 
!
      IF(myrank.eq.0) THEN
         ALLOCATE(celem(nn))
         ALLOCATE(num_neigh_tmp(nn))
         ALLOCATE(xloc_tmp(nn,ndim))
!........Local alloc
         ALLOCATE(npb_tmp(nn))
!
         ALLOCATE(nzone_tmp(nn),nmaterial_tmp(nn))
         ALLOCATE(icore_tmp(nn))
         ALLOCATE(vol_tmp(nn),poro_tmp(nn))
         ALLOCATE(hydraulicd_tmp(nn),sl_tmp(nn,ndim),sgap_tmp(nn,ndim))             !LSJ 161122 porous
      ELSE
         ALLOCATE(xloc_tmp(1,ndim))
!........Local alloc
         ALLOCATE(npb_tmp(1))
!
         ALLOCATE(nzone_tmp(1),nmaterial_tmp(1))
         ALLOCATE(icore_tmp(1))
         ALLOCATE(vol_tmp(1),poro_tmp(1))
         ALLOCATE(hydraulicd_tmp(1),sl_tmp(1,ndim),sgap_tmp(1,ndim))             !LSJ 161122 porous
      ENDIF
!
!     ncell_zone(20))
!     ncell_zone(:)=0
!.....Initialize local variables 
!
      IF(myrank.eq.0) THEN
         icore_tmp(:)=0
         hydraulicd_tmp(:)=0.d0      !LSJ 161122 porous
         sl_tmp(:,:)=0.d0            !LSJ 161122 porous 
         sgap_tmp(:,:)=0.d0            !LSJ 161122 porous 
      ENDIF
!
!.....Set ncell to value of nn
!
      ncell=nn 
      IF(myrank.eq.0) ALLOCATE(i_neigh_tmp(ncell+1))
!      if ( MG_solver ) CALL read_input_mg
      if ( MG_solver == .true. .and. np == 1 ) ifrink = 1
!            
!.....Read cell data only myrank=0 do work
!
      IF(myrank.eq.0) THEN
         CALL read_data(npb_tmp,nzone_tmp,nmaterial_tmp,celem,vol_tmp,poro_tmp, &
                        n_face,nb_cell_node)
         IF(.not.(mesh_openfoam.eq.1.and.ifrink.ge.1)) ALLOCATE(i_cell_node_tmp(1))
      ELSE
         ALLOCATE(perm_tmp1(1))
      ENDIF ! myrank
!
!.....Broadcast needed variables and add ALLOCATE for myrank # 0
!
      IF(np.gt.1) THEN
         IF(myrank.eq.0) THEN
            itmp(1)=ncell_fluid_all
            itmp(2)=ncell_cond_all
            itmp(3)=maxmt_nfluid
            itmp(4)=maxmt_cell
         ENDIF
         CALL broadcast_i(itmp,4)
         ncell_fluid_all=itmp(1)
         ncell_cond_all =itmp(2)
         maxmt_nfluid   =itmp(3)
         maxmt_cell     =itmp(4)
         IF(mesh_openfoam.eq.1.and.ifrink.ge.1) CALL broadcast_i1(nb_cell_node)
      ENDIF
      IF(myrank.ne.0) THEN
         ALLOCATE(sv_tmp1(1,ndim),xfc_tmp1(1,ndim))
         ALLOCATE(neigh_face_tmp1(1))
         IF(mesh_openfoam.eq.1.and.ifrink.ge.1) THEN
            ALLOCATE(i_cell_node_tmp(nn+1))
         ELSE 
            ALLOCATE(i_cell_node_tmp(1))
         ENDIF
      ELSE
         IF(mesh_openfoam.ne.1) ALLOCATE(neigh_face_tmp1(1))
      ENDIF
!
!.....Modify the geometric data for specific problem
!
      IF(udfl_grid_user) CALL udfn_grid_user(ncell,xloc_tmp,npb_tmp)
!      IF(myrank.eq.0) THEN
!         IF(udfl_permeability) CALL udfn_permeability
!      ENDIF
      IF(udfl_porous_user) CALL udfn_porous_user(ncell,                                                   &
                                                 xloc_tmp,vol_tmp,poro_tmp,                               &
                                                 nmaterial_tmp,nzone_tmp,sl_tmp,hydraulicd_tmp,icore_tmp)

!
!.....Subchannel geometry setting
!     For subchannel analysis, geometry settings are defined here 
!      - channel info, porosity, Hydradulic diameter, gap distance, etc.
      
      IF(l_subchannel)then
         CALL udfn_subchannel_geo(ncell,vol_tmp,poro_tmp,hydraulicd_tmp,sl_tmp,sgap_tmp)
      ENDIF    
!
!.....Consistency of two permeabilities at common face of two cells 
!.....whose nzones are different from each other
!
!
      IF(myrank.eq.0) THEN
         DO i=1,ncell
            DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
               k=j_neigh_tmp(j)
               IF(k.gt.0) THEN
                  IF(nzone_tmp(i).ne.nzone_tmp(k))THEN
                     DO jj=i_neigh_tmp(k),i_neigh_tmp(k+1)-1
                        kk=j_neigh_tmp(jj)
                        IF(kk.eq.i) THEN
                           IF(perm_tmp1(j).ne.perm_tmp1(jj))then
                              permi=MIN(perm_tmp1(j),perm_tmp1(jj))
                              perm_tmp1(j)=permi
                              perm_tmp1(jj)=permi
                           ENDIF
                           GOTO 100
                        ENDIF
                     ENDDO
100                  CONTINUE
                  ENDIF   
               ENDIF   
            ENDDO  
         ENDDO       
      ENDIF
!
!.....Call the function to divide sub-domains
!
      celem_max=0
      IF(metis.eq.2.and.np.gt.1) THEN
         IF(myrank.eq.0) THEN
            DO
               INQUIRE(unit=wUnit,opened=check)
               IF(check==.false.)EXIT
               wUnit=wUnit+1
            ENDDO
            OPEN(unit=wUnit,file='celem',status='old',iostat=colorin)
            IF(colorin.eq.0)THEN
               DO i=1,nn
                  READ(wUnit,*)celem(i)
               ENDDO
               CLOSE(wUnit)
               IF(myrank.eq.0)PRINT *,'          METIS domain decomposition was read.'
            ELSE
               IF(myrank.eq.0)PRINT *,'          celem file is needed when metis is 2.'
               PAUSE
               STOP
            ENDIF
            celem_max=maxval(celem(1:nn))
            WRITE(*,*) '          The number of subdomains is', celem_max
            WRITE(*,*) '          The number of processors is', np
            IF(np.ne.celem_max) THEN
               WRITE(*,*) '          The number of subdomains differs from the number of processors !!!'
               PAUSE
               STOP
            ENDIF
            WRITE(unit_log,*) '          The number of subdomains is', celem_max
            WRITE(unit_log,*) '          The number of processors is', np         
         ENDIF
      ENDIF
      IF(mesh_openfoam.eq.1.and.ifrink.ge.1)THEN
         IF(np.gt.1) CALL broadcast_i(i_cell_node_tmp,ncell_fluid_all+1)
      ENDIF
      IF(np.gt.1) THEN
         CALL subdomain_info(ncell_fluid_all,i_cell_node_tmp) 
      ELSE
         CALL subdomain_info_ser(ncell_fluid_all) 
      ENDIF
!
! PMG
      IF(MG_solver) THEN
         IF(myrankt .EQ. 0) THEN
           CALL subdomain_infor_MG
         END IF
      END IF
!
      IF(mesh_openfoam.eq.1.and.ifrink.ge.1)THEN
         IF(np.gt.1) THEN
            IF(myrank.eq.0) THEN
               itmp(1)=nd_max
               itmp(2)=nmax_vertex
            ENDIF
            CALL broadcast_i(itmp,2)
            nd_max=itmp(1)
            nmax_vertex=itmp(2)
         ENDIF
!
         IF(myrank.eq.0) THEN
            ALLOCATE(nd_tmp1(maxmt_nfluid))
            DO i=1,ncell_fluid_all
               j0=i_neigh_tmp(i)-1
               DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
                  ii=neigh_face_tmp1(j)
                  nd_tmp1(j)=nd(ii)
               ENDDO
            ENDDO
            ALLOCATE(cell_node_tmp1(nb_cell_node))
            ALLOCATE(rwcn_tmp1(nb_cell_node))
            DO i=1,ncell_fluid_all
               j0=i_cell_node_tmp(i)-1
               DO j=1,num_cell_node_tmp(i)
                  cell_node_tmp1(j0+j)=cell_node_tmp(j,i)
                  rwcn_tmp1(j0+j)=rwcn_tmp(j,i)
               ENDDO
            ENDDO
            ALLOCATE(dxr_tmp1(nb_cell_node,ndim))
            DO ix=1,ndim
               DO i=1,ncell_fluid_all
                  j0=i_cell_node_tmp(i)-1
                  DO j=1,num_cell_node_tmp(i)
                     k=cell_node_tmp(j,i)
                     dxr_tmp1(j0+j,ix)=xnode(k,ix)-xloc_tmp(i,ix)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            ALLOCATE(num_cell_node_tmp(1))
            ALLOCATE(cell_node_tmp1(1))
            ALLOCATE(rwcn_tmp1(1))
            ALLOCATE(dxr_tmp1(1,ndim))
            ALLOCATE(nd_tmp1(1))
         ENDIF
      ELSE
            ALLOCATE(cell_node_tmp1(1))
            ALLOCATE(rwcn_tmp1(1))
            ALLOCATE(dxr_tmp1(1,ndim))
            ALLOCATE(nd_tmp1(1))
      ENDIF
!  
!.....Count the number of zones
!
      IF(myrank.eq.0) THEN
         num_fzone=0
         num_pzone=0
         num_szone=0
!
!.....Only fluid zone      
!
         DO i=1,ncell
            IF(nmaterial_tmp(i).eq.0)num_fzone=MAX(num_fzone,nzone_tmp(i))  
         ENDDO
!
!.....Only porous zone         
!
         DO i=1,ncell
            IF(nmaterial_tmp(i).lt.0)num_pzone=MAX(num_pzone,nzone_tmp(i)-num_fzone) 
         ENDDO
!
!.....Only solid zone
!
         DO i=1,ncell
            num_szone=MAX(num_szone,nzone_tmp(i))  
         ENDDO
         num_szone=num_szone-num_fzone-num_pzone 
!
!.....The number of fluid zone and porous zone      
!
         num_fzone=num_fzone+num_pzone 
      ENDIF
      IF(np.gt.1) THEN
        IF(myrank.eq.0) THEN
         itmp(1)=num_fzone
         itmp(2)=num_pzone
         itmp(3)=num_szone
        ENDIF
         CALL broadcast_i(itmp,3)
         num_fzone=itmp(1)
         num_pzone=itmp(2)
         num_szone=itmp(3)
      ENDIF
!
!     DO j=1,num_fzone+num_rzone+num_pzone+num_szone
!        ncell_zone(j)=0
!     ENDDO
!     DO i=1,ncell
!        DO j=1,num_fzone+num_rzone+num_pzone+num_szone
!..............The number of cells of jth zone            
!            IF(nzone_tmp(i).eq.j) ncell_zone(j)=ncell_zone(j)+1
!        ENDDO
!     ENDDO
!
!.....Allocate preliminary cell info variables 
!
      CALL ALLOCATE_preliminary_var(ncell_fluid,ncell_fp)
!     CALL ALLOCATE_geo_var(ncell_fp)
!
      CALL scatter_data(npb_tmp,nzone_tmp,sl_tmp,sgap_tmp,nmaterial_tmp,icore_tmp,vol_tmp, &
                        poro_tmp,hydraulicd_tmp,                                  &
                        nd_tmp1,cell_node_tmp1,rwcn_tmp1,dxr_tmp1,                &
                        n_face,nb_cell_node)
!
!.....Allocate variables 
!
      CALL ALLOCATE_geo_var(ncell_fp,ncell_fluid)
!
!.....Communicate xloc / porosity - real(8)
!
      IF(np.gt.1) THEN
         CALL communicate_1d(porosity)
         CALL communicate_2d(xloc)
      ENDIF
!
!.....Calculate local variables
!
      CALL calc_geo
!
      ncell_cond=0
      ncell_ps=0
      IF(ncell_cond_all.gt.0) CALL calc_geo_solid(nmaterial_tmp,             &
                                                  vol_tmp,poro_tmp,xloc_tmp, &
                                                  nzone_tmp)
!     DEALLOCATE(i_neigh_tmp)
!     DEALLOCATE(j_neigh_tmp)
!     DEALLOCATE(j_nbcon_tmp)
!
      CALL rv_calc_geo_2d
      CALL rv_calc_geo_1d
!
!.....DeALLOCATE tmp variables
!
!     used in write_fieldview.f90
!     DEALLOCATE(celem)
      IF(myrank.eq.0) THEN
!     used in user_def_output.f90
!        DEALLOCATE(num_neigh_tmp)
         DEALLOCATE(npb_tmp)
         DEALLOCATE(nzone_tmp,nmaterial_tmp)
         DEALLOCATE(icore_tmp)
         DEALLOCATE(vol_tmp,poro_tmp)
         DEALLOCATE(hydraulicd_tmp,sl_tmp,sgap_tmp)
         IF(mesh_openfoam.eq.1)THEN
            DEALLOCATE(xnode)
         ENDIF
         IF(mesh_openfoam.eq.1.and.ifrink.ge.1)THEN
            DEALLOCATE(nd_tmp1,cell_node_tmp1)
            DEALLOCATE(num_cell_node_tmp)
            DEALLOCATE(cell_node_tmp)
            DEALLOCATE(rwcn_tmp1)
            DEALLOCATE(rwcn_tmp)
         ENDIF
      ENDIF
      IF(mesh_openfoam.eq.1.and.ifrink.ge.1)THEN
         DEALLOCATE(dxr_tmp1)
      ENDIF
      DEALLOCATE(perm_tmp1)
      DEALLOCATE(sv_tmp1,xfc_tmp1)
!     IF(ALLOCATED(neigh_tmp))     DEALLOCATE(neigh_tmp)
!     IF(ALLOCATED(xloc_tmp))      DEALLOCATE(xloc_tmp)
!     IF(ALLOCATED(ncell_zone))    DEALLOCATE(ncell_zone)
!
      CALL find_icelltype
!
!.....Allocate variables 
!
      CALL ALLOCATE_physical_var
      CALL ALLOCATE_model_var
      CALL ALLOCATE_user_var
!
!.....Calculate walln
!
      CALL calc_walln
!
!.....Find the most close wall and distance      
!      
      IF(iturb.gt.0.or.Nwlf.ne.-1.0d0)THEN
         CALL find_walls
      ENDIF
!
!..... Compute iface_wall1 needed in heat_partition
!      Assume cell i has only one face with nbcon -2,-3,-4,-5,-6 
!
      IF(iheatpart.gt.0.or.rv_htmodel_forCFD.gt.0) THEN
         DO i=1,ncell_fluid
            iface_wall1(i)=0
            j0=i_neigh(i)-1
            DO j=i_neigh(i),i_neigh(i+1)-1
               IF(nbcon(j).le.-2)THEN
                  iface_wall1(i)=j-j0
                  sa_walll(i)=sa(j)
                  EXIT
               ENDIF
            ENDDO
         ENDDO
      ENDIF      
!
!.....Save the cell number with wall, pressure and velocity boundary conditions
!
      !nwb(:)=0
      cellvin(:)=0
      cellpin(:)=0    
      DO i=1,ncell_fluid
         j0=i_neigh(i)-1
         DO j=i_neigh(i),i_neigh(i+1)-1
            !IF(nbcon(j).le.-1) nwb(i)=1
            IF(nbcon(j).ge.1.and.nbcon(j).le.nin_max) cellvin(nbcon(j))=i
            IF(nbcon(j).gt.nin_max.and.nbcon(j).le.nb_max) cellpin(nbcon(j)-nin_max)=i
         ENDDO
      ENDDO   
!
!     get the maxmt to ALLOCATE ja_a for compute elements only
!
      ALLOCATE(ia_a(ncell_fluid+1),ju_a(ncell_fluid),iend(ncell_fluid))
      ia_a(1)=1
      maxmt=0
      max_neigh=0
      DO i=1,ncell_fluid
         j0=i_neigh(i)-1 
!
!........Get neighbors count for compute elements only
!
         DO j=i_neigh(i),i_neigh(i+1)-1
            IF(nbcon(j).ne.0) exit
            maxmt=maxmt+1 
         ENDDO 
!
!........Add 1 for new diagonal entry
!
         maxmt=maxmt+1 
         ia_a(i+1)=maxmt+1
         max_neigh=max(max_neigh,ia_a(i+1)-ia_a(i))
      ENDDO 
!
      ALLOCATE(ja_a(maxmt))
      CALL csr_neigh(maxmt,ncell_fluid,ia_a,ja_a,ju_a,iend,i_neigh,neigh,nbcon)
!
!.....Process ilup 
!
      write(*,*)'test,pre ilu'

      CALL reorder_ilup

      write(*,*)'test,after ilu'
! Rearrange A in blocks of equal number of neighbors to vectorize
      CALL gener_vect_size(ncell_fluid,max_neigh,ia_a,ngroup)
!     write(*,*) 'ngroup',ncell_fluid,ngroup,max_neigh
      ALLOCATE(au(maxmt))
      ALLOCATE(iap(2,ngroup+1))
      ALLOCATE(iaa(2,ngroup+1),jaa(ncell_fluid),jaar(ncell_fluid))
      ALLOCATE(nbgroup(3,ngroup))
      CALL gener_vect_u(ncell_fluid,maxmt_pad,ncell_fluid_pad,ncell_fluid_padv,max_neigh, &
                        ia_a,iaa,iap,                                                     &
                        jaa,jaar,ngroup,nbgroup)
      ALLOCATE(jap(maxmt_pad))
      ALLOCATE(ap(maxmt_pad))
      CALL copy_ja_vector(ncell_fluid,maxmt,maxmt_pad, &
                          ja_a,ia_a,jap,iap,           &
                          jaa,iaa,ngroup,nbgroup)
!
!     DO i=1,ngroup
!          write(*,*) '--',i,nbgroup(1,i),nbgroup(2,i),nbgroup(3,i)
!     ENDDO
!
!.....Set boundary condition for MARS coupling
!
      n_marsbc=0
      DO i=1,ncell_fluid
         j0=i_neigh(i)-1
         DO j=i_neigh(i),i_neigh(i+1)-1
            IF(nbcon(j).ge.nb_mars)THEN
               n_marsbc=n_marsbc+1
               WRITE(*,"(a,1i5,1i8,3e17.10)")'nbcon,i,sa,len,vol=',nbcon(j),jperm(i),sa(j),vol(i)/sa(j),vol(i)
            ENDIF
         ENDDO
      ENDDO
      IF(np.gt.1) CALL allreducei_i1(n_marsbc)
!      
!.....Kill code if mismatche occurs in cupvols    
!  
      IF(cupid_mars.eq.1.and.n_marsbc.eq.0)then
         WRITE(*,*)'Grid file should have more than one cupvol!'
         WRITE(unit_log,*)'Grid file should have more than one cupvol!'
         PAUSE
         STOP
      ELSEIF(cupid_mars.eq.0.and.n_marsbc.ne.0)then
         WRITE(*,*)'Grid file should not have any cupvols!'
         WRITE(unit_log,*)'Grid file should not have any cupvols!'
         PAUSE
         STOP
      ENDIF
!      
      IF(l_subchannel .or. cplmaster.gt.0)then
         CALL udfn_subchannel_info
      ENDIF

!...Preprocessing for MG(Multi-Grid) pressure solver...!

         IF (MG_solver) THEN

            CALL read_mesh_MPI
      
            CALL Prep_fine_P
            
            IF (myrankt .eq. 0) THEN
               WRITE (*, *) 'pre-MG'
            END IF

            CALL Prep_MG_GarL
            IF (myrankt .eq. 0) THEN
               WRITE (*, *) 'pre-MG-finish'
            END IF

         END IF      

      IF(myrank.eq.0.and.mesh_openfoam.eq.1) THEN
         DEALLOCATE(nd, node_face)
      ENDIF

      END SUBROUTINE read_grid
!
      SUBROUTINE searchx_list(n,list,m,a,indexo)
!
      IMPLICIT NONE
!
!     input
      INTEGER n,m
      INTEGER list(n)
      INTEGER a(m)
!     output
      INTEGER indexo(m)
!     local variables
      INTEGER i,ip,ip1,ip2
      INTEGER x
!
      DO i=1,m
         ip1=1
         ip2=n
110      CONTINUE
         if(ip1.gt.ip2) goto 100
         x=a(i)
         ip=(ip1+ip2)/2
         IF(x.eq.list(ip)) then
            indexo(i)=ip
            ip1=ip+1
            GOTO 100
         ELSEIF(x.lt.list(ip)) then
            ip2=ip-1
            GOTO 110
         ELSEIF(x.gt.list(ip)) then
            ip1=ip+1
            GOTO 110
         ENDIF
100      CONTINUE
      ENDDO
!
      END SUBROUTINE searchx_list
