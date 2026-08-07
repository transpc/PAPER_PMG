!
      SUBROUTINE calc_geo_solid(nmaterial_tmp,vol_tmp,poro_tmp,xloc_tmp, &
                                nzone_tmp)
!
!.....This routine compute solid data, proceeds with mi decomposition
!.....and distributes global variables across mpi tasks then
!.....then maps against fluid cells
!
      USE Zinterface
      USE Zmpi        , ONLY: ncell_ps,celem,celem_c,jperm,maxmt_ncond, &
                              maxmt_nncond,                             &
                              maxmt_fluid
      USE Zzone       , ONLY: ncell_fluid,nmaterial,nmaterial_c,ncell_cond_all,ncell_cond,nzone_c
      USE Zcore       , ONLY: np,myrank
      USE Zparam      , ONLY: nn,ns,ndim
      USE Zmodel      , ONLY: use_porous
      USE Znum_cell   , ONLY: neigh_c,n_fluid,                                     &
                              i_neigh_tmp,j_neigh_tmp,j_nbcon_tmp,index_sort_tmp1, &
                              i_neigh_c_tmp,j_neigh_c_tmp,j_nbcon_c_tmp,           &
                              i_neigh,neigh,i_neigh_c,                             &
                              perm_tmp1,sv_tmp1,xfc_tmp1,                          &
                              index_sort_c
      USE Zbc_index   , ONLY: nbcon_c,nbcon
      USE Znormal     , ONLY: nji_c
      USE Zcoord1     , ONLY: xloc_c,xloc
      USE Zcoord2     , ONLY: fac_c,fac1_c,xfc,fac,fac1,xfc_c
      USE Zcoord3     , ONLY: volp_c
      USE Zcoord4     , ONLY: sap_c,dji_a_c,dji,dji_x,sap,dji_a
!      USE Zuserdefined, ONLY: udfl_perm_c_user
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: nmaterial_tmp(nn),nzone_tmp(nn)
      REAL(8) :: vol_tmp(nn),poro_tmp(nn)
      REAL(8) :: xloc_tmp(nn,ndim)
!.....Local variables
      INTEGER :: izone=1
      INTEGER :: i,j,k,k0,k1,m,ix,ii,jj,j0,j1,j2
      INTEGER :: ncells
      INTEGER :: ncellf0,ncells0
      INTEGER :: nb
      REAL(8) :: dr,dr0,dr1,dr2,dr3,dr4,xloc_k
      REAL(8) :: eps
!.....Local arrays
      REAL(8) :: tmp(ns)
      INTEGER :: i_neigh_fsw(ncell_fluid+1),j_neigh_fsw(maxmt_fluid)
      INTEGER,DIMENSION(:),ALLOCATABLE :: n_fluid1
      INTEGER,DIMENSION(:),ALLOCATABLE :: ns_o,ns_n,n_fluid_tmp
!.....Local domain arrays
      REAL(8) :: dji_x_c(ndim)
      REAL(8),DIMENSION(:),ALLOCATABLE :: perm_c
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: xn_c,sv_c     !,xfc_c   !ST2-CT-01
!.....Global array 
      INTEGER,DIMENSION(:),ALLOCATABLE :: nzone_c_tmp
      INTEGER,DIMENSION(:),ALLOCATABLE :: nmaterial_c_tmp
      REAL(8),DIMENSION(:),ALLOCATABLE :: vol_c_tmp,poro_c_tmp
!.....Global array ndim
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: xloc_c_tmp
!.....Global arrays csr
      REAL(8),DIMENSION(:),ALLOCATABLE :: perm_c_tmp1
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: sv_c_tmp1,xfc_c_tmp1
!.....Comm buffers 
      REAL(8),DIMENSION(:),ALLOCATABLE :: buff_comm_r
!
      ncells=ncell_cond_all
      use_porous=0
!
      IF(myrank.eq.0) THEN
         ALLOCATE(celem_c(ncells))
         ALLOCATE(ns_o(ncells),ns_n(nn),n_fluid_tmp(ncells))
         n_fluid_tmp(:)=0
         ncellf0=0
         ncells0=0
         DO i=1,nn
            IF(nmaterial_tmp(i).le.0)THEN
               ncellf0=ncellf0+1
               IF(nmaterial_tmp(i).lt.0)THEN
                  ncells0=ncells0+1
                  ns_o(ncells0)=i
                  ns_n(i)=ncells0
                  n_fluid_tmp(ncells0)=ncellf0
               ENDIF
            ELSE
               ncells0=ncells0+1
               ns_o(ncells0)=i
               ns_n(i)=ncells0
            ENDIF
         ENDDO
      ELSE
         ALLOCATE(n_fluid_tmp(1))
         ALLOCATE(ns_n(1))
      ENDIF
!  
      IF(myrank.eq.0) THEN
         ALLOCATE(i_neigh_c_tmp(ncell_cond_all+1))
         i_neigh_c_tmp(1)=1
         DO i=1,ncell_cond_all
            k=ns_o(i)
            i_neigh_c_tmp(i+1)=i_neigh_c_tmp(i)+(i_neigh_tmp(k+1)-i_neigh_tmp(k))
         ENDDO
         nb=i_neigh_c_tmp(ncell_cond_all+1)-1
      ENDIF
      IF(np.gt.1) CALL broadcast_i1(nb) 
      maxmt_nncond=nb
!
      IF(myrank.eq.0) THEN
         ALLOCATE(j_neigh_c_tmp(nb),j_nbcon_c_tmp(nb))
         DO i=1,ncell_cond_all 
            k=ns_o(i)
            celem_c(i)=celem(k)
            j0=i_neigh_tmp(k)-1
            j1=i_neigh_c_tmp(i)-1
            DO j=i_neigh_tmp(k),i_neigh_tmp(k+1)-1
               IF(j_nbcon_tmp(j).eq.0)THEN
                  m=j_neigh_tmp(j)
                  IF(nmaterial_tmp(m).ne.0)THEN
                     j_neigh_c_tmp(j-j0+j1)=ns_n(m)
                     j_nbcon_c_tmp(j-j0+j1)=j_nbcon_tmp(j)
                  ELSE
                     j_neigh_c_tmp(j-j0+j1)=0
                     j_nbcon_c_tmp(j-j0+j1)=-1
                  ENDIF
               ELSE
                  j_neigh_c_tmp(j-j0+j1)=0
                  j_nbcon_c_tmp(j-j0+j1)=j_nbcon_tmp(j)
               ENDIF
            ENDDO
         ENDDO
      ENDIF
!
      IF(np.gt.1) THEN
         CALL subdomain_info_solid(ncell_cond_all)
      ELSE
         CALL subdomain_info_solid_ser(ncell_cond_all)
      ENDIF
!
      DEALLOCATE(j_neigh_c_tmp,j_nbcon_c_tmp)
!
      eps=1.0d-8
!
      IF(myrank.eq.0) THEN
         ALLOCATE(nmaterial_c_tmp(ncells),nzone_c_tmp(ncells))
         ALLOCATE(vol_c_tmp(ncells),poro_c_tmp(ncells))
         ALLOCATE(xloc_c_tmp(ncells,ndim))
!
!........Collect global xloc_c,sv_c,xfc_c from csr global arrays
!........Convert to csr form immediatly
!
         ALLOCATE(perm_c_tmp1(nb))
         ALLOCATE(sv_c_tmp1(nb,ndim))
         ALLOCATE(xfc_c_tmp1(nb,ndim))
         DO i=1,ncell_cond_all 
            k=ns_o(i)
            j0=i_neigh_c_tmp(i)-1
            j1=i_neigh_tmp(k)-1
            IF(n_fluid_tmp(i).eq.0) then
               DO j=i_neigh_tmp(k),i_neigh_tmp(k+1)-1
                  IF(nmaterial_tmp(k).lt.0) THEN
                     perm_c_tmp1(j0+j-j1)=1.0d0-perm_tmp1(j)
                  ELSE
                     perm_c_tmp1(j0+j-j1)=perm_tmp1(j)
                  ENDIF
               ENDDO
            ELSE
               DO j=i_neigh_tmp(k),i_neigh_tmp(k+1)-1
                  k1=index_sort_tmp1(j)
                  IF(nmaterial_tmp(k).lt.0) THEN
                     perm_c_tmp1(j0+j-j1)=1.0d0-perm_tmp1(k1+j1)
                  ELSE
                     perm_c_tmp1(j0+j-j1)=perm_tmp1(k1+j1)
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
!
         DO i=1,ncell_cond_all 
            k=ns_o(i)
            nmaterial_c_tmp(i)=nmaterial_tmp(k)
            nzone_c_tmp(i)=nzone_tmp(k)
            vol_c_tmp(i)=vol_tmp(k)
            IF(nmaterial_tmp(k).lt.0) THEN
               poro_c_tmp(i)=1.0d0-poro_tmp(k)
            ELSE
               poro_c_tmp(i)=poro_tmp(k)
            ENDIF
         ENDDO
         DO ix=1,ndim
            DO i=1,ncell_cond_all 
               k=ns_o(i)
               xloc_c_tmp(i,ix)=xloc_tmp(k,ix)
            ENDDO
         ENDDO
         DO ix=1,ndim
            DO i=1,ncell_cond_all 
               k=ns_o(i)
               j0=i_neigh_c_tmp(i)-1
               j1=i_neigh_tmp(k)-1
               IF(n_fluid_tmp(i).eq.0) then
                  DO j=i_neigh_tmp(k),i_neigh_tmp(k+1)-1
                     sv_c_tmp1(j0+j-j1,ix)=sv_tmp1(j,ix)
                     xfc_c_tmp1(j0+j-j1,ix)=xfc_tmp1(j,ix)
                  ENDDO
               ELSE
                  DO j=i_neigh_tmp(k),i_neigh_tmp(k+1)-1
                     k1=index_sort_tmp1(j)
                     sv_c_tmp1(j0+j-j1,ix)=sv_tmp1(j1+k1,ix)
                     xfc_c_tmp1(j0+j-j1,ix)=xfc_tmp1(j1+k1,ix)
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ELSE
         ALLOCATE(nmaterial_c_tmp(1),nzone_c_tmp(1))
         ALLOCATE(perm_c_tmp1(1))
         ALLOCATE(sv_c_tmp1(1,ndim))
         ALLOCATE(xfc_c_tmp1(1,ndim))
         ALLOCATE(xloc_c_tmp(1,ndim))
      ENDIF
!
!.....read input file for radiation model between zones  
!  
      IF(myrank.eq.0)THEN
         CALL radiation_component_init(ns_n,nn,nzone_tmp,xloc_tmp)
      ELSE
         CALL radiation_component_init(ns_n,1,nzone_tmp,xloc_tmp)
      ENDIF   
!      
      IF(myrank.eq.0) THEN
         DEALLOCATE(ns_o,ns_n)
      ENDIF
!
!.....Allocate global variables 
!
      CALL allocate_var_solid(ncell_cond,ncell_ps)
!
!.....Allocate local variables
!
      ALLOCATE(n_fluid1(ncell_cond))
      if(maxmt_ncond.gt.0) then
         ALLOCATE(xfc_c(maxmt_ncond,ndim))
         ALLOCATE(sv_c(maxmt_ncond,ndim),xn_c(maxmt_ncond,ndim),perm_c(maxmt_ncond))
      ELSE
!     when maxmt_ncond=0 on some procs avoid troubles with check bounds
         ALLOCATE(xfc_c(1,ndim))
         ALLOCATE(sv_c(1,ndim),xn_c(1,ndim),perm_c(1))
      ENDIF
!
      xn_c(:,:)=0.0d0
      perm_c(:)=0.0d0
      sv_c(:,:)=0.0d0
      xfc_c(:,:)=0.0d0
!
!.....Save global variables into local vaiables
!
      CALL scatterv_i(nzone_c    ,nzone_c_tmp    ,ncell_cond,ncell_cond_all,izone)
      CALL scatterv_i(nmaterial_c,nmaterial_c_tmp,ncell_cond,ncell_cond_all,izone)
      CALL scatterv_i(n_fluid1   ,n_fluid_tmp    ,ncell_cond,ncell_cond_all,izone)
!
      IF(myrank.eq.0) THEN
         ALLOCATE(buff_comm_r(ncell_cond_all))
         DO i=1,ncell_cond_all
            buff_comm_r(i)=vol_c_tmp(i)*poro_c_tmp(i)
         ENDDO
      ELSE
         ALLOCATE(buff_comm_r(1))
      ENDIF
      CALL  scatterv_r(volp_c,buff_comm_r,ncell_cond,ncell_cond_all,izone)
      DEALLOCATE(buff_comm_r)
!
      CALL scatterv_ndim_fp_r(xloc_c,xloc_c_tmp,ncell_cond,ncell_ps,ncells,izone)
!
      CALL scatterv_csr_r(perm_c,maxmt_ncond,perm_c_tmp1,nb,ncell_cond_all,nb,i_neigh_c_tmp,izone)
!
!.....Get sorted permeability
!
      DO i=1,ncell_cond
         j0=i_neigh_c(i)-1
         DO j=i_neigh_c(i),i_neigh_c(i+1)-1
            k=index_sort_c(j)
            tmp(j-j0)=perm_c(k+j0)
         ENDDO
         DO j=i_neigh_c(i),i_neigh_c(i+1)-1
            perm_c(j)=tmp(j-j0)
         ENDDO
      ENDDO
!
      DO ix=1,ndim
         CALL scatterv_csr_r(sv_c(1,ix),maxmt_ncond,sv_c_tmp1(1,ix),nb,ncell_cond_all,nb,i_neigh_c_tmp,izone)
         DO i=1,ncell_cond
            j0=i_neigh_c(i)-1
            DO j=i_neigh_c(i),i_neigh_c(i+1)-1
               k=index_sort_c(j)
               tmp(j-j0)=sv_c(k+j0,ix)
            ENDDO
            DO j=i_neigh_c(i),i_neigh_c(i+1)-1
               sv_c(j,ix)=tmp(j-j0)
            ENDDO
         ENDDO
      ENDDO
!
      DO ix=1,ndim
         CALL scatterv_csr_r(xfc_c(1,ix),maxmt_ncond,xfc_c_tmp1(1,ix),nb,ncell_cond_all,nb,i_neigh_c_tmp,izone)
         DO i=1,ncell_cond
            j0=i_neigh_c(i)-1
            DO j=i_neigh_c(i),i_neigh_c(i+1)-1
               k=index_sort_c(j)
               tmp(j-j0)=xfc_c(k+j0,ix)
            ENDDO
            DO j=i_neigh_c(i),i_neigh_c(i+1)-1
               xfc_c(j,ix)=tmp(j-j0)
            ENDDO
         ENDDO
      ENDDO
!
      IF(np.gt.1) CALL communicate_2d_c(xloc_c)
!
      IF(myrank.eq.0) THEN
         DEALLOCATE(n_fluid_tmp)
         DEALLOCATE(nzone_c_tmp,nmaterial_c_tmp)
         DEALLOCATE(vol_c_tmp,poro_c_tmp)
         DEALLOCATE(xloc_c_tmp)
         DEALLOCATE(perm_c_tmp1)
      ENDIF
      DEALLOCATE(sv_c_tmp1,xfc_c_tmp1)
      DEALLOCATE(i_neigh_c_tmp)
!
      DO i=1,ncell_cond
!
!........Primary geo data
!
         j0=i_neigh_c(i)-1
         DO j=i_neigh_c(i),i_neigh_c(i+1)-1
            sap_c(j)=0.0d0
            DO ix=1,ndim
               sap_c(j)=sap_c(j)+sv_c(j,ix)**2
            ENDDO
            sap_c(j)=SQRT(sap_c(j))
            DO ix=1,ndim
               xn_c(j,ix)=sv_c(j,ix)/sap_c(j)
            ENDDO
            dr1=0.0d0
            dr2=0.0d0
            DO ix=1,ndim
               IF(nbcon_c(j).eq.0)THEN
                  k=neigh_c(j)
                  xloc_k=xloc_c(k,ix)
               ELSE
                  xloc_k=xfc_c(j,ix)
               ENDIF
               dji_x_c(ix)=xloc_k-xloc_c(i,ix)
               dr1=dr1+(xfc_c(j,ix)-xloc_c(i,ix))**2
               dr2=dr2+(xfc_c(j,ix)-xloc_k)**2
            ENDDO
            dr1=SQRT(dr1)
            dr2=SQRT(dr2)
            fac1_c(j)=dr2/(dr1+dr2)
            fac_c(j) =dr1/(dr1+dr2)
!           fac_c(j)=1.d0-fac1_c(j)
!
!...........Secondary variables
!
            dr=0.0d0
            DO ix=1,ndim
               dr=dr+dji_x_c(ix)*xn_c(j,ix)
            ENDDO
            dji_a_c(j)=dr
            sap_c(j)=sap_c(j)*perm_c(j)/dji_a_c(j)
!
!...........Modified cell center coordinates for non-orthogonal grid
!
!           dr=0.0d0
!           DO ix=1,ndim
!             dr=dr+(xfc_c(j,i,ix)-xloc_c(i,ix))*xn_c(j,i,ix)
!           ENDDO
!   xloc_m_c never used
!           DO ix=1,ndim
!              xloc_m_c(j,i,ix)=xfc_c(j,i,ix)-dr*xn_c(j,i,ix)-xloc_c(i,ix)
!           ENDDO
!
         ENDDO
!
      ENDDO
!
!.....Find fluid cells which have interfaces with solid cells
!
!.....Get all the fluid fsw
      i_neigh_fsw(1)=1
      DO i=1,ncell_fluid
         j0=0
         j1=i_neigh_fsw(i)-1
         j2=i_neigh(i)-1
         DO j=i_neigh(i),i_neigh(i+1)-1
            IF(nbcon(j).eq.-2)THEN
               j0=j0+1
               j_neigh_fsw(j0+j1)=j
            ENDIF
         ENDDO
         i_neigh_fsw(i+1)=i_neigh_fsw(i)+j0
      ENDDO
!
      DO i=1,ncell_cond
!
!........Fluid cells having porous solid
!
         IF(nmaterial_c(i).lt.0)THEN
            use_porous=1
            k=0
            DO ii=1,ncell_fluid
               IF(nmaterial(ii).lt.0)THEN
                  IF(n_fluid1(i).eq.jperm(ii))THEN
                     n_fluid(i)=ii
                     k=1
                  ENDIF
               ENDIF
            ENDDO
            IF(k.eq.0)THEN
               STOP '### Fluid cell has not match with a porous solid cell ###'
            ENDIF
         ENDIF
!
!........Case solid cell is the same as fluid cell avoid dr0=0.d0 divide
!
         IF(nmaterial_c(i).lt.0) CYCLE
!
         DO j=i_neigh_c(i),i_neigh_c(i+1)-1
!
!...........Fluid-solid interface cells
!
            IF(nbcon_c(j).eq.-2)THEN
               neigh_c(j)=0
               DO ii=1,ncell_fluid
                  DO j0=i_neigh_fsw(ii),i_neigh_fsw(ii+1)-1
                     jj=j_neigh_fsw(j0)
                     dr=0.0d0
                     DO ix=1,ndim
                        dr=dr+(xfc(jj,ix)-xfc_c(j,ix))**2
                     ENDDO
                     IF(dr.lt.eps)THEN
                        neigh_c(j)=ii
                        neigh(jj)=i
!
                        dr1=0.0d0
                        dr2=0.0d0
                        dr3=0.0d0
                        dr4=0.0d0
                        sap_c(j)=sap_c(j)*dji_a_c(j)
!bug                    sap(jj,ii)=sap_c(jj,ii)*dji_a_c(jj,ii)
                        DO ix=1,ndim
                           dr0=xloc(ii,ix)-xloc_c(i,ix)
                           dr4=dr4+dr0**2
                           dr1=dr1+(xfc_c(j,ix)-xloc_c(i,ix))**2
                           dr2=dr2+(xfc_c(j,ix)-xloc(ii,ix))**2
                           dr3=dr3+dr0*xn_c(j,ix)
                           dji_x(jj,ix)=-dr0
                        ENDDO
                        dji_a_c(j)=dr3 
                        dr1=SQRT(dr1)
                        dr2=SQRT(dr2)
                        fac1_c(j)=dr2/(dr1+dr2)
                        fac_c(j) =dr1/(dr1+dr2)
!                       fac_c(j)=1.d0-fac1_c(j)
                        sap_c(j)=sap_c(j)/dji_a_c(j)
                        dji(jj)=SQRT(dr4)
                        dji_a(jj)=dji_a_c(j)
                        fac1(jj)=fac_c(j)
                        fac(jj) =fac1_c(j)
                        sap(jj)=sap_c(j)
                        GOTO 100
                     ENDIF
                  ENDDO
               ENDDO
100            CONTINUE
               IF(neigh_c(j).eq.0)THEN
                  PAUSE '### Solid-fluid Interface is not found ###'
                  STOP '### Solid-fluid Interface is not found ###'
               ENDIF
            ENDIF
!
         ENDDO
      ENDDO
!
      DEALLOCATE(n_fluid1)
!
!.....Make fac_c and fac1_c symmetric
!
!
      DO i=1,ncell_cond
         DO j=i_neigh_c(i),i_neigh_c(i+1)-1
            IF(nbcon_c(j).eq.0)THEN
               k=neigh_c(j)
               IF(k.lt.i)THEN
                  k0=i_neigh_c(k)-1
                  fac_c(nji_c(j)+k0)=fac1_c(j)
                  fac1_c(nji_c(j)+k0)=fac_c(j)
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      IF(np.gt.1) THEN
         CALL communicate_1d_c_csr(fac_c,i_neigh_c)
         CALL communicate_1d_c_csr(fac1_c,i_neigh_c)
      ENDIF
!
      DEALLOCATE(celem_c)
      DEALLOCATE(xn_c)
      DEALLOCATE(perm_c,sv_c)     !,xfc_c   !!ST2-CT-01
!
      END SUBROUTINE calc_geo_solid
