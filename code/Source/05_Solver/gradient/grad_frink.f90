!
      SUBROUTINE grad_frink(s,dsdx,ig)
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell            
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: np,myrank
      USE Zparam       , ONLY: ndim,mesh_openfoam
      USE Zvec_param   , ONLY: nf_nonk,nf_non,nf_mcc,nf_inl,nf_out,nf_adw,nf_fsw,nf_ctw,nf_chw,nf_sym,nf_totk
      USE Znum_cell    , ONLY: i_neigh,istart_nf, &
                               nf_number_nb,lens, &
                               right_nb_k,istart_nf,nf_number_id,istart_nfs
      USE Zvec_index   , ONLY: left_nf,right_non,jneigh_nf,kneigh_non
      USE Zcoord3      , ONLY: volr
      USE Zgradoption  , ONLY: ifrink
      USE Zconst2      , ONLY: grav,gfactor
      USE Znode        , ONLY: n_node,num_cell_node,cell_node,num_nd,node_face_cell,rwcn,swcn
      USE Zbc_index    , ONLY: ngrad
      USE Zvec_geo     , ONLY: sv_nf,dxfc_nf,          &
                               dxfc_non_k,             &
                               fac_non,fac1_non
! 
      IMPLICIT NONE
!
!.....Input
      INTEGER :: ig
      REAL(8),DIMENSION(ncell_fp) :: s
!.....Output
      REAL(8),DIMENSION(ncell_fp,ndim) :: dsdx
!.....Local variables
      INTEGER :: i,j,k,jd,k0,k1,i2
      INTEGER :: ii,jj,kk,jk
      INTEGER :: nv,nf_number,istart0,istart,len,i0,i1
      REAL(8) :: fie,fie_i
      REAL(8) :: p_i,p_k
!.....Local arrays      
      REAL(8),DIMENSION(n_node) :: sn
      REAL(8),DIMENSION(nf_totk,ndim) :: fie_nf
!
!.....Check correct environment
!
      IF(mesh_openfoam.ne.1 .or. ifrink.ne.1) THEN
         IF(myrank.eq.0) THEN
            write(*,*) 'grad_frink should not be called in this environment'
            write(*,*) 'mesh_openfoam=',mesh_openfoam
            write(*,*) 'ifrink       =',ifrink
         ENDIF
         CALL finalize_mpi
         STOP 
      ENDIF  
!
!.....Define boundary nodes and cells
!.....Interpolate node values using inverse distance weighting
!
      sn(:)=0.0d0
      DO i=1,ncell_fluid
         DO j=1,num_cell_node(i)
            k=cell_node(j,i)
            sn(k)=sn(k)+rwcn(j,i)*s(i)
         ENDDO
      ENDDO
      DO i=1,n_node
         sn(i)=sn(i)/swcn(i)
      ENDDO
!
!......Allreduce sn
!
      IF(np.gt.1) CALL allreducei_r(sn,n_node)
!
      nf_number_nb=8
      nf_number_id(-1)=-1
      nf_number_id(0)=0
      nf_number_id(1)=1
      nf_number_id(2)=2
      nf_number_id(3)=3
      nf_number_id(4)=4
      nf_number_id(5)=5
      nf_number_id(6)=6
      nf_number_id(7)=7
      nf_number_id(8)=8
      istart_nfs(-1)=0
      istart_nfs(0)=istart_nfs(-1)+nf_nonk
      istart_nfs(1)=istart_nfs(0)+nf_non
      istart_nfs(2)=istart_nfs(1)+nf_mcc
      istart_nfs(3)=istart_nfs(2)+nf_inl
      istart_nfs(4)=istart_nfs(3)+nf_out
      istart_nfs(5)=istart_nfs(4)+nf_adw
      istart_nfs(6)=istart_nfs(5)+nf_fsw
      istart_nfs(7)=istart_nfs(6)+nf_ctw
      istart_nfs(8)=istart_nfs(7)+nf_chw
      lens         =istart_nfs(8)+nf_sym
!
      if(ndim.eq.2) THEN
         nv=0
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(nv)
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len
            i0=istart0+i
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            i2=i_neigh(ii)-1
            jj=jneigh_nf(i1)
!
            IF(ngrad(ii).eq.0)THEN
               fie_i=0.d0
               DO jd=1,num_nd(jj+i2)
                  k=node_face_cell(jd,jj+i2)
                  fie_i=fie_i+sn(k)
               ENDDO
               fie=fie_i/DFLOAT(num_nd(jj+i2))
!
               IF(ig.eq.3)THEN
                  p_i=grav(1)*dxfc_nf(i1,1)  +grav(2)*dxfc_nf(i1,2)
                  p_k=grav(1)*dxfc_non_k(i,1)+grav(2)*dxfc_non_k(i,2)
                  fie=fie+0.5*(cell%rhom(ii)*p_i+cell%rhom(kk)*p_k)*gfactor(ii)
               ENDIF
               fie_nf(i0,1)=fie*sv_nf(i1,1)
               fie_nf(i0,2)=fie*sv_nf(i1,2)
            ELSE
               IF(ig.eq.2)THEN
                  p_i=s(ii)
                  p_k=s(kk)
                  fie=fac1_non(i)*p_i+fac_non(i)*p_k
               ELSE
                  p_i=grav(1)*dxfc_nf(i1,1)  +grav(2)*dxfc_nf(i1,2)
                  p_k=grav(1)*dxfc_non_k(i,1)+grav(2)*dxfc_non_k(i,2)
                  fie=(        fac1_non(i)*p_i+fac_non(i)   *p_k               &
                       +0.5*(cell%rhom(ii)*p_i+cell%rhom(kk)*p_k))*gfactor(ii)
               ENDIF
               fie_nf(i0,1)=fie*sv_nf(i1,1)
               fie_nf(i0,2)=fie*sv_nf(i1,2)
            ENDIF
         ENDDO
!
         nv=-1
         nf_number=nf_number_id(nv)
         len   =istart_nf(2,nf_number)
         DO i=1,len
            k=right_nb_k(i)
            ii=right_non(k)
            kk=left_nf(k)
            k0=i_neigh(ii)-1
            jk=kneigh_non(k)
!
            IF(ngrad(ii).eq.0)THEN
               fie_i=0.d0
               DO jd=1,num_nd(jk+k0)
                  k1=node_face_cell(jd,jk+k0)
                  fie_i=fie_i+sn(k1)
               ENDDO
               fie=fie_i/DFLOAT(num_nd(jk+k0))
!
               IF(ig.eq.3)THEN
                  p_i=grav(1)*dxfc_nf(k,1)   +grav(2)*dxfc_nf(k,2)
                  p_k=grav(1)*dxfc_non_k(k,1)+grav(2)*dxfc_non_k(k,2)
                  fie=fie+0.5*(cell%rhom(ii)*p_i+cell%rhom(kk)*p_k)*gfactor(ii)
               ENDIF
               fie_nf(i,1)=-fie*sv_nf(k,1)
               fie_nf(i,2)=-fie*sv_nf(k,2)
            ELSE
               IF(ig.eq.2)THEN
                  p_i=s(ii)
                  p_k=s(kk)
                  fie=fac1_non(k)*p_i+fac_non(k)*p_k
               ELSE
                  p_i=grav(1)*dxfc_nf(k,1)   +grav(2)*dxfc_nf(k,2)
                  p_k=grav(1)*dxfc_non_k(k,1)+grav(2)*dxfc_non_k(k,2)
                  fie=(        fac1_non(k)*p_i+fac_non(k)   *p_k               &
                       +0.5*(cell%rhom(ii)*p_i+cell%rhom(kk)*p_k))*gfactor(ii)
               ENDIF
               fie_nf(i,1)=-fie*sv_nf(k,1)
               fie_nf(i,2)=-fie*sv_nf(k,2)
            ENDIF
         ENDDO
!
         DO nv=1,8
            nf_number=nf_number_id(nv)
            istart0=istart_nfs(nv)
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            IF(nv.ge.4 .and. nv.le.7) THEN
               DO i=1,len
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  i2=i_neigh(ii)-1
                  jj=jneigh_nf(i1)
!
                  IF(ngrad(ii).eq.0)THEN
                     fie_i=0.d0
                     DO jd=1,num_nd(jj+i2)
                        k=node_face_cell(jd,jj+i2)
                        fie_i=fie_i+sn(k)
                     ENDDO
                     fie=fie_i/DFLOAT(num_nd(jj+i2))
!
                     IF(ig.eq.3)THEN
                        p_i=grav(1)*dxfc_nf(i1,1)+grav(2)*dxfc_nf(i1,2)
                        fie=fie+cell%rhom(ii)*p_i*gfactor(ii)
                     ENDIF
                     fie_nf(i0,1)=fie*sv_nf(i1,1)
                     fie_nf(i0,2)=fie*sv_nf(i1,2)
                  ENDIF
               ENDDO
            ELSE
               DO i=1,len
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
!
                  IF(ngrad(ii).eq.0)THEN
                     fie_nf(i0,1)=0.d0
                     fie_nf(i0,2)=0.d0
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
!
         DO nv=1,8
            nf_number=nf_number_id(nv)
            istart0=istart_nfs(nv)
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            IF(nv.eq.1 .or. nv.eq.3 .or. nv.eq.8) THEN
               DO i=1,len
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  kk=right_non(i)
                  i2=i_neigh(ii)-1
                  jj=jneigh_nf(i1)
!
                  IF(ngrad(ii).ne.0)THEN
                    fie=s(ii)
                    fie_nf(i0,1)=fie*sv_nf(i1,1)
                    fie_nf(i0,2)=fie*sv_nf(i1,2)
                  ENDIF
               ENDDO
            ELSE
               DO i=1,len
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  kk=right_non(i)
                  jj=jneigh_nf(i1)
!
                  IF(ngrad(ii).ne.0)THEN
                     IF(ig.eq.2) THEN
                       fie=s(ii)
                     ELSE
                        p_i=grav(1)*dxfc_nf(i1,1)+grav(2)*dxfc_nf(i1,2)
                        fie=s(ii)+cell%rhom(ii)*p_i*gfactor(ii)
                     ENDIF
                     fie_nf(i0,1)=fie*sv_nf(i1,1)
                     fie_nf(i0,2)=fie*sv_nf(i1,2)
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ELSE
         nv=0
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(nv)
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len
            i0=istart0+i
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            i2=i_neigh(ii)-1
            jj=jneigh_nf(i1)
!
            IF(ngrad(ii).eq.0)THEN
               fie_i=0.d0
               DO jd=1,num_nd(jj+i2)
                  k=node_face_cell(jd,jj+i2)
                  fie_i=fie_i+sn(k)
               ENDDO
               fie=fie_i/DFLOAT(num_nd(jj+i2))
!
               IF(ig.eq.3)THEN
                  p_i=grav(1)*dxfc_nf(i1,1)  +grav(2)*dxfc_nf(i1,2)  +grav(3)*dxfc_nf(i1,3)
                  p_k=grav(1)*dxfc_non_k(i,1)+grav(2)*dxfc_non_k(i,2)+grav(3)*dxfc_non_k(i,3)
                  fie=fie+0.5*(cell%rhom(ii)*p_i+cell%rhom(kk)*p_k)*gfactor(ii)
               ENDIF
               fie_nf(i0,1)=fie*sv_nf(i1,1)
               fie_nf(i0,2)=fie*sv_nf(i1,2)
               fie_nf(i0,3)=fie*sv_nf(i1,3)
            ELSE
               IF(ig.eq.2)THEN
                  p_i=s(ii)
                  p_k=s(kk)
                  fie=fac1_non(i)*p_i+fac_non(i)*p_k
               ELSE
                  p_i=grav(1)*dxfc_nf(i1,1)  +grav(2)*dxfc_nf(i1,2)  +grav(3)*dxfc_nf(i1,3)
                  p_k=grav(1)*dxfc_non_k(i,1)+grav(2)*dxfc_non_k(i,2)+grav(3)*dxfc_non_k(i,3)
                  fie=(        fac1_non(i)*p_i+fac_non(i)   *p_k               &
                       +0.5*(cell%rhom(ii)*p_i+cell%rhom(kk)*p_k))*gfactor(ii)
               ENDIF
               fie_nf(i0,1)=fie*sv_nf(i1,1)
               fie_nf(i0,2)=fie*sv_nf(i1,2)
               fie_nf(i0,3)=fie*sv_nf(i1,3)
            ENDIF
         ENDDO
!
         nv=-1
         nf_number=nf_number_id(nv)
         len   =istart_nf(2,nf_number)
         DO i=1,len
            k=right_nb_k(i)
            ii=right_non(k)
            kk=left_nf(k)
            k0=i_neigh(ii)-1
            jk=kneigh_non(k)
!
            IF(ngrad(ii).eq.0)THEN
               fie_i=0.d0
               DO jd=1,num_nd(jk+k0)
                  k1=node_face_cell(jd,jk+k0)
                  fie_i=fie_i+sn(k1)
               ENDDO
               fie=fie_i/DFLOAT(num_nd(jk+k0))
!
               IF(ig.eq.3)THEN
                  p_i=grav(1)*dxfc_nf(k,1)   +grav(2)*dxfc_nf(k,2)   +grav(3)*dxfc_nf(k,3)
                  p_k=grav(1)*dxfc_non_k(k,1)+grav(2)*dxfc_non_k(k,2)+grav(3)*dxfc_non_k(k,3)
                  fie=fie+0.5*(cell%rhom(ii)*p_i+cell%rhom(kk)*p_k)*gfactor(ii)
               ENDIF
               fie_nf(i,1)=-fie*sv_nf(k,1)
               fie_nf(i,2)=-fie*sv_nf(k,2)
               fie_nf(i,3)=-fie*sv_nf(k,3)
            ELSE
               IF(ig.eq.2)THEN
                  p_i=s(ii)
                  p_k=s(kk)
                  fie=fac1_non(k)*p_i+fac_non(k)*p_k
               ELSE
                  p_i=grav(1)*dxfc_nf(k,1)   +grav(2)*dxfc_nf(k,2)   +grav(3)*dxfc_nf(k,3)
                  p_k=grav(1)*dxfc_non_k(k,1)+grav(2)*dxfc_non_k(k,2)+grav(3)*dxfc_non_k(k,3)
                  fie=(        fac1_non(k)*p_i+fac_non(k)   *p_k               &
                       +0.5*(cell%rhom(ii)*p_i+cell%rhom(kk)*p_k))*gfactor(ii)
               ENDIF
               fie_nf(i,1)=-fie*sv_nf(k,1)
               fie_nf(i,2)=-fie*sv_nf(k,2)
               fie_nf(i,3)=-fie*sv_nf(k,3)
            ENDIF
         ENDDO
!
         DO nv=1,8
            nf_number=nf_number_id(nv)
            istart0=istart_nfs(nv)
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            IF(nv.ge.4 .and. nv.le.7) THEN
               DO i=1,len
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  i2=i_neigh(ii)-1
                  jj=jneigh_nf(i1)
!
                  IF(ngrad(ii).eq.0)THEN
                     fie_i=0.d0
                     DO jd=1,num_nd(jj+i2)
                        k=node_face_cell(jd,jj+i2)
                        fie_i=fie_i+sn(k)
                     ENDDO
                     fie=fie_i/DFLOAT(num_nd(jj+i2))
!
                     IF(ig.eq.3)THEN
                        p_i=grav(1)*dxfc_nf(i1,1)+grav(2)*dxfc_nf(i1,2)+grav(3)*dxfc_nf(i1,3)
                        fie=fie+cell%rhom(ii)*p_i*gfactor(ii)
                     ENDIF
                     fie_nf(i0,1)=fie*sv_nf(i1,1)
                     fie_nf(i0,2)=fie*sv_nf(i1,2)
                     fie_nf(i0,3)=fie*sv_nf(i1,3)
                  ENDIF
               ENDDO
            ELSE
               DO i=1,len
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
!
                  IF(ngrad(ii).eq.0)THEN
                     fie_nf(i0,1)=0.d0
                     fie_nf(i0,2)=0.d0
                     fie_nf(i0,3)=0.d0
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
!
         DO nv=1,8
            nf_number=nf_number_id(nv)
            istart0=istart_nfs(nv)
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            IF(nv.eq.1 .or. nv.eq.3 .or. nv.eq.8) THEN
               DO i=1,len
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  kk=right_non(i)
                  i2=i_neigh(ii)-1
                  jj=jneigh_nf(i1)
!
                  IF(ngrad(ii).ne.0)THEN
                    fie=s(ii)
                    fie_nf(i0,1)=fie*sv_nf(i1,1)
                    fie_nf(i0,2)=fie*sv_nf(i1,2)
                    fie_nf(i0,3)=fie*sv_nf(i1,3)
                  ENDIF
               ENDDO
            ELSE
               DO i=1,len
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  kk=right_non(i)
                  jj=jneigh_nf(i1)
!
                  IF(ngrad(ii).ne.0)THEN
                     IF(ig.eq.2) THEN
                       fie=s(ii)
                     ELSE
                        p_i=grav(1)*dxfc_nf(i1,1)+grav(2)*dxfc_nf(i1,2)+grav(3)*dxfc_nf(i1,3)
                        fie=s(ii)+cell%rhom(ii)*p_i*gfactor(ii)
                     ENDIF
                     fie_nf(i0,1)=fie*sv_nf(i1,1)
                     fie_nf(i0,2)=fie*sv_nf(i1,2)
                     fie_nf(i0,3)=fie*sv_nf(i1,3)
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDIF
!
      CALL sum_nf_ndim(0,0,ncell_fp, &
                       fie_nf,dsdx)
!
      IF(ndim.eq.2) THEN
         DO i=1,ncell_fluid
            dsdx(i,1)=dsdx(i,1)*volr(i)
            dsdx(i,2)=dsdx(i,2)*volr(i)
         ENDDO
      ELSE
         DO i=1,ncell_fluid
            dsdx(i,1)=dsdx(i,1)*volr(i)
            dsdx(i,2)=dsdx(i,2)*volr(i)
            dsdx(i,3)=dsdx(i,3)*volr(i)
         ENDDO
      ENDIF
!
      END SUBROUTINE grad_frink
!
      SUBROUTINE grad_frink_vel(s,dsdx,vb,vin)
!
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: np
      USE Zparam       , ONLY: ndim,nb_max
      USE Zvec_param   , ONLY: nf_nonk,nf_non,nf_mcc,nf_inl,nf_out,nf_adw,nf_fsw,nf_ctw,nf_chw,nf_sym,nf_totk
      USE Zvec_index   , ONLY: left_nf,right_non,jneigh_nf,kneigh_non,nbcon_nf
      USE Znum_cell    , ONLY: i_neigh,istart_nf,istart_nbcon_nf, &
                               nf_number_nb,lens, &
                               right_nb_k,istart_nf,nf_number_id,istart_nfs
      USE Znode        , ONLY: n_node,num_cell_node,cell_node,num_nd,node_face_cell,rwcn,swcn
      USE Zuserdefined , ONLY: vel_bc_profile_inl
      USE c3com_cupid  , ONLY: i3invtbl
      USE Zbc_index    , ONLY: vin_norm,ngrad
      USE Zvec_geo     , ONLY: xn_nf,sv_nf,     &
                               fac_non,fac1_non
! 
      IMPLICIT NONE
      INCLUDE '../../10_LinkToMARS/c3com.h' 
!
!.....Input
      REAL(8),DIMENSION(ncell_fp,ndim) :: s
      REAL(8),DIMENSION(nb_max,ndim) :: vb
      REAL(8),DIMENSION(nb_max) :: vin
!.....Output
      REAL(8),DIMENSION(ncell_fp,ndim,ndim) :: dsdx
!.....Local variables
      INTEGER :: i,j,k,jd,idx,ix,k0
      INTEGER :: ii,jj,kk,jk
      INTEGER :: nv,nf_number,istart,len,istart0,istart2,i0,i1,i2
      REAL(8) f_profile
      REAL(8) fie1,fie2,fie3
!.....Local arrays
      REAL(8),DIMENSION(n_node,ndim) :: sn
!.....Local vector arrays      
      REAL(8),DIMENSION(nf_totk,ndim,ndim) :: fie_nf
!
!.....Define boundary nodes and cells
!.....Interpolate node values using inverse distance weighting
!
      sn(:,:)=0.0d0
      if(ndim.eq.2) THEN
         DO i=1,ncell_fluid
            DO j=1,num_cell_node(i)
               k=cell_node(j,i)
               sn(k,1)=sn(k,1)+rwcn(j,i)*s(i,1)
               sn(k,2)=sn(k,2)+rwcn(j,i)*s(i,2)
            ENDDO
         ENDDO
         DO i=1,n_node
            sn(i,1)=sn(i,1)/swcn(i)
            sn(i,2)=sn(i,2)/swcn(i)
         ENDDO
      ELSE
         DO i=1,ncell_fluid
            DO j=1,num_cell_node(i)
               k=cell_node(j,i)
               sn(k,1)=sn(k,1)+rwcn(j,i)*s(i,1)
               sn(k,2)=sn(k,2)+rwcn(j,i)*s(i,2)
               sn(k,3)=sn(k,3)+rwcn(j,i)*s(i,3)
            ENDDO
         ENDDO
         DO i=1,n_node
            sn(i,1)=sn(i,1)/swcn(i)
            sn(i,2)=sn(i,2)/swcn(i)
            sn(i,3)=sn(i,3)/swcn(i)
         ENDDO
      ENDIF
!
!......Allreduce sn
!
      IF(np.gt.1) THEN
         DO ix=1,ndim
            CALL allreducei_r(sn(1,ix),n_node)
         ENDdo
      ENDIF
!
!.....Build summation info for non,inl
!
      nf_number_nb=8
      nf_number_id(-1)=-1
      nf_number_id(0)=0
      nf_number_id(1)=1
      nf_number_id(2)=2
      nf_number_id(3)=3
      nf_number_id(4)=4
      nf_number_id(5)=5
      nf_number_id(6)=6
      nf_number_id(7)=7
      nf_number_id(8)=8
      istart_nfs(-1)=0
      istart_nfs(0)=istart_nfs(-1)+nf_nonk
      istart_nfs(1)=istart_nfs(0) +nf_non
      istart_nfs(2)=istart_nfs(1) +nf_mcc
      istart_nfs(3)=istart_nfs(2) +nf_inl
      istart_nfs(4)=istart_nfs(3) +nf_out
      istart_nfs(5)=istart_nfs(4) +nf_adw
      istart_nfs(6)=istart_nfs(5) +nf_fsw
      istart_nfs(7)=istart_nfs(6) +nf_ctw
      istart_nfs(8)=istart_nfs(7) +nf_chw
      lens         =istart_nfs(8) +nf_sym
!
!....Cell gradient by Frink's method
!
!
      IF(ndim.eq.2) THEN
         DO nv=0,8
            nf_number=nf_number_id(nv)
            istart0=istart_nfs(nv)
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               i2=i_neigh(ii)-1
               jj=jneigh_nf(i1)
               IF(ngrad(ii).eq.0)THEN
                  fie1=0.d0
                  fie2=0.d0
                  DO jd=1,num_nd(jj+i2)
                     k=node_face_cell(jd,jj+i2)
                     fie1=fie1+sn(k,1)
                     fie2=fie2+sn(k,2)
                 ENDDO
                 fie1=fie1/DFLOAT(num_nd(jj+i0))
                 fie2=fie2/DFLOAT(num_nd(jj+i0))
               ELSE
                 IF(nv.eq.0) THEN
                    kk=right_non(i)
                    fie1=fac1_non(i)*s(ii,1)+fac_non(i)*s(kk,1)
                    fie2=fac1_non(i)*s(ii,2)+fac_non(i)*s(kk,2)
                 ELSEIF(nv.eq.1) THEN
                    idx=i3invtbl(i)
                    fie1=c3vl(1,idx)*xn_nf(i1,1)
                    fie2=c3vl(1,idx)*xn_nf(i1,2)
                 ELSEIF(nv.eq.2) THEN
                    istart2=istart_nbcon_nf(nf_number)
                    i2=istart2+i
                    k=nbcon_nf(i2)
                    f_profile=vel_bc_profile_inl(i)
                    IF(vin_norm(k).eq.0)THEN
                       fie1=vb(k,1)*f_profile
                       fie2=vb(k,2)*f_profile
                    ELSE
                       fie1=vin(k)*xn_nf(i1,1)*f_profile
                       fie2=vin(k)*xn_nf(i1,2)*f_profile
                    ENDIF
                 ELSEIF(nv.eq.3 .or. nv.eq.8) THEN
                    fie1=s(ii,1)
                    fie2=s(ii,2)
                 ELSE
                    fie1=0.d0
                    fie2=0.d0
                 ENDIF
               ENDIF
               fie_nf(i0,1,1)=fie1*sv_nf(i1,1)
               fie_nf(i0,2,1)=fie1*sv_nf(i1,2)
               fie_nf(i0,1,2)=fie2*sv_nf(i1,1)
               fie_nf(i0,2,2)=fie2*sv_nf(i1,2)
            ENDDO
         ENDDO
!
         nv=-1
         nf_number=nf_number_id(nv)
         len   =istart_nf(2,nf_number)
         DO i=1,len
            k=right_nb_k(i)
            ii=right_non(k)
            kk=left_nf(k)
            k0=i_neigh(ii)-1
            jk=kneigh_non(k)
!
            IF(ngrad(ii).eq.0)THEN
               fie1=0.d0
               fie2=0.d0
               DO jd=1,num_nd(jk+k0)
                  k=node_face_cell(jd,jk+k0)
                  fie1=fie1-sn(k,1)
                  fie2=fie2-sn(k,2)
               ENDDO
               fie1=fie1/DFLOAT(num_nd(jk+k0))
               fie2=fie2/DFLOAT(num_nd(jk+k0))
            ELSE
               fie1=-fac1_non(k)*s(ii,1)-fac_non(k)*s(kk,1)
               fie2=-fac1_non(k)*s(ii,2)-fac_non(k)*s(kk,2)
            ENDIF
            fie_nf(i,1,1)=fie1*sv_nf(k,1)
            fie_nf(i,2,1)=fie1*sv_nf(k,2)
            fie_nf(i,1,2)=fie2*sv_nf(k,1)
            fie_nf(i,2,2)=fie2*sv_nf(k,2)
         ENDDO
      ELSE
         DO nv=0,8
            nf_number=nf_number_id(nv)
            istart0=istart_nfs(nv)
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               i2=i_neigh(ii)-1
               jj=jneigh_nf(i1)
               IF(ngrad(ii).eq.0)THEN
                  fie1=0.d0
                  fie2=0.d0
                  fie3=0.d0
                  DO jd=1,num_nd(jj+i2)
                     k=node_face_cell(jd,jj+i2)
                     fie1=fie1+sn(k,1)
                     fie2=fie2+sn(k,2)
                     fie3=fie3+sn(k,3)
                 ENDDO
                 fie1=fie1/DFLOAT(num_nd(jj+i0))
                 fie2=fie2/DFLOAT(num_nd(jj+i0))
                 fie3=fie3/DFLOAT(num_nd(jj+i0))
               ELSE
                 IF(nv.eq.0) THEN
                    kk=right_non(i)
                    fie1=fac1_non(i)*s(ii,1)+fac_non(i)*s(kk,1)
                    fie2=fac1_non(i)*s(ii,2)+fac_non(i)*s(kk,2)
                    fie3=fac1_non(i)*s(ii,3)+fac_non(i)*s(kk,3)
                 ELSEIF(nv.eq.1) THEN
                    idx=i3invtbl(i)
                    fie1=c3vl(1,idx)*xn_nf(i1,1)
                    fie2=c3vl(1,idx)*xn_nf(i1,2)
                    fie3=c3vl(1,idx)*xn_nf(i1,3)
                 ELSEIF(nv.eq.2) THEN
                    istart2=istart_nbcon_nf(nf_number)
                    i2=istart2+i
                    k=nbcon_nf(i2)
                    f_profile=vel_bc_profile_inl(i)
                    IF(vin_norm(k).eq.0)THEN
                       fie1=vb(k,1)*f_profile
                       fie2=vb(k,2)*f_profile
                       fie3=vb(k,3)*f_profile
                    ELSE
                       fie1=vin(k)*xn_nf(i1,1)*f_profile
                       fie2=vin(k)*xn_nf(i1,2)*f_profile
                       fie3=vin(k)*xn_nf(i1,3)*f_profile
                    ENDIF
                 ELSEIF(nv.eq.3 .or. nv.eq.8) THEN
                    fie1=s(ii,1)
                    fie2=s(ii,2)
                    fie3=s(ii,3)
                 ELSE
                    fie1=0.d0
                    fie2=0.d0
                    fie3=0.d0
                 ENDIF
               ENDIF
               fie_nf(i0,1,1)=fie1*sv_nf(i1,1)
               fie_nf(i0,2,1)=fie1*sv_nf(i1,2)
               fie_nf(i0,3,1)=fie1*sv_nf(i1,3)
               fie_nf(i0,1,2)=fie2*sv_nf(i1,1)
               fie_nf(i0,2,2)=fie2*sv_nf(i1,2)
               fie_nf(i0,3,2)=fie2*sv_nf(i1,3)
               fie_nf(i0,1,3)=fie3*sv_nf(i1,1)
               fie_nf(i0,2,3)=fie3*sv_nf(i1,2)
               fie_nf(i0,3,3)=fie3*sv_nf(i1,3)
            ENDDO
         ENDDO
!
         nv=-1
         nf_number=nf_number_id(nv)
         len   =istart_nf(2,nf_number)
         DO i=1,len
            k=right_nb_k(i)
            ii=right_non(k)
            kk=left_nf(k)
            k0=i_neigh(ii)-1
            jk=kneigh_non(k)
!
            IF(ngrad(ii).eq.0)THEN
               fie1=0.d0
               fie2=0.d0
               fie3=0.d0
               DO jd=1,num_nd(jk+k0)
                  k=node_face_cell(jd,jk+k0)
                  fie1=fie1-sn(k,1)
                  fie2=fie2-sn(k,2)
                  fie3=fie3-sn(k,3)
               ENDDO
               fie1=fie1/DFLOAT(num_nd(jk+k0))
               fie2=fie2/DFLOAT(num_nd(jk+k0))
               fie3=fie3/DFLOAT(num_nd(jk+k0))
            ELSE
               fie1=-fac1_non(k)*s(ii,1)-fac_non(k)*s(kk,1)
               fie2=-fac1_non(k)*s(ii,2)-fac_non(k)*s(kk,2)
               fie3=-fac1_non(k)*s(ii,3)-fac_non(k)*s(kk,3)
            ENDIF
            fie_nf(i,1,1)=fie1*sv_nf(k,1)
            fie_nf(i,2,1)=fie1*sv_nf(k,2)
            fie_nf(i,3,1)=fie1*sv_nf(k,3)
            fie_nf(i,1,2)=fie2*sv_nf(k,1)
            fie_nf(i,2,2)=fie2*sv_nf(k,2)
            fie_nf(i,3,2)=fie2*sv_nf(k,3)
            fie_nf(i,1,3)=fie3*sv_nf(k,1)
            fie_nf(i,2,3)=fie3*sv_nf(k,2)
            fie_nf(i,3,3)=fie3*sv_nf(k,3)
         ENDDO
      ENDIF
!
      CALL sum_nf_ndim2(0,0,ncell_fp, &
                        fie_nf,dsdx)
!
      END SUBROUTINE grad_frink_vel
!
