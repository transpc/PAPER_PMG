!
      SUBROUTINE subdomain_info(nelem,i_cell_node0) 
!
!     This routine defines local arrays for each subdomain.  All the 
!     arrays can be automatically defined based on the coloring (manually
!     or METIS). It includes Window METIS (Ver.5.0) with having been 
!     activated in preprocessor option (metis=.true.)
!
      USE Zinterface
      USE Zmpi        , ONLY: celem,ncell_fp,niut,iut,ri,si,rintf,sintf,    &
                              maxmt_fp,maxmt_fluid,maxmt_cell,maxmt_nfluid, &
                              jperm,jjperm,                                 &
                              ncell_fluid1,ncell_fluid1_dsp,                &
                              ncell_fluid1_2d,ncell_fluid1_2d_dsp,          &
                              ncell_csr_sz,ncell_csr_dsp,                   &
                              ncell_ndim_sz,ncell_ndim_dsp,                 &
                              ncell_node_csr_sz,ncell_node_csr_dsp,         &
                              ncell_csr_sz_nbcon0,ncell_csr_dsp_nbcon0,     &
                              i3perm,                                       &
                              metis,mapping_ext
      USE Zzone       , ONLY: ncell_fluid,ncell_fluid_all
      USE Zcore       , ONLY: np,myrank
      USE Zparam      , ONLY: nn,ns,ndim,nin_max,nb_max,mesh_openfoam
      USE Zbc_index   , ONLY: nbcon
      USE Znum_cell   , ONLY: i_neigh_tmp,j_neigh_tmp,j_nbcon_tmp,index_sort_tmp1, &
                              i_neigh,neigh,index_sort,indexr_sort
      USE Znormal     , ONLY: nji,i_neigh_nbcon0,maxmt_fluid_nbcon0
      USE Zgradoption , ONLY: ifrink
      USE Zio_unit    , ONLY: unit_log
! PMG 
      USE MD_parameter, ONLY: nf_max
      USE MD_geometry, ONLY: nelem_mg,num_neigh_mg,neigh_mg
!
      IMPLICIT NONE      
!
!.....Input
      INTEGER :: nelem
      INTEGER :: i_cell_node0(nn+1)
!.....Local variables
      INTEGER :: i,j,k,jj,ie,ne
      INTEGER :: j0,j1,j2,k1
      INTEGER :: ip,jp,iptr
      INTEGER :: nelem_sub
      INTEGER :: celem_max
      INTEGER :: cinter,cext0
      INTEGER :: maxmt_fluid_nbcon00
!.....Local arrays
      INTEGER :: neighl(ns)
      INTEGER :: temp_sort(ns)
      INTEGER :: iutjp(np)
      INTEGER :: cext(np)
      INTEGER :: irecv_cnt(np),jsend_cnt(np)
      INTEGER :: icount(np)
      INTEGER :: flag(nelem),flagt(nelem)
      INTEGER :: ia(np+1),ja(nelem)
      INTEGER :: jperms(nelem)
      INTEGER :: itmp(2)
      INTEGER :: i_neigh_tmp_nbcon0(nelem+1)
!.....Local allocatable arrays
      INTEGER,ALLOCATABLE :: num_neigh(:)
      INTEGER,ALLOCATABLE :: irecv(:),jsend(:)
      INTEGER,DIMENSION(:),ALLOCATABLE :: ia_sub,ja_sub
      INTEGER,DIMENSION(:),ALLOCATABLE :: nji_cell
!
!.....Delete zero entries.
!
      IF(myrank.eq.0) THEN
         ALLOCATE(ia_sub(nelem+1))
         ia_sub(1)=1
         k=1
         DO i=1,nelem
            DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
               IF(j_nbcon_tmp(j).eq.0)THEN
                  k=k+1
               ENDIF
            ENDDO
            ia_sub(i+1)=k
         ENDDO
         nelem_sub=ia_sub(nelem+1)-1 
         ALLOCATE(ja_sub(nelem_sub))
         celem_max=0
         k=0
         DO i=1,nelem
            DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
               IF(j_nbcon_tmp(j).eq.0)THEN
                  k=k+1
                  ja_sub(k)=j_neigh_tmp(j)
               ENDIF
            ENDDO
            celem_max=max(celem(i),celem_max)
         ENDDO
         itmp(1)=celem_max
         itmp(2)=nelem_sub
      ENDIF
      CALL broadcast_i(itmp,2)
      celem_max=itmp(1)
      nelem_sub=itmp(2)
!
!.....Window METIS Ver.5.0 API
!
      IF(metis.eq.1) then
         IF(myrank.eq.0) THEN
            CALL domain_decomposition(nelem,nelem_sub,ia_sub,ja_sub,celem,celem_max)
            WRITE(*,*) '          1The number of subdomains is', celem_max
            WRITE(*,*) '          1The number of processors is', np
            WRITE(unit_log,*) '          1The number of subdomains is', celem_max
            WRITE(unit_log,*) '          1The number of processors is', np         
         ENDIF   
      ELSEIF(metis.eq.0) then
         IF(np.ne.celem_max) THEN
            IF(myrank.eq.0) THEN
               WRITE(*,*) '          2The number of subdomains is', celem_max
               WRITE(*,*) '          2The number of processors is', np
               PRINT *,'          2The number of subdomains differs from the number of processors !!!'
               WRITE(unit_log,*) '          2The number of subdomains is', celem_max
               WRITE(unit_log,*) '          2The number of processors is', np
            ENDIF
            CALL finalize_mpi
            STOP
         ENDIF
      ENDIF
      IF(myrank.ne.0) THEN
         ALLOCATE(celem(nelem))
         ALLOCATE(i_neigh_tmp(nelem+1))
         ALLOCATE(ia_sub(nelem+1))
         ALLOCATE(ja_sub(nelem_sub))
      ENDIF
      CALL broadcast_i(celem,nelem)
      CALL broadcast_i(i_neigh_tmp,nelem+1)
      CALL broadcast_i(ia_sub,nelem+1) 
      CALL broadcast_i(ja_sub,nelem_sub)
! PMG copy to PMG solver
!      ndom = np
      nelem_mg = nelem
      nf_max = ns
      allocate(num_neigh_mg(nelem_mg),neigh_mg(nf_max,nelem_mg))
      DO i=1,nelem
         num_neigh_mg(i)=ia_sub(i+1)-ia_sub(i)
         k=0
         DO j=ia_sub(i),ia_sub(i+1)-1
            k=k+1 
            neigh_mg(k,i)=ja_sub(j)
         ENDDO
      ENDDO
!
      ALLOCATE(ncell_fluid1(np),ncell_fluid1_dsp(np))
      ALLOCATE(ncell_fluid1_2d(np),ncell_fluid1_2d_dsp(np))
      ALLOCATE(ncell_csr_sz(np),ncell_csr_dsp(np))
      ALLOCATE(ncell_ndim_sz(np),ncell_ndim_dsp(np))
      ALLOCATE(ncell_csr_sz_nbcon0(np),ncell_csr_dsp_nbcon0(np))
!
!.....Build index based on proc id to access global cell array directly
!.....Get ncell_fluid for all procs
!
      ia(1)=1
      DO ip=1,np
         ia(ip+1)=0
         icount(ip)=0
      ENDDO
      DO ie=1,nelem
         ip=celem(ie)
         ia(ip+1)=ia(ip+1)+1
      ENDDO
      DO ip=1,np
         ncell_fluid1(ip)=ia(ip+1)
         ncell_fluid1_2d(ip)=ia(ip+1)*ndim
         ia(ip+1)=ia(ip+1)+ia(ip)
      ENDDO
      DO ie=1,nelem
         ip=celem(ie)
         ja(ia(ip)+icount(ip))=ie
         icount(ip)=icount(ip)+1
      ENDDO
!
!.....Count interior cells, external=receive cells
!.....Count cell only once on external cells
!
      ip=myrank+1
      ncell_fluid=ncell_fluid1(ip)
      ALLOCATE(jperm(ncell_fluid))
      cinter=0
      DO i=1,nelem
         flag(i)=0
         jperms(i)=0
      ENDDO
      DO jp=1,np
         irecv_cnt(jp)=0
      ENDDO
      iptr=0
      DO jj=ia(ip),ia(ip+1)-1 
         ie=ja(jj)
         DO j=ia_sub(ie),ia_sub(ie+1)-1
            k=ja_sub(j)
            jp=celem(k)
            IF(jp.ne.ip) THEN
!
!..............Cext cells
!
               DO j1=j,ia_sub(ie+1)-1
                  k=ja_sub(j1)
                  jp=celem(k)
                  IF(jp.ne.ip) THEN
                     IF(flag(k).eq.0) then
                        flag(k)=1
                        iptr=iptr+1
                        flagt(iptr)=k
                        irecv_cnt(jp)=irecv_cnt(jp)+1
                     ENDIF
                  ENDIF
               ENDDO
               GOTO 100
            ENDIF
         ENDDO
!........Cintr cells
!........Get the mapping for interior cell
         cinter=cinter+1
         jperm(cinter)=ie
         jperms(ie)=cinter
100      CONTINUE 
      ENDDO
      cext0=0
      DO jp=1,np
         cext0=cext0+irecv_cnt(jp)
      ENDDO
      ncell_fp=ncell_fluid+cext0
!
!.....Get all cext needed to build ALLGATHERV parameters 
!
      CALL allgather_i(cext0,cext)
!
!.....Receive cells
!.....Remove zero count
      niut=0
      DO jp=1,np
         IF(irecv_cnt(jp).ne.0) niut=niut+1
      ENDDO
      ALLOCATE(iut(niut),ri(niut+1),si(niut+1))
      ri(1)=1
      j=0
      DO jp=1,np
         IF(irecv_cnt(jp).ne.0) THEN
            j=j+1
            iut(j)=jp
            iutjp(jp)=j
            ri(j+1)=irecv_cnt(jp)
         ENDIF
      ENDDO
      DO j=1,niut
         ri(j+1)=ri(j+1)+ri(j)
      ENDDO
      ALLOCATE(irecv(ri(niut+1)-1))
      DO jp=1,np
         icount(jp)=0
      ENDDO
!
!.....Minimize the zeroing of flag via flagt
!
      DO i=1,cext0
         k=flagt(i)
         jp=celem(k)
         j=iutjp(jp)
         irecv(ri(j)+icount(jp))=k
         icount(jp)=icount(jp)+1
         flag(k)=0
      ENDDO
!
!.....Count remote=send cells
!.....Count cell only once on remote cells for each ip
!
      jp=myrank+1
      iptr=0
      DO ip=1,np
         jsend_cnt(ip)=0
         IF(jp.ne.ip) THEN
            DO jj=ia(ip),ia(ip+1)-1 
               ie=ja(jj)
!..............Bypass scan of interior points
               IF(jperms(ie).eq.0) THEN
                  DO j=ia_sub(ie),ia_sub(ie+1)-1
                     k=ja_sub(j)
                     j0=celem(k)
                     IF(j0.eq.jp) THEN
                        IF(flag(k).eq.0) then
                           flag(k)=1
                           jsend_cnt(ip)=jsend_cnt(ip)+1
                           flagt(iptr+jsend_cnt(ip))=k
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
!...........Minimize the zeroing of flag via flagt
            DO i=1,jsend_cnt(ip)
               k=flagt(iptr+i)
               flag(k)=0
            ENDDO
            iptr=iptr+jsend_cnt(ip)
         ENDIF
      ENDDO
      DEALLOCATE(ia_sub,ja_sub)
!
!.....Send cells
!
      si(1)=1
      j=0
      DO jp=1,np
         IF(jsend_cnt(jp).ne.0) THEN
            j=j+1
!           iut(j)=jp
!           iutjp(jp)=j
            si(j+1)=jsend_cnt(jp)
         ENDIF
      ENDDO
      DO j=1,niut
         si(j+1)=si(j+1)+si(j)
      ENDDO
!
      iut=iut-1
!
!.....Get all neighbors counted once jsend
!
      ALLOCATE(jsend(si(niut+1)-1))
      DO i=1,si(niut+1)-1
         jsend(i)=flagt(i)
      ENDDO
!
!.....Get the mapping for remote cells counted once.
!
      k=cinter
      DO i=1,si(niut+1)-1
         k1=jsend(i)
         IF(flag(k1).eq.0) then
            k=k+1
            flag(k1)=1
            jperm(k)=k1
            jperms(k1)=k
         ENDIF
      ENDDO
!
!.....Get the mapping for ext cells
!
      DO i=1,ri(niut+1)-1
         k1=irecv(i)
         jperms(k1)=ncell_fluid+i
      ENDDO
!
!.....Convert global numbering to local numbering for send,receive cells
!
      ALLOCATE(rintf(ri(niut+1)-1),sintf(si(niut+1)-1))
      DO i=1,ri(niut+1)-1
         rintf(i)=ncell_fluid+i
      ENDDO
      DO i=1,si(niut+1)-1
         k1=jsend(i)
         sintf(i)=jperms(k1)
      ENDDO
!
!.....Get size and displacement to call ALLGATHERV
!
      ip=1
      ncell_fluid1_dsp(ip)=0
      ncell_fluid1_2d_dsp(ip)=0
      DO ip=2,np
         ncell_fluid1_dsp(ip)=ncell_fluid1_dsp(ip-1)+ncell_fluid1(ip-1)
         ncell_fluid1_2d_dsp(ip)=ncell_fluid1_2d_dsp(ip-1)+ncell_fluid1_2d(ip-1)
      ENDDO
!
!.....Build jjperm out of local jperm
!
      ALLOCATE(jjperm(nelem))
      CALL allgather_vec_i(jperm,ncell_fluid,jjperm,nelem,ncell_fluid1,ncell_fluid1_dsp)
!
!.....Get size and displacement to call ndim type ALLGATHERV
!
      DO ip=1,np
         ncell_ndim_sz(ip)=ncell_fluid1(ip)*ndim
      ENDDO
      ip=1
      ncell_ndim_dsp(ip)=0
      DO ip=2,np
         ncell_ndim_dsp(ip)=ncell_ndim_dsp(ip-1)+ncell_ndim_sz(ip-1)
      ENDDO
!
!.....Get size and displacement to call csr type ALLGATHERV
!
      DO ip=1,np
         k=0
         DO jj=ia(ip),ia(ip+1)-1 
            ie=ja(jj)
            k=k+(i_neigh_tmp(ie+1)-i_neigh_tmp(ie))
         ENDDO
         ncell_csr_sz(ip)=k
      ENDDO
      ip=1
      ncell_csr_dsp(ip)=0
      DO ip=2,np
         ncell_csr_dsp(ip)=ncell_csr_dsp(ip-1)+ncell_csr_sz(ip-1)
      ENDDO
!
!.....Get size and displacement to call csr type ALLGATHERV for num_cell_node_tmp
!
      IF(mesh_openfoam.eq.1.and.ifrink.ge.1)THEN
         ALLOCATE(ncell_node_csr_sz(np),ncell_node_csr_dsp(np))
         DO ip=1,np
            k=0
            DO jj=ia(ip),ia(ip+1)-1 
               ie=ja(jj)
               k=k+(i_cell_node0(ie+1)-i_cell_node0(ie))
            ENDDO
            ncell_node_csr_sz(ip)=k
         ENDDO
         ip=1
         ncell_node_csr_dsp(ip)=0
         DO ip=2,np
            ncell_node_csr_dsp(ip)=ncell_node_csr_dsp(ip-1)+ncell_node_csr_sz(ip-1)
         ENDDO
      ENDIF
!
!.....Define local neigh & num_neigh for fluid region
!
      ALLOCATE(num_neigh(ncell_fp))
!.....Need to get hold of maxmt_fluid here move it from read_grid
      ALLOCATE(i_neigh(ncell_fp+1))
      i_neigh(1)=1
      DO i=1,ncell_fluid
         ne=jperm(i)
         num_neigh(i)=i_neigh_tmp(ne+1)-i_neigh_tmp(ne)
         i_neigh(i+1)=i_neigh(i)+num_neigh(i)
      ENDDO
      maxmt_fluid=i_neigh(ncell_fluid+1)-1
      ALLOCATE(index_sort(maxmt_fluid),indexr_sort(maxmt_fluid))
      ALLOCATE(neigh(maxmt_fluid))
      ALLOCATE(nbcon(maxmt_fluid))
!
      IF(myrank.ne.0) ALLOCATE(j_neigh_tmp(1),j_nbcon_tmp(1))
      CALL scatterv_csr_i(neigh,maxmt_fluid,j_neigh_tmp,maxmt_cell,  &
                          ncell_fluid_all,maxmt_nfluid,i_neigh_tmp,0)
      CALL scatterv_csr_i(nbcon,maxmt_fluid,j_nbcon_tmp,maxmt_cell,  &
                          ncell_fluid_all,maxmt_nfluid,i_neigh_tmp,0)
      DO i=1,ncell_fluid
         j1=i_neigh(i)-1
!
!........Sort in the vector order
!
         DO j=1,num_neigh(i)
            index_sort(j+j1)=j
            IF(nbcon(j1+j).eq.0)THEN
               k=neigh(j1+j)
               neighl(j)=jperms(k)
            ELSE
               IF    (nbcon(j1+j).ge.201)THEN
                  neighl(j)=ncell_fp+1
               ELSEIF(nbcon(j1+j).gt.0.and.nbcon(j1+j).le.nin_max)THEN
                  neighl(j)=ncell_fp+2
               ELSEIF(nbcon(j1+j).gt.nin_max.and.nbcon(j1+j).le.nb_max)THEN
                  neighl(j)=ncell_fp+3
               ELSEIF(nbcon(j1+j).eq.-1)THEN
                  neighl(j)=ncell_fp+4
               ELSEIF(nbcon(j1+j).eq.-2)THEN
                  neighl(j)=ncell_fp+5
               ELSEIF(nbcon(j1+j).eq.-3.or.nbcon(j1+j).eq.-4)THEN
                  neighl(j)=ncell_fp+6
               ELSEIF(nbcon(j1+j).eq.-5.or.nbcon(j1+j).eq.-6)THEN
                  neighl(j)=ncell_fp+7
               ELSEIF(nbcon(j1+j).eq.101)THEN
                  neighl(j)=ncell_fp+8
               ELSEIF(nbcon(j1+j).ge.-39 .and. nbcon(j1+j).le.-31)THEN
                  neighl(j)=ncell_fp+9
               ELSEIF(nbcon(j1+j).ge.-49 .and. nbcon(j1+j).le.-41)THEN
                  neighl(j)=ncell_fp+10
               ELSEIF(nbcon(j1+j).ge.-59 .and. nbcon(j1+j).le.-51)THEN
                  neighl(j)=ncell_fp+11
               ENDIF
            ENDIF
         ENDDO
         CALL sortx_i(neighl,index_sort(j1+1),num_neigh(i))
!
!........Restore neigh to 0 for non compute
!
         DO j=num_neigh(i),1,-1
            if(neighl(j).le.ncell_fp) exit
            neighl(j)=0
         ENDDO
!........Get sorted nbcon
         DO j=1,num_neigh(i)
            k=index_sort(j+j1)
            indexr_sort(k+j1)=j
            temp_sort(j)=nbcon(j1+k)
         ENDDO
         DO j=1,num_neigh(i)
            neigh(j+j1)=neighl(j)
            nbcon(j+j1)=temp_sort(j)
         ENDDO
      ENDDO
!
!.....Get index_sort global
!
      IF(myrank.eq.0) THEN
         ALLOCATE(index_sort_tmp1(maxmt_cell))
      ELSE
         ALLOCATE(index_sort_tmp1(1))
      ENDIF
      CALL gatherv_csr_i(index_sort,maxmt_fluid,index_sort_tmp1,maxmt_cell, &
                         ncell_fluid_all,maxmt_nfluid,i_neigh_tmp,0)
!
!.....Get sorted global j_neigh_tmp to get nji
!.....Get sorted global j_nbcon_tmp to be used in solid
!
      IF(myrank.eq.0) THEN
         DO ne=1,nelem
            j0=i_neigh_tmp(ne)-1
            DO j=i_neigh_tmp(ne),i_neigh_tmp(ne+1)-1
               k=index_sort_tmp1(j)
               temp_sort(j-j0)=j_neigh_tmp(j0+k)
            ENDDO
            DO j=i_neigh_tmp(ne),i_neigh_tmp(ne+1)-1
               j_neigh_tmp(j)=temp_sort(j-j0)
            ENDDO
            DO j=i_neigh_tmp(ne),i_neigh_tmp(ne+1)-1
               k=index_sort_tmp1(j)
               temp_sort(j-j0)=j_nbcon_tmp(j0+k)
            ENDDO
            DO j=i_neigh_tmp(ne),i_neigh_tmp(ne+1)-1
               j_nbcon_tmp(j)=temp_sort(j-j0)
            ENDDO
         ENDDO
      ENDIF
!
      IF(myrank.eq.0) THEN
         i_neigh_tmp_nbcon0(1)=1
         j0=1
         DO ne=1,ncell_fluid_all
            DO j=i_neigh_tmp(ne),i_neigh_tmp(ne+1)-1
               IF(j_nbcon_tmp(j).ne.0) EXIT
               j0=j0+1
            ENDDO
            i_neigh_tmp_nbcon0(ne+1)=j0
         ENDDO
         maxmt_fluid_nbcon00=i_neigh_tmp_nbcon0(ncell_fluid_all+1)-1
         ALLOCATE(nji_cell(maxmt_fluid_nbcon00))
!
!........Compute nji in the global
!
         DO ne=1,nelem
            j2=i_neigh_tmp_nbcon0(ne)
            DO j=i_neigh_tmp(ne),i_neigh_tmp(ne+1)-1
               IF(j_nbcon_tmp(j).ne.0) CYCLE
                  k=j_neigh_tmp(j)
                  j0=i_neigh_tmp(k)-1
                  DO j1=i_neigh_tmp(k),i_neigh_tmp(k+1)-1
                     IF(j_neigh_tmp(j1).eq.ne) then
                         nji_cell(j2)=j1-j0
                         j2=j2+1
                         EXIT
                     ENDIF
                  ENDDO
            ENDDO
         ENDDO
      ELSE
         ALLOCATE(nji_cell(1))
      ENDIF
!
      ALLOCATE(i_neigh_nbcon0(ncell_fluid+1))
      i_neigh_nbcon0(1)=1
      j0=1
      DO i=1,ncell_fluid
         DO j=i_neigh(i),i_neigh(i+1)-1
            IF(nbcon(j).ne.0) EXIT
            j0=j0+1
         ENDDO
         i_neigh_nbcon0(i+1)=j0
      ENDDO
      maxmt_fluid_nbcon0=i_neigh_nbcon0(ncell_fluid+1)-1
      CALL allgather_i(maxmt_fluid_nbcon0,ncell_csr_sz_nbcon0)
      ip=1
      ncell_csr_dsp_nbcon0(ip)=0
      DO ip=2,np
         ncell_csr_dsp_nbcon0(ip)=ncell_csr_dsp_nbcon0(ip-1)+ncell_csr_sz_nbcon0(ip-1)
      ENDDO
!
      ALLOCATE(nji(maxmt_fluid_nbcon0))
      CALL scatterv_csr_i_nbcon0(nji,maxmt_fluid_nbcon0,nji_cell,maxmt_fluid_nbcon00, &
                                 ncell_fluid_all,maxmt_nfluid,i_neigh_tmp_nbcon0,0)
      DEALLOCATE(nji_cell)
!
      CALL communicate_1d_int(num_neigh)
      DO i=ncell_fluid+1,ncell_fp
         i_neigh(i+1)=i_neigh(i)+num_neigh(i)
      ENDDO
      maxmt_fp=i_neigh(ncell_fp+1)-1
!
!.....Global cell to local cell for CUPID_RV
!
      ALLOCATE(i3perm(nelem))
      i3perm(:)=0
      DO i=1,ncell_fluid
         i3perm(jperm(i))=i
      ENDDO        
!      
      ALLOCATE(mapping_ext(ncell_fp))
      CALL communicate_1d_int0(jsend,mapping_ext)      
!
!.....Deallocate the local variables for making subdomain array
!
      DEALLOCATE(num_neigh)
      DEALLOCATE(irecv,jsend)
!
      END SUBROUTINE subdomain_info
