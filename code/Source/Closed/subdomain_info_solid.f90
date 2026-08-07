!
      SUBROUTINE subdomain_info_solid(nelem)
!
!     This routine defines local arrays for each subdomain.  All the 
!     arrays can be automatically defined based on the coloring (manually
!     or METIS). It includes Window METIS (Ver.5.0) with having been 
!     activated in preprocessor option (metis=.true.)
!
      USE Zinterface
      USE Zmpi      , ONLY: celem_c,ncell_ps,niut_c,iut_c,ri_c,si_c,rintf_c,sintf_c, &
                            maxmt_ncond,maxmt_ps,                                    &
                            maxmt_nncond,                                            &
                            jperm_c,jjperm_c,                                        &
                            ncell_fluid1_c,ncell_fluid1_dsp_c,                       &
                            ncell_ndim_sz_c,ncell_ndim_dsp_c,                        &
                            ncell_csr_sz_c,ncell_csr_dsp_c,                          &
                            metis
      USE Zzone     , ONLY: ncell_cond,ncell_cond_all
      USE Zcore     , ONLY: np,myrank
      USE Zparam    , ONLY: ns ,ndim
      USE Zbc_index , ONLY: nbcon_c
      USE Znum_cell , ONLY: neigh_c,i_neigh_c_tmp,j_neigh_c_tmp,j_nbcon_c_tmp, &
                            i_neigh_c,index_sort_c
      USE Znormal   , ONLY: nji_c
      USE Zio_unit  , ONLY: unit_log
!
      IMPLICIT NONE      
!
!.....Input
      INTEGER :: nelem
!     local variables
      INTEGER :: i,j,k,jj,ie,ne,proc
      INTEGER :: j0,j1,j2,j3,k1
      INTEGER :: ip,jp,iptr
      INTEGER :: nelem_sub
      INTEGER :: celem_max
      INTEGER :: cinter,cext0
!.....Local arrays
      INTEGER temp_sort(ns)
      INTEGER :: iutjp(np)
      INTEGER :: cext(np)
      INTEGER :: irecv_cnt(np),jsend_cnt(np)
      INTEGER :: icount(np)
      INTEGER :: flag(nelem),flagt(nelem)
      INTEGER :: ia(np+1),ja(nelem)
      INTEGER :: jperms(nelem)
      INTEGER :: index_sort0(maxmt_nncond)
!.....Local allocatable arrays
      INTEGER,DIMENSION(:),ALLOCATABLE :: num_neigh_c,ia_sub,ja_sub
      INTEGER,DIMENSION(:),ALLOCATABLE :: irecv,jsend
!
      IF(myrank.eq.0) THEN
!
!.....Delete zero entries.
!
         ALLOCATE(ia_sub(nelem+1))
         ia_sub(1)=1
         k=1
         DO i=1,nelem
            DO j=i_neigh_c_tmp(i),i_neigh_c_tmp(i+1)-1
               IF(j_nbcon_c_tmp(j).eq.0)THEN
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
            DO j=i_neigh_c_tmp(i),i_neigh_c_tmp(i+1)-1
               IF(j_nbcon_c_tmp(j).eq.0)THEN
                  k=k+1
                  ja_sub(k)=j_neigh_c_tmp(j)
               ENDIF
            ENDDO
            celem_max=max(celem_c(i),celem_max)
         ENDDO
      ENDIF
      IF(np.gt.1) CALL broadcast_i1(celem_max)
!
      IF(np.eq.1) celem_c=1
!
!.....Window METIS Ver.5.0 API
!
      IF(np.gt.1) then
         IF(metis.eq.1) THEN
            IF(myrank.eq.0) THEN
               CALL domain_decomposition(nelem,nelem_sub,ia_sub,ja_sub,celem_c,celem_max)
               WRITE(*,*) '          3The number of subdomains is', celem_max
               WRITE(*,*) '          3The number of processors is', np
               WRITE(unit_log,*) '          3The number of subdomains is', celem_max
               WRITE(unit_log,*) '          3The number of processors is', np         
            ENDIF   
         ELSEIF(metis.eq.0) THEN
            IF(np.ne.celem_max) THEN
               IF(myrank.eq.0) THEN
                  WRITE(*,*) '          4The number of subdomains is', celem_max
                  WRITE(*,*) '          4The number of processors is', np
                  PRINT *,'          4The number of subdomains differs from the number of processors !!!'
                  WRITE(unit_log,*) '          4The number of subdomains is', celem_max
                  WRITE(unit_log,*) '          4The number of processors is', np
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      IF(np.gt.1) CALL broadcast_i1(nelem_sub)
      IF(myrank.ne.0) THEN
         ALLOCATE(celem_c(nelem))
         ALLOCATE(i_neigh_c_tmp(ncell_cond_all+1))
         ALLOCATE(j_neigh_c_tmp(maxmt_nncond))
         ALLOCATE(j_nbcon_c_tmp(maxmt_nncond))
         ALLOCATE(ia_sub(nelem+1))
         ALLOCATE(ja_sub(nelem_sub))
      ENDIF
      IF(np.gt.1) THEN
         CALL broadcast_i(celem_c,nelem)
         CALL broadcast_i(i_neigh_c_tmp,ncell_cond_all+1)
         CALL broadcast_i(j_neigh_c_tmp,maxmt_nncond)
         CALL broadcast_i(j_nbcon_c_tmp,maxmt_nncond)
         CALL broadcast_i(ia_sub,nelem+1) 
         CALL broadcast_i(ja_sub,nelem_sub)
      ENDIF
!
      ALLOCATE(ncell_fluid1_dsp_c(np),ncell_fluid1_c(np))
      ALLOCATE(ncell_ndim_sz_c(np),ncell_ndim_dsp_c(np))
      ALLOCATE(ncell_csr_sz_c(np),ncell_csr_dsp_c(np))
!
!.....Build index based on proc id to access global cell array directly
!.....Get ncell_cond for all procs
!
      ia(1)=1
      DO ip=1,np
         ia(ip+1)=0
         icount(ip)=0
      ENDDO
      DO ie=1,nelem
         ip=celem_c(ie)
         ia(ip+1)=ia(ip+1)+1
      ENDDO
      DO ip=1,np
         ncell_fluid1_c(ip)=ia(ip+1)
         ia(ip+1)=ia(ip+1)+ia(ip)
      ENDDO
      DO ie=1,nelem
         ip=celem_c(ie)
         ja(ia(ip)+icount(ip))=ie
         icount(ip)=icount(ip)+1
      ENDDO
!
!.....Count interior cells, external=receive cells
!.....Count cell only once on external cells
!
      ip=myrank+1
      ncell_cond=ncell_fluid1_c(ip)
      ALLOCATE(jperm_c(ncell_cond))
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
            jp=celem_c(k)
            IF(jp.ne.ip) THEN
!..............Cext cells
               DO j1=j,ia_sub(ie+1)-1
                  k=ja_sub(j1)
                  jp=celem_c(k)
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
         jperm_c(cinter)=ie
         jperms(ie)=cinter
100      CONTINUE 
      ENDDO
      cext0=0
      DO jp=1,np
         cext0=cext0+irecv_cnt(jp)
      ENDDO
      ncell_ps=ncell_cond+cext0
!
!.....Get all cext needed to build ALLGATHERV parameters 
!
      CALL allgather_i(cext0,cext)
!
!.....Receive cells
!
!.....Remove zero count
      niut_c=0
      DO jp=1,np
         IF(irecv_cnt(jp).ne.0) niut_c=niut_c+1
      ENDDO
      ALLOCATE(iut_c(niut_c),ri_c(niut_c+1),si_c(niut_c+1))
      ri_c(1)=1
      j=0
      DO jp=1,np
         IF(irecv_cnt(jp).ne.0) THEN
            j=j+1
            iut_c(j)=jp
            iutjp(jp)=j
            ri_c(j+1)=irecv_cnt(jp)
         ENDIF
      ENDDO
      DO j=1,niut_c
         ri_c(j+1)=ri_c(j+1)+ri_c(j)
      ENDDO
      ALLOCATE(irecv(ri_c(niut_c+1)-1))
      DO jp=1,np
         icount(jp)=0
      ENDDO
!
!.....Minimize the zeroing of flag via flagt
!
      DO i=1,cext0
         k=flagt(i)
         jp=celem_c(k)
         j=iutjp(jp)
         irecv(ri_c(j)+icount(jp))=k
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
                     j0=celem_c(k)
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
!
!...........Minimize the zeroing of flag via flagt
!
            DO i=1,jsend_cnt(ip)
               k=flagt(iptr+i)
               flag(k)=0
            ENDDO
            iptr=iptr+jsend_cnt(ip)
         ENDIF
      ENDDO
!
!.....Send cells
!
      si_c(1)=1
      j=0
      DO jp=1,np
         IF(jsend_cnt(jp).ne.0) THEN
            j=j+1
            si_c(j+1)=jsend_cnt(jp)
         ENDIF
      ENDDO
      DO j=1,niut_c
         si_c(j+1)=si_c(j+1)+si_c(j)
      ENDDO
!
      iut_c=iut_c-1
!
!.....Get all neighbors counted once jsend
!
      ALLOCATE(jsend(si_c(niut_c+1)-1))
      DO i=1,si_c(niut_c+1)-1
         jsend(i)=flagt(i)
      ENDDO
!
!.....Get the mapping for remote cells counted once.
!
      k=cinter
      DO i=1,si_c(niut_c+1)-1
         k1=jsend(i)
         IF(flag(k1).eq.0) then
            k=k+1
            flag(k1)=1
            jperm_c(k)=k1
            jperms(k1)=k
         ENDIF
      ENDDO
!
!.....Get the mapping for ext cells
!
      DO i=1,ri_c(niut_c+1)-1
         k1=irecv(i)
         jperms(k1)=ncell_cond+i
      ENDDO
!
!.....Convert global numbering to local numbering for send,receive cells
!
      ALLOCATE(rintf_c(ri_c(niut_c+1)-1),sintf_c(si_c(niut_c+1)-1))
      DO i=1,ri_c(niut_c+1)-1
         rintf_c(i)=ncell_cond+i
      ENDDO
      DO i=1,si_c(niut_c+1)-1
         k1=jsend(i)
         sintf_c(i)=jperms(k1)
      ENDDO
!
!.....Get size and displacement to call ALLGATHERV
!
         proc=1
         ncell_fluid1_dsp_c(proc)=0
      DO proc=2,np
         ncell_fluid1_dsp_c(proc)=ncell_fluid1_dsp_c(proc-1)+ncell_fluid1_c(proc-1)
      ENDDO
!
!.....Build jjperm out of local jperm
!
      ALLOCATE(jjperm_c(nelem))
      IF(np.gt.1) THEN
         CALL allgather_vec_i(jperm_c,ncell_cond,jjperm_c,nelem,ncell_fluid1_c,ncell_fluid1_dsp_c)
      ELSE
         DO i=1,nelem
            jjperm_c(i)=jperm_c(i)
         ENDDO
      ENDIF
!
!.....Get size and displacement to call ndim type ALLGATHERV
!
      DO ip=1,np
         ncell_ndim_sz_c(ip)=ncell_fluid1_c(ip)*ndim
      ENDDO
         proc=1
         ncell_ndim_dsp_c(proc)=0
      DO proc=2,np
         ncell_ndim_dsp_c(proc)=ncell_ndim_dsp_c(proc-1)+ncell_ndim_sz_c(proc-1)
      ENDDO
!
!.....Get size and displacement to call csr type ALLGATHERV
!
      DO ip=1,np
         k=0
         DO jj=ia(ip),ia(ip+1)-1 
            ie=ja(jj)
            k=k+(i_neigh_c_tmp(ie+1)-i_neigh_c_tmp(ie))
         ENDDO
         ncell_csr_sz_c(ip)=k
      ENDDO
         proc=1
         ncell_csr_dsp_c(proc)=0
      DO proc=2,np
         ncell_csr_dsp_c(proc)=ncell_csr_dsp_c(proc-1)+ncell_csr_sz_c(proc-1)
      ENDDO
!
!
!.....Define local neigh & num_neigh for fluid region
!.....Define local nji
!
      ALLOCATE(num_neigh_c(ncell_ps))
      ALLOCATE(i_neigh_c(ncell_ps+1))
      i_neigh_c(1)=1
      DO i=1,ncell_cond
         ne=jperm_c(i)
         num_neigh_c(i)=i_neigh_c_tmp(ne+1)-i_neigh_c_tmp(ne)
         i_neigh_c(i+1)=i_neigh_c(i)+(i_neigh_c_tmp(ne+1)-i_neigh_c_tmp(ne))
      ENDDO
      maxmt_ncond=i_neigh_c(ncell_cond+1)-1
      ALLOCATE(index_sort_c(maxmt_ncond))
      ALLOCATE(nji_c(maxmt_ncond))
      ALLOCATE(neigh_c(maxmt_ncond))
      ALLOCATE(nbcon_c(maxmt_ncond))
!
      DO i=1,ncell_cond
         ne=jperm_c(i)
         j0=i_neigh_c_tmp(ne)-1
         j1=i_neigh_c(i)-1
!
!........Sort in the vector order
!
         DO j=i_neigh_c(i),i_neigh_c(i+1)-1
            index_sort_c(j)=j-j1
            IF(j_nbcon_c_tmp(j-j1+j0).eq.0)THEN
               k=j_neigh_c_tmp(j-j1+j0)
               neigh_c(j)=jperms(k)
            ELSE
               IF    (j_nbcon_c_tmp(j-j1+j0).eq.-2)THEN
                  neigh_c(j)=ncell_ps+1
               ELSEIF(j_nbcon_c_tmp(j-j1+j0).eq.-3.or.j_nbcon_c_tmp(j-j1+j0).eq.-4)THEN
                  neigh_c(j)=ncell_ps+2
               ELSEIF(j_nbcon_c_tmp(j-j1+j0).eq.-5.or.j_nbcon_c_tmp(j-j1+j0).eq.-6)THEN
                  neigh_c(j)=ncell_ps+3
               ELSE
                  neigh_c(j)=ncell_ps+4
               ENDIF
            ENDIF
         ENDDO
         CALL sortx_i(neigh_c(j1+1),index_sort_c(j1+1),i_neigh_c(i+1)-i_neigh_c(i))
!
!........Restore neigh to 0 for non compute
!
         DO j=i_neigh_c(i+1)-1,i_neigh_c(i),-1
            if(neigh_c(j).le.ncell_ps) exit
            neigh_c(j)=0
         ENDDO
!........Get sorted nbcon
         DO j=i_neigh_c(i),i_neigh_c(i+1)-1
            k=index_sort_c(j)
            temp_sort(j-j1)=j_nbcon_c_tmp(j0+k)
         ENDDO
         DO j=i_neigh_c(i),i_neigh_c(i+1)-1
            nbcon_c(j)=temp_sort(j-j1)
         ENDDO
      ENDDO
!
!.....Get index_sort global
!
      IF(np.gt.1) THEN
         CALL allgatherv_csr_i(index_sort_c,maxmt_ncond,index_sort0,maxmt_nncond, &
                               ncell_cond_all,maxmt_nncond,i_neigh_c_tmp,1)
      ELSEIF(np.eq.1) THEN
         DO i=1,maxmt_nncond
            index_sort0(i)=index_sort_c(i)
         ENDDO
      ENDIF
!
!.....Get sorted global j_neigh_c_tmp to get nji
!
      DO ne=1,nelem
         j0=i_neigh_c_tmp(ne)-1
         DO j=i_neigh_c_tmp(ne),i_neigh_c_tmp(ne+1)-1
            k=index_sort0(j)
            temp_sort(j-j0)=j_neigh_c_tmp(j0+k)
         ENDDO
         DO j=i_neigh_c_tmp(ne),i_neigh_c_tmp(ne+1)-1
            j_neigh_c_tmp(j)=temp_sort(j-j0)
         ENDDO
         DO j=i_neigh_c_tmp(ne),i_neigh_c_tmp(ne+1)-1
            k=index_sort0(j)
            temp_sort(j-j0)=j_nbcon_c_tmp(j0+k)
         ENDDO
         DO j=i_neigh_c_tmp(ne),i_neigh_c_tmp(ne+1)-1
            j_nbcon_c_tmp(j)=temp_sort(j-j0)
         ENDDO
      ENDDO
!
      DO i=1,ncell_cond
         ne=jperm_c(i)
         j0=i_neigh_c_tmp(ne)-1
         j3=i_neigh_c(i)-1
         DO j=i_neigh_c(i),i_neigh_c(i+1)-1
            IF(nbcon_c(j).eq.0)THEN
               k=j_neigh_c_tmp(j-j3+j0)
               j2=i_neigh_c_tmp(k)-1
               DO j1=i_neigh_c_tmp(k),i_neigh_c_tmp(k+1)-1
                  IF (j_neigh_c_tmp(j1).eq.ne) then
                      nji_c(j)=j1-j2
                      GOTO 200
                  ENDIF
               ENDDO
200            CONTINUE
            ELSE
!              nji_c(j,i)=0
            ENDIF
         ENDDO
      ENDDO
!
      IF(np.gt.1) CALL communicate_1d_c_int(num_neigh_c)
      DO i=ncell_cond+1,ncell_ps
         i_neigh_c(i+1)=i_neigh_c(i)+num_neigh_c(i)
      ENDDO
      maxmt_ps=i_neigh_c(ncell_ps+1)-1
!
!.....Deallocate the local variables for making subdomain array
!
      DEALLOCATE(num_neigh_c)
      DEALLOCATE(irecv,jsend)
      DEALLOCATE(ia_sub,ja_sub)
!
      END SUBROUTINE subdomain_info_solid
!
