!
      SUBROUTINE rv_subdomain_info_fuel_rod(nelem,neigh0)
!
!     This routine defines local arrays for each subdomain.  All the 
!     arrays can be automatically defined based on the coloring (manually
!     or METIS). It includes Window METIS (Ver.5.0) with having been 
!     activated in preprocessor option (metis=.true.)
!
      USE Zrv_mpi   , ONLY: celem_fuel_rod,ncell_fuel_rod_p,niut_fuel_rod,iut_fuel_rod,ri_fuel_rod,si_fuel_rod, &
                            rintf_fuel_rod,sintf_fuel_rod,                                                      &
                            jperm_fuel_rod,jjperm_fuel_rod,                                                     &
                            ncell_fuel1_rod,ncell_fuel1_rod_dsp,                                                &
                             ncell_fuel1_rod_2d,ncell_fuel1_rod_2d_dsp
      USE Zcore     , ONLY: np,myrank
      USE Zrv_ncell , ONLY: ncell_fuel_rod,neigh_fuel_rod
      USE Zrv_hts_2d , ONLY: nr_2d
!
      IMPLICIT NONE      
!
!     input
      INTEGER nelem
      INTEGER neigh0(2,nelem)
!     output
!     local variables
      INTEGER i,j,k,jj,ie,ne,proc
      INTEGER :: j0,j1,k1
      INTEGER ip,jp,cell,iptr
      INTEGER :: nelem_sub
      INTEGER :: cinter,cext0
!     local arrays
      INTEGER :: iutjp(np)
      INTEGER :: irecv_cnt(np),jsend_cnt(np)
      INTEGER :: icount(np)
      INTEGER :: flag(nelem),flagt(nelem)
      INTEGER :: ia(np+1),ja(nelem)
      INTEGER :: jperms(nelem)
      INTEGER :: ia_sub(nelem+1)
!
      INTEGER,ALLOCATABLE :: irecv(:),jsend(:)
      INTEGER,ALLOCATABLE :: ja_sub(:)
!
!.....Delete zero entries.
!
      ia_sub(1)=1
      j=1
      DO cell=1,nelem
         DO ne=1,2
            IF(neigh0(ne,cell).gt.0)THEN
               j=j+1
            ENDIF
         ENDDO
         ia_sub(cell+1)=j
      ENDDO
      nelem_sub=ia_sub(nelem+1)-1
      ALLOCATE(ja_sub(ia_sub(nelem+1)-1))
      j=0
      DO cell=1,nelem
         DO ne=1,2
            IF(neigh0(ne,cell).gt.0)THEN
               j=j+1
               ja_sub(j)=neigh0(ne,cell)
            ENDIF
         ENDDO
      ENDDO
!
      IF(np.eq.1) celem_fuel_rod(:)=1
!
      ALLOCATE(ncell_fuel1_rod(np),ncell_fuel1_rod_dsp(np))
      ALLOCATE(ncell_fuel1_rod_2d(np),ncell_fuel1_rod_2d_dsp(np))
!
!     build index based on proc id to access global cell array directly
!     get ncell_fluid for all procs
!
      ia(1)=1
      DO ip=1,np
         ia(ip+1)=0
         icount(ip)=0
      ENDDO
      DO ie=1,nelem
         ip=celem_fuel_rod(ie)
         ia(ip+1)=ia(ip+1)+1
      ENDDO
      DO ip=1,np
         ncell_fuel1_rod(ip)=ia(ip+1)
         ncell_fuel1_rod_2d(ip)=ia(ip+1)*nr_2d
         ia(ip+1)=ia(ip+1)+ia(ip)
      ENDDO
      DO ie=1,nelem
         ip=celem_fuel_rod(ie)
         ja(ia(ip)+icount(ip))=ie
         icount(ip)=icount(ip)+1
      ENDDO
!
!     count interior cells, external=receive cells
!     count cell only once on external cells
!
      ip=myrank+1
      ncell_fuel_rod=ia(ip+1)-ia(ip)
      ALLOCATE(jperm_fuel_rod(ncell_fuel_rod))
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
            jp=celem_fuel_rod(k)
            IF(jp.ne.ip) THEN
!     cext cells
               DO j1=j,ia_sub(ie+1)-1
                  k=ja_sub(j1)
                  jp=celem_fuel_rod(k)
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
!     cintr cells
!     get the mapping for interior cell
         cinter=cinter+1
         jperm_fuel_rod(cinter)=ie
         jperms(ie)=cinter
100      CONTINUE 
      ENDDO
      cext0=0
      DO jp=1,np
         cext0=cext0+irecv_cnt(jp)
      ENDDO
      ncell_fuel_rod_p=ncell_fuel_rod+cext0
!
!     receive cells
!
!     remove zero count
      niut_fuel_rod=0
      DO jp=1,np
         IF(irecv_cnt(jp).ne.0) niut_fuel_rod=niut_fuel_rod+1
      ENDDO
      ALLOCATE(iut_fuel_rod(niut_fuel_rod),ri_fuel_rod(niut_fuel_rod+1),si_fuel_rod(niut_fuel_rod+1))
      ri_fuel_rod(1)=1
      j=0
      DO jp=1,np
         IF(irecv_cnt(jp).ne.0) THEN
            j=j+1
            iut_fuel_rod(j)=jp
            iutjp(jp)=j
            ri_fuel_rod(j+1)=irecv_cnt(jp)
         ENDIF
      ENDDO
      DO j=1,niut_fuel_rod
         ri_fuel_rod(j+1)=ri_fuel_rod(j+1)+ri_fuel_rod(j)
      ENDDO
      ALLOCATE(irecv(ri_fuel_rod(niut_fuel_rod+1)-1))
      DO jp=1,np
         icount(jp)=0
      ENDDO
!     minimize the zeroing of flag via flagt
      DO i=1,cext0
         k=flagt(i)
         jp=celem_fuel_rod(k)
         j=iutjp(jp)
         irecv(ri_fuel_rod(j)+icount(jp))=k
         icount(jp)=icount(jp)+1
         flag(k)=0
      ENDDO
!
!     count remote=send cells
!     count cell only once on remote cells for each ip
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
                     j0=celem_fuel_rod(k)
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
!
!     send cells
!
      si_fuel_rod(1)=1
      j=0
      DO jp=1,np
         IF(jsend_cnt(jp).ne.0) THEN
            j=j+1
            si_fuel_rod(j+1)=jsend_cnt(jp)
         ENDIF
      ENDDO
      DO j=1,niut_fuel_rod
         si_fuel_rod(j+1)=si_fuel_rod(j+1)+si_fuel_rod(j)
      ENDDO
!
      iut_fuel_rod=iut_fuel_rod-1
!
!     get all neighbors counted once jsend
!
      ALLOCATE(jsend(si_fuel_rod(niut_fuel_rod+1)-1))
      DO i=1,si_fuel_rod(niut_fuel_rod+1)-1
         jsend(i)=flagt(i)
      ENDDO
!
!     get the mapping for remote cells counted once.
!
      k=cinter
      DO i=1,si_fuel_rod(niut_fuel_rod+1)-1
         k1=jsend(i)
         IF(flag(k1).eq.0) then
            k=k+1
            flag(k1)=1
            jperm_fuel_rod(k)=k1
            jperms(k1)=k
         ENDIF
      ENDDO
!
!     get the mapping for ext cells
!
      DO i=1,ri_fuel_rod(niut_fuel_rod+1)-1
         k1=irecv(i)
         jperms(k1)=ncell_fuel_rod+i
      ENDDO
!
!     convert global numbering to local numbering for send,receive cells
!
      ALLOCATE(rintf_fuel_rod(ri_fuel_rod(niut_fuel_rod+1)-1),sintf_fuel_rod(si_fuel_rod(niut_fuel_rod+1)-1))
      DO i=1,ri_fuel_rod(niut_fuel_rod+1)-1
         rintf_fuel_rod(i)=ncell_fuel_rod+i
      ENDDO
      DO i=1,si_fuel_rod(niut_fuel_rod+1)-1
         k1=jsend(i)
         sintf_fuel_rod(i)=jperms(k1)
      ENDDO
!
!     get size and displacement to call ALLGATHERV
!
         proc=1
         ncell_fuel1_rod_dsp(proc)=0
         ncell_fuel1_rod_2d_dsp(proc)=0
      DO proc=2,np
         ncell_fuel1_rod_dsp(proc)=ncell_fuel1_rod_dsp(proc-1)+ncell_fuel1_rod(proc-1)
         ncell_fuel1_rod_2d_dsp(proc)=ncell_fuel1_rod_2d_dsp(proc-1)+ncell_fuel1_rod_2d(proc-1)
      ENDDO
!
!     build jjperm out of local jperm
!
      ALLOCATE(jjperm_fuel_rod(nelem))
      IF(np.gt.1) THEN
         CALL allgather_vec_i(jperm_fuel_rod,ncell_fuel_rod,jjperm_fuel_rod,nelem,ncell_fuel1_rod,ncell_fuel1_rod_dsp)
      ELSE
         DO i=1,nelem
            jjperm_fuel_rod(i)=jperm_fuel_rod(i)
         ENDDO
      ENDIF

!
!.....Writing array for CUPID calculation
!
      ALLOCATE(neigh_fuel_rod(ncell_fuel_rod_p,2))
!
!.....Define local neigh
!
      DO i=1,ncell_fuel_rod
         ne=jperm_fuel_rod(i)
         DO j=1,2
            k=neigh0(j,ne) 
            IF(neigh0(j,ne).gt.0)THEN
               neigh_fuel_rod(i,j)=jperms(k)
            ELSE
               neigh_fuel_rod(i,j)=0
            ENDIF
         ENDDO
      ENDDO
!
!.....Deallocate the local variables for making subdomain array
!
      DEALLOCATE(irecv,jsend)
      DEALLOCATE(ja_sub)
!
      END SUBROUTINE rv_subdomain_info_fuel_rod
!
