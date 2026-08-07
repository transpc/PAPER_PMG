!      
      SUBROUTINE vectorize_index_solid
!
!     This routine prepares indeces for vectorization
!
      USE Zmpi             , ONLY: maxmt_ncond
      USE Zbc_index        , ONLY: nbcon_c
      USE Zcore            , ONLY: myrank
      USE Zparam           , ONLY: nb_max
      USE Zvec_param       , ONLY: nfc_tot,nfc_tot1,nfc_nonk,nfc_non,nfc_fsw,nfc_ctw,nfc_chw, &
                                   nfc_tot_nb1,nfc_chtcw
      USE Znormal          , ONLY: nji_c      
      USE Znum_cell        , ONLY: i_neigh_c,neigh_c,                       &
                                   istartc_nf,istartc_nb1,                  &
                                   iac_nb,rightc_nb_k,icellc_nb,iptrc_nb_k, &
                                   nbc_nf
      USE Zzone            , ONLY: ncell_cond
      USE Zvec_index_solid , ONLY: nbcon_solid_ctw,nbcon_solid_chw,nbcon_solid_chtcw
      USE Zvec_index_solid , ONLY: left_solid_nf,jneigh_solid_nf,right_solid_non,right_solid_fsw
      USE Zvec_index_solid , ONLY: left_solid_k,right_solid_k
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,j,k,ii,j0,k0,jk,nb
      INTEGER :: jj
      INTEGER :: k_non
      INTEGER :: nf_number,istart,len,i1,idisp
      INTEGER :: istart1
      INTEGER :: icount_kt
!.....Local arrays
      INTEGER :: icount(-1:nbc_nf),icount_nb(-1:nbc_nf)
      INTEGER :: neighji_c(maxmt_ncond),neighji_c_nf(maxmt_ncond)
!
!.....Caculate vector size and index
!
      nfc_non=0 
      nfc_fsw=0 
      nfc_ctw=0 
      nfc_chw=0 
      nfc_chtcw=0      
      nfc_nonk=0 
      DO nf_number=-1,nbc_nf
         istartc_nb1(2,nf_number)=0
      ENDDO
!      
      DO i=1,ncell_cond
         DO nf_number=-1,nbc_nf
            icount(nf_number)=0
         ENDDO
         DO j=i_neigh_c(i),i_neigh_c(i+1)-1
            neighji_c_nf(j)=99
            IF(nbcon_c(j).eq.0)THEN
               nf_number=0
               k=neigh_c(j)
               IF(k.gt.i) THEN
                  nfc_non=nfc_non+1
                  neighji_c(j)=nfc_non
                  icount(nf_number)=icount(nf_number)+1
               ELSE
                  nfc_nonk=nfc_nonk+1
                  k0=i_neigh_c(k)-1
                  jk=nji_c(j)
                  neighji_c(j)=-neighji_c(jk+k0)
                  icount(-1)=icount(-1)+1
               ENDIF
               neighji_c_nf(j)=nf_number
            ELSEIF(nbcon_c(j).eq.-2)THEN
               nfc_fsw=nfc_fsw+1
               neighji_c(j)=nfc_fsw
               nf_number=1
               icount(nf_number)=icount(nf_number)+1
               neighji_c_nf(j)=nf_number
            ELSEIF(nbcon_c(j).eq.-3.or.nbcon_c(j).eq.-4)THEN
               nfc_ctw=nfc_ctw+1
               neighji_c(j)=nfc_ctw
               nf_number=2
               icount(nf_number)=icount(nf_number)+1
               neighji_c_nf(j)=nf_number
            ELSEIF(nbcon_c(j).eq.-5.or.nbcon_c(j).eq.-6)THEN
               nfc_chw=nfc_chw+1
               neighji_c(j)=nfc_chw
               nf_number=3
               icount(nf_number)=icount(nf_number)+1
               neighji_c_nf(j)=nf_number
            ELSEIF(nbcon_c(j).ge.-39 .and. nbcon_c(j).le.-31)THEN
               nfc_chtcw=nfc_chtcw+1
               neighji_c(j)=nfc_chtcw
               nf_number=4
               icount(nf_number)=icount(nf_number)+1
               neighji_c_nf(j)=nf_number
            ELSEIF(nbcon_c(j).eq.-1)THEN
            ELSEIF(nbcon_c(j).ge.1.and.nbcon_c(j).le.nb_max)THEN 
            ELSE
                WRITE(*,*)'Error in solid face identification during vectorization!',nbcon_c(j)
                STOP 99
            ENDIF
         ENDDO
         DO nf_number=-1,nbc_nf
            if(icount(nf_number).gt.0) istartc_nb1(2,nf_number)=istartc_nb1(2,nf_number)+1
         ENDDO
      ENDDO
      nfc_tot=nfc_non+nfc_fsw+nfc_ctw+nfc_chw+nfc_chtcw
      nfc_tot1=nfc_non+nfc_fsw
      goto 100
      IF(myrank.eq.0) THEN
      write(*,*)'ncell_cond',ncell_cond
      if(nfc_non.ne.0) write(*,*)'nfc_non,nfc_nonk ',nfc_non,nfc_nonk
      if(nfc_fsw.ne.0) write(*,*)'nfc_fsw ',nfc_fsw
      if(nfc_ctw.ne.0) write(*,*)'nfc_ctw ',nfc_ctw
      if(nfc_chw.ne.0) write(*,*)'nfc_chw ',nfc_chw
      if(nfc_chtcw.ne.0) write(*,*)'nfc_chtcw ',nfc_chtcw      
      write(*,*)'nfc_tot',nfc_tot
      ENDIF
100   continue
!
      nfc_tot_nb1=0
      DO nf_number=-1,nbc_nf
         IF(istartc_nb1(2,nf_number).ne.0) THEN
            nfc_tot_nb1=nfc_tot_nb1+istartc_nb1(2,nf_number)+1 ! add 1 per entry csr style format
         ENDIF
      ENDDO
      istartc_nb1(1,-1)=0
      DO nf_number=0,nbc_nf
         IF(istartc_nb1(2,nf_number-1).ne.0) THEN
            istartc_nb1(1,nf_number)=istartc_nb1(1,nf_number-1)+istartc_nb1(2,nf_number-1)+1 ! add 1 per entry csr style format
         ELSE
            istartc_nb1(1,nf_number)=istartc_nb1(1,nf_number-1)
         ENDIF
      ENDDO
      goto 105
      DO nf_number=-1,nbc_nf
         if(istartc_nb1(2,nf_number).ne.0) THEN
            write(*,*) 'nb1',nf_number,istartc_nb1(1,nf_number),istartc_nb1(2,nf_number)
         ENDIF
      ENDDO
      write(*,*)'nf_tot_nb',nfc_tot_nb1
105   continue
      ALLOCATE(icellc_nb(nfc_tot_nb1))
      ALLOCATE(iptrc_nb_k(ncell_cond))
      ALLOCATE(iac_nb(nfc_tot_nb1))
      ALLOCATE(rightc_nb_k(nfc_nonk))
!
      DO nf_number=-1,nbc_nf
         icount_nb(nf_number)=0
         IF(istartc_nb1(2,nf_number).ne.0) THEN
            istart=istartc_nb1(1,nf_number)
            i1=istart+1
            iac_nb(i1)=1
         ENDIF
      ENDDO
!
      icount_kt=0
      DO i=1,ncell_cond
         DO nf_number=-1,nbc_nf
            icount(nf_number)=0
         ENDDO
         DO j=i_neigh_c(i),i_neigh_c(i+1)-1
            nf_number=neighji_c_nf(j)
            IF(nf_number.ne.99) THEN
            idisp=neighji_c(j)
            IF(idisp.gt.0) THEN
               icount(nf_number)=icount(nf_number)+1
            ELSEIF(idisp.lt.0) THEN
               icount(-1)=icount(-1)+1
               icount_kt=icount_kt+1
               rightc_nb_k(icount_kt)=-idisp
            endif
            endif
         ENDDO
         DO nf_number=-1,nbc_nf
            istart1=istartc_nb1(1,nf_number)
            IF(icount(nf_number).gt.0) THEN
              ii=icount_nb(nf_number)
              ii=ii+1
              IF(nf_number.eq.-1) then
                 iptrc_nb_k(i)=ii
              ENDIF
              icount_nb(nf_number)=ii
              i1=istart1+ii
              iac_nb(i1+1)=icount(nf_number)
              icellc_nb(i1)=i
            ELSE
              IF(nf_number.eq.-1) then
                 iptrc_nb_k(i)=0
              ENDIF
            ENDIF
         ENDDO
      ENDDO
!
!     DO i=1,20
!        write(*,200) i,(nbcon(j,i),j=1,num_neigh(i))
!     ENDDO
200    format(8i6)
      DO nf_number=-1,nbc_nf
         istart=istartc_nb1(1,nf_number)
         len   =istartc_nb1(2,nf_number)
         DO i=2,len+1
            i1=istart+i
            iac_nb(i1)=iac_nb(i1)+iac_nb(i1-1)
         ENDDO
      ENDDO
!.....Nf arrray pointers
      istartc_nf(2,-1)=nfc_nonk
      istartc_nf(2,0)=nfc_non
      istartc_nf(2,1)=nfc_fsw
      istartc_nf(2,2)=nfc_ctw
      istartc_nf(2,3)=nfc_chw
      istartc_nf(2,4)=nfc_chtcw !pik-chtcw      
      istartc_nf(1,0)=0
      DO i=1,nbc_nf
         istartc_nf(1,i)=istartc_nf(1,i-1)+istartc_nf(2,i-1)
      ENDDO
!
      ALLOCATE(left_solid_k(nfc_nonk))
      ALLOCATE(right_solid_k(nfc_nonk))
      ALLOCATE(left_solid_nf(nfc_tot))
      ALLOCATE(jneigh_solid_nf(nfc_tot))
      ALLOCATE(right_solid_non(nfc_non))
      ALLOCATE(right_solid_fsw(nfc_fsw))
!
      k=0
      DO i=1,ncell_cond
         j0=i_neigh_c(i)-1
         DO j=i_neigh_c(i),i_neigh_c(i+1)-1
            nf_number=neighji_c_nf(j)
            IF(nf_number.ne.99) THEN
               istart=istartc_nf(1,nf_number) 
               ii=neighji_c(j)
               IF(ii.gt.0) THEN
                  i1=istart+ii
                  left_solid_nf(i1)=i
                  jneigh_solid_nf(i1)=j-j0
               ELSE
                  k=k+1
                  left_solid_k(k)=i
               ENDIF
            ENDIF
         ENDDO
      ENDDO
!
      k_non=0
      DO i=1,ncell_cond
         DO j=i_neigh_c(i),i_neigh_c(i+1)-1
            nf_number=neighji_c_nf(j)
            IF(nf_number.ne.99) THEN
               istart=istartc_nf(1,nf_number) 
               ii=neighji_c(j)
               k=neigh_c(j)
               IF(nf_number.eq.0)THEN
                  IF(ii.gt.0)THEN
                     right_solid_non(ii)=k
                  ELSE
                     k_non=k_non+1
                     right_solid_k(k_non)=k
                  ENDIF               
               ELSEIF(nf_number.eq.1)THEN
                  right_solid_fsw(ii)=k
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      DO i=1,ncell_cond
!           write(*,'(20(i5,2x))') i,(neigh_c(j),j=i_neigh_c(i),i_neigh_c(i+1)-1)
!           write(*,'(20(i5,2x))') i,(nbcon_c(j),j=i_neigh_c(i),i_neigh_c(i+1)-1)
      ENDDO
         nf_number=-1
         len   =istartc_nb1(2,nf_number)
         DO nb=1,len
            ii=icellc_nb(nb)
            DO i=iac_nb(nb),iac_nb(nb+1)-1
!           write(*,'(20(i5,2x))') nb,i,ii,left_solid_k(i),right_solid_k(i)
            ENDDO
         ENDDO
!      
!.....Ncon
!
      ALLOCATE(nbcon_solid_ctw(nfc_ctw))
      ALLOCATE(nbcon_solid_chw(nfc_chw))
      ALLOCATE(nbcon_solid_chtcw(nfc_chtcw))      
!      
      nf_number=2
      istart=istartc_nf(1,nf_number)
      len   =istartc_nf(2,nf_number)
      DO i=1,len
         i1=istart+i
         ii=left_solid_nf(i1)
         jj=jneigh_solid_nf(i1)
         j0=i_neigh_c(ii)-1
         nbcon_solid_ctw(i)=nbcon_c(jj+j0)
      ENDDO           
!
      nf_number=3
      istart=istartc_nf(1,nf_number)
      len   =istartc_nf(2,nf_number)
      DO i=1,len
         i1=istart+i
         ii=left_solid_nf(i1)
         jj=jneigh_solid_nf(i1)
         j0=i_neigh_c(ii)-1
         nbcon_solid_chw(i)=nbcon_c(jj+j0)
      ENDDO      
!pik-chtcw
      nf_number=4 
      istart=istartc_nf(1,nf_number)
      len   =istartc_nf(2,nf_number)
      DO i=1,len
         i1=istart+i
         ii=left_solid_nf(i1)
         jj=jneigh_solid_nf(i1)
         j0=i_neigh_c(ii)-1
         nbcon_solid_chtcw(i)=nbcon_c(jj+j0)
      ENDDO      
!      
      END SUBROUTINE vectorize_index_solid      
      
