      SUBROUTINE vectorize_index
!
!     This routine prepares indeces for vectorization
!
! Summation
! The structure of the arrays is identical (ns,ncell_fluid) like nbcon
! with the differrence we have  from 0:11 arrays and add the k array pointed
! with nf_number=-1
! rather than using the 2 dimension array with ns,and since the number
! of neighbors is variable we are using CSR style format via the buildup
! of ia_nb, then for each nf_number except the k we do  have the
! elements to be summed nicely sitting stride 1 from
!    ia_nb(i)=>ia(nb(i+1)-1   and icell_nb holds the cell number
!    equivalent to left_nf
!    once we get the partial sum for each cell and for each nf_number
!    then we have to do a global sum on that cell number.   
!    nf_totnb1 size of ia_nb have 1 additional element for each nf_number as
!    compared to nf_tot_nb that addresses icell_nb
!
!
!
!     nf_number -1:11                     icell_nb   
!     istart_nb(2,nf_number)
!     ia_nb ---->   x0                             
!           --|     x0                      i    
!             |->   x1
!
!
!     nf_totn_b1                           nf_tot_nb 
!
!     icell_nb_indx contains cell number where elements summation > 1 from ia_nb(i) to ia_nb(i+1)-1
!
      USE Zmpi          , ONLY: ncell_fp,maxmt_fluid,ia_nrhs
      USE Zzone         , ONLY: ncell_fluid
      USE Zcore         , ONLY: myrank
      USE Znormal       , ONLY: nji,i_neigh_nbcon0
      USE Zbc_index     , ONLY: nbcon,npb
      USE Znum_cell     , ONLY: nb_nf,istart_nf,istart_nbcon_nf,istart_svp_nf, &
                                ia_nb,icell_nb,iptr_nb_k,                      &
                                right_nb_k,                                    &
                                istart_nb1,                                    &
                                i_neigh,neigh
      USE Zparam        , ONLY: ndim,nin_max,nb_max
      USE Zvec_param    , ONLY: nf_nonk,nf_non,nf_mcc,nf_inl,nf_out,nf_adw,nf_fsw,nf_ctw,nf_chw,nf_sym, &
                                nf_hconv,nf_hcond,nf_hvad,                                              &
                                nf_totk,nf_tot,nf_tot_nbcon,nf_tot_svp,nf_nbcon_change,                 &
                                nf_flux,nf_fluxk,nf_flux1,nf_fluxk1,nf_flux2,nf_fluxk2,                 &
                                nf_tot_nb1
      USE Zvec_index    , ONLY: left_nf,jneigh_nf,nbcon_nf, &
                                right_non,right_fsw,        &
                                kneigh_non,                 &
                                delta_npb
      USE Zvector       , ONLY: vl_f_non,vg_f_non
      USE Zscalar_coeff , ONLY: sfg_nf,sfl_nf,sfd_nf,           &
                                sfg_non_k,sfl_non_k,sfd_non_k,  &
                                sfg6_nf,sfl6_nf,sfd6_nf,        &
                                sfg6_non_k,sfl6_non_k,sfd6_non_k
      USE Zbc_index    , ONLY: iface_wall,iface_wall0
      USE Zuserdefined , ONLY: vel_bc_profile_inl
      USE Zmars        , ONLY: mass_nf_mcc      
!
      IMPLICIT NONE
!
      INTEGER :: nf_number,istart,len,i1,idisp
      INTEGER :: istart1,istart2,i2
      INTEGER :: i,j,k,ii,j0,k0
      INTEGER :: jj,jk,nf_face
      INTEGER :: i_non,i_fsw
      INTEGER :: icount_kt
      INTEGER :: icount(-1:nb_nf),icount_nb(-1:nb_nf)
      INTEGER,DIMENSION(maxmt_fluid) :: neighji,neighji_nf
!
!.....Build neighji mapping of cell j,i to face loop.  nf_non ii>0 kk<0
!     nf_non=0,nf_mcc=1,nf_inl=2,nf_out=3,nf_adw=4,nf_fsw=5,nf_ctw=6,nf_chw=7,nf_sym=8
!
      nf_nonk=0
      nf_non=0
      nf_inl=0
      nf_out=0
      nf_adw=0
      nf_fsw=0
      nf_ctw=0
      nf_chw=0
      nf_sym=0
      nf_mcc=0
      nf_hconv=0
      nf_hcond=0
      nf_hvad=0
      nf_face=0
      DO nf_number=-1,nb_nf
         istart_nb1(2,nf_number)=0
      ENDDO
      DO i=1,ncell_fluid
         DO nf_number=-1,nb_nf
            icount(nf_number)=0
         ENDDO
         nf_face=nf_face+(i_neigh(i+1)-i_neigh(i))
         j0=i_neigh_nbcon0(i)
         DO j=i_neigh(i),i_neigh(i+1)-1
            IF(nbcon(j).eq.0)THEN
               k=neigh(j)
               k0=i_neigh(k)-1
               IF(k.gt.i) then
                  nf_non=nf_non+1
                  neighji(j)=nf_non
               ELSE
                  nf_nonk=nf_nonk+1
                  jk=nji(j0)
                  j0=j0+1
                  neighji(j)=-neighji(jk+k0)
               ENDIF
               nf_number=0
            ELSEIF(nbcon(j).ge.201)THEN
               nf_mcc=nf_mcc+1
               neighji(j)=nf_mcc
               nf_number=1
            ELSEIF(nbcon(j).gt.0.and.nbcon(j).le.nin_max)THEN
               nf_inl=nf_inl+1
               neighji(j)=nf_inl
               nf_number=2
            ELSEIF(nbcon(j).gt.nin_max.and.nbcon(j).le.nb_max)THEN
               nf_out=nf_out+1
               neighji(j)=nf_out
               nf_number=3
            ELSEIF(nbcon(j).eq.-1)THEN
               nf_adw=nf_adw+1
               neighji(j)=nf_adw
               nf_number=4
            ELSEIF(nbcon(j).eq.-2)THEN
               nf_fsw=nf_fsw+1
               neighji(j)=nf_fsw
               nf_number=5
            ELSEIF(nbcon(j).eq.-3.or.nbcon(j).eq.-4)THEN
               nf_ctw=nf_ctw+1
               neighji(j)=nf_ctw
               nf_number=6
            ELSEIF(nbcon(j).eq.-5.or.nbcon(j).eq.-6)THEN
               nf_chw=nf_chw+1
               neighji(j)=nf_chw
               nf_number=7
            ELSEIF(nbcon(j).eq.101)THEN
               nf_sym=nf_sym+1
               neighji(j)=nf_sym
               nf_number=8
            ELSEIF(nbcon(j).ge.-39 .and. nbcon(j).le.-31)THEN
               nf_hconv=nf_hconv+1
               neighji(j)=nf_hconv
               nf_number=9
            ELSEIF(nbcon(j).ge.-49 .and. nbcon(j).le.-41)THEN
               nf_hcond=nf_hcond+1
               neighji(j)=nf_hcond
               nf_number=10
            ELSEIF(nbcon(j).ge.-59 .and. nbcon(j).le.-51)THEN
               nf_hvad=nf_hvad+1
               neighji(j)=nf_hvad
               nf_number=11
            ELSE
               WRITE(*,*)'Error during vectorization=>unknown nbcon!' ,k,i,nbcon(j)
               STOP 99
            ENDIF
            neighji_nf(j)=nf_number
            IF(neighji(j).gt.0) THEN
               icount(nf_number)=icount(nf_number)+1
            ELSE
               icount(-1)=icount(-1)+1
            ENDIF
         ENDDO
         DO nf_number=-1,nb_nf
            if(icount(nf_number).gt.0) istart_nb1(2,nf_number)=istart_nb1(2,nf_number)+1
         ENDDO
      ENDDO
      nf_tot=nf_non+nf_mcc+nf_inl+nf_out+nf_adw+nf_fsw+nf_ctw+nf_chw+nf_sym+nf_hconv+nf_hcond+nf_hvad
      nf_totk=nf_nonk+nf_tot
      nf_flux=nf_non+nf_mcc+nf_inl+nf_out
      nf_fluxk=nf_nonk+nf_flux
      nf_flux1=nf_non+nf_mcc+nf_inl
      nf_fluxk1=nf_nonk+nf_flux1
      nf_flux2=nf_non+nf_mcc
      nf_fluxk2=nf_nonk+nf_flux2
      nf_nbcon_change=nf_inl+nf_out+nf_adw
      nf_tot_nbcon=nf_mcc+nf_inl+nf_out+nf_adw+nf_fsw+nf_ctw+nf_chw+nf_hconv+nf_hcond+nf_hvad
      nf_tot_svp=nf_non+nf_mcc+nf_inl+nf_out+nf_adw+nf_sym
!     nf arrray pointers to access nbcon , index_nbcon
      istart_nf(2,0)=0
      istart_nf(2,1)=nf_mcc
      istart_nf(2,2)=nf_inl
      istart_nf(2,3)=nf_out
      istart_nf(2,4)=nf_adw
      istart_nf(2,5)=nf_fsw
      istart_nf(2,6)=nf_ctw
      istart_nf(2,7)=nf_chw
      istart_nf(2,8)=nf_hconv
      istart_nf(2,9)=nf_hcond
      istart_nf(2,10)=nf_hvad
      istart_nf(2,11)=0
      istart_nbcon_nf(0)=0
      DO i=1,nb_nf
         istart_nbcon_nf(i)=istart_nbcon_nf(i-1)+istart_nf(2,i-1)
      ENDDO
!        write(*,200) myrank,(istart_nbcon_nf(i),i=0,8)
!     nf arrray pointers to access svp 
      istart_nf(2,0)=nf_non
      istart_nf(2,1)=nf_mcc
      istart_nf(2,2)=nf_inl
      istart_nf(2,3)=nf_out
      istart_nf(2,4)=nf_adw
      istart_nf(2,5)=0
      istart_nf(2,6)=0
      istart_nf(2,7)=0
      istart_nf(2,8)=nf_sym
      istart_nf(2,9)=nf_hconv
      istart_nf(2,10)=nf_hcond
      istart_nf(2,11)=nf_hvad
      istart_svp_nf(0)=0
      DO i=1,nb_nf
         istart_svp_nf(i)=istart_svp_nf(i-1)+istart_nf(2,i-1)
      ENDDO
!.....nf arrray pointers
      istart_nf(2,-1)=nf_nonk
      istart_nf(2,0)=nf_non
      istart_nf(2,1)=nf_mcc
      istart_nf(2,2)=nf_inl
      istart_nf(2,3)=nf_out
      istart_nf(2,4)=nf_adw
      istart_nf(2,5)=nf_fsw
      istart_nf(2,6)=nf_ctw
      istart_nf(2,7)=nf_chw
      istart_nf(2,8)=nf_sym
      istart_nf(2,9)=nf_hconv
      istart_nf(2,10)=nf_hcond
      istart_nf(2,11)=nf_hvad
      istart_nf(1,0)=0
      DO i=1,nb_nf
         istart_nf(1,i)=istart_nf(1,i-1)+istart_nf(2,i-1)
      ENDDO
      goto 100
!     IF(myrank.eq.0) THEN
      write(*,*) myrank,'ncell_fluid',ncell_fluid
      if(nf_non  .ne.0) write(*,*) myrank,'nf_non ',nf_non,nf_nonk,istart_nf(1,0)
      if(nf_mcc  .ne.0) write(*,*) myrank,'nf_mcc ',nf_mcc,istart_nf(1,1)
      if(nf_inl  .ne.0) write(*,*) myrank,'nf_inl ',nf_inl,istart_nf(1,2)
      if(nf_out  .ne.0) write(*,*) myrank,'nf_out ',nf_out,istart_nf(1,3)
      if(nf_adw  .ne.0) write(*,*) myrank,'nf_adw ',nf_adw,istart_nf(1,4)
      if(nf_fsw  .ne.0) write(*,*) myrank,'nf_fsw ',nf_fsw,istart_nf(1,5)
      if(nf_ctw  .ne.0) write(*,*) myrank,'nf_ctw ',nf_ctw,istart_nf(1,6)
      if(nf_chw  .ne.0) write(*,*) myrank,'nf_chw ',nf_chw,istart_nf(1,7)
      if(nf_sym  .ne.0) write(*,*) myrank,'nf_sym ',nf_sym,istart_nf(1,8)
      if(nf_hconv.ne.0) write(*,*) myrank,'nf_hconv ',nf_hconv,istart_nf(1,9)
      if(nf_hcond.ne.0) write(*,*) myrank,'nf_hcond ',nf_hcond,istart_nf(1,10)
      if(nf_hvad .ne.0)  write(*,*) myrank,'nf_hvad ',nf_hvad,istart_nf(1,11)
      write(*,*)'nf_tot',nf_tot,nf_face,nf_flux,nf_nbcon_change,nf_tot_nbcon
!     DO i=1,ncell_fluid
!        if(i.eq.39) then
!        write(*,110) i,(neigh(j,i),j=i_neigh(i),i_neigh(i+1)-1)
!        write(*,110) i,(nbcon(j),j=i_neigh(i),i_neigh(i+1)-1)
!        write(*,110) i,(neighji_nf(j),j=i_neigh(i),i_neigh(i+1)-1)
!        endif
!     ENDDO
!     ENDIF
110   format(20(i6,2x))
100   continue
!.....Build csr pointer for scalar matrix to access non,mcc,inl
      ALLOCATE(ia_nrhs(ncell_fluid))
      DO i=1,ncell_fluid
         i2=0
         DO j=i_neigh(i),i_neigh(i+1)-1
            nf_number=neighji_nf(j)
            IF(nf_number.gt.2) exit
            i2=i2+1
         ENDDO
         ia_nrhs(i)=i2
      ENDDO
!
      nf_tot_nb1=0
      DO nf_number=-1,nb_nf
         IF(istart_nb1(2,nf_number).ne.0) THEN
            nf_tot_nb1=nf_tot_nb1+istart_nb1(2,nf_number)+1 ! add 1 per entry csr style format
         ENDIF
      ENDDO
      istart_nb1(1,-1)=0
      DO nf_number=0,nb_nf
         IF(istart_nb1(2,nf_number-1).ne.0) THEN
            istart_nb1(1,nf_number)=istart_nb1(1,nf_number-1)+istart_nb1(2,nf_number-1)+1 ! add 1 per entry csr style format
         ELSE
            istart_nb1(1,nf_number)=istart_nb1(1,nf_number-1)
         ENDIF
      ENDDO
      goto 105
      IF(myrank.eq.0) THEN
      DO nf_number=-1,nb_nf
         if(istart_nb1(2,nf_number).ne.0) THEN
            write(*,*) 'nb1',nf_number,istart_nb1(1,nf_number),istart_nb1(2,nf_number)
         ENDIF
      ENDDO
      write(*,*)'nf_tot_nb',nf_tot_nb1
      ENDIF
105   continue
      ALLOCATE(icell_nb(nf_tot_nb1))
      ALLOCATE(iptr_nb_k(ncell_fluid))
      ALLOCATE(ia_nb(nf_tot_nb1))
      ALLOCATE(right_nb_k(nf_nonk))
!
      DO nf_number=-1,nb_nf
         icount_nb(nf_number)=0
         IF(istart_nb1(2,nf_number).ne.0) THEN
            istart=istart_nb1(1,nf_number)
            i1=istart+1
            ia_nb(i1)=1
         ENDIF
      ENDDO
      icount_kt=0
      DO i=1,ncell_fluid
         DO nf_number=-1,nb_nf
            icount(nf_number)=0
         ENDDO
         DO j=i_neigh(i),i_neigh(i+1)-1
            nf_number=neighji_nf(j)
            idisp=neighji(j)
            IF(idisp.gt.0) THEN
               icount(nf_number)=icount(nf_number)+1
            ELSEIF(idisp.lt.0) THEN
               icount(-1)=icount(-1)+1
               icount_kt=icount_kt+1
               right_nb_k(icount_kt)=-idisp
            endif
         ENDDO
         DO nf_number=-1,nb_nf
            istart1=istart_nb1(1,nf_number)
            IF(icount(nf_number).gt.0) THEN
              ii=icount_nb(nf_number)
              ii=ii+1
              IF(nf_number.eq.-1) then
                 iptr_nb_k(i)=ii
              ENDIF
              icount_nb(nf_number)=ii
              i1=istart1+ii
              ia_nb(i1+1)=icount(nf_number)
              icell_nb(i1)=i
            ELSE
              IF(nf_number.eq.-1) then
                 iptr_nb_k(i)=0
              ENDIF
            ENDIF
         ENDDO
      ENDDO
      DO nf_number=-1,nb_nf
         istart1=istart_nb1(1,nf_number)
         len    =istart_nb1(2,nf_number)
         IF(len.ne.0) THEN
            i1=istart1+len+1
            icell_nb(i1)=0
         ENDIF
      ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!
!     DO i=1,20
!        write(*,200) i,(nbcon(j),j=i_neigh(i),i_neigh(i+1)-1)
!     ENDDO
200    format(8i6)
      DO nf_number=-1, nb_nf
         istart=istart_nb1(1,nf_number)
         len   =istart_nb1(2,nf_number)
         DO i=2,len+1
            i1=istart+i
            ia_nb(i1)=ia_nb(i1)+ia_nb(i1-1)
         ENDDO
      ENDDO
!
      ALLOCATE(left_nf(nf_tot),jneigh_nf(nf_tot),nbcon_nf(nf_tot_nbcon))
      DO i=1,ncell_fluid
         j0=i_neigh(i)-1
         DO j=i_neigh(i),i_neigh(i+1)-1
            nf_number=neighji_nf(j)
            istart=istart_nf(1,nf_number)
            idisp=neighji(j)
            IF(idisp.gt.0) THEN
                i1=istart+idisp
                left_nf(i1)=i
                jneigh_nf(i1)=j-j0
            ENDIF
         ENDDO
!
!.....Wall
!
         j=iface_wall(i)
         IF(j.eq.0) THEN
            i1=0
         ELSE
            nf_number=neighji_nf(j+j0)
            istart=istart_nf(1,nf_number)
            ii=neighji(j+j0)
            i1=istart+ii
         ENDIF
         iface_wall0(i)=i1
      ENDDO
      DO nf_number=1,7
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            istart2=istart_nbcon_nf(nf_number)
            DO i=1,len  
               i1=istart+i
               i2=istart2+i
               ii=left_nf(i1)
               j0=i_neigh(ii)-1
               jj=jneigh_nf(i1)
               nbcon_nf(i2)=nbcon(jj+j0)
            ENDDO  
      ENDDO  
!
      ALLOCATE(right_non(nf_non),kneigh_non(nf_non))
      ALLOCATE(right_fsw(nf_fsw))
!      
      i_non=0 
      i_fsw=0 
!
      DO i=1,ncell_fluid
         j0=i_neigh_nbcon0(i)
         DO j=i_neigh(i),i_neigh(i+1)-1
            IF(nbcon(j).eq.0)THEN
               k=neigh(j)
               IF(k.gt.i)THEN
                  i_non=i_non+1
                  right_non(i_non)=k
                  jk=nji(j0)
                  kneigh_non(i_non)=jk
               ENDIF
               j0=j0+1
            ELSEIF(nbcon(j).eq.-2)THEN
               k=neigh(j)
               i_fsw=i_fsw+1
               right_fsw(i_fsw)=k             
            ENDIF
         ENDDO
      ENDDO
!
!.....vg_f_non,vl_f_non
      ALLOCATE(vg_f_non(nf_non,ndim))
      ALLOCATE(vl_f_non(nf_non,ndim))
      IF(ndim.eq.2) THEN
         DO i=1,nf_non
            vg_f_non(i,1)=0.d0
            vg_f_non(i,2)=0.d0
            vl_f_non(i,1)=0.d0
            vl_f_non(i,2)=0.d0
         ENDDO
      ELSE
         DO i=1,nf_non
            vg_f_non(i,1)=0.d0
            vg_f_non(i,2)=0.d0
            vg_f_non(i,3)=0.d0
            vl_f_non(i,1)=0.d0
            vl_f_non(i,2)=0.d0
            vl_f_non(i,3)=0.d0
         ENDDO
      ENDIF
!.....sfl,sfg,sfd
      ALLOCATE(sfl_nf(nf_flux,5),sfg_nf(nf_flux,5),sfd_nf(nf_flux,5))
      ALLOCATE(sfl_non_k(nf_nonk,5),sfg_non_k(nf_nonk,5),sfd_non_k(nf_nonk,5))
      ALLOCATE(sfl6_nf(nf_flux),sfg6_nf(nf_flux),sfd6_nf(nf_flux))
      ALLOCATE(sfl6_non_k(nf_nonk),sfg6_non_k(nf_nonk),sfd6_non_k(nf_nonk))
!.....npb
!      
      ALLOCATE(delta_npb(ncell_fp))
!          
      DO i=1,ncell_fp
         IF(npb(i).ne.0)THEN
            delta_npb(i)=0.0d0
         ELSE  
            delta_npb(i)=1.0d0
         ENDIF  
      ENDDO
!
!.....Initialize vel_bc_profile_inl
!
      IF(nf_inl.gt.0) ALLOCATE(vel_bc_profile_inl(nf_inl))
      CALL udfn_vel_bc_profile_inl
!
!.....mass_nf_mcc
!
      IF(nf_mcc.gt.0) ALLOCATE(mass_nf_mcc(nf_mcc))
!
      END SUBROUTINE vectorize_index
