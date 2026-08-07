      SUBROUTINE vectorize_geo_var
!
!     This routine vectorizes geometrical variables.
!
      USE Zzone        , ONLY: ncell_fluid,ncell_cond
      USE Zparam       , ONLY: ndim
      USE Zvec_param   , ONLY: nf_tot,nf_nonk,nf_non,nf_out,nf_fsw,nf_nbcon_change, &
                               nfc_tot,nfc_tot1
      USE Zvec_index   , ONLY: left_nf,right_non,jneigh_nf,kneigh_non
      USE Znum_cell    , ONLY: istart_nf,istartc_nf,              &
                               i_neigh,i_neigh_c,                 &
                               right_nb_k,neigh,&
                               nbc_nf
      USE Zvec_index_solid , ONLY: left_solid_nf,jneigh_solid_nf
      USE Zgradoption  , ONLY: iavgtype
      USE Zwall_HTC    , ONLY: f_direc
      USE Znormal      , ONLY: xn,nji,i_neigh_nbcon0
      USE Zbc_index     , ONLY: nbcon
      USE Zcoord1      , ONLY: xloc,xloc_m
      USE Zcoord2      , ONLY: fac,fac1,xfc,fac_c,fac1_c
      USE Zcoord3      , ONLY: sv,svp,permeability
      USE Zcoord4      , ONLY: sap,sa,saa,sad,dji_a,dji_x,dji,dnj, &
                               sap_c,dji_a_c
      USE Zvec_geo     , ONLY: xn_nf,sv_nf,svp_nf,xfc_nf,dxfc_nf, &
                               sap_nf,sa_nf,saa_nf,               &
                               dji_nf,djia_nf,dji_x_nf,           &
                               fac_nf,fac1_nf,                    &
                               djir_non,                          &
                               dnj_non,xloc_m_non_i,              &
                               dxfc_non_k,                        &
                               xloc_m_non_k,                      &
                               fac_non,fac_fsw,                   &
                               fac1_non,fac1_fsw,                 &
                               f0,f1,                             &
                               perm_non,perm_out,perm_nf,         &
                               sad_non,                           &
                               fac_c_nf,fac1_c_nf,                &
                               sap_c_nf,dji_a_c_nf,               &
                               xn_non_k,djia_non_k,dji_x_non_k,   &
                               sv_non_k,fac1_non_k
      USE Zcore        , ONLY: myrank
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,j,k
      INTEGER :: ix,ii,kk,jk,j0,k0,j1
      INTEGER :: nf_number,istart,len,i1
!
!.....1D arrays
!
      ALLOCATE(sap_nf(nf_tot),sa_nf(nf_tot),saa_nf(nf_tot))
      ALLOCATE(dji_nf(nf_tot))
      ALLOCATE(djia_nf(nf_tot))
      ALLOCATE(djia_non_k(nf_nonk))
      ALLOCATE(fac_nf(nf_tot),fac1_nf(nf_tot))
      ALLOCATE(fac1_non_k(nf_nonk))
      DO nf_number=0,8
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len
            i1=istart+i
            ii=left_nf(i1)
            j0=i_neigh(ii)-1
            j=jneigh_nf(i1)+j0
            sap_nf(i1) =sap(j)
            sa_nf(i1)  =sa(j)
            saa_nf(i1) =saa(j)
            dji_nf(i1) =dji(j)
            djia_nf(i1)=dji_a(j)
            fac_nf(i1) =fac(j)
            fac1_nf(i1)=fac1(j)
         ENDDO
         DO i=1,nf_nonk
            k=right_nb_k(i)
            djia_non_k(i)=djia_nf(k)
            fac1_non_k(i)=fac1_nf(k)
         ENDDO
      ENDDO
!
!.....Special case needed for nbcon change to extract perm_out
!
      ALLOCATE(perm_nf(nf_nbcon_change))
      j1=istart_nf(1,2)
      DO nf_number=2,4
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len
            i1=istart+i
            ii=left_nf(i1)
            j0=i_neigh(ii)-1
            j=jneigh_nf(i1)+j0
            perm_nf(i1-j1)=permeability(j)
         ENDDO
      ENDDO
!
!.....2D arrays
!
      ALLOCATE(xn_nf(nf_tot,ndim))
      ALLOCATE(xn_non_k(nf_nonk,ndim))
      ALLOCATE(sv_nf(nf_tot,ndim))
      ALLOCATE(sv_non_k(nf_nonk,ndim))
      ALLOCATE(svp_nf(nf_tot,ndim))
      ALLOCATE(xfc_nf(nf_tot,ndim))
      ALLOCATE(dxfc_nf(nf_tot,ndim))
      ALLOCATE(dji_x_nf(nf_tot,ndim))
      ALLOCATE(dji_x_non_k(nf_nonk,ndim))
      DO nf_number=0,8
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO ix=1,ndim
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               j0=i_neigh(ii)-1
               j=jneigh_nf(i1)+j0
               xn_nf(i1,ix)  =xn(j,ix)
               xn_nf(i1,ix)  =xn(j,ix)
               sv_nf(i1,ix)  =sv(j,ix)
               svp_nf(i1,ix) =svp(j,ix)
               xfc_nf(i1,ix) =xfc(j,ix)
               dxfc_nf(i1,ix)=xfc(j,ix)-xloc(ii,ix)
               dji_x_nf(i1,ix)=dji_x(j,ix)
            ENDDO
         ENDDO
      ENDDO
      DO ix=1,ndim
         DO i=1,nf_nonk
            k=right_nb_k(i)
            xn_non_k(i,ix)=xn_nf(k,ix)
            dji_x_non_k(i,ix)=dji_x_nf(k,ix)
            sv_non_k(i,ix)=sv_nf(k,ix)
         ENDDO
      ENDDO
!.....Check symetry
      GOTO 100
      i=1
      IF(myrank.eq.0) write(*,'(20(i4,2x),2(f15.7))') i,(neigh(j),j=i_neigh(i),i_neigh(i+1)-1)
      IF(myrank.eq.0) write(*,'(20(f15.7,2x))') (sv(j,1),j=i_neigh(i),i_neigh(i+1)-1)
      i=2
      IF(myrank.eq.0) write(*,'(20(i4,2x),2(f15.7))') i,(neigh(j),j=i_neigh(i),i_neigh(i+1)-1)
      IF(myrank.eq.0) write(*,'(20(f15.7,2x))') (sv(j,1),j=i_neigh(i),i_neigh(i+1)-1)
      DO i=1,ncell_fluid
         j0=i_neigh_nbcon0(i)
         j1=i_neigh(i)-1
         DO j=i_neigh(i),i_neigh(i+1)-1
            IF(nbcon(j).eq.0)THEN
               k=neigh(j)
               k0=i_neigh(k)-1
               IF(k.lt.i) THEN
                  jk=nji(j0)
                  j0=j0+1
                  IF(sv(j,1).eq.0.d0) CYCLE 
                  IF(sv(j,1).eq.-sv(jk+k0,1)) THEN
                     IF(myrank.eq.0) write(*,'(a,4(i4,2x),2(f15.7))') 'sv NOT symetric',i,j-j1,k,jk,sv(j,1),svp(jk+k0,1)
                  ELSE
                     IF(myrank.eq.0) write(*,'(a,4(i4,2x),2(f15.7))') 'sv symetric',i,j-j1,k,jk,sv(j,1),sv(jk+k0,1)
                  ENDIF
                  CALL finalize_mpi
                  STOP
               ENDIF
            ENDIF
         ENDDO
      ENDDO
100   CONTINUE
!
!.....1D non arrays
!
      ALLOCATE(djir_non(nf_non))
      ALLOCATE(fac_non(nf_non),fac_fsw(nf_fsw))
      ALLOCATE(fac1_non(nf_non),fac1_fsw(nf_fsw))
      ALLOCATE(f0(nf_non),f1(nf_non))
      ALLOCATE(perm_non(nf_non),perm_out(nf_out))
      ALLOCATE(sad_non(nf_non))
      nf_number=0
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         i1=istart+i
         ii=left_nf(i1)
         j0=i_neigh(ii)-1
         j=jneigh_nf(i1)+j0
         djir_non(i)=1.d0/dji(j)
         fac_non(i)= fac(j)
         fac1_non(i)=fac1(j)
         perm_non(i)=permeability(j)
         sad_non(i)=sad(j)
      ENDDO
!
!.....2D non arrays
!
      ALLOCATE(dnj_non(nf_non,ndim))
      ALLOCATE(xloc_m_non_i(nf_non,ndim))
      ALLOCATE(xloc_m_non_k(nf_non,ndim))
      ALLOCATE(dxfc_non_k(nf_non,ndim))
      DO ix=1,ndim
         DO i=1,len
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            j0=i_neigh(ii)-1
            j=jneigh_nf(i1)+j0
            dnj_non(i,ix)     =dnj(j,ix)
            xloc_m_non_i(i,ix)=xloc_m(j,ix)
            k0=i_neigh(kk)-1
            jk=kneigh_non(i)+k0
            xloc_m_non_k(i,ix)=xloc_m(jk,ix)
            dxfc_non_k(i,ix)  =xfc(j,ix)-xloc(kk,ix)
         ENDDO
      ENDDO
      IF(iavgtype.eq.1)THEN
         DO i=1,nf_non
            f0(i)=0.5d0
            f1(i)=0.5d0
         ENDDO
      ELSE
         DO i=1,nf_non
            f0(i)=fac_non(i)
            f1(i)=fac1_non(i)
         ENDDO
      ENDIF
!
!.....Arrays out
!
      nf_number=3
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i1=istart+i
         ii=left_nf(i1)
         j0=i_neigh(ii)-1
         j=jneigh_nf(i1)+j0
         perm_out(i)=permeability(j)
      ENDDO
!
!.....Arrays fsw
!
      nf_number=5
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i1=istart+i
         ii=left_nf(i1)
         j0=i_neigh(ii)-1
         j=jneigh_nf(i1)+j0
         fac_fsw(i)= fac(j)
         fac1_fsw(i)=fac1(j)
      ENDDO
!
!.....Solid arrays
!
      IF(ncell_cond.gt.0) THEN
         ALLOCATE(fac_c_nf(nfc_tot1),fac1_c_nf(nfc_tot1))
         ALLOCATE(sap_c_nf(nfc_tot),dji_a_c_nf(nfc_tot))
         DO nf_number=0,nbc_nf
            istart=istartc_nf(1,nf_number)
            len   =istartc_nf(2,nf_number)
            IF(nf_number.eq.0 .or. nf_number.eq.1) THEN
               DO i=1,len
                  i1=istart+i
                  ii=left_solid_nf(i1)
                  j0=i_neigh_c(ii)-1
                  j=jneigh_solid_nf(i1)+j0
                  fac_c_nf(i1)=fac_c(j)
                  fac1_c_nf(i1)=fac1_c(j)
               ENDDO
            ENDIF
            DO i=1,len
               i1=istart+i
               ii=left_solid_nf(i1)
               j0=i_neigh_c(ii)-1
               j=jneigh_solid_nf(i1)+j0
               sap_c_nf(i1)=sap_c(j)
            ENDDO
            IF(nf_number.eq.1 .or. nf_number.eq.3) THEN
               DO i=1,len
                  i1=istart+i
                  ii=left_solid_nf(i1)
                  j0=i_neigh_c(ii)-1
                  j=jneigh_solid_nf(i1)+j0
                  dji_a_c_nf(i1)=dji_a_c(j)
               ENDDO
            ENDIF
         ENDDO
      ENDIF
!         
!.....!!!cyj: Calculate dZ of each cell for QF calc. when Reflood option is ON.
!       
!     IF(reflood)THEN
!     IF(initial)THEN
!        IF(ALLOCATED(dz_cell)) DEALLOCATE(dz_cell)
!        ALLOCATE(dz_cell(ncell_fluid))       
!     ENDIF          
!        dz_cell(:)=0.0d0  
!        DO i=1, ncell_fluid 
!           h_min=1.0d6
!           h_max=-1.0d6
!           DO j=1, num_neigh(i)
!              h_min=DMIN1(h_min,xfc(f_direc,j,i))
!              h_max=DMAX1(h_max,xfc(f_direc,j,i))
!           ENDDO
!           dz_cell(i)=h_max-h_min
!        ENDDO
!     ENDIF          
!
      END SUBROUTINE vectorize_geo_var
