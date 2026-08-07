!
      SUBROUTINE GeoData_user_opr1000_rod01_ser(nn)
!
!     Define user-defined geometry data for OPR1000 Full vessel (rod-scale)
!
      USE Zcore       , ONLY: myrank,np
      USE Zporous     , ONLY: chn_type_tmp
      USE Zrv_ncell   , ONLY: asm_nx,asm_ny,asm_nz,asm_ni,asm_ni2,chn_nx,chn_ny
      USE Zparam      , ONLY: ndim
      USE Znum_cell   , ONLY: i_neigh_tmp,j_neigh_tmp,j_nbcon_tmp,sv_tmp1
      USE Zmpi        , ONLY: maxmt_cell

      USE Zrv_hts_2d  , ONLY: dz_fuel00
      USE MASTER4     , ONLY: NXY_TH,NZ_TH,NXYF,NPINX,nchn
      USE Zporous     , ONLY: nz_nk,dz_nk,nz_th0,hz_nk
      USE Zporous     , ONLY: n_sg,h_sg,sg_loc
      USE Zporous     , ONLY: l_spacer_grid
      USE Zio_unit    , ONLY: unit_log
!
      IMPLICIT NONE

      !input
      INTEGER nn

      !open file index
      INTEGER scin      
      
      INTEGER i,ii,j
      INTEGER i0,i1,i2,i3,i4,i5,i6,i7,i8
      INTEGER j0,m,k,k0
      REAL(8) xn_tmp(maxmt_cell,ndim)
      REAL(8) sa

      INTEGER z0,nz0

       DATA scin/1/
!   
!======================================================================
! READ OPR1000 GeoData
!======================================================================

      IF(myrank.eq.0) THEN
         OPEN(unit=1550,file='subchannel_info.in',status='old', iostat=scin)
      ENDIF ! myrank
      IF(np.gt.1) CALL broadcast_i1(scin)
      IF(scin.ne.0)then
         IF(myrank.eq.0) WRITE(* ,*)'Program was terminated due to lack of <subchannel_info.in>!'
         IF(myrank.eq.0) WRITE(unit_log,*)'Program was terminated due to lack of <subchannel_info.in>!'
         CALL finalize_mpi
         STOP
      ENDIF

      ALLOCATE(asm_nx(nn),asm_ny(nn),asm_nz(nn),asm_ni(nn),asm_ni2(nn))
      ALLOCATE(chn_nx(nn),chn_ny(nn))
      ALLOCATE(chn_type_tmp(nn))

      IF(myrank.eq.0) then
         READ(1550,*)i0,z0
         nz0=z0
         ! Vessel application (number of nz is different)
         !IF(vv_prob.eq.'OPR1000_fullcore_modmesh02_rv_vessel')then
         !   nz0=z0-2
         !   nz0=z0
         !ENDIF
      ENDIF

      IF(np.gt.1) CALL broadcast_i1(i0)
      IF(np.gt.1) CALL broadcast_i1(nz0)

      DO i=1,nn
         asm_nx(i)=0
         asm_ny(i)=0
         asm_nz(i)=0
         asm_ni(i)=0
         asm_ni2(i)=0
         chn_nx(i)=0
         chn_ny(i)=0
         chn_type_tmp(i)=0
      ENDDO

if(0)then
else

      IF(myrank.eq.0)then
         DO ii=1,i0
            READ(1550,*)i,i1,i2,i3,i4,i5,i6,i7,i8
            asm_nx(i)=i1
            asm_ny(i)=i2
            asm_nz(i)=i3
            asm_ni(i)=i4
            asm_ni2(i)=i5
            chn_nx(i)=i6
            chn_ny(i)=i7
            IF(i7.eq.17)asm_ni2(i)=0
            chn_type_tmp(i)=i8
         ENDDO
      ENDIF
      
      IF(np.gt.1) THEN
         CALL broadcast_i(asm_nx,nn)
         CALL broadcast_i(asm_ny,nn)
         CALL broadcast_i(asm_nz,nn)
         CALL broadcast_i(asm_ni,nn)
         CALL broadcast_i(asm_ni2,nn)
         CALL broadcast_i(chn_nx,nn)
         CALL broadcast_i(chn_ny,nn)
         CALL broadcast_i(chn_type_tmp,nn)
      ENDIF
      
endif

if(0)then
      IF(myrank.eq.0) THEN
         DO i=1,nn
            j0=i_neigh_tmp(i)-1
            DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
               sa=sv_tmp1(j,1)**2+sv_tmp1(j,2)**2+sv_tmp1(j,3)**2
               sa=1.d0/DSQRT(sa)
               xn_tmp(j,1)=sv_tmp1(j,1)*sa
               xn_tmp(j,2)=sv_tmp1(j,2)*sa
               xn_tmp(j,3)=sv_tmp1(j,3)*sa
            ENDDO
         ENDDO

         DO i=1,nn
            IF(chn_type_tmp(i).eq.2 .or. chn_type_tmp(i).eq.3)then
               j0=i_neigh_tmp(i)-1
               DO j=1,6 
                  IF(dabs(xn_tmp(j+j0,1)).gt.0.5 .or. dabs(xn_tmp(j+j0,2)).gt.0.5)then
                     k=j_neigh_tmp(j+j0)
                     IF(k.ne.0 .and. chn_type_tmp(k).eq.0)then
                        j_nbcon_tmp(j+j0)=-1
                       !j_neigh_tmp(j+j0)=0
                        k0=i_neigh_tmp(k)-1
                        DO m=1,6
                           IF(j_neigh_tmp(m+k0).eq.i)then
                              j_nbcon_tmp(m+k0)=-1
                             !j_neigh_tmp(m+k0)=0
                             !j_neigh_tmp(j+j0)=0
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDIF
if(myrank.eq.0)write(*,*)'DONE nji to wall wow'
endif
!.....Axial 
!
      ALLOCATE(dz_fuel00(nz0))
      IF(myrank.eq.0) THEN
         DO j=1,nz0
            READ(1550,*)dz_fuel00(j)
         ENDDO
      ENDIF
      IF(np.gt.1) CALL broadcast_r(dz_fuel00,nz0)
!
!.....MASTER variables => Fixed size
      nz_th0=nz0
      NXY_TH=(177+64)*4
      NZ_TH =26
      NXYF=177
      NPINX=16
      nchn=maxval(asm_ni(:))
      nchn=int(dsqrt(dble(i0/nz0/nchn)))
      IF(myrank.eq.0)write(*,*)'nchn is ',nchn
!
!.....Power distribution for Non-uniform mesh between CUPID and NK code
!
      IF(myrank.eq.0)then
         READ(1550,*)nz_nk
         IF(nz_nk.eq.0)then
            WRITE(unit_log,*)'axial cell height of NK code is not defined in subchannel_info.in'
           !WRITE(97,*)'axial cell height of NK code is not defined in subchannel_info.in'
            WRITE(* ,*)'axial cell height of NK code is not defined in subchannel_info.in'
            nz_nk=nz0
            WRITE(unit_log,*)'axial cell height of NK code is to be same as CUPID cell height'
           !WRITE(97,*)'axial cell height of NK code is to be same as CUPID cell height'
            WRITE(* ,*)'axial cell height of NK code is to be same as CUPID cell height'
            ALLOCATE(dz_nk(nz_nk))
            ALLOCATE(hz_nk(nz_nk+1))
            dz_nk=dz_fuel00
            hz_nk(1)=0.0d0
            DO i=2,nz_nk+1
               hz_nk(i)=hz_nk(i-1)+dz_fuel00(i-1)
            ENDDO
         ELSE
            WRITE(unit_log,*)'axial cell number of NK code is',nz_nk
           !WRITE(97,*)'axial cell number of NK code is ',nz_nk
            WRITE(* ,*)'axial cell number of NK code is ',nz_nk
            ALLOCATE(dz_nk(nz_nk))
            ALLOCATE(hz_nk(nz_nk+1))
            dz_nk=0.0d0
            hz_nk=0.0d0
            DO i=1,nz_nk+1
               READ(1550,*)hz_nk(i)
            ENDDO
            DO i=1,nz_nk
               dz_nk(i)=hz_nk(i+1)-hz_nk(i)
            ENDDO
         ENDIF
      ENDIF
      IF(np.gt.1) CALL broadcast_i1(nz_nk)
      IF(myrank.ne.0) then
         ALLOCATE(dz_nk(nz_nk))
         ALLOCATE(hz_nk(nz_nk+1))
         dz_nk=0.0d0
         hz_nk=0.0d0
      ENDIF
      IF(np.gt.1) CALL broadcast_r(dz_nk,nz_nk)
      IF(np.gt.1) CALL broadcast_r(hz_nk,nz_nk+1)
!
!........Spacer Grid
!
      IF(l_spacer_grid)then
         IF(myrank.eq.0) then
            READ(1550,*) n_sg
            IF(n_sg.eq.0)then
               WRITE(unit_log,*)'Spacer Grid is not defined in subchannel_info.in'
              !WRITE(97,*)'Spacer Grid is not defined in subchannel_info.in'
               WRITE(* ,*)'Spacer Grid is not defined in subchannel_info.in'
               goto 1551
            ENDIF
            ALLOCATE(sg_loc(n_sg))
            READ(1550,*) h_sg
            DO i=1,n_sg
               READ(1550,*)sg_loc(i)
            ENDDO
         ENDIF
         IF(np.gt.1) THEN
            CALL broadcast_i1(n_sg)
            CALL broadcast_r1(h_sg(1))
         ENDIF
         IF(myrank.ne.0) ALLOCATE(sg_loc(n_sg))
         IF(np.gt.1) CALL broadcast_r(sg_loc,n_sg)
         IF(myrank.eq.0) then
               WRITE(unit_log,*)'Spacer Grid is successfully read from subchannel_info.in'
              !WRITE(97,*)'Spacer Grid is successfully read from subchannel_info.in'
               WRITE(* ,*)'Spacer Grid is successfully read from subchannel_info.in'
         ENDIF
1551     continue
      ELSE
         IF(myrank.eq.0)then
            WRITE(unit_log,*)'Spacer Grid model is not applied in somaFlow.in'
           !WRITE(97,*)'Spacer Grid model is not applied in somaFlow.in'
            WRITE(* ,*)'Spacer Grid model is not applied in somaFlow.in'
         ENDIF
      ENDIF
!
      END SUBROUTINE GeoData_user_opr1000_rod01_ser
!
