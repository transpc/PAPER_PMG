!
      SUBROUTINE wall_condensation_RBLA
!
!     This routine calculates gamma_wall(i) using RBLA model.
!
      USE Zinterface
      USE VOL_DATA         , ONLY: cell
      USE Wall_DATA        , ONLY: face
      USE SOLID_DATA       , ONLY: solid
      USE Zvec_param       , ONLY: nf_ctw,nf_fsw
      USE Zvec_index       , ONLY: left_nf,nbcon_nf,right_fsw
      USE Zncg             , ONLY: wmole_gas
      USE Znum_cell        , ONLY: istart_nf,istart_nbcon_nf,                &
                                   nf_number_nb,lens,nf_number_id,istart_nfs
      USE STM_TBL_cupid    , ONLY: st_tbl,   &
                                   nt,ndxstd
      USE Zb_condition     , ONLY: twall
      USE Zcoord3          , ONLY: volr
      USE Zturb            , ONLY: wallnr
      USE Zqvol            , ONLY: gamma_wall
      USE Zmodel           , ONLY: qconden,coef_diff
      USE Zvec_geo         , ONLY: saa_nf
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER i,k,ii,kk
      INTEGER :: nv,nf_number,isize,istart0,istart,istart2,i0,i1,i2
      LOGICAL err
      REAL(8) Wv,Wv1r,dWvdx
      REAL(8) ppsTi,Xvi,Wvi,Molv
!.....Local arrays
      REAL(8) s(36)
!.....Local vector arrays
      REAL(8) :: gamma_wall_nf(nf_fsw+nf_ctw),qconden_nf(nf_fsw+nf_ctw)
!
      DATA Molv/18.02d0/   ! molecular weights of water vapor & air
!
!.....Build summation info for fsw,ctw
!
      nf_number_nb=1
      nf_number_id(0)=5
      nf_number_id(1)=6
      istart_nfs(0)=0
      istart_nfs(1)=istart_nfs(0)+nf_fsw
      lens         =istart_nfs(1)+nf_ctw
!
      nv=0
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
!
      DO i=1,isize
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)
!
         kk=right_fsw(i)
         face%twall_partition(ii)=solid%tsol_o(kk)
         IF(face%twall_partition(ii).lt.cell%ts(ii))THEN
            Wv=1.0d0-cell%quala(ii)
            Wv1r=1.0d0/(Wv-1.0d0)
!
            s(1)=face%twall_partition(ii)
            s(9)=1.0d0
            CALL sth2x1_cupid(st_tbl(ndxstd),st_tbl(ndxstd+nt),s,err)
            ppsTi=s(10)
            Xvi=ppsTi/cell%p(ii)
            Wvi=(Xvi*Molv)/((Xvi*Molv)+(1.0d0-Xvi)*wmole_gas(ii))
            dWvdx=(Wvi-Wv)*wallnr(ii)
!
            qconden_nf(i0)=-Wv1r*cell%rhog(ii)*coef_diff(ii)*dWvdx*(cell%hg(ii)-cell%hlsat(ii))
            gamma_wall_nf(i0)=-Wv1r*cell%rhog(ii)*coef_diff(ii)*dWvdx*saa_nf(i1)*volr(ii)
         ELSE
            qconden_nf(i0)=0.d0
            gamma_wall_nf(i0)=0.d0
         ENDIF
      ENDDO
!
      nv=1
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
      istart2=istart_nbcon_nf(nf_number)
!
      DO i=1,isize
         i0=istart0+i
         i1=istart+i
         i2=istart2+i
         ii=left_nf(i1)
!
         k=-nbcon_nf(i2)
         face%twall_partition(ii)=twall(k)
         IF(face%twall_partition(ii).lt.cell%ts(ii))THEN
            Wv=1.0d0-cell%quala(ii)
            Wv1r=1.0d0/(Wv-1.0d0)
!
            s(1)=face%twall_partition(ii)
            s(9)=1.0d0
            CALL sth2x1_cupid(st_tbl(ndxstd),st_tbl(ndxstd+nt),s,err)
            ppsTi=s(10)
            Xvi=ppsTi/cell%p(ii)
            Wvi=(Xvi*Molv)/((Xvi*Molv)+(1.0d0-Xvi)*wmole_gas(ii))
            dWvdx=(Wvi-Wv)*wallnr(ii)
!
            qconden_nf(i0)=-Wv1r*cell%rhog(ii)*coef_diff(ii)*dWvdx*(cell%hg(ii)-cell%hlsat(ii))
            gamma_wall_nf(i0)=-Wv1r*cell%rhog(ii)*coef_diff(ii)*dWvdx*saa_nf(i1)*volr(ii)
         ELSE
            qconden_nf(i0)=0.d0
            gamma_wall_nf(i0)=0.d0
         ENDIF
      ENDDO
!      
      CALL sum_nf(0,-1   ,                  &
                  gamma_wall_nf,gamma_wall, &
                  qconden_nf   ,qconden)
!
      END SUBROUTINE wall_condensation_RBLA
!
!------------------------------------------------------------------------------
!
      SUBROUTINE wall_condensation_RBLA_porous
!
      USE VOL_DATA        , ONLY: cell
      USE SOLID_DATA      , ONLY: solid
      USE STM_TBL_cupid   , ONLY: st_tbl,   &
                                  nt,ndxstd
      USE Zcoord3         , ONLY: volp,aporous
      USE Znum_cell       , ONLY: n_fluid
      USE Zqvol           , ONLY: qporous_gas,qporous_gamma,gamma_wall
      USE Zturb           , ONLY: wallnr
      USE Zmodel          , ONLY: coef_diff,qconden
      USE Zzone           , ONLY: ncell_cond,nmaterial_c
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER i,ii
      LOGICAL err
      REAL(8) twallp
      REAL(8) Wv,Wv1r,dWvdx
      REAL(8) ppsTi,Xvi,Wvi,Molv,Mola
!.....Local arrays
      REAL(8) s(36)
!
      DATA Molv,Mola/18.02d0,28.97d0/   ! molecular weights of water vapor & air
!
      qconden(:)=0.0d0
!
      DO i=1,ncell_cond
         IF(nmaterial_c(i).lt.0)THEN
            ii=n_fluid(i)
!
            twallp=solid%tsol(i)
            IF(twallp.lt.cell%ts(ii))THEN
               Wv=1.0d0-cell%quala(ii)
               Wv1r=1.0d0/(Wv-1.0d0)
!
               s(1)=twallp
               s(9)=1.0d0
               CALL sth2x1_cupid(st_tbl(ndxstd),st_tbl(ndxstd+nt),s,err)
               ppsTi=s(10)
               Xvi=ppsTi/cell%p(ii)
               Wvi=(Xvi*Molv)/((Xvi*Molv)+(1.0d0-Xvi)*Mola)
               dWvdx=(Wvi-Wv)*wallnr(ii)
!
               gamma_wall(ii)=-Wv1r*cell%rhog(ii)*coef_diff(ii)*dWvdx*aporous(ii)/volp(ii)
               qconden(ii)=-Wv1r*cell%rhog(ii)*coef_diff(ii)*dWvdx*(cell%hg(ii)-cell%hlsat(ii)) !W/m2
               qporous_gas(ii)=0.0d0
               qporous_gamma(ii)=qconden(ii)*aporous(ii)
            ENDIF
         ENDIF
      ENDDO
!
      END SUBROUTINE wall_condensation_RBLA_porous
