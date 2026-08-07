!
      SUBROUTINE wall_condensation_HMTA
!
!     This routine calculates gamma_wall(i) using HMTA model.
!
      USE Zinterface
      USE VOL_DATA         , ONLY: cell
      USE SOLID_DATA       , ONLY: solid
      USE Wall_DATA        , ONLY: face
      USE Zvec_param       , ONLY: nf_ctw,nf_fsw
      USE Zvec_index       , ONLY: left_nf,nbcon_nf,right_fsw
      USE Zncg             , ONLY: wmole_gas
      USE Znum_cell        , ONLY: istart_nf,istart_nbcon_nf,                &
                                   nf_number_nb,lens,nf_number_id,istart_nfs
      USE STM_TBL_cupid    , ONLY: st_tbl,   &
                                   nt,ndxstd
      USE Zb_condition     , ONLY: twall
      USE Zcoord3          , ONLY: volr
      USE Zqvol            , ONLY: gamma_wall
      USE Zturb            , ONLY: utaug,yplusg
      USE Zmodel           , ONLY: coef_diff,qconden
      USE Zvec_geo         , ONLY: saa_nf
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER i,k,ii,kk
      INTEGER :: nv,nf_number,isize,istart0,istart,istart2,i0,i1,i2
      LOGICAL err
      REAL(8) hm,Wvi,Wvb
      REAL(8) h,h1,h2,Pr,Prt,PrPrt,P,Tref,Econst,vK
      REAL(8) Sc
      REAL(8) ppsTi,Xvi,Molv
!.....Local arrays
      REAL(8) s(36)
!.....Local vector arrays
      REAL(8) :: gamma_wall_nf(nf_fsw+nf_ctw),qconden_nf(nf_fsw+nf_ctw)
!
      DATA Molv/18.02d0/
!
!.....Build summation info for fsw,ctw
!
      nf_number_nb=1
      nf_number_id(0)=5
      nf_number_id(1)=6
      istart_nfs(0)=0
      istart_nfs(1)=istart_nfs(0)+nf_fsw
      lens=         istart_nfs(1)+nf_ctw
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
            h1=(face%twall_partition(ii)-cell%tg(ii))*cell%rhog(ii)*cell%cpg(ii)*utaug(ii)
!
            Pr=cell%lviscosg(ii)*cell%cpg(ii)/cell%lcondg(ii)
            Prt=0.9d0
            PrPrt=Pr/Prt
            P=9.24d0*(-1.d0+PrPrt**0.75d0)*(1.d0+0.28d0*EXP(-0.007d0*PrPrt))
!
            Tref=cell%tg(ii)
            Econst=9.793d0
            vK=0.42d0
            h2=(face%twall_partition(ii)-Tref)*Prt*(LOG(Econst*yplusg(ii))/vK+P)
            IF(yplusg(ii).lt.11.225d0) h2=(face%twall_partition(ii)-Tref)*Pr*yplusg(ii)
            h=h1/h2
!
            Sc=cell%lviscosg(ii)/(cell%rhog(ii)*coef_diff(ii))
            hm=h*(cell%rhog(ii)*coef_diff(ii)/cell%lcondg(ii))*(Sc/Pr)**0.333333
!               
            s(1)=face%twall_partition(ii)
            s(9)=1.d0
            CALL sth2x1_cupid(st_tbl(ndxstd),st_tbl(ndxstd+nt),s,err)
            ppsTi=s(10)
            Xvi=ppsTi/cell%p(ii)
            Wvi=(Xvi*Molv)/((Xvi*Molv)+(1.d0-Xvi)*wmole_gas(ii))
            Wvb=1.d0-cell%quala(ii)
!
            gamma_wall_nf(i0)=-hm*LOG((1.d0-Wvi)/(1.d0-Wvb))*saa_nf(i1)*volr(ii)
            qconden_nf(i0)=-hm*LOG((1.d0-Wvi)/(1.d0-Wvb))*(cell%hg(ii)-cell%hlsat(ii))
         ELSE
            gamma_wall_nf(i0)=0.d0
            qconden_nf(i0)   =0.d0
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
            h1=(face%twall_partition(ii)-cell%tg(ii))*cell%rhog(ii)*cell%cpg(ii)*utaug(ii)
!
            Pr=cell%lviscosg(ii)*cell%cpg(ii)/cell%lcondg(ii)
            Prt=0.9d0
            PrPrt=Pr/Prt
            P=9.24d0*(-1.d0+PrPrt**0.75d0)*(1.d0+0.28d0*EXP(-0.007d0*PrPrt))
!
            Tref=cell%tg(ii)
            Econst=9.793d0
            vK=0.42d0
            h2=(face%twall_partition(ii)-Tref)*Prt*(LOG(Econst*yplusg(ii))/vK+P)
            IF(yplusg(ii).lt.11.225d0) h2=(face%twall_partition(ii)-Tref)*Pr*yplusg(ii)
            h=h1/h2
!
            Sc=cell%lviscosg(ii)/(cell%rhog(ii)*coef_diff(ii))
            hm=h*(cell%rhog(ii)*coef_diff(ii)/cell%lcondg(ii))*(Sc/Pr)**0.333333
!               
            s(1)=face%twall_partition(ii)
            s(9)=1.d0
            CALL sth2x1_cupid(st_tbl(ndxstd),st_tbl(ndxstd+nt),s,err)
            ppsTi=s(10)
            Xvi=ppsTi/cell%p(ii)
            Wvi=(Xvi*Molv)/((Xvi*Molv)+(1.d0-Xvi)*wmole_gas(ii))
            Wvb=1.d0-cell%quala(ii)
!
            gamma_wall_nf(i0)=-hm*LOG((1.d0-Wvi)/(1.d0-Wvb))*saa_nf(i1)*volr(ii)
            qconden_nf(i0)=-hm*LOG((1.d0-Wvi)/(1.d0-Wvb))*(cell%hg(ii)-cell%hlsat(ii))
         ELSE
            gamma_wall_nf(i0)=0.d0
            qconden_nf(i0)=0.d0
         ENDIF
      ENDDO
!      
      CALL sum_nf(0,-1,                     &
                  gamma_wall_nf,gamma_wall, &
                  qconden_nf   ,qconden)
!
      END SUBROUTINE wall_condensation_HMTA
!
!------------------------------------------------------------------------------
!
      SUBROUTINE wall_condensation_HMTA_porous
!
      USE VOL_DATA        , ONLY: cell
      USE SOLID_DATA      , ONLY: solid
      USE STM_TBL_cupid   , ONLY: st_tbl,   &
                                  nt,ndxstd
      USE Zcoord3         , ONLY: volp,aporous
      USE Znum_cell       , ONLY: n_fluid
      USE Zqvol           , ONLY: qporous_gas,qporous_gamma,gamma_wall
      USE Zturb           , ONLY: utaug,yplusg
      USE Zmodel          , ONLY: coef_diff,qconden
      USE Zzone           , ONLY: ncell_cond,nmaterial_c
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER i,ii
      LOGICAL err
      REAL(8) twallp
      REAL(8) hm,Wvi,Wvb
      REAL(8) h,h1,h2,Pr,Prt,PrPrt,P,Tref,Econst,vK
      REAL(8) Sc
      REAL(8) s(36),ppsTi,Xvi,Molv,Mola
!
      DATA Molv,Mola/18.02d0,28.97d0/   ! molecular weights of water vapor & air
!
      qconden(:)=0.d0
!
      DO i=1,ncell_cond
         IF(nmaterial_c(i).lt.0)THEN
            ii=n_fluid(i)
!
            twallp=solid%tsol(i)
            IF(twallp.lt.cell%ts(ii))THEN
               h1=(twallp-cell%tg(ii))*cell%rhog(ii)*cell%cpg(ii)*utaug(ii)
!
               Pr=cell%lviscosg(ii)*cell%cpg(ii)/cell%lcondg(ii)
               Prt=0.9d0
               PrPrt=Pr/Prt
               P=9.24d0*(-1.d0+PrPrt**0.75d0)*(1.d0+0.28d0*DEXP(-0.007d0*PrPrt))
!
               Tref=cell%tg(ii)
               Econst=9.793d0
               vK=0.42d0
               h2=(twallp-Tref)*Prt*(DLOG(Econst*yplusg(ii))/vK+P)
               IF(yplusg(ii).lt.11.225d0) h2=(twallp-Tref)*Pr*yplusg(ii)
               h=h1/h2
!
               Sc=cell%lviscosg(ii)/(cell%rhog(ii)*coef_diff(ii))
               hm=h*(cell%rhog(ii)*coef_diff(ii)/cell%lcondg(ii))*(Sc/Pr)**0.333333
!               
               s(1)=twallp
               s(9)=1.d0
               CALL sth2x1_cupid(st_tbl(ndxstd),st_tbl(ndxstd+nt),s,err)
               ppsTi=s(10)
               Xvi=ppsTi/cell%p(ii)
               Wvi=(Xvi*Molv)/((Xvi*Molv)+(1.d0-Xvi)*Mola)
               Wvb=1.d0-cell%quala(ii)
!
               gamma_wall(ii)=-hm*DLOG((1.d0-Wvi)/(1.d0-Wvb))*aporous(ii)/volp(ii)
               qconden(ii)=-hm*DLOG((1.d0-Wvi)/(1.d0-Wvb))*(cell%hg(ii)-cell%hlsat(ii)) !W/m2
               qporous_gas(ii)=0.d0
               qporous_gamma(ii)=qconden(ii)*aporous(ii)
            ENDIF
         ENDIF
      ENDDO
!
      END SUBROUTINE wall_condensation_HMTA_porous
