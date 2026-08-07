!
      SUBROUTINE wall_condensation_uchida
!
!     This routine calculates gamma_wall(i) using Uchida correlation model.
!
      USE Zinterface
      USE VOL_DATA         , ONLY: cell
      USE Wall_DATA        , ONLY: face
      USE SOLID_DATA       , ONLY: solid
      USE Zvec_param       , ONLY: nf_ctw,nf_fsw
      USE Zvec_index       , ONLY: left_nf,nbcon_nf,right_fsw
      USE Znum_cell        , ONLY: istart_nf,istart_nbcon_nf, &
                                   nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zb_condition     , ONLY: twall
      USE Zcoord3          , ONLY: volr,volp
      USE Zconst2          , ONLY: dt
      USE Zqvol            , ONLY: gamma_wall
      USE Zmodel           , ONLY: qconden
      USE Zvec_geo         , ONLY: saa_nf
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,k,ii,kk
      INTEGER :: nv,nf_number,isize,istart0,istart,istart2,i0,i1,i2
      REAL(8) :: huchida,qx,gamma_wall_min
!.....Local vector arrays
      REAL(8) :: gamma_wall_nf(nf_fsw+nf_ctw),qconden_nf(nf_fsw+nf_ctw)
      LOGICAL, SAVE:: initial=.true.
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
!.....Fluid-Solid Interface
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
         kk=right_fsw(i)
!         
!Nuscale-03Pool
         IF(initial) THEN
            face%twall_partition(ii)=solid%tsol_o(kk)            !!! tsol is not the temp. at face. It should be considered.
            initial=.false.         
         ENDIF
!
         IF(face%twall_partition(ii).lt.cell%ts(ii))THEN
            qx=DMAX1(0.1d0,cell%quala(ii))
            huchida=(1.0d0-qx)/qx  !Nuscale-03Pool
            huchida=380.0d0*huchida**0.7d0
!
            qconden_nf(i0)=-huchida*(cell%ts(ii)-face%twall_partition(ii))
!            gamma_wall_nf(i0)=qconden_nf(i0)*saa_nf(i1)*volr(ii)/(cell%hg(ii)-cell%hlsat(ii))
!
            gamma_wall_nf(i0)=qconden_nf(i0)*saa_nf(i1)/(cell%hg(ii)-cell%hlsat(ii))
            gamma_wall_min=-cell%alphag(ii)*volp(ii)*cell%rhol(ii)/dt
            gamma_wall_nf(i0)=DMAX1(gamma_wall_nf(i0),gamma_wall_min)
            gamma_wall_nf(i0)=gamma_wall_nf(i0)*volr(ii)
            IF(cell%alphag(i).le.0.1d0) gamma_wall_nf(i0)=0.0d0
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
            qx=DMAX1(0.1d0,cell%quala(ii))
            huchida=(1.0d0-qx)/qx  !Nuscale-03Pool
            huchida=380.0d0*huchida**0.7d0
!
            qconden_nf(i0)=-huchida*(cell%ts(ii)-face%twall_partition(ii))
!            gamma_wall_nf(i0)=qconden_nf(i0)*saa_nf(i1)*volr(ii)/(cell%hg(ii)-cell%hlsat(ii))
!
            gamma_wall_nf(i0)=qconden_nf(i0)*saa_nf(i1)/(cell%hg(ii)-cell%hlsat(ii))
            gamma_wall_min=-cell%alphag(ii)*volp(ii)*cell%rhol(ii)/dt
            gamma_wall_nf(i0)=DMAX1(gamma_wall_nf(i0),gamma_wall_min)
            gamma_wall_nf(i0)=gamma_wall_nf(i0)*volr(ii)
            IF(cell%alphag(i).le.0.1d0) gamma_wall_nf(i0)=0.0d0
         ELSE
            qconden_nf(i0)=0.d0
            gamma_wall_nf(i0)=0.d0 
         ENDIF
      ENDDO
!
      CALL sum_nf(0,-1,                     &
                  gamma_wall_nf,gamma_wall, &
                  qconden_nf   ,qconden)
!
      END SUBROUTINE wall_condensation_uchida
!
!---------------------------------------------------------------------------------
!
      SUBROUTINE wall_condensation_uchida_porous
!
      USE VOL_DATA        , ONLY: cell
      USE SOLID_DATA      , ONLY: solid
      USE Znum_cell       , ONLY: n_fluid
      USE Zzone           , ONLY: ncell_cond,nmaterial_c
      USE Zcoord3         , ONLY: volp,aporous
      USE Zqvol           , ONLY: qporous_gas,qporous_gamma, &
                                  gamma_wall
      USE Zmodel          , ONLY: qconden
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER i,ii
      REAL(8) twallp,huchida
!
      qconden(:)=0.0d0
!
      DO i=1,ncell_cond
         IF(nmaterial_c(i).lt.0)THEN
            ii=n_fluid(i)
!
            twallp=solid%tsol(i)
            IF(twallp.lt.cell%ts(ii))THEN
               huchida=(1.0d0-cell%quala(ii))/cell%quala(ii)
               huchida=380.0d0*huchida**0.7d0
!
               qconden(ii)=-huchida*(cell%ts(ii)-twallp) !W/m2
               gamma_wall(ii)=qconden(ii)*aporous(ii)/volp(ii)/(cell%hg(ii)-cell%hlsat(ii)) 
               qporous_gas(ii)=0.0d0
               qporous_gamma(ii)=qconden(ii)*aporous(ii)
            ENDIF
         ENDIF
      ENDDO
!
      END SUBROUTINE wall_condensation_uchida_porous
