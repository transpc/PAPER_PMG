!
      SUBROUTINE scalar_matrix_solid(poiss_diag_solid,poiss_solid_non_i,poiss_solid_non_k,src_solid)
!
!     This routine defines matrix elements for the solid conduction equation
!
      USE Zinterface
      USE VOL_DATA         , ONLY: cell    
      USE WALL_DATA        , ONLY: face
      USE SOLID_DATA       , ONLY: solid
      USE Zzone            , ONLY: ncell_cond,nmaterial_c
      USE Zcore            , ONLY: np
      USE Zconst2          , ONLY: dt
      USE Zvec_param       , ONLY: nfc_nonk,nfc_non,nfc_fsw,nfc_ctw,nfc_chw
      USE Znum_cell        , ONLY: n_fluid,istartc_nf,rightc_nb_k,           &
                                   lens,nf_number_nb,nf_number_id,istart_nfs
      USE Zvec_index_solid , ONLY: left_solid_k,right_solid_k,left_solid_nf,right_solid_non,right_solid_fsw, &
                                   nbcon_solid_ctw,nbcon_solid_chw,f_fluid_fsw,flux_fsw
      USE Zb_condition     , ONLY: twall
      USE Zconst1          , ONLY: iheatpart
      USE Zcoord3          , ONLY: volp_c
      USE Zqvol            , ONLY: qwall_solid,qvol_ice_solid,qporous_liq,qporous_gas,qporous_gamma
      USE Zuserdefined     , ONLY: hflux_bc_profile_chw_c
      USE Zmodel           , ONLY: i_fs_temp_intpol,qconden,qrad
      USE Zvec_geo         , ONLY: fac_c_nf,fac1_c_nf,sap_c_nf,dji_a_c_nf
      USE Zrad_comp        , ONLY: qrad_sol
      USE Zqvol            , ONLY: qconv_sol 
!
      IMPLICIT NONE
!
!.....Output
      REAL(8),DIMENSION(ncell_cond) :: poiss_diag_solid,src_solid
      REAL(8),DIMENSION(nfc_non) :: poiss_solid_non_i
      REAL(8),DIMENSION(nfc_nonk) :: poiss_solid_non_k
!.....Local variables
      INTEGER :: i,k,ii,kk
      INTEGER :: nv,nf_number,len,istart0,istart,istart2,i0,i1,i2
      REAL(8) :: a,b
      REAL(8) :: tmp,cf
      REAL(8) :: condsi,h_profile
      REAL(8) :: partition,heat_partition_ratio
!.....Local vector arrays
      REAL(8),DIMENSION(nfc_nonk+nfc_non+nfc_fsw+nfc_ctw) :: poiss_diag_solid_nf
      REAL(8),DIMENSION(nfc_fsw+nfc_ctw+nfc_chw) :: src_solid_nf
!
      IF(np.gt.1) CALL communicate_1d_c(solid%conds)
!
!.....The gradient has to be considered here for non-orthogonal grid
!
!.....Build summation info for non,fsw,ctw
!
      nf_number_nb=2
      nf_number_id(-1)=-1
      nf_number_id(0)=0
      nf_number_id(1)=1
      nf_number_id(2)=2
      nf_number_id(3)=3                     ! used to address src
      istart_nfs(-1)=0
      istart_nfs(0)=istart_nfs(-1)+nfc_nonk
      istart_nfs(1)=istart_nfs(0) +nfc_non
      istart_nfs(2)=istart_nfs(1) +nfc_fsw
      istart_nfs(3)=istart_nfs(2) +nfc_ctw  ! used to address src
      lens         =istart_nfs(2) +nfc_ctw
!
!.....Set the matrix for the solid conduction equation
!
      nv=0
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istartc_nf(1,nf_number)
      len   =istartc_nf(2,nf_number)
      DO i=1,len
         i0=istart0+i
         i1=istart+i
         ii=left_solid_nf(i1)
         kk=right_solid_non(i)
         tmp=dt*sap_c_nf(i1)
         condsi=1.d0/(fac1_c_nf(i1)/solid%conds(ii)+fac_c_nf(i1)/solid%conds(kk)) 
         cf=(tmp*condsi)/(solid%rhocps(ii)*volp_c(ii))
         poiss_solid_non_i(i)=-cf
         poiss_diag_solid_nf(i0)=cf
      ENDDO
!
      nv=-1
      nf_number=nf_number_id(nv)
      len   =istartc_nf(2,nf_number)
      DO i=1,len
         ii=left_solid_k(i)
         kk=right_solid_k(i)
         k=rightc_nb_k(i)
         tmp=dt*sap_c_nf(k)
         condsi=1.d0/(fac_c_nf(k)/solid%conds(ii)+fac1_c_nf(k)/solid%conds(kk)) 
         cf=(tmp*condsi)/(solid%rhocps(ii)*volp_c(ii))
         poiss_solid_non_k(i)=-cf
         poiss_diag_solid_nf(i)=cf
      ENDDO
!
      nv=1                            !Fluid-solid interface
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart2=istart0-istart_nfs(1)
      istart=istartc_nf(1,nf_number)
      len   =istartc_nf(2,nf_number)
      IF(i_fs_temp_intpol.ge.1)THEN         
         DO i=1,len
            i0=istart0+i
            i1=istart+i
            i2=istart2+i
            ii=left_solid_nf(i1)
            kk=right_solid_fsw(i)     !Right: Fluid Side
            tmp=dt*sap_c_nf(i1)/fac_c_nf(i1)
            condsi=solid%conds(ii)
            cf=tmp*condsi/(solid%rhocps(ii)*volp_c(ii))
            IF(nmaterial_c(ii).ge.1)THEN !fluid-solid 
               partition=heat_partition_ratio(cell%alphal(kk))
               IF(iheatpart.eq.0)THEN
                  a= fac_c_nf(i1) *(partition*cell%condl(kk)*cell%tl_o(kk)+(1.d0-partition)*cell%condg(kk)*cell%tg_o(kk))  &
                    +fac1_c_nf(i1)*solid%conds(ii)*solid%tsol_o(ii)
                  b= fac_c_nf(i1) *(partition*cell%condl(kk)              +(1.d0-partition)*cell%condg(kk)              )  &
                    +fac1_c_nf(i1)*solid%conds(ii)
                  face%twall_partition(kk)=a/b
               ELSE 
                  a= fac_c_nf(i1) *(partition*cell%condl(kk)*cell%tl_o(kk)+(1.d0-partition)*cell%condg(kk)*cell%tg_o(kk))  &
                    +fac1_c_nf(i1)*solid%conds(ii)*solid%tsol_o(ii)
                  b= fac_c_nf(i1) *(partition*cell%condl(kk)              +(1.d0-partition)*cell%condg(kk)              )  & 
                    +fac1_c_nf(i1)*solid%conds(ii)
                  face%twall_partition(kk)=a/b
               ENDIF
            ENDIF
            src_solid_nf(i2)=cf*face%twall_partition(kk)
            poiss_diag_solid_nf(i0)=cf
         ENDDO
      ELSE
         DO i=1,len
            i0=istart0+i
            i1=istart+i
            i2=istart2+i
            ii=left_solid_nf(i1)
            kk=right_solid_fsw(i)
            cf=dt/(solid%rhocps(ii)*volp_c(ii))
            src_solid_nf(i2)=-cf*flux_fsw(f_fluid_fsw(i))                                       &  !flux_fsw=(qliq_fsw+qgas_gsw)*saa from scalar_energy_diffusion
                             -cf*sap_c_nf(i1)/fac_c_nf(i1)*dji_a_c_nf(i1)*(qconden(kk)+qrad(kk))  !qconden from condensation model
            poiss_diag_solid_nf(i0)=0.d0
         ENDDO
      ENDIF
!      
      nv=2                            !constant wall temp. BC
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart2=istart0-istart_nfs(1)
      istart=istartc_nf(1,nf_number)
      len   =istartc_nf(2,nf_number)
      DO i=1,len
         i0=istart0+i
         i1=istart+i
         i2=istart2+i
         ii=left_solid_nf(i1)
         k=-nbcon_solid_ctw(i)
         tmp=dt*sap_c_nf(i1)
         condsi=solid%conds(ii)
         cf=tmp*condsi/(solid%rhocps(ii)*volp_c(ii))
         src_solid_nf(i2)=cf*twall(k)
         poiss_diag_solid_nf(i0)=cf                           
      ENDDO
!      
      nv=3                            !constant heat flux BC
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart2=istart0-istart_nfs(1)
      istart=istartc_nf(1,nf_number)
      len   =istartc_nf(2,nf_number)
      IF(len.gt.0) CALL udfn_hflux_bc_profile_chw_c
      DO i=1,len
         i1=istart+i
         i2=istart2+i
         ii=left_solid_nf(i1)
         k=-nbcon_solid_chw(i)
         tmp=dt*sap_c_nf(i1)
         condsi=solid%conds(ii)
         cf=tmp*condsi/(solid%rhocps(ii)*volp_c(ii))
         h_profile=hflux_bc_profile_chw_c(i)
         src_solid_nf(i2)=cf*h_profile*qwall_solid(k)*dji_a_c_nf(i1)/condsi
      ENDDO
!
      CALL sum_nf_solid(1,0,0,                                &
                        poiss_diag_solid_nf,poiss_diag_solid)
!
!.....Build summation info for fsw,ctw,chw
!
      nf_number_nb=3
      nf_number_id(1)=1
      nf_number_id(2)=2
      nf_number_id(3)=3
      istart_nfs(1)=0
      istart_nfs(2)=istart_nfs(1)+nfc_fsw
      istart_nfs(3)=istart_nfs(2)+nfc_ctw
      lens         =istart_nfs(3)+nfc_chw
!
      CALL sum_nf_solid(1,0,1,                  &
                        src_solid_nf,src_solid)
!
!.....Complete Matrix: poiss_diag_solid,src_solid
!
      DO i=1,ncell_cond
         poiss_diag_solid(i)=poiss_diag_solid(i)+1.0d0
         src_solid(i)=src_solid(i)+(qvol_ice_solid(i)/solid%rhocps(i)*dt+solid%tsol_o(i)) &
                     +qrad_sol(i)/(solid%rhocps(i)*volp_c(i))*dt & 
                     +qconv_sol(i)/(solid%rhocps(i)*volp_c(i))*dt      
         IF(nmaterial_c(i).lt.0)THEN
            ii=n_fluid(i)
            src_solid(i)=src_solid(i)- (qporous_liq(ii)+qporous_gas(ii)+qporous_gamma(ii)) &
                                      /(solid%rhocps(i)*volp_c(i))*dt
         ENDIF
      ENDDO
!
!.....Clean the off-diagonal elements when porosity=0
!
!!!!!!!!!assuming volp_c(i) never be 0.d0 if not divide by zero
!     DO i=1,ncell_cond
!        IF(volp_c(i).eq.0.d0)THEN
!           src_solid(i)=0.d0
!           poiss_diag_solid(i)=1.d0
!           poiss_solid(:,i)=0.d0
!        ENDIF
!     ENDDO
!
      END SUBROUTINE scalar_matrix_solid
