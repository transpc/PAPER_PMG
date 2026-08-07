!
      SUBROUTINE continuity_for_smac3_imp
!
!     This routine solves non-conservative form of continuity equation.
!     This routine is activated when smac=3.
!
      USE Zinterface
      USE VOL_DATA    , ONLY: cell
      USE Zzone       , ONLY: ncell_fluid
      USE Zcore       , ONLY: np
      USE Zconst2     , ONLY: dt
      USE Zvec_param  , ONLY: nf_nonk,nf_non,nf_mcc,nf_inl,nf_out,nf_flux,nf_fluxk
      USE Znum_cell   , ONLY: istart_nf,istart_nbcon_nf,        &
                              nf_number_nb,lens,                &
                              right_nb_k,                       &
                              istart_nfs,nf_number_id,istart_nf
      USE Zvec_index  , ONLY: left_nf,right_non,nbcon_nf
      USE Zbc_index   , ONLY: npb
      USE Z2nd_order  , ONLY: mass_conv_2nd
      USE Zare        , ONLY: ar_liq,ar_gas
      USE Zb_condition, ONLY: alphab_liq,alphab_gas,rhob_liq,rhob_gas
      USE Zcoord3     , ONLY: volp,volpr
      USE Zdecoupled  , ONLY: al_min_c,ag_min_c
      USE Zdel_scalar , ONLY: del_el,del_eg,del_x
      USE Zimplicit   , ONLY: eps_imp_alpha,max_iter_alpha
      USE Zpress      , ONLY: pp
      USE Zqvol       , ONLY: gamma,gamma_wall
      USE Ztimecon    , ONLY: iso_thermal
      USE Zvec_major  , ONLY: flux_l_nf,flux_g_nf,liq_conv_nf,vap_conv_nf
      USE Zvec_scalar , ONLY: arli_nf,argi_nf
      ! OPR1000 rod-scale (EVVD)                              
      USE Zporous     , ONLY: tm_mas_l,tm_mas_g, &
                              vd_mas_l,vd_mas_g, &
                              mixing_vane_l
      USE Zporous     , ONLY: l_subchannel,l_mixing_vane
!
      IMPLICIT NONE

!
!.....Local variables
      INTEGER :: i,k
      INTEGER :: ii,kk
      INTEGER :: iswitch
      INTEGER :: nv,nf_number,len,istart0,istart,istart2,i0,i1,i2
      REAL(8) :: vt
      REAL(8) :: fl_1,fl_2,fg_1,fg_2
      REAL(8) :: ag,al
!.....Local arrays
      REAL(8),DIMENSION(ncell_fluid) :: alg,diag,src,rhog,rhol
!
      !OPR1000 rod-scale (EVVD)
      REAL(8),DIMENSION(ncell_fluid) :: evvdg,evvdl
!.....Local vector arrays
      REAL(8),DIMENSION(nf_non) :: off_diag_non_i
      REAL(8),DIMENSION(nf_nonk) :: off_diag_non_k
      REAL(8),DIMENSION(nf_non) :: fl_1_non,fl_2_non,fg_1_non,fg_2_non
      REAL(8),DIMENSION(nf_nonk+nf_non+nf_out) :: diag_nf
      REAL(8),DIMENSION(nf_fluxk) :: src_nf,alg_nf
      REAL(8),DIMENSION(nf_flux) :: al_nf,ag_nf
!
!.....Estimated density according to pp in n+1 step
!
      DO i=1,ncell_fluid
         rhog(i)=cell%rhog(i)+cell%drhogdp(i)*pp(i)
         rhol(i)=cell%rhol(i)+cell%drholdp(i)*pp(i)
      ENDDO
      IF(iso_thermal.eq.0)THEN
         DO i=1,ncell_fluid
            rhog(i)=rhog(i)+cell%drhogde(i)*del_eg(i)+cell%drhogdx(i)*del_x(i)
            rhol(i)=rhol(i)+cell%drholde(i)*del_el(i)
         ENDDO
      ENDIF
!
!.....EVVD
      IF(l_subchannel)then
         DO i=1,ncell_fluid
            evvdg(i)=tm_mas_g(i)+vd_mas_g(i)
            evvdl(i)=tm_mas_l(i)+vd_mas_l(i)
         ENDDO
      ELSE
         DO i=1,ncell_fluid
            evvdg(i)=0.0d0
            evvdl(i)=0.0d0
         ENDDO
      ENDIF   
!
!.....Write a matrix for continuity
!
      IF(l_subchannel) THEN
         IF(l_mixing_vane) THEN
            DO i=1,ncell_fluid
               evvdl(i)=evvdl(i)+mixing_vane_l(2,i)*volpr(i)
            ENDDO
         ENDIF   
         DO i=1,ncell_fluid
            vt=volp(i)/dt
            IF(cell%alphag_o(i).lt.ag_min_c) THEN
               diag(i)=-cell%rhog(i)*vt
               src(i) =-cell%rhog(i)*cell%alphal_o(i)*vt+(gamma(i)+gamma_wall(i)+evvdg(i))*dt !EVVD
            ELSEIF(cell%alphal_o(i).lt.al_min_c) THEN
               diag(i)=cell%rhol(i)*vt
               src(i) =cell%rhol(i)*cell%alphal_o(i)*vt-(gamma(i)+gamma_wall(i)-evvdl(i))*dt !EVVD
            ELSE
               diag(i)=(rhol(i)-rhog(i))*vt
               src(i)=(ar_liq(i)+ar_gas(i)-rhog(i))*vt+(evvdg(i)+evvdl(i))*vt*dt !EVVD            
            ENDIF
         ENDDO
      ELSE
         DO i=1,ncell_fluid
            vt=volp(i)/dt
            IF(cell%alphag_o(i).lt.ag_min_c) THEN
               diag(i)=-cell%rhog(i)*vt
               src(i) =-cell%rhog(i)*cell%alphal_o(i)*vt+(gamma(i)+gamma_wall(i))*dt
            ELSEIF(cell%alphal_o(i).lt.al_min_c) THEN
               diag(i)=cell%rhol(i)*vt
               src(i) =cell%rhol(i)*cell%alphal_o(i)*vt-(gamma(i)+gamma_wall(i))*dt
            ELSE
               diag(i)=(rhol(i)-rhog(i))*vt
               src(i)=(ar_liq(i)+ar_gas(i)-rhog(i))*vt
            ENDIF
         ENDDO
      ENDIF
!
!.....Convection
!
!
!.....Mass convection between computing cells
!
      nf_number=0
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)
         fl_1_non(i)=MAX(flux_l_nf(i1),0.d0)*cell%rhol(ii)
         fg_1_non(i)=MAX(flux_g_nf(i1),0.d0)*cell%rhog(ii)
         fl_2_non(i)=MIN(flux_l_nf(i1),0.d0)*cell%rhol(kk)
         fg_2_non(i)=MIN(flux_g_nf(i1),0.d0)*cell%rhog(kk)
      ENDDO
!
!.....2nd order convection
!
      IF(mass_conv_2nd.gt.0) CALL mass_2nd_conv_imp(fl_1_non,fl_2_non,fg_1_non,fg_2_non,src)
!
!.....Build summation info for non,mcc,inl,out
!
      nf_number_nb=3
      nf_number_id(-1)=-1
      nf_number_id(0)=0
      nf_number_id(1)=1
      nf_number_id(2)=2
      nf_number_id(3)=3
      istart_nfs(-1)=0
      istart_nfs(0)=istart_nfs(-1)+nf_nonk
      istart_nfs(1)=istart_nfs(0) +nf_non
      istart_nfs(2)=istart_nfs(1) +nf_mcc
      istart_nfs(3)=istart_nfs(2) +nf_inl
      lens         =istart_nfs(3) +nf_out
!
      nv=0
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)
         fl_1=fl_1_non(i)
         fg_1=fg_1_non(i)
         fl_2=fl_2_non(i)
         fg_2=fg_2_non(i)
         IF    (cell%alphag_o(ii).lt.ag_min_c)THEN
            diag_nf(i0)=-fg_1
            src_nf(i0)=-fg_1-fg_2
            off_diag_non_i(i)=-fg_2
         ELSEIF(cell%alphal_o(ii).lt.al_min_c)THEN
            diag_nf(i0)=fl_1
            src_nf(i0)=0.d0
            off_diag_non_i(i)=fl_2
         ELSE
            diag_nf(i0)=fl_1-fg_1
            src_nf(i0)=-fg_1-fg_2
            off_diag_non_i(i)=fl_2-fg_2
         ENDIF
      ENDDO
!
      nv=-1
      nf_number=nf_number_id(nv)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         k=right_nb_k(i)
         ii=right_non(k)
         fl_1=fl_1_non(k)
         fg_1=fg_1_non(k)
         fl_2=fl_2_non(k)
         fg_2=fg_2_non(k)
         IF    (cell%alphag_o(ii).lt.ag_min_c) THEN
            diag_nf(i)=fg_2
            src_nf(i)=fg_2+fg_1
            off_diag_non_k(i)=fg_1
         ELSEIF(cell%alphal_o(ii).lt.al_min_c) THEN
            diag_nf(i)=-fl_2
            src_nf(i)=0.d0
            off_diag_non_k(i)=-fl_1
         ELSE
            diag_nf(i)=fg_2-fl_2
            src_nf(i)=fg_2+fg_1
            off_diag_non_k(i)=-fl_1+fg_1
         ENDIF
      ENDDO
!
!.....Source for MARS interface
!
      nv=1
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)
         IF    (cell%alphag_o(ii).lt.ag_min_c) THEN
            src_nf(i0)=-argi_nf(i1)*flux_g_nf(i1)
         ELSEIF(cell%alphal_o(ii).lt.al_min_c) THEN
            src_nf(i0)=-arli_nf(i1)*flux_l_nf(i1)
         ELSE
            src_nf(i0)=-arli_nf(i1)*flux_l_nf(i1)-argi_nf(i1)*flux_g_nf(i1)
         ENDIF
      ENDDO
!
!.....Source term contribution for inlet boundary
!
      nv=2
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      DO i=1,len  
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)
         i2=istart2+i
         k=nbcon_nf(i2)
         IF    (cell%alphag_o(ii).lt.ag_min_c)THEN
            src_nf(i0)=-alphab_gas(k)*rhob_gas(k)*flux_g_nf(i1)
         ELSEIF(cell%alphal_o(ii).lt.al_min_c)THEN
            src_nf(i0)=-alphab_liq(k)*rhob_liq(k)*flux_l_nf(i1)
         ELSE
            src_nf(i0)=-alphab_liq(k)*rhob_liq(k)*flux_l_nf(i1)-alphab_gas(k)*rhob_gas(k)*flux_g_nf(i1)
         ENDIF
      ENDDO
!
!.....Source term contribution for outlet boundary
!
      nv=3
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)
         IF    (cell%alphag_o(ii).lt.ag_min_c)THEN
            src_nf(i0)=-flux_g_nf(i1)*cell%rhog(ii)
         ELSEIF(cell%alphal_o(ii).lt.al_min_c)THEN
            src_nf(i0)=0.d0
         ELSE
            src_nf(i0)=-flux_g_nf(i1)*cell%rhog(ii)
         ENDIF
      ENDDO
!
      CALL sum_nf(1,0,        &
                  src_nf,src)
!
!.....Build summation info for non,mcc,inl,out
!
      nf_number_nb=1
      nf_number_id(-1)=-1
      nf_number_id(0)=0
      nf_number_id(1)=3
      istart_nfs(-1)=0
      istart_nfs(0)=istart_nfs(-1)+nf_nonk
      istart_nfs(1)=istart_nfs(0) +nf_non
      lens         =istart_nfs(1) +nf_out
!
      nv=1
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)
         IF    (cell%alphag_o(ii).lt.ag_min_c)THEN
            diag_nf(i0)=-flux_g_nf(i1)*cell%rhog(ii)
         ELSEIF(cell%alphal_o(ii).lt.al_min_c)THEN
            diag_nf(i0)=flux_l_nf(i1)*cell%rhol(ii)
         ELSE
            diag_nf(i0)=flux_l_nf(i1)*cell%rhol(ii)-flux_g_nf(i1)*cell%rhog(ii)
         ENDIF
      ENDDO
      CALL sum_nf(1,0,          &
                  diag_nf,diag)
!
!.....Matrix Calculation for liquid fraction (cell%alphal(:))
!
!
! build directly solverCSR  array here
      CALL csr_build_a(diag,off_diag_non_i,off_diag_non_k)
!
      CALL csr_cg_solvers_scalar(diag,src,cell%alphal,eps_imp_alpha,max_iter_alpha)
!
!.....Void fraction update based on liquid fraction
!
!.....Following code is to enhance numerical stability when cell%alphag,l_o are small
!.....We check first if such situation is met and bypass that code if not  
!
      iswitch=0
      DO i=1,ncell_fluid
         IF    (cell%alphag_o(i).lt.ag_min_c) THEN
           iswitch=1
           exit
         ELSEIF(cell%alphal_o(i).lt.al_min_c) THEN
           iswitch=1
           exit
         ENDIF
      ENDDO
      IF(iswitch.eq.0) GOTO 200
!
!.....Build summation info for non,mcc,inl,out
!
      nf_number_nb=3
      nf_number_id(-1)=-1
      nf_number_id(0)=0
      nf_number_id(1)=1
      nf_number_id(2)=2
      nf_number_id(3)=3
      istart_nfs(-1)=0
      istart_nfs(0)=istart_nfs(-1)+nf_nonk
      istart_nfs(1)=istart_nfs(0) +nf_non
      istart_nfs(2)=istart_nfs(1) +nf_mcc
      istart_nfs(3)=istart_nfs(2) +nf_inl
      lens         =istart_nfs(3) +nf_out
!
!......Cell non,mcc,inl,out
!
      DO nf_number=0,nf_number_nb
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len  
            i1=istart+i
            ii=left_nf(i1)
            al_nf(i1)=-liq_conv_nf(i1)*flux_l_nf(i1)
            ag_nf(i1)=-vap_conv_nf(i1)*flux_g_nf(i1)
         ENDDO
      ENDDO
!
      nv=0
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)
         IF    (cell%alphag_o(ii).lt.ag_min_c) THEN
            alg_nf(i0)=ag_nf(i1)
         ELSEIF(cell%alphal_o(ii).lt.al_min_c) THEN
            alg_nf(i0)=al_nf(i1)
         ELSE
            alg_nf(i0)=0.d0
         ENDIF
      ENDDO
!
      nv=-1
      nf_number=nf_number_id(nv)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         k=right_nb_k(i)
         ii=right_non(k)
         IF    (cell%alphag_o(ii).lt.ag_min_c) THEN
            alg_nf(i)=-ag_nf(k)
         ELSEIF(cell%alphal_o(ii).lt.al_min_c) THEN
            alg_nf(i)=-al_nf(k)
         ELSE
            alg_nf(i)=0.d0
         ENDIF
      ENDDO
!
      DO nv=1,nf_number_nb
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(nv)
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len
            i0=istart0+i
            i1=istart+i
            ii=left_nf(i1)
            IF    (cell%alphag_o(ii).lt.ag_min_c) THEN
               alg_nf(i0)=ag_nf(i1)
            ELSEIF(cell%alphal_o(ii).lt.al_min_c) THEN
               alg_nf(i0)=al_nf(i1)
            ELSE
               alg_nf(i0)=0.d0
            ENDIF
         ENDDO
      ENDDO
!
      CALL sum_nf(0,0,        &
                  alg_nf,alg)
!
      DO i=1,ncell_fluid
         vt=dt*volpr(i)
         IF    (cell%alphag_o(i).lt.ag_min_c) THEN
           !ag=(ag+(gamma(i)+gamma_wall(i))*dt)/cell%rhog(i)
            ag=alg(i)*vt
            ag=(ag+(gamma(i)+gamma_wall(i)+evvdg(i))*dt)/cell%rhog(i) !EVVD               
            cell%alphag(i)=cell%alphag_o(i)+ag
            cell%alphag(i)=MIN(MAX(0.0d0,cell%alphag(i)),1.0d0)
            cell%alphal(i)=1.d0-cell%alphag(i)
         ELSEIF(cell%alphal_o(i).lt.al_min_c) THEN
           !al=(al-(gamma(i)+gamma_wall(i))*dt)/cell%rhol(i)
            al=alg(i)*vt
            al=(al-(gamma(i)+gamma_wall(i)-evvdl(i))*dt)/cell%rhol(i) !EVVD              
            cell%alphal(i)=cell%alphal_o(i)+al
            cell%alphal(i)=MIN(MAX(0.0d0,cell%alphal(i)),1.0d0)
            cell%alphag(i)=1.d0-cell%alphal(i)
         ENDIF
      ENDDO
200   CONTINUE
!         
      DO i=1,ncell_fluid
         IF(npb(i).ne.0) THEN
            cell%alphal(i)=cell%alphal_o(i)
         ELSE
            cell%alphal(i)=MIN(MAX(0.0d0,cell%alphal(i)),1.0d0)
         ENDIF  
         cell%alphag(i)=1.d0-cell%alphal(i)
         cell%alphad(i)=0.0d0         
      ENDDO
!
!.....Isothermal option to update temperatures
!
      IF(iso_thermal.eq.0) THEN
         DO i=1,ncell_fluid
            cell%ts(i)=cell%ts_o(i)+cell%dtsdp(i)*pp(i)
            cell%tg(i)=cell%tg_o(i)+cell%dtgdp(i)*pp(i)
            cell%tl(i)=cell%tl_o(i)+cell%dtldp(i)*pp(i)
         ENDDO
      ENDIF
!
!.....MPI communication
!
      IF(np.gt.1) CALL communicate_1d(cell%alphag, &
                                      cell%alphal)
!
      END SUBROUTINE continuity_for_smac3_imp
