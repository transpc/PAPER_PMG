!
      SUBROUTINE scalar_energy_diffusion
!
!     This routine calculates energy diffusive fluxes through the cell face
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell
      USE Wall_DATA    , ONLY: face
      USE SOLID_DATA   , ONLY: solid
      USE Zmpi         , ONLY: ncell_fp
      USE Zcore        , ONLY: np
      USE Zparam       , ONLY: ndim
      USE Zvec_param   , ONLY: nf_non,nf_inl,nf_fsw,nf_ctw,nf_chw
      USE Znum_cell    , ONLY: istart_nf,istart_nbcon_nf, &
                               nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zvec_index   , ONLY: left_nf,right_non,right_fsw,nbcon_nf
      USE Zconst1      , ONLY: iheatpart,wconden,rv_htmodel_forCFD
      USE Zb_condition , ONLY: alphab_liq,alphab_gas,tb_liq,tb_gas,twall
      USE Zenergy_diff , ONLY: ediff_liq,ediff_gas
      USE Zface        , ONLY: twall_model,laminar
      USE Zgradoption  , ONLY: non_orth_eng
      USE Zqvol        , ONLY: qwall_solid
      USE Zturb        , ONLY: wcd_liq,wcd_gas,wallnr
      USE Zuserdefined , ONLY: udfl_erg_diff,hflux_bc_profile_chw
      USE Zvec_geo     , ONLY: f1,f0,fac1_fsw,fac_fsw,     &
                               dnj_non,saa_nf,sap_nf,sa_nf
      USE Zvec_index_solid , ONLY: flux_fsw
      USE Zrv_model    , ONLY: rv_valve      
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,k
      integer :: ii,kk
      INTEGER :: nv,nf_number,len,istart0,istart,istart2,i0,i1,i2
!
      REAL(8) :: tl1,tg1,tl2,tg2,h_profile
      REAL(8) :: akli,akgi
      REAL(8) :: dtli1,dtli2,dtli3
      REAL(8) :: dtgi1,dtgi2,dtgi3
      REAL(8) :: dtlj,dtgj
      REAL(8) :: sv1
      REAL(8) :: partition,heat_partition_ratio
!     REAL(8) :: twallbc1,twallbc2 !only for vv_prob='nat_conv_ampofo'
!
      REAL(8) :: ali_tmp,alk_tmp
      REAL(8) :: agi_tmp,agk_tmp
!.....Local arrays
      REAL(8),DIMENSION(ncell_fp,ndim) :: dtldx,dtgdx
!.....Local vector arrays
      REAL(8),DIMENSION(nf_non) :: akli_non,akgi_non
      REAL(8),DIMENSION(nf_non+nf_inl+nf_fsw+nf_ctw+nf_chw) :: fluxl_diff_nf,fluxg_diff_nf
!
      fluxl_diff_nf(:)=0.0d0
      fluxg_diff_nf(:)=0.0d0
!
!.....Calculate temperature gradient at cell center for non-orthogonal grid
!
      IF(non_orth_eng.eq.1) THEN
         CALL grad_temp(cell%tl_o,dtldx,tb_liq, &
                        cell%tg_o,dtgdx,tb_gas)
      ELSEIF(non_orth_eng.eq.2) THEN
         CALL grad_temp(cell%tl_o,dtldx,tb_liq, &  ! will be replaced with Frink method
                        cell%tg_o,dtgdx,tb_gas)
      ENDIF
      IF(non_orth_eng.gt.0)THEN
         IF(np.gt.1) CALL communicate_2d(dtldx, &
                                         dtgdx)
      ENDIF
!
!.....Build summation info for non,inl,fsw,ctw,chw
!
      nf_number_nb=4
      nf_number_id(0)=0
      nf_number_id(1)=2
      nf_number_id(2)=5
      nf_number_id(3)=6
      nf_number_id(4)=7
      istart_nfs(0)=0
      istart_nfs(1)=istart_nfs(0)+nf_non
      istart_nfs(2)=istart_nfs(1)+nf_inl
      istart_nfs(3)=istart_nfs(2)+nf_fsw
      istart_nfs(4)=istart_nfs(3)+nf_ctw
      lens         =istart_nfs(4)+nf_chw
!
!.....Computing cell
!
      nv=0
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)
         ali_tmp=cell%alphal(ii)
         alk_tmp=cell%alphal(kk)
         agi_tmp=cell%alphag(ii)
         agk_tmp=cell%alphag(kk)
         akli_non(i)=f1(i)*ali_tmp*cell%condl(ii)+f0(i)*alk_tmp*cell%condl(kk)
         akgi_non(i)=f1(i)*agi_tmp*cell%condg(ii)+f0(i)*agk_tmp*cell%condg(kk)
      ENDDO
      DO i=1,len  
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)
         akli=akli_non(i)
         akgi=akgi_non(i)
         tl1=cell%tl_o(ii)
         tg1=cell%tg_o(ii)
         tl2=cell%tl_o(kk)
         tg2=cell%tg_o(kk)
         fluxl_diff_nf(i0)=akli*(tl2-tl1)*sap_nf(i1)
         fluxg_diff_nf(i0)=akgi*(tg2-tg1)*sap_nf(i1)
      ENDDO
!
!........Non-orthogonal grid contribution
!
      IF(non_orth_eng.gt.0)THEN
         IF(ndim.eq.2)THEN
            DO i=1,len  
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
               akli=akli_non(i)
               akgi=akgi_non(i)
               dtli1=f1(i)*dtldx(ii,1)+f0(i)*dtldx(kk,1)
               dtli2=f1(i)*dtldx(ii,2)+f0(i)*dtldx(kk,2)
               dtgi1=f1(i)*dtgdx(ii,1)+f0(i)*dtgdx(kk,1)
               dtgi2=f1(i)*dtgdx(ii,2)+f0(i)*dtgdx(kk,2)
               dtlj=dtli1*dnj_non(i,1)+dtli2*dnj_non(i,2)
               dtgj=dtgi1*dnj_non(i,1)+dtgi2*dnj_non(i,2)
               fluxl_diff_nf(i0)=fluxl_diff_nf(i0)+akli*dtlj*sa_nf(i1)
               fluxg_diff_nf(i0)=fluxg_diff_nf(i0)+akgi*dtgj*sa_nf(i1)
            ENDDO
         ELSE
            DO i=1,len  
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
               akli=akli_non(i)
               akgi=akgi_non(i)
               dtli1=f1(i)*dtldx(ii,1)+f0(i)*dtldx(kk,1)
               dtli2=f1(i)*dtldx(ii,2)+f0(i)*dtldx(kk,2)
               dtli3=f1(i)*dtldx(ii,3)+f0(i)*dtldx(kk,3)
               dtgi1=f1(i)*dtgdx(ii,1)+f0(i)*dtgdx(kk,1)
               dtgi2=f1(i)*dtgdx(ii,2)+f0(i)*dtgdx(kk,2)
               dtgi3=f1(i)*dtgdx(ii,3)+f0(i)*dtgdx(kk,3)
               dtlj=dtli1*dnj_non(i,1)+dtli2*dnj_non(i,2)+dtli3*dnj_non(i,3)
               dtgj=dtgi1*dnj_non(i,1)+dtgi2*dnj_non(i,2)+dtgi3*dnj_non(i,3)
               fluxl_diff_nf(i0)=fluxl_diff_nf(i0)+akli*dtlj*sa_nf(i1)
               fluxg_diff_nf(i0)=fluxg_diff_nf(i0)+akgi*dtgj*sa_nf(i1)
            ENDDO
         ENDIF
      ENDIF
!
!.....valve model
!      
      IF(rv_valve.eq.1) CALL valve_model_scalar_energy_diffusion(fluxl_diff_nf,fluxg_diff_nf)
!
!.....outlet, adwall, symetric condition: flux=0
!
      !Do nothing
!      
!.....Inlet
!
      nv=1
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart2=istart_nbcon_nf(nf_number)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i0=istart0+i
         i1=istart+i
         i2=istart2+i
         ii=left_nf(i1)
         k=nbcon_nf(i2)
         akli=alphab_liq(k)*cell%condl(ii)
         akgi=alphab_gas(k)*cell%condg(ii)
         tl1=cell%tl_o(ii)
         tg1=cell%tg_o(ii)
         tl2=tb_liq(k)
         tg2=tb_gas(k)
         fluxl_diff_nf(i0)=akli*(tl2-tl1)*sap_nf(i1)
         fluxg_diff_nf(i0)=akgi*(tg2-tg1)*sap_nf(i1)
      ENDDO
!
      IF(rv_htmodel_forCFD.gt.0)GOTO 1
!
!.....Fluid-Solid interface
!
      nv=2
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)
         kk=right_fsw(i)
         ali_tmp=cell%alphal(ii)
         agi_tmp=cell%alphag(ii)
!         
         akli=ali_tmp*cell%condl(ii)
         akgi=agi_tmp*cell%condg(ii)
         partition=heat_partition_ratio(ali_tmp)
         tl1=cell%tl_o(ii)
         tg1=cell%tg_o(ii)
         tl2=((partition*cell%condl(ii)*tl1                                           &
             +(1.0d0-partition)*cell%condg(ii)*tg1)*fac1_fsw(i)                       &
             +fac_fsw(i)*solid%conds(kk)*solid%tsol_o(kk))/((partition*cell%condl(ii) &
             +(1.0d0-partition)*cell%condg(ii))*fac1_fsw(i)+fac_fsw(i)*solid%conds(kk))
         tg2=tl2
!
         IF(Twall_Model.eq.Laminar)THEN
            fluxl_diff_nf(i0)=akli*(tl2-tl1)*sap_nf(i1)/fac_fsw(i)
            fluxg_diff_nf(i0)=akgi*(tg2-tg1)*sap_nf(i1)/fac_fsw(i)
            flux_fsw(i)=(fluxl_diff_nf(i0)+fluxg_diff_nf(i0))
         ELSE
            sv1=sa_nf(i1)*wallnr(ii)
            fluxl_diff_nf(i0)=ali_tmp*wcd_liq(ii)*(tl2-tl1)*sv1
            fluxg_diff_nf(i0)=agi_tmp*wcd_gas(ii)*(tg2-tg1)*sv1            
            flux_fsw(i)=(fluxl_diff_nf(i0)+fluxg_diff_nf(i0))
         ENDIF
!==>next
         IF(iheatpart.gt.0)THEN
            fluxl_diff_nf(i0)=face%wall_fluxl_diff(ii)
            fluxg_diff_nf(i0)=face%wall_fluxg_diff(ii)
            flux_fsw(i)=fluxl_diff_nf(i0)+fluxg_diff_nf(i0)
         !ELSEIF(wconden.ne.0)THEN   
            !!!qliq_fsw: partitioned heat flux from fluid to the wall (total - latent). Not considered in the current condensation model.
            !fluxl_diff_nf(i0)=qliq_fsw(i)*sap_nf(i1)/fac_fsw(i)*djia_nf(i1)
            !fluxg_diff_nf(i0)=qgas_fsw(i)*sap_nf(i1)/fac_fsw(i)*djia_nf(i1) 
            !flux_fsw(i)=fluxl_diff_nf(i0)+fluxg_diff_nf(i0)
         ENDIF 
         face%twall_partition(ii)=tl2 !solid%tsol_o(kk) !<==See udfl_input_Twall_temp
         !tsol_fsw(i)=solid%tsol_o(kk)
         !tliq_fsw(i)=cell%tl_o(ii)
      ENDDO
!
1     CONTINUE       
!
!.....Constant wall temperature
!
      nv=3
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i0=istart0+i
         i1=istart+i
         i2=istart2+i
         ii=left_nf(i1)
         k=-nbcon_nf(i2)
         ali_tmp=cell%alphal(ii)
         agi_tmp=cell%alphag(ii)
!         
         tl1=cell%tl_o(ii)
         tg1=cell%tg_o(ii)
         IF(Twall_Model.eq.Laminar)THEN
            fluxl_diff_nf(i0)=ali_tmp*cell%condl(ii)*(twall(k)-tl1)*sap_nf(i1)
            fluxg_diff_nf(i0)=agi_tmp*cell%condg(ii)*(twall(k)-tg1)*sap_nf(i1)
         ELSE
            sv1=sa_nf(i1)*wallnr(ii)
            fluxl_diff_nf(i0)=ali_tmp*wcd_liq(ii)*(twall(k)-tl1)*sv1
            fluxg_diff_nf(i0)=agi_tmp*wcd_gas(ii)*(twall(k)-tg1)*sv1            
         ENDIF
!==>next
         IF(iheatpart.gt.0)THEN
            fluxl_diff_nf(i0)=face%wall_fluxl_diff(ii)
            fluxg_diff_nf(i0)=face%wall_fluxg_diff(ii)
         ELSEIF(wconden.ne.0)THEN  
!            fluxl_diff_nf(i1)=qliq_ctw(i)*sap_nf(i1)/fac_nf(i1)*djia_nf(i1)
!            fluxg_diff_nf(i1)=qgas_ctw(i)*sap_nf(i1)/fac_nf(i1)*djia_nf(i1)           
            fluxl_diff_nf(i0)=0.0d0
            fluxg_diff_nf(i0)=0.0d0 
         ENDIF
!         
      ENDDO
!
!.....Constant heat flux
!
      nv=4
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      len   =istart_nf(2,nf_number)
      IF(len.gt.0) CALL udfn_hflux_bc_profile_chw
      DO i=1,len  
         i0=istart0+i
         i1=istart+i
         i2=istart2+i
         ii=left_nf(i1)
         k=-nbcon_nf(i2)
         ali_tmp=cell%alphal(ii)
!             
         partition=heat_partition_ratio(ali_tmp)
         h_profile=hflux_bc_profile_chw(i) 
         fluxl_diff_nf(i0)=partition*h_profile*qwall_solid(k)*saa_nf(i1)
         fluxg_diff_nf(i0)=(1.0d0-partition)*h_profile*qwall_solid(k)*saa_nf(i1)
!==>next
         IF(iheatpart.gt.0)THEN
            fluxl_diff_nf(i0)=face%wall_fluxl_diff(ii)
            fluxg_diff_nf(i0)=face%wall_fluxg_diff(ii)
         ELSEIF(wconden.ne.0)THEN  
!            fluxl_diff_nf(i1)=qliq_chw(i)*sap_nf(i1)/fac_nf(i1)*djia_nf(i1)
!            fluxg_diff_nf(i1)=qgas_chw(i)*sap_nf(i1)/fac_nf(i1)*djia_nf(i1)           
            fluxl_diff_nf(i0)=0.0d0
            fluxg_diff_nf(i0)=0.0d0 
         ENDIF
!
      ENDDO
!
      CALL sum_nf(0,-1,                    &
                  fluxl_diff_nf,ediff_liq, &
                  fluxg_diff_nf,ediff_gas)
!
      IF(udfl_erg_diff) CALL udfn_erg_diff
!
      END SUBROUTINE scalar_energy_diffusion
!
!---------------------------------------------------------------------------------------
!
      FUNCTION heat_partition_ratio(alphal)
!
!     This function defines ratio of heat partition into each phases
!
!
      IMPLICIT NONE
!
!.....Input
      REAL(8) :: alphal
!.....Output
      REAL(8) :: heat_partition_ratio
!.....Local variables
      REAL(8) :: alphag
!
      SELECT CASE(2)
         CASE(1)
            heat_partition_ratio=alphal
         CASE(2)
            IF((1-alphal).lt.0.9d0)heat_partition_ratio=1.0d0
            IF((1-alphal).ge.0.9d0)heat_partition_ratio=1.0d0-10.0d0*(0.1d0-alphal)
         CASE(3)
            IF((1-alphal).lt.0.10d0)heat_partition_ratio=1.0d0
            IF((1-alphal).ge.0.10d0)heat_partition_ratio=1.0d0-1.0d0/0.8d0*(0.9d0-alphal)
            IF((1-alphal).ge.0.90d0)heat_partition_ratio=0.0d0
         CASE(4)
            alphag=1.0d0-alphal
            IF(alphag.ge.0.95d0)THEN
               heat_partition_ratio=0.0d0
            ELSEIF(alphag.ge.0.90d0)THEN
               heat_partition_ratio=1.0d0-(alphag-0.9d0)/0.05d0 
            ELSE
               heat_partition_ratio=1.0d0
            ENDIF
         END SELECT         
!
      END FUNCTION heat_partition_ratio
