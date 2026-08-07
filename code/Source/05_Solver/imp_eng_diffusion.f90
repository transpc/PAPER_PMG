!
      SUBROUTINE imp_eng_diffusion(diag_g,diag_l,diag_x,src_g,src_l,src_x,             &
                                   off_diag_g_non_i,off_diag_l_non_i,off_diag_x_non_i, &
                                   off_diag_g_non_k,off_diag_l_non_k,off_diag_x_non_k, &
                                   iter)      
!
!     Implicit energy diffusion      
!
      USE Zinterface
      USE VOL_DATA      , ONLY: cell
      USE SOLID_DATA    , ONLY: solid
      USE Wall_DATA     , ONLY: face
      USE Zmpi          , ONLY: ncell_fp
      USE Zzone         , ONLY: ncell_fluid
      USE Zcore         , ONLY: np
      USE Zparam        , ONLY: ndim
      USE Zvec_param    , ONLY: nf_nonk,nf_non,nf_inl,nf_fsw,nf_ctw,nf_chw
      USE Znum_cell     , ONLY: istart_nf,istart_nbcon_nf,                 &
                                nf_number_nb,lens,nf_number_id,istart_nfs, &
                                right_nb_k
      USE Zconst1       , ONLY: iheatpart,rv_htmodel_forCFD
      USE Zb_condition  , ONLY: alphab_gas,alphab_liq,tb_gas,tb_liq,qualab,rhob_gas
      USE Zface         , ONLY: twall_model,laminar
      USE Zgradoption   , ONLY: non_orth_eng
      USE Zncg          , ONLY: ncg_diff
      USE Zqvol         , ONLY: qwall_solid
      USE Zrv_model     , ONLY: rv_ht_w
      USE Zturb         , ONLY: wcd_liq,wcd_gas,wallnr
      USE Zuserdefined  , ONLY: hflux_bc_profile_chw
      USE Zvec_index    , ONLY: left_nf,right_non,nbcon_nf,right_fsw
      USE Zvec_geo      , ONLY: dnj_non,fac1_non,fac_non,fac1_fsw,fac_fsw, &
                                saa_nf,dji_nf,sap_nf,sa_nf
      USE Zvec_index_solid , ONLY: flux_fsw
      USE Zrv_model     , ONLY: rv_valve
!
      IMPLICIT NONE
!.....Input
      INTEGER :: iter
!.....Output
      REAL(8),DIMENSION(nf_non) :: off_diag_g_non_i,off_diag_l_non_i,off_diag_x_non_i
      REAL(8),DIMENSION(nf_nonk) :: off_diag_g_non_k,off_diag_l_non_k,off_diag_x_non_k
      REAL(8),DIMENSION(ncell_fluid) :: diag_g,diag_l,diag_x, &
                                        src_g,src_l,src_x
!.....Local variables
      INTEGER :: i,k
      INTEGER :: ii,kk
      INTEGER :: nv,nf_number,len,istart0,istart,istart2,i0,i1,i2
      REAL(8) :: dtli1,dtgi1
      REAL(8) :: dtli2,dtgi2
      REAL(8) :: dtli3,dtgi3
      REAL(8) :: cf_g,cf_l,sv1,tg1,tg2,tl1,tl2,part,part1,heat_partition_ratio,denom,cl,cg,cs
      REAL(8) :: cf_x,cf_xe,h_profile
      REAL(8) :: dtlj,dtgj
!.....Local arrays
      REAL(8),DIMENSION(ncell_fluid) :: diag_gg,diag_ll,diag_xx
      REAL(8),DIMENSION(ncell_fp,ndim) :: dtldx,dtgdx
!.....Local vector arrays
      REAL(8),DIMENSION(nf_non) :: cf_g_non,cf_l_non
      REAL(8),DIMENSION(nf_inl) :: cf_g_inl,cf_l_inl
      REAL(8),DIMENSION(nf_ctw) :: cf_g_ctw,cf_l_ctw,twall_ctw
      REAL(8),DIMENSION(nf_fsw) :: cf_g_fsw,cf_l_fsw
      REAL(8),DIMENSION(nf_non+nf_inl) :: diag_x_nf,src_x_nf
      REAL(8),DIMENSION(nf_non+nf_inl+nf_fsw+nf_ctw) :: diag_g_nf,diag_l_nf
      REAL(8),DIMENSION(nf_non+nf_inl+nf_fsw+nf_ctw+nf_chw) :: src_g_nf,src_l_nf
!
!.....Communicate needed variables
!
      IF(np.gt.1) CALL communicate_1d(cell%dtgde, &
                                      cell%dtlde)
!
!.....Calculate temperature gradient at cell center for non-orthogonal grid
!
      IF(iter.eq.1)THEN
         IF(non_orth_eng.eq.1) THEN
            CALL grad_temp(cell%tl_o,dtldx,tb_liq, &
                           cell%tg_o,dtgdx,tb_gas)
         ELSEIF(non_orth_eng.eq.2) THEN
            CALL grad_temp(cell%tl_o,dtldx,tb_liq, &  ! will be replaced with Frinl method
                           cell%tg_o,dtgdx,tb_gas)
         ENDIF
      ENDIF
      IF(non_orth_eng.gt.0)THEN
         IF(np.gt.1) CALL communicate_2d(dtldx, &
                                         dtgdx)
      ENDIF
!
!.....Build summation info for non,inl to sum for diag_x_nf,src_x
!
      nf_number_nb=1
      nf_number_id(-1)=-1
      nf_number_id(0)=0
      nf_number_id(1)=2
      istart_nfs(0)=0
      istart_nfs(1)=istart_nfs(0)+nf_non
      lens         =istart_nfs(1)+nf_inl
!
!.....Computing cells
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
         kk=right_non(i)
!
         cf_g_non(i)=( fac1_non(i)*cell%alphag_o(ii)*cell%condg(ii) &
                      +fac_non(i) *cell%alphag_o(kk)*cell%condg(kk) &
                     )*sa_nf(i1)
         cf_l_non(i)=( fac1_non(i)*cell%alphal_o(ii)*cell%condl(ii) &
                      +fac_non(i) *cell%alphal_o(kk)*cell%condl(kk) &
                     )*sa_nf(i1)
         src_g_nf(i0)=0.d0
         src_l_nf(i0)=0.d0
      ENDDO
!
!.....valve model
!       
      IF(rv_valve.eq.1) CALL valve_model_imp_end_diffusion1(cf_g_non,cf_l_non)
!
!.....Non-orthogonal grid contribution
!
      IF(non_orth_eng.eq.1.and.iter.eq.1)THEN
         IF(ndim.eq.2)THEN
            DO i=1,len  
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
!
               dtli1=fac1_non(i)*dtldx(ii,1)+fac_non(i)*dtldx(kk,1)
               dtgi1=fac1_non(i)*dtgdx(ii,1)+fac_non(i)*dtgdx(kk,1)
               dtli2=fac1_non(i)*dtldx(ii,2)+fac_non(i)*dtldx(kk,2)
               dtgi2=fac1_non(i)*dtgdx(ii,2)+fac_non(i)*dtgdx(kk,2)
               dtlj=dtli1*dnj_non(i,1)+dtli2*dnj_non(i,2)
               dtgj=dtgi1*dnj_non(i,1)+dtgi2*dnj_non(i,2)
               src_g_nf(i0)=src_g_nf(i0)+cf_g_non(i)*dtgj
               src_l_nf(i0)=src_l_nf(i0)+cf_l_non(i)*dtlj
            ENDDO
         ELSEIF(ndim.eq.3)THEN
            DO i=1,len  
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
!
               dtli1=fac1_non(i)*dtldx(ii,1)+fac_non(i)*dtldx(kk,1)
               dtgi1=fac1_non(i)*dtgdx(ii,1)+fac_non(i)*dtgdx(kk,1)
               dtli2=fac1_non(i)*dtldx(ii,2)+fac_non(i)*dtldx(kk,2)
               dtgi2=fac1_non(i)*dtgdx(ii,2)+fac_non(i)*dtgdx(kk,2)
               dtli3=fac1_non(i)*dtldx(ii,3)+fac_non(i)*dtldx(kk,3)
               dtgi3=fac1_non(i)*dtgdx(ii,3)+fac_non(i)*dtgdx(kk,3)
               dtlj=dtli1*dnj_non(i,1)+dtli2*dnj_non(i,2)+dtli3*dnj_non(i,3)
               dtgj=dtgi1*dnj_non(i,1)+dtgi2*dnj_non(i,2)+dtgi3*dnj_non(i,3)
               src_g_nf(i0)=src_g_nf(i0)+cf_g_non(i)*dtgj
               src_l_nf(i0)=src_l_nf(i0)+cf_l_non(i)*dtlj
            ENDDO
         ENDIF
      ENDIF
!
      DO i=1,len  
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)
!
         cf_g=cf_g_non(i)/dji_nf(i1)
         cf_l=cf_l_non(i)/dji_nf(i1)
!
         diag_g_nf(i0)=cf_g
         diag_l_nf(i0)=cf_l
         off_diag_g_non_i(i)=off_diag_g_non_i(i)-diag_g_nf(i0)*cell%dtgde(kk)
         off_diag_l_non_i(i)=off_diag_l_non_i(i)-diag_l_nf(i0)*cell%dtlde(kk)
      ENDDO
!
      IF(iter.eq.1) THEN
         DO i=1,len  
            i0=istart0+i
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
!
            cf_g=cf_g_non(i)/dji_nf(i1)
            cf_l=cf_l_non(i)/dji_nf(i1)
!
            src_g_nf(i0)=src_g_nf(i0)+cf_g*(cell%tg_o(kk)-cell%tg_o(ii))
            src_l_nf(i0)=src_l_nf(i0)+cf_l*(cell%tl_o(kk)-cell%tl_o(ii))
         ENDDO
      ENDIF         
!
      IF(ncg_diff.gt.0) THEN
         DO i=1,len  
            i0=istart0+i
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
!
            cf_x =( fac1_non(i)*cell%alphag_o(ii)*cell%rhog(ii)*cell%mdiff(ii) &
                   +fac_non(i) *cell%alphag_o(kk)*cell%rhog(kk)*cell%mdiff(kk) &
                  )*sap_nf(i1)
!            
            diag_x_nf(i0)=cf_x
            off_diag_x_non_i(i)=off_diag_x_non_i(i)-diag_x_nf(i0)
         ENDDO
!
!.........valve model
!         
         IF(rv_valve.eq.1) CALL valve_model_imp_end_diffusion2(diag_x_nf,off_diag_x_non_i)         
!         
         IF(iter.eq.1) THEN
            DO i=1,len  
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
!
               cf_x =( fac1_non(i)*cell%alphag_o(ii)*cell%rhog(ii)*cell%mdiff(ii) &
                      +fac_non(i) *cell%alphag_o(kk)*cell%rhog(kk)*cell%mdiff(kk) &
                     )*sap_nf(i1) 
               cf_xe= fac1_non(i)*cell%alphag_o(ii)*cell%rhog(ii)*cell%mdiff(ii)*(cell%ha(ii)-cell%hg(ii)) &
                     +fac_non(i) *cell%alphag_o(kk)*cell%rhog(kk)*cell%mdiff(kk)*(cell%ha(kk)-cell%hg(kk))
               cf_xe=cf_xe*sap_nf(i1)
!                 
               src_x_nf(i0)=cf_x*(cell%quala_o(kk)-cell%quala_o(ii))
               src_g_nf(i0)=src_g_nf(i0)+cf_xe*(cell%quala_o(kk)-cell%quala_o(ii))
            ENDDO
!
!.........valve model
!         
            IF(rv_valve.eq.1) CALL valve_model_imp_end_diffusion3(src_x_nf,src_g_nf)         
!            
         ENDIF               
      ELSE
         DO i=1,len
            i0=istart0+i
            diag_x_nf(i0)=0.d0
            src_x_nf(i0)=0.d0
         ENDDO
      ENDIF
!
      nv=-1
      nf_number=nf_number_id(nv)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         k=right_nb_k(i)
         kk=left_nf(k)
         off_diag_g_non_k(i)=off_diag_g_non_k(i)-diag_g_nf(k)*cell%dtgde(kk)
         off_diag_l_non_k(i)=off_diag_l_non_k(i)-diag_l_nf(k)*cell%dtlde(kk)
         off_diag_x_non_k(i)=off_diag_x_non_k(i)-diag_x_nf(k)
      ENDDO
!
!.....Inlet
!
      nv=1
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
         k=nbcon_nf(i2)
!
         cf_g_inl(i)=alphab_gas(k)*cell%condg(ii)*sap_nf(i1)
         cf_l_inl(i)=alphab_liq(k)*cell%condl(ii)*sap_nf(i1)
         diag_g_nf(i0)=cf_g_inl(i)
         diag_l_nf(i0)=cf_l_inl(i)
      ENDDO
!
      IF(iter.eq.1) THEN
         DO i=1,len  
            i0=istart0+i
            i1=istart+i
            i2=istart2+i
            ii=left_nf(i1)
            k=nbcon_nf(i2)
            src_g_nf(i0)=cf_g_inl(i)*(tb_gas(k)-cell%tg_o(ii))
            src_l_nf(i0)=cf_l_inl(i)*(tb_liq(k)-cell%tl_o(ii))
         ENDDO
      ENDIF         
!
      IF(ncg_diff.gt.0) THEN
         DO i=1,len  
            i0=istart0+i
            i1=istart+i
            i2=istart2+i
            ii=left_nf(i1)
            k=nbcon_nf(i2)
            cf_x=alphab_gas(k)*rhob_gas(k)*cell%mdiff(ii)*sap_nf(i1)
            cf_xe=cf_x*(cell%ha(ii)-cell%hg(ii))
            diag_x_nf(i0)=cf_x
         ENDDO
         IF(iter.eq.1) THEN
            DO i=1,len  
               i0=istart0+i
               i1=istart+i
               i2=istart2+i
               ii=left_nf(i1)
               k=nbcon_nf(i2)
               cf_x=alphab_gas(k)*rhob_gas(k)*cell%mdiff(ii)*sap_nf(i1)
               cf_xe=cf_x*(cell%ha(ii)-cell%hg(ii))
               src_x_nf(i0)=cf_x*(qualab(k)-cell%quala_o(ii))
               src_g_nf(i0)=src_g_nf(i0)+cf_xe*(qualab(k)-cell%quala_o(ii))
            ENDDO
         ENDIF            
      ELSE
         DO i=1,len  
            i0=istart0+i
            diag_x_nf(i0)=0.d0
            src_x_nf(i0)=0.d0
         ENDDO
      ENDIF
!
      IF(iter.eq.1) CALL sum_nf(1,-1,           &
                                src_x_nf,src_x)
      CALL sum_nf(0,1,               &
                  diag_x_nf,diag_xx) 
!
!.....Build summation info for fsw,ctw for diag_g,diag_l
!
      nf_number_nb=3
      nf_number_id(2)=5
      nf_number_id(3)=6
      istart_nfs(2)=istart_nfs(1)+nf_inl
      istart_nfs(3)=istart_nfs(2)+nf_fsw
      lens         =istart_nfs(3)+nf_ctw
!
!.....Fluid-Solid wall
!
      nv=2
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      IF(iheatpart.gt.0.or.rv_ht_w.eq.1)THEN
         DO i=1,len  
            i0=istart0+i
            diag_g_nf(i0)=0.d0
            diag_l_nf(i0)=0.d0
         ENDDO         
         IF(iter.eq.1)THEN
            DO i=1,len  
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               src_g_nf(i0)=face%wall_fluxg_diff(ii)
               src_l_nf(i0)=face%wall_fluxl_diff(ii)
            ENDDO
         ELSE  
            DO i=1,len  
               i0=istart0+i
               src_g_nf(i0)=0.0d0
               src_l_nf(i0)=0.0d0
            ENDDO            
         ENDIF
      ELSEIF(rv_htmodel_forCFD.gt.0)THEN
         DO i=1,len  
            i0=istart0+i
            diag_g_nf(i0)=0.d0
            diag_l_nf(i0)=0.d0
         ENDDO 
         DO i=1,len  
            i0=istart0+i
            src_g_nf(i0)=0.0d0
            src_l_nf(i0)=0.0d0
         ENDDO
      ELSE
         IF(Twall_Model.eq.Laminar)THEN
            DO i=1,len  
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
!
               sv1=sap_nf(i1)/fac_fsw(i)
               cf_g_fsw(i)=cell%alphag_o(ii)*cell%condg(ii)*sv1
               cf_l_fsw(i)=cell%alphal_o(ii)*cell%condl(ii)*sv1
            ENDDO
         ELSE
            DO i=1,len  
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
!
               sv1=sa_nf(i1)*wallnr(ii)
               cf_g_fsw(i)=cell%alphag_o(ii)*wcd_gas(ii)*sv1
               cf_l_fsw(i)=cell%alphal_o(ii)*wcd_liq(ii)*sv1
            ENDDO
         ENDIF
         DO i=1,len  
            i0=istart0+i
            i1=istart+i
            ii=left_nf(i1)
            kk=right_fsw(i)
            part=heat_partition_ratio(cell%alphal_o(ii))
            part1=1.d0-part
            denom=1.d0/( fac1_fsw(i)*(part*cell%condl(ii)+part1*cell%condg(ii)) &
                         +fac_fsw(i) * solid%conds(kk))
            cl=fac1_fsw(i)*part*cell%condl(ii)*denom
            cg=fac1_fsw(i)*part1*cell%condg(ii)*denom
            cs=fac_fsw(i)*solid%conds(kk)*solid%tsol_o(kk)*denom
!
            diag_g_nf(i0)=cf_g_fsw(i)*(1.d0-cg)
            diag_l_nf(i0)=cf_l_fsw(i)*(1.d0-cl)
         ENDDO
!
         IF(iter.eq.1)THEN
            DO i=1,len  
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               kk=right_fsw(i)
               part=heat_partition_ratio(cell%alphal_o(ii))
               part1=1.d0-part
               denom=1.d0/( fac1_fsw(i)*(part*cell%condl(ii)+part1*cell%condg(ii)) &
                            +fac_fsw(i) * solid%conds(kk))
               cl=fac1_fsw(i)*part*cell%condl(ii)*denom
               cg=fac1_fsw(i)*part1*cell%condg(ii)*denom
               cs=fac_fsw(i)*solid%conds(kk)*solid%tsol_o(kk)*denom
!
               tg1=cell%tg_o(ii)
               tl1=cell%tl_o(ii)
               tl2=cg*tg1+cs
               tg2=cl*tl1+cs
               src_g_nf(i0)=cf_g_fsw(i)*(tg2+cg*tg1-tg1)
               src_l_nf(i0)=cf_l_fsw(i)*(tl2+cl*tl1-tl1)
!
!Nuscale-03Pool
               flux_fsw(i)=src_g_nf(i0)+src_l_nf(i0)
               face%twall_partition(ii)= tg2+cg*tg1
            ENDDO
         ELSE  
            DO i=1,len  
               i0=istart0+i
               src_g_nf(i0)=0.0d0
               src_l_nf(i0)=0.0d0
            ENDDO            
         ENDIF            
      ENDIF
!
!.....Constant temp. wall
!
      CALL udfn_tw_profile(twall_ctw)
      nv=3
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      IF(iheatpart.gt.0)THEN
         DO i=1,len  
            i0=istart0+i
            diag_g_nf(i0)=0.d0
            diag_l_nf(i0)=0.d0
         ENDDO
         IF(iter.eq.1)THEN
            DO i=1,len  
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               src_g_nf(i0)=face%wall_fluxg_diff(ii)
               src_l_nf(i0)=face%wall_fluxl_diff(ii)
            ENDDO
         ENDIF
      ELSE
         IF(Twall_Model.eq.Laminar)THEN
            DO i=1,len  
               i1=istart+i
               ii=left_nf(i1)
               cf_g_ctw(i)=cell%alphag_o(ii)*cell%condg(ii)*sap_nf(i1)
               cf_l_ctw(i)=cell%alphal_o(ii)*cell%condl(ii)*sap_nf(i1)
            ENDDO
         ELSE
            DO i=1,len  
               i1=istart+i
               ii=left_nf(i1)
               sv1=sa_nf(i1)*wallnr(ii)
               cf_g_ctw(i)=cell%alphag_o(ii)*wcd_gas(ii)*sv1
               cf_l_ctw(i)=cell%alphal_o(ii)*wcd_liq(ii)*sv1
            ENDDO
         ENDIF
         DO i=1,len
            i0=istart0+i
            diag_g_nf(i0)=cf_g_ctw(i)
            diag_l_nf(i0)=cf_l_ctw(i)
         ENDDO
         IF(iter.eq.1)THEN
            DO i=1,len
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               src_g_nf(i0)=cf_g_ctw(i)*(twall_ctw(i)-cell%tg_o(ii))
               src_l_nf(i0)=cf_l_ctw(i)*(twall_ctw(i)-cell%tl_o(ii))
            ENDDO
         ENDIF            
      ENDIF
      CALL sum_nf(0,1,               &
                  diag_g_nf,diag_gg, &
                  diag_l_nf,diag_ll) 
!
!.....Build summation info for fsw,ctw,chw for src_g,src_l
!
      nf_number_nb=4
      nf_number_id(4)=7
      istart_nfs(4)=istart_nfs(3)+nf_ctw
      lens         =istart_nfs(4)+nf_chw
!
!.....Constant heat flux wall
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
!
         src_g_nf(i0)=0.d0
         src_l_nf(i0)=0.d0
         IF(iter.eq.1)THEN
            IF(iheatpart.gt.0)THEN
               src_g_nf(i0)=face%wall_fluxg_diff(ii)
               src_l_nf(i0)=face%wall_fluxl_diff(ii)
            ELSE
               part=heat_partition_ratio(cell%alphal_o(ii))
               h_profile=hflux_bc_profile_chw(i)
               src_g_nf(i0)=(1.d0-part)*h_profile*qwall_solid(k)*saa_nf(i1)
               src_l_nf(i0)=part*h_profile*qwall_solid(k)*saa_nf(i1)
            ENDIF
         ENDIF
      ENDDO
      IF(iter.eq.1) CALL sum_nf(1,-1,           &
                                src_g_nf,src_g, &
                                src_l_nf,src_l)
!
      DO i=1,ncell_fluid
         diag_g(i)=diag_g(i)+diag_gg(i)*cell%dtgde(i)
         diag_l(i)=diag_l(i)+diag_ll(i)*cell%dtlde(i)
         diag_x(i)=diag_x(i)+diag_xx(i)
      ENDDO
!
      END SUBROUTINE imp_eng_diffusion
