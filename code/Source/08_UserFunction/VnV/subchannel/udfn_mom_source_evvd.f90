      SUBROUTINE udfn_mom_source_evvd

!     Modifies the momentum source terms at the free surface cells

      USE Zinterface
      USE VOL_DATA     , ONLY: cell
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: ndim
      USE Zvec_param   , ONLY: nf_nonk,nf_non
      USE Znum_cell    , ONLY: right_nb_k,istart_nf,                     &
                               nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Zconst2      , ONLY: hydraulicd
      USE Zporous      , ONLY: tm_mas_l,tm_mas_g,tm_eng_l,tm_eng_g,vd_mas_l,vd_mas_g,vd_eng_l,vd_eng_g,cell_area, &
                               s_ij_non_i,s_ij_non_k
      USE Zporous      , ONLY: l_2p_multiplier_evvd
      USE Zporous      , ONLY: ka,beta
      USE Zvector      , ONLY: vl_o,vg_o
      USE Zm_src       , ONLY: src_gas,src_liq
      ! OPR1000 rod-scale
      USE Zconst1      , ONLY: vv_prob
      USE Zporous      , ONLY: chn_type
      USE Zvec_geo     , ONLY: xn_nf
!
      IMPLICIT NONE      
!.....Local variables
      INTEGER :: i,k
      INTEGER :: ii,kk
      INTEGER :: nv,nf_number,istart0,istart,len,i0,i1
      REAL(8) :: tm_mom_l,vd_mom_l                          !LSJ porous   
      REAL(8) :: tm_mom_g,vd_mom_g                          !LSJ porous 
      REAL(8) :: arj_g,ari_g,aruj_g,arui_g,arej_g,arei_g
      REAL(8) :: arj_l,ari_l,aruj_l,arui_l,arej_l,arei_l
      REAL(8) :: arj_lvd,ari_lvd,aruj_lvd,arui_lvd,arej_lvd,arei_lvd
      REAL(8) :: g_i,g_k,g_avg,rho_avg
      REAL(8) :: rho_avgi,rho_avgk
      REAL(8) :: TPM_M,hd,x_M,x_0,x_x,Rel,velo
      REAL(8) :: beta_i,beta_k
      REAL(8) xn3
      REAL(8) :: vt
!.....Local arrays
      REAL(8),DIMENSION(ncell_fp) :: g_ik,rho_avgik,arij_l,arij_lvd,arij_g
      REAL(8),DIMENSION(ncell_fluid) :: TPM
      REAL(8),DIMENSION(ncell_fluid,ndim) :: turb_mixing_l_mom,turb_mixing_g_mom, &
                                             void_drift_l_mom,void_drift_g_mom
!.....Local vector arrays
      REAL(8),DIMENSION(nf_non) :: turb_mixing_l_mas_non,turb_mixing_g_mas_non, &
                                   void_drift_l_mas_non,void_drift_g_mas_non,   &
                                   turb_mixing_l_eng_non,turb_mixing_g_eng_non, &
                                   void_drift_l_eng_non,void_drift_g_eng_non
      REAL(8),DIMENSION(nf_non,ndim) :: turb_mixing_l_mom_non,turb_mixing_g_mom_non, &
                                        void_drift_l_mom_non,void_drift_g_mom_non
      REAL(8),DIMENSION(nf_non) :: vt_non_i,vt_non_k
      REAL(8),DIMENSION(nf_nonk+nf_non) :: turb_mixing_l_mas_nf, &
                                           turb_mixing_g_mas_nf, &
                                           turb_mixing_l_eng_nf, &
                                           turb_mixing_g_eng_nf, &
                                           void_drift_l_mas_nf,  &
                                           void_drift_g_mas_nf,  &
                                           void_drift_l_eng_nf,  &
                                           void_drift_g_eng_nf
      REAL(8),DIMENSION(nf_nonk+nf_non,ndim) :: turb_mixing_l_mom_nf, &
                                                turb_mixing_g_mom_nf, &
                                                void_drift_l_mom_nf,  &
                                                void_drift_g_mom_nf
!
      IF(vv_prob.eq.'OPR1000_fullcore_modmesh02_rv'        .or. &
         vv_prob.eq.'OPR1000_fullcore_modmesh02_rv_vessel' .or. &
         vv_prob.eq.'OPR1000_single_assem'                 .or. &
         vv_prob.eq.'KSMR'                                        ) THEN   !PSH
         CALL udfn_mom_source_evvd_rod
         RETURN
      ENDIF
        
!
!.....Two-phase multiplier
!
      IF(l_2p_multiplier_evvd) THEN
         IF(ndim.eq.2) THEN
            DO i=1,ncell_fluid
               hd=hydraulicd(i)
               g_i=cell%alphag(i)*cell%rhog(i)*vg_o(i,ndim)+cell%alphal(i)*cell%rhol(i)*vl_o(i,ndim)
               velo=SQRT(vl_o(i,1)**2+vl_o(i,2)**2)
               Rel=MAX(1.d0,cell%rhol(i)*velo*hd/cell%lviscosl(i))
               x_x=cell%quals(i)
               x_M= ((0.4d0*SQRT(cell%rhol(i)*(cell%rhol(i)-cell%rhog(i))*9.8d0*hd))/g_i+0.6d0)  &
                   /(SQRT(cell%rhol(i)/cell%rhog(i))+0.6d0) !Wallis Model
               TPM_M=5.d0 !Faya Model 
               IF(x_x.lt.x_M) THEN
                  TPM(i)=1.d0+(TPM_M-1.d0)*x_x/x_M
               ELSE
                  x_0=0.75d0*x_M*Rel**0.0417
                  TPM(i)=1.d0+(TPM_M-1.d0)*(x_M-x_0)/(x_x-x_0)
               ENDIF               
            ENDDO
         ELSE
            DO i=1,ncell_fluid
               hd=hydraulicd(i)
               g_i=cell%alphag(i)*cell%rhog(i)*vg_o(i,ndim)+cell%alphal(i)*cell%rhol(i)*vl_o(i,ndim)
               velo=SQRT(vl_o(i,1)**2+vl_o(i,2)**2+vl_o(i,3)**2)
               Rel=MAX(1.d0,cell%rhol(i)*velo*hd/cell%lviscosl(i))
               x_x=cell%quals(i)
               x_M= ((0.4d0*SQRT(cell%rhol(i)*(cell%rhol(i)-cell%rhog(i))*9.8d0*hd))/g_i+0.6d0)  &
                   /(SQRT(cell%rhol(i)/cell%rhog(i))+0.6d0) !Wallis Model
               TPM_M=5.d0 !Faya Model 
               IF(x_x.lt.x_M) THEN
                  TPM(i)=1.d0+(TPM_M-1.d0)*x_x/x_M
               ELSE
                  x_0=0.75d0*x_M*Rel**0.0417
                  TPM(i)=1.d0+(TPM_M-1.d0)*(x_M-x_0)/(x_x-x_0)
               ENDIF               
            ENDDO
         ENDIF
      ELSE
         DO i=1,ncell_fluid
            TPM(i)=1.d0
         ENDDO
      ENDIF
!
!.....Build summation info for non
!
      nf_number_nb=0
      nf_number_id(-1)=-1
      nf_number_id(0)=0
      istart_nfs(-1)=0
      istart_nfs(0)=istart_nfs(-1)+nf_nonk
      lens         =istart_nfs(0) +nf_non
!
!.....> check that all those cell%..(kk) are communicated
      DO i=1,ncell_fp
         g_ik(i)=cell%alphag(i)*cell%rhog(i)*vg_o(i,ndim)+cell%alphal(i)*cell%rhol(i)*vl_o(i,ndim)
         rho_avgik(i)=(cell%alphag(i)*cell%rhog(i)+cell%alphal(i)*cell%rhol(i))*cell_area(i)
         arij_l(i)  =(1.d0-cell%alphag(i))*cell%rhol(i)
         arij_lvd(i)=      cell%alphag(i) *cell%rhol(i)
         arij_g(i)=        cell%alphag(i) *cell%rhog(i)
      ENDDO
!
!.....Cells non
!
      nf_number=0
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
!
      DO i=1,len
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)
         xn3=xn_nf(i1,3)
            beta_i=beta
            beta_k=beta
!...........g_avg
            g_i=g_ik(ii)
            g_k=g_ik(kk)
            g_avg=(g_i*cell_area(ii)+g_k*cell_area(kk))/(cell_area(ii)+cell_area(kk))
!...........rho_avg
            rho_avgi=rho_avgik(ii)
            rho_avgk=rho_avgik(kk)
            rho_avg=(rho_avgi+rho_avgk)/(cell_area(ii)+cell_area(kk))
!...........vt with sgap
            vt_non_i(i)=s_ij_non_i(i)*(beta_i*g_avg/rho_avg)
            vt_non_k(i)=s_ij_non_k(i)*(beta_k*g_avg/rho_avg)
!
            ari_g  =arij_g(ii)  !PSH
            ari_l  =arij_l(ii)
            ari_lvd=arij_lvd(ii)
            arei_g=ari_g*cell%hg(ii)
            arei_l=ari_l*cell%hl(ii)
            arei_lvd=ari_lvd*cell%hl(ii)
            
            arj_g  =arij_g(kk)
            arj_l  =arij_l(kk)
            arj_lvd=arij_lvd(kk)
            arej_g=arj_g*cell%hg(kk)
            arej_l=arj_l*cell%hl(kk)
            arej_lvd=arj_lvd*cell%hl(kk)         
!            
            turb_mixing_l_mas_non(i)=    (arj_l  -ari_l )
            turb_mixing_g_mas_non(i)=    (arj_g  -ari_g )
            turb_mixing_l_eng_non(i)=    (arej_l -arei_l)
            turb_mixing_g_eng_non(i)=    (arej_g -arei_g)

            void_drift_l_mas_non(i) =-ka*(arj_lvd +ari_lvd )*(g_i-g_k)/(g_i+g_k)
            void_drift_g_mas_non(i) = ka*(arj_g   +ari_g   )*(g_i-g_k)/(g_i+g_k)
            void_drift_l_eng_non(i) =-ka*(arej_lvd+arei_lvd)*(g_i-g_k)/(g_i+g_k)
            void_drift_g_eng_non(i) = ka*(arej_g  +arei_g  )*(g_i-g_k)/(g_i+g_k)
!
            arui_g=ari_g*vg_o(ii,1)
            aruj_g=arj_g*vg_o(kk,1)
            arui_l=ari_l*vl_o(ii,1)
            aruj_l=arj_l*vl_o(kk,1)
            arui_lvd=ari_lvd*vl_o(ii,1)
            aruj_lvd=arj_lvd*vl_o(kk,1)         
            turb_mixing_l_mom_non(i,1)=    (aruj_l  -arui_l  )
            turb_mixing_g_mom_non(i,1)=    (aruj_g  -arui_g  )
            void_drift_l_mom_non(i,1) =-ka*(aruj_lvd+arui_lvd)*(g_i-g_k)/(g_i+g_k)
            void_drift_g_mom_non(i,1) = ka*(aruj_g  +arui_g  )*(g_i-g_k)/(g_i+g_k)
!
            arui_g=ari_g*vg_o(ii,2)
            aruj_g=arj_g*vg_o(kk,2)
            arui_l=ari_l*vl_o(ii,2)
            aruj_l=arj_l*vl_o(kk,2)
            arui_lvd=ari_lvd*vl_o(ii,2)
            aruj_lvd=arj_lvd*vl_o(kk,2)         
            turb_mixing_l_mom_non(i,2)=    (aruj_l  -arui_l  )
            turb_mixing_g_mom_non(i,2)=    (aruj_g  -arui_g  )
            void_drift_l_mom_non(i,2) =-ka*(aruj_lvd+arui_lvd)*(g_i-g_k)/(g_i+g_k)
            void_drift_g_mom_non(i,2) = ka*(aruj_g  +arui_g  )*(g_i-g_k)/(g_i+g_k)
!
            arui_g=ari_g*vg_o(ii,3)
            aruj_g=arj_g*vg_o(kk,3)
            arui_l=ari_l*vl_o(ii,3)
            aruj_l=arj_l*vl_o(kk,3)
            arui_lvd=ari_lvd*vl_o(ii,3)
            aruj_lvd=arj_lvd*vl_o(kk,3)         
            turb_mixing_l_mom_non(i,3)=    (aruj_l  -arui_l  )
            turb_mixing_g_mom_non(i,3)=    (aruj_g  -arui_g  )
            void_drift_l_mom_non(i,3) =-ka*(aruj_lvd+arui_lvd)*(g_i-g_k)/(g_i+g_k)
            void_drift_g_mom_non(i,3) = ka*(aruj_g  +arui_g  )*(g_i-g_k)/(g_i+g_k)
      ENDDO
!
      nv=0
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart =istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)
         xn3=xn_nf(i1,3)
         IF(chn_type(ii).ne.0 .and. ABS(xn3).le.0.5d0) THEN
            vt=vt_non_i(i)
            turb_mixing_l_mas_nf(i0)=vt*turb_mixing_l_mas_non(i)
            turb_mixing_g_mas_nf(i0)=vt*turb_mixing_g_mas_non(i)
            turb_mixing_l_eng_nf(i0)=vt*turb_mixing_l_eng_non(i)
            turb_mixing_g_eng_nf(i0)=vt*turb_mixing_g_eng_non(i)
            void_drift_l_mas_nf(i0)=vt*void_drift_l_mas_non(i)
            void_drift_g_mas_nf(i0)=vt*void_drift_g_mas_non(i)
            void_drift_l_eng_nf(i0)=vt*void_drift_l_eng_non(i)
            void_drift_g_eng_nf(i0)=vt*void_drift_g_eng_non(i)
!
            turb_mixing_l_mom_nf(i0,1)=vt*turb_mixing_l_mom_non(i,1)
            turb_mixing_g_mom_nf(i0,1)=vt*turb_mixing_g_mom_non(i,1)
            void_drift_l_mom_nf(i0,1) =vt*void_drift_l_mom_non(i,1)
            void_drift_g_mom_nf(i0,1) =vt*void_drift_g_mom_non(i,1)
!
            turb_mixing_l_mom_nf(i0,2)=vt*turb_mixing_l_mom_non(i,2)
            turb_mixing_g_mom_nf(i0,2)=vt*turb_mixing_g_mom_non(i,2)
            void_drift_l_mom_nf(i0,2) =vt*void_drift_l_mom_non(i,2)
            void_drift_g_mom_nf(i0,2) =vt*void_drift_g_mom_non(i,2)
!
            turb_mixing_l_mom_nf(i0,3)=vt*turb_mixing_l_mom_non(i,3)
            turb_mixing_g_mom_nf(i0,3)=vt*turb_mixing_g_mom_non(i,3)
            void_drift_l_mom_nf(i0,3) =vt*void_drift_l_mom_non(i,3)
            void_drift_g_mom_nf(i0,3) =vt*void_drift_g_mom_non(i,3)
         ELSE
            turb_mixing_l_mas_nf(i0)=0.d0
            turb_mixing_g_mas_nf(i0)=0.d0
            turb_mixing_l_eng_nf(i0)=0.d0
            turb_mixing_g_eng_nf(i0)=0.d0
            void_drift_l_mas_nf(i0)=0.d0
            void_drift_g_mas_nf(i0)=0.d0
            void_drift_l_eng_nf(i0)=0.d0
            void_drift_g_eng_nf(i0)=0.d0
!
            turb_mixing_l_mom_nf(i0,1)=0.d0
            turb_mixing_g_mom_nf(i0,1)=0.d0
            void_drift_l_mom_nf(i0,1)=0.d0
            void_drift_g_mom_nf(i0,1)=0.d0
!
            turb_mixing_l_mom_nf(i0,2)=0.d0
            turb_mixing_g_mom_nf(i0,2)=0.d0
            void_drift_l_mom_nf(i0,2)=0.d0
            void_drift_g_mom_nf(i0,2)=0.d0
!
            turb_mixing_l_mom_nf(i0,3)=0.d0
            turb_mixing_g_mom_nf(i0,3)=0.d0
            void_drift_l_mom_nf(i0,3)=0.d0
            void_drift_g_mom_nf(i0,3)=0.d0
         ENDIF
      ENDDO
!
      nv=-1
      nf_number=nf_number_id(nv)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         k=right_nb_k(i)
         ii=right_non(k)
         xn3=xn_nf(k,3)
         IF(chn_type(ii).ne.0 .and. ABS(xn3).le.0.5d0) THEN
            vt=-vt_non_k(k)
            turb_mixing_l_mas_nf(i)=vt*turb_mixing_l_mas_non(k)
            turb_mixing_g_mas_nf(i)=vt*turb_mixing_g_mas_non(k)
            turb_mixing_l_eng_nf(i)=vt*turb_mixing_l_eng_non(k)
            turb_mixing_g_eng_nf(i)=vt*turb_mixing_g_eng_non(k)
            void_drift_l_mas_nf(i)=vt*void_drift_l_mas_non(k)
            void_drift_g_mas_nf(i)=vt*void_drift_g_mas_non(k)
            void_drift_l_eng_nf(i)=vt*void_drift_l_eng_non(k)
            void_drift_g_eng_nf(i)=vt*void_drift_g_eng_non(k)
!
            turb_mixing_l_mom_nf(i,1)=vt*turb_mixing_l_mom_non(k,1)
            turb_mixing_g_mom_nf(i,1)=vt*turb_mixing_g_mom_non(k,1)
            void_drift_l_mom_nf(i,1)=vt*void_drift_l_mom_non(k,1)
            void_drift_g_mom_nf(i,1)=vt*void_drift_g_mom_non(k,1)
!
            turb_mixing_l_mom_nf(i,2)=vt*turb_mixing_l_mom_non(k,2)
            turb_mixing_g_mom_nf(i,2)=vt*turb_mixing_g_mom_non(k,2)
            void_drift_l_mom_nf(i,2)=vt*void_drift_l_mom_non(k,2)
            void_drift_g_mom_nf(i,2)=vt*void_drift_g_mom_non(k,2)
!
            turb_mixing_l_mom_nf(i,3)=vt*turb_mixing_l_mom_non(k,3)
            turb_mixing_g_mom_nf(i,3)=vt*turb_mixing_g_mom_non(k,3)
            void_drift_l_mom_nf(i,3)=vt*void_drift_l_mom_non(k,3)
            void_drift_g_mom_nf(i,3)=vt*void_drift_g_mom_non(k,3)
         ELSE
            turb_mixing_l_mas_nf(i)=0.d0
            turb_mixing_g_mas_nf(i)=0.d0
            turb_mixing_l_eng_nf(i)=0.d0
            turb_mixing_g_eng_nf(i)=0.d0
            void_drift_l_mas_nf(i)=0.d0
            void_drift_g_mas_nf(i)=0.d0
            void_drift_l_eng_nf(i)=0.d0
            void_drift_g_eng_nf(i)=0.d0
!
            turb_mixing_l_mom_nf(i,1)=0.d0
            turb_mixing_g_mom_nf(i,1)=0.d0
            void_drift_l_mom_nf(i,1)=0.d0
            void_drift_g_mom_nf(i,1)=0.d0
!
            turb_mixing_l_mom_nf(i,2)=0.d0
            turb_mixing_g_mom_nf(i,2)=0.d0
            void_drift_l_mom_nf(i,2)=0.d0
            void_drift_g_mom_nf(i,2)=0.d0
!
            turb_mixing_l_mom_nf(i,3)=0.d0
            turb_mixing_g_mom_nf(i,3)=0.d0
            void_drift_l_mom_nf(i,3)=0.d0
            void_drift_g_mom_nf(i,3)=0.d0
         ENDIF
      ENDDO
!
      CALL sum_nf(0,0,                           &
                  turb_mixing_l_mas_nf,tm_mas_l, &
                  turb_mixing_g_mas_nf,tm_mas_g, &
                  turb_mixing_l_eng_nf,tm_eng_l, &
                  turb_mixing_g_eng_nf,tm_eng_g, &
                  void_drift_l_mas_nf, vd_mas_l, &
                  void_drift_g_mas_nf, vd_mas_g, &
                  void_drift_l_eng_nf, vd_eng_l, &
                  void_drift_g_eng_nf, vd_eng_g)
!
      CALL sum_nf_ndim(0,0,ncell_fluid,                        &
                       turb_mixing_l_mom_nf,turb_mixing_l_mom, &
                       turb_mixing_g_mom_nf,turb_mixing_g_mom, &
                       void_drift_l_mom_nf,void_drift_l_mom,   &
                       void_drift_g_mom_nf,void_drift_g_mom)
!
      DO i=1,ncell_fluid
         tm_mas_l(i)=tm_mas_l(i)*TPM(i)/cell_area(i)
         tm_mas_g(i)=tm_mas_g(i)*TPM(i)/cell_area(i)
         tm_eng_l(i)=tm_eng_l(i)*TPM(i)/cell_area(i)
         tm_eng_g(i)=tm_eng_g(i)*TPM(i)/cell_area(i)
         vd_mas_l(i)=vd_mas_l(i)*TPM(i)/cell_area(i)
         vd_mas_g(i)=vd_mas_g(i)*TPM(i)/cell_area(i)
         vd_eng_l(i)=vd_eng_l(i)*TPM(i)/cell_area(i)
         vd_eng_g(i)=vd_eng_g(i)*TPM(i)/cell_area(i)
      ENDDO

      DO i=1,ncell_fluid
         tm_mom_l=turb_mixing_l_mom(i,1)*TPM(i)/cell_area(i)
         tm_mom_g=turb_mixing_g_mom(i,1)*TPM(i)/cell_area(i)
         vd_mom_l=void_drift_l_mom(i,1)*TPM(i)/cell_area(i)
         vd_mom_g=void_drift_g_mom(i,1)*TPM(i)/cell_area(i)
         src_liq(i,1)=src_liq(i,1)+(tm_mom_l+vd_mom_l)
         src_gas(i,1)=src_gas(i,1)+(tm_mom_g+vd_mom_g)
!
         tm_mom_l=turb_mixing_l_mom(i,2)*TPM(i)/cell_area(i)
         tm_mom_g=turb_mixing_g_mom(i,2)*TPM(i)/cell_area(i)
         vd_mom_l=void_drift_l_mom(i,2)*TPM(i)/cell_area(i)
         vd_mom_g=void_drift_g_mom(i,2)*TPM(i)/cell_area(i)
         src_liq(i,2)=src_liq(i,2)+(tm_mom_l+vd_mom_l)
         src_gas(i,2)=src_gas(i,2)+(tm_mom_g+vd_mom_g)
!
         tm_mom_l=turb_mixing_l_mom(i,3)*TPM(i)/cell_area(i)
         tm_mom_g=turb_mixing_g_mom(i,3)*TPM(i)/cell_area(i)
         vd_mom_l=void_drift_l_mom(i,3)*TPM(i)/cell_area(i)
         vd_mom_g=void_drift_g_mom(i,3)*TPM(i)/cell_area(i)
         src_liq(i,3)=src_liq(i,3)+(tm_mom_l+vd_mom_l)
         src_gas(i,3)=src_gas(i,3)+(tm_mom_g+vd_mom_g)
      ENDDO
!
      END SUBROUTINE udfn_mom_source_evvd
!
      SUBROUTINE udfn_mom_source_evvd_old

!     Modifies the momentum source terms at the free surface cells

      USE VOL_DATA     , ONLY: cell
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: ndim
      USE Zvec_param   , ONLY: nf_non
      USE Znum_cell    , ONLY: istart_nf
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Zconst2      , ONLY: hydraulicd
      USE Zporous      , ONLY: tm_mas_l,tm_mas_g,tm_eng_l,tm_eng_g,vd_mas_l,vd_mas_g,vd_eng_l,vd_eng_g,cell_area, &
                               s_ij_non_i,s_ij_non_k
      USE Zporous      , ONLY: l_2p_multiplier_evvd
      USE Zporous      , ONLY: ka,beta
      USE Zvector      , ONLY: vl_o,vg_o
      USE Zm_src       , ONLY: src_gas,src_liq
      ! OPR1000 rod-scale
      USE Zconst1      , ONLY: vv_prob
      USE Zporous      , ONLY: chn_type
      USE Zvec_geo     , ONLY: xn_nf
!
      IMPLICIT NONE      
!.....Local variables
      INTEGER :: i,ix
      INTEGER :: ii,kk
      INTEGER :: nf_number,istart,len,i1
      REAL(8) :: tm_mom_l,vd_mom_l                          !LSJ porous   
      REAL(8) :: tm_mom_g,vd_mom_g                          !LSJ porous 
      REAL(8) :: arj_g,ari_g,aruj_g,arui_g,arej_g,arei_g
      REAL(8) :: arj_l,ari_l,aruj_l,arui_l,arej_l,arei_l
      REAL(8) :: arj_lvd,ari_lvd,aruj_lvd,arui_lvd,arej_lvd,arei_lvd
      REAL(8) :: g_i,g_k,g_avg,rho_avg
      REAL(8) :: rho_avgi,rho_avgk
      REAL(8) :: TPM_M,hd,x_M,x_0,x_x,Rel,velo
!.....Local arrays
      REAL(8) :: TPM(ncell_fluid)
      REAL(8) :: turb_mixing_l_mom(ncell_fluid,ndim),turb_mixing_g_mom(ncell_fluid,ndim)
      REAL(8) :: void_drift_l_mom(ncell_fluid,ndim),void_drift_g_mom(ncell_fluid,ndim)
!.....Local vector arrays
      REAL(8) :: beta_i(nf_non),beta_k(nf_non)
      REAL(8) :: turb_mixing_l_mas_non(nf_non),turb_mixing_g_mas_non(nf_non), &
                 void_drift_l_mas_non(nf_non),void_drift_g_mas_non(nf_non)
      REAL(8) :: turb_mixing_l_eng_non(nf_non),turb_mixing_g_eng_non(nf_non), &
                 void_drift_l_eng_non(nf_non),void_drift_g_eng_non(nf_non)
      REAL(8) :: turb_mixing_l_mom_non(nf_non,ndim),turb_mixing_g_mom_non(nf_non,ndim), &
                 void_drift_l_mom_non(nf_non,ndim),void_drift_g_mom_non(nf_non,ndim)
      REAL(8) :: vt_non_i(nf_non),vt_non_k(nf_non)
      REAL(8) xn3
!
      IF(vv_prob.eq.'OPR1000_fullcore_modmesh02_rv'        .or. &
         vv_prob.eq.'OPR1000_fullcore_modmesh02_rv_vessel' .or. &
         vv_prob.eq.'OPR1000_single_assem'                       ) THEN
         CALL udfn_mom_source_evvd_rod
         RETURN
      ENDIF
        
!
!.....Two-phase multiplier
!
      IF(l_2p_multiplier_evvd) THEN
         IF(ndim.eq.2) THEN
            DO i=1,ncell_fluid
               hd=hydraulicd(i)
               g_i=cell%alphag(i)*cell%rhog(i)*vg_o(i,ndim)+cell%alphal(i)*cell%rhol(i)*vl_o(i,ndim)
               velo=SQRT(vl_o(i,1)**2+vl_o(i,2)**2)
               Rel=MAX(1.d0,cell%rhol(i)*velo*hd/cell%lviscosl(i))
               x_x=cell%quals(i)
               x_M= ((0.4d0*SQRT(cell%rhol(i)*(cell%rhol(i)-cell%rhog(i))*9.8d0*hd))/g_i+0.6d0)  &
                   /(SQRT(cell%rhol(i)/cell%rhog(i))+0.6d0) !Wallis Model
               TPM_M=5.d0 !Faya Model 
               IF(x_x.lt.x_M) THEN
                  TPM(i)=1.d0+(TPM_M-1.d0)*x_x/x_M
               ELSE
                  x_0=0.75d0*x_M*Rel**0.0417
                  TPM(i)=1.d0+(TPM_M-1.d0)*(x_M-x_0)/(x_x-x_0)
               ENDIF               
            ENDDO
         ELSE
            DO i=1,ncell_fluid
               hd=hydraulicd(i)
               g_i=cell%alphag(i)*cell%rhog(i)*vg_o(i,ndim)+cell%alphal(i)*cell%rhol(i)*vl_o(i,ndim)
               velo=SQRT(vl_o(i,1)**2+vl_o(i,2)**2+vl_o(i,3)**2)
               Rel=MAX(1.d0,cell%rhol(i)*velo*hd/cell%lviscosl(i))
               x_x=cell%quals(i)
               x_M= ((0.4d0*SQRT(cell%rhol(i)*(cell%rhol(i)-cell%rhog(i))*9.8d0*hd))/g_i+0.6d0)  &
                   /(SQRT(cell%rhol(i)/cell%rhog(i))+0.6d0) !Wallis Model
               TPM_M=5.d0 !Faya Model 
               IF(x_x.lt.x_M) THEN
                  TPM(i)=1.d0+(TPM_M-1.d0)*x_x/x_M
               ELSE
                  x_0=0.75d0*x_M*Rel**0.0417
                  TPM(i)=1.d0+(TPM_M-1.d0)*(x_M-x_0)/(x_x-x_0)
               ENDIF               
            ENDDO
         ENDIF
      ELSE
         DO i=1,ncell_fluid
            TPM(i)=1.d0
         ENDDO
      ENDIF
!
!.....Cells non
!
      nf_number=0
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
!
      DO i=1,len
         beta_i(i)=beta
         beta_k(i)=beta
      ENDDO
!
      DO i=1,len
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)
         xn3=xn_nf(i1,3)
!========> wrong IF not valid for kk
         IF(chn_type(ii).eq.0)CYCLE
         IF(DABS(xn3).gt.0.5d0)CYCLE
!        g_avg
         g_i=cell%alphag(ii)*cell%rhog(ii)*vg_o(ii,ndim)+cell%alphal(ii)*cell%rhol(ii)*vl_o(ii,ndim)
         g_k=cell%alphag(kk)*cell%rhog(kk)*vg_o(kk,ndim)+cell%alphal(kk)*cell%rhol(kk)*vl_o(kk,ndim)
         g_avg=(g_i*cell_area(ii)+g_k*cell_area(kk))/(cell_area(ii)+cell_area(kk))
!        rho_avg
         rho_avgi=(cell%alphag(ii)*cell%rhog(ii)+cell%alphal(ii)*cell%rhol(ii))*cell_area(ii)
         rho_avgk=(cell%alphag(kk)*cell%rhog(kk)+cell%alphal(kk)*cell%rhol(kk))*cell_area(kk)
         rho_avg=(rho_avgi+rho_avgk)/(cell_area(ii)+cell_area(kk))
!        vt with sgap
         vt_non_i(i)=s_ij_non_i(i)*(beta_i(i)*g_avg/rho_avg)
         vt_non_k(i)=s_ij_non_k(i)*(beta_k(i)*g_avg/rho_avg)
!
         ari_g=      cell%alphag(ii) *cell%rhog(ii)            
         ari_l=(1.d0-cell%alphag(ii))*cell%rhol(ii)
         ari_lvd=    cell%alphag(ii) *cell%rhol(ii)
         arei_g=ari_g*cell%hg(ii)
         arei_l=ari_l*cell%hl(ii)
         arei_lvd=ari_lvd*cell%hl(ii)
         arj_g=      cell%alphag(kk) *cell%rhog(kk)
         arj_l=(1.d0-cell%alphag(kk))*cell%rhol(kk)
         arj_lvd=    cell%alphag(kk) *cell%rhol(kk)         
         arej_g=arj_g*cell%hg(kk)
         arej_l=arj_l*cell%hl(kk)
         arej_lvd=arj_lvd*cell%hl(kk)         
!            
         turb_mixing_l_mas_non(i)=    (arj_l  -ari_l )
         turb_mixing_g_mas_non(i)=    (arj_g  -ari_g )
         turb_mixing_l_eng_non(i)=    (arej_l -arei_l)
         turb_mixing_g_eng_non(i)=    (arej_g -arei_g)

         void_drift_l_mas_non(i) =-ka*(arj_lvd +ari_lvd )*(g_i-g_k)/(g_i+g_k)
         void_drift_g_mas_non(i) = ka*(arj_g   +ari_g   )*(g_i-g_k)/(g_i+g_k)
         void_drift_l_eng_non(i) =-ka*(arej_lvd+arei_lvd)*(g_i-g_k)/(g_i+g_k)
         void_drift_g_eng_non(i) = ka*(arej_g  +arei_g  )*(g_i-g_k)/(g_i+g_k)
!
         DO ix=1,ndim
            arui_g=ari_g*vg_o(ii,ix)
            aruj_g=arj_g*vg_o(kk,ix)
            arui_l=ari_l*vl_o(ii,ix)
            aruj_l=arj_l*vl_o(kk,ix)
            arui_lvd=ari_lvd*vl_o(ii,ix)
            aruj_lvd=arj_lvd*vl_o(kk,ix)         
!
            turb_mixing_l_mom_non(i,ix)=    (aruj_l  -arui_l  )
            turb_mixing_g_mom_non(i,ix)=    (aruj_g  -arui_g  )
            void_drift_l_mom_non(i,ix) =-ka*(aruj_lvd+arui_lvd)*(g_i-g_k)/(g_i+g_k)
            void_drift_g_mom_non(i,ix) = ka*(aruj_g  +arui_g  )*(g_i-g_k)/(g_i+g_k)
         ENDDO
      ENDDO
!
!     CALL sum_nf0_8v_ndim_4v(turb_mixing_l_mas_non,turb_mixing_g_mas_non,turb_mixing_l_eng_non,turb_mixing_g_eng_non, &
!                             void_drift_l_mas_non,void_drift_g_mas_non,void_drift_l_eng_non,void_drift_g_eng_non, &
!                             turb_mixing_l_mom_non,turb_mixing_g_mom_non,void_drift_l_mom_non,void_drift_g_mom_non, &
!                             vt_non_i,vt_non_k,                                                                     &
!                             tm_mas_l,tm_mas_g,tm_eng_l,tm_eng_g,vd_mas_l,vd_mas_g,vd_eng_l,vd_eng_g,               &
!                             turb_mixing_l_mom,turb_mixing_g_mom,void_drift_l_mom,void_drift_g_mom)
!
      DO i=1,ncell_fluid
         tm_mas_l(i)=tm_mas_l(i)*TPM(i)/cell_area(i)
         tm_mas_g(i)=tm_mas_g(i)*TPM(i)/cell_area(i)
         tm_eng_l(i)=tm_eng_l(i)*TPM(i)/cell_area(i)
         tm_eng_g(i)=tm_eng_g(i)*TPM(i)/cell_area(i)
         vd_mas_l(i)=vd_mas_l(i)*TPM(i)/cell_area(i)
         vd_mas_g(i)=vd_mas_g(i)*TPM(i)/cell_area(i)
         vd_eng_l(i)=vd_eng_l(i)*TPM(i)/cell_area(i)
         vd_eng_g(i)=vd_eng_g(i)*TPM(i)/cell_area(i)
      ENDDO

      DO ix=1,ndim
         DO i=1,ncell_fluid
            tm_mom_l=turb_mixing_l_mom(i,ix)*TPM(i)/cell_area(i)
            tm_mom_g=turb_mixing_g_mom(i,ix)*TPM(i)/cell_area(i)
            vd_mom_l=void_drift_l_mom(i,ix)*TPM(i)/cell_area(i)
            vd_mom_g=void_drift_g_mom(i,ix)*TPM(i)/cell_area(i)
            src_liq(i,ix)=src_liq(i,ix)+(tm_mom_l+vd_mom_l)
            src_gas(i,ix)=src_gas(i,ix)+(tm_mom_g+vd_mom_g)
         ENDDO
      ENDDO
!
      END SUBROUTINE udfn_mom_source_evvd_old

