!
      SUBROUTINE scalar_update_smac3
!
!     This routine solves non-conservative energy equations
!
      USE Zinterface
      USE VOL_DATA      , ONLY: cell
      USE Zzone         , ONLY: ncell_fluid
      USE Znum_cell     , ONLY: istart_nf,istart_nb1,icell_nb,            &
                                nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zvec_param    , ONLY: nf_non,nf_mcc,nf_inl,nf_out,nf_flux
      USE Zare          , ONLY: ar_gas,ar_liq
      USE Zbc_index     , ONLY: npb,icell_type
      USE Zconst1       , ONLY: p_work_face,wconden,lwconden_alphal0
      USE Zconst2       , ONLY: dt
      USE Zcoord3       , ONLY: volpr
      USE Zcore         , ONLY: np,myrank      
      USE Zdel_scalar   , ONLY: del_eg,del_el,del_ag,del_x,limit_eng_src_opt,  &
                                limit_iht_opt,suspend_iht_opt,prn_div_eng,dsrc,err_stm, &
                                max_ihtc_opt,relax_interface_dtemp,max_ihtc_opt_coeff
      USE Zenergy_diff  , ONLY: ediff_liq,ediff_gas
      USE Zimplicit     , ONLY: imp_scalar_diff,imp_scalar_conv,iter_scalar
      USE Zmass_diff    , ONLY: mdiff_gas,mediff_gas
      USE Zncg          , ONLY: ncg_diff
      USE Zpress        , ONLY: pp
      USE Zqvol         , ONLY: h_ig,h_il,h_gf,qvol_gas,qvol_liq,gamma,gamma_wall,qporous_gas,qporous_liq
      USE Zmodel        , ONLY: rad_source,rad_model
      USE Ztimecon      , ONLY: repeat_smac
      USE Zvec_major    , ONLY: flux_l_nf,flux_g_nf,liq_conv_nf,vap_conv_nf,                &
                                ecnvc_l_nf,ecnvc_g_nf,al_conv_nf,void_conv_nf,quala_conv_nf
      USE Zrv_model     , ONLY: free_model,rv_ht_i      
      USE Zporous       , ONLY: mixing_vane_l      
      USE Zporous       , ONLY: l_subchannel,l_mixing_vane
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,iter,nb
      INTEGER :: nf_number,istart,len,i1
      INTEGER :: is
      REAL(8) :: eg,el,qu
      LOGICAL :: repeat
      REAL(8) :: ca,hmax,eg_min,el_min,deg,del
!.....Local array
      REAL(8),DIMENSION(ncell_fluid) :: dtvol
!      
      REAL(8),DIMENSION(ncell_fluid) :: arevg,arvg,avg,arxvg,     &
                                        arevl,arvl,avl,           &
                                        arvg_p,avg_p,arvl_p,avl_p
!      
      REAL(8),DIMENSION(ncell_fluid) :: eg_1,el_1,xn_1, &
                                        eg_2,el_2,xn_2
!      
      REAL(8),DIMENSION(ncell_fluid) :: parg,parl,      &
                                        dtsg,dtsl,dtgl, &
                                        xndt,gg,rgp,rlp
!
      REAL(8),DIMENSION(ncell_fluid) :: Hig_1,Hil_1,Hig_2,Hil_2,    &
                                        Hig_3,Hil_3,Hgl_1,          &
                                        hi_gas,hi_liq,del_hi,       &
                                        h0_gas,h0_liq,PsP,          &
                                        hi_gas_w,hi_liq_w
!.....Local vector arrays
      REAL(8),DIMENSION(nf_flux) :: arevg_nf,arvg_nf,avg_nf,arxvg_nf, &
                                    arevl_nf,arvl_nf,avl_nf
!.....Local allocatable arrays
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: ip
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: bm,a
      REAL(8),DIMENSION(:,:,:),ALLOCATABLE :: sm
!
      ALLOCATE(a(ncell_fluid,3),bm(ncell_fluid,3),ip(ncell_fluid,2),sm(ncell_fluid,3,3))
!
      IF(free_model)CALL int_htc
      CALL int_swap(2)
      IF(rv_ht_i.gt.0)CALL rv_int_ht
      CALL int_swap(22)       
!
      DO i=1,ncell_fluid
         H_il(i)=MAX(1.d0,H_il(i))
         H_ig(i)=MAX(1.d0,H_ig(i))
      ENDDO
!
      DO nf_number=0,3
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len  
            i1=istart+i
            arevg_nf(i1)=ecnvc_g_nf(i1)   *flux_g_nf(i1)
            arvg_nf(i1) =vap_conv_nf(i1)  *flux_g_nf(i1)
            avg_nf(i1)  =void_conv_nf(i1) *flux_g_nf(i1)
            arxvg_nf(i1)=quala_conv_nf(i1)*flux_g_nf(i1)
!
            arevl_nf(i1)=ecnvc_l_nf(i1)   *flux_l_nf(i1)
            arvl_nf(i1) =liq_conv_nf(i1)  *flux_l_nf(i1)           
            avl_nf(i1)  =al_conv_nf(i1)   *flux_l_nf(i1)           
         ENDDO
      ENDDO
!
!.....Surface option for the pressure work
!
      IF(p_work_face.gt.0) CALL press_work_face(arvg_nf,avg_nf,arvl_nf,avl_nf,arvg_p,avg_p,arvl_p,avl_p)
!
!.....Build summation info for non,mcc,inl,out
!
      nf_number_nb=3
      nf_number_id(0)=0
      nf_number_id(1)=1
      nf_number_id(2)=2
      nf_number_id(3)=3
      istart_nfs(0)=0
      istart_nfs(1)=istart_nfs(0)+nf_non
      istart_nfs(2)=istart_nfs(1)+nf_mcc
      istart_nfs(3)=istart_nfs(2)+nf_inl
      lens         =istart_nfs(3)+nf_out
!
      CALL sum_nf(0,-1,           &
                  arevg_nf,arevg, & 
                  arvg_nf ,arvg,  &
                  avg_nf  ,avg,   &
                  arxvg_nf,arxvg, &
                  arevl_nf,arevl, &
                  arvl_nf ,arvl,  &
                  avl_nf  ,avl)  
!
      DO i=1,ncell_fluid
         IF(gamma(i).ge.0.d0)THEN
            hi_gas(i)=cell%hgsat(i)
            hi_liq(i)=cell%hl(i)
         ELSE
            hi_gas(i)=cell%hg(i)
            hi_liq(i)=cell%hlsat(i)
         ENDIF
         IF(gamma_wall(i).ge.0.d0)THEN
            hi_gas_w(i)=cell%hgsat(i)
            hi_liq_w(i)=cell%hl(i)
         ELSE
            hi_gas_w(i)=cell%hg(i)
            hi_liq_w(i)=cell%hlsat(i)
         ENDIF         
!
         PsP(i)=cell%pps_o(i)/cell%p_o(i)
         h0_gas(i)=cell%eg_o(i)+cell%p_o(i)/cell%rhog(i)
         h0_liq(i)=cell%el_o(i)+cell%p_o(i)/cell%rhol(i)
         del_hi(i)=1.d0/(hi_gas(i)-hi_liq(i))
!
         parg(i)=cell%p_o(i)*cell%alphag_o(i)/cell%rhog(i)
         parl(i)=cell%p_o(i)*cell%alphal_o(i)/cell%rhol(i)
      ENDDO 
!
      IF(repeat_smac)THEN
         DO i=1,ncell_fluid
            dtsg(i)=0.0d0
            dtsl(i)=0.0d0
            dtgl(i)=0.0d0
         ENDDO
      ELSE
         DO i=1,ncell_fluid
            dtsg(i)=cell%ts_o(i)-cell%tg_o(i)+(cell%dtsdp(i)-cell%dtgdp(i))*pp(i)
            dtsl(i)=cell%ts_o(i)-cell%tl_o(i)+(cell%dtsdp(i)-cell%dtldp(i))*pp(i)
            dtgl(i)=cell%tg_o(i)-cell%tl_o(i)+(cell%dtgdp(i)-cell%dtldp(i))*pp(i)
         ENDDO
      ENDIF
!
      repeat=.false.      
!
100   CONTINUE
!    
      IF(max_ihtc_opt.eq.1) THEN !somaflow.in
         ca=max_ihtc_opt_coeff !somaflow.in
         DO i=1,ncell_fluid
!
!...........Gas interface
!
            IF(dtsg(i).lt.0.0d0)THEN
               IF(cell%alphal_o(i).le.1.0D-8)THEN
                  H_ig(i)=0.0d0
               ELSE
                  hmax=-ca*ar_liq(i)*(cell%hgsat(i)-cell%hl(i))/(dt*DMAX1(PsP(i),1.0d-5)*dtsg(i))
                  hmax=DMAX1(0.0d0,hmax)
                  H_ig(i)=DMIN1(hmax,H_ig(i))
               ENDIF
            ELSEIF(dtsg(i).gt.0.0d0)THEN
               hmax=ca*(1.0d0-cell%quala_o(i))*ar_gas(i)*(cell%hg(i)-cell%hlsat(i))/(dt*DMAX1(PsP(i),1.0d-5)*dtsg(i))
               hmax=DMAX1(0.0d0,hmax)
               H_ig(i)=DMIN1(hmax,H_ig(i))
            ENDIF
!
!..........Liquid interface
!
            IF(dtsl(i).gt.0.0d0)THEN
               IF(cell%alphag_o(i).le.1.0D-8)THEN
                  H_il(i)=0.0d0
               ELSE
                  hmax=ca*(1.0d0-cell%quala_o(i))*ar_gas(i)*(cell%hg(i)-cell%hlsat(i))/(dt*dtsl(i))
                  hmax=DMAX1(0.0d0,hmax)
                  H_il(i)=DMIN1(hmax,H_il(i))
               ENDIF
            ELSEIF(dtsl(i).lt.0.0d0)THEN
               hmax=-ca*ar_liq(i)*(cell%hgsat(i)-cell%hl(i))/(dt*dtsl(i))
               hmax=DMAX1(0.0d0,hmax)
               H_il(i)=DMIN1(hmax,H_il(i))
            ENDIF
!
         ENDDO
      ENDIF
!
      DO i=1,ncell_fluid
         Hig_1(i)=(hi_liq(i)-h0_gas(i))*del_hi(i)*PsP(i)*H_ig(i)
         Hig_2(i)=(hi_liq(i)-h0_liq(i))*del_hi(i)*PsP(i)*H_ig(i)
         Hig_3(i)=del_hi(i)*PsP(i)*H_ig(i)
         Hil_1(i)=(hi_gas(i)-h0_gas(i))*del_hi(i)*H_il(i)
         Hil_2(i)=(hi_gas(i)-h0_liq(i))*del_hi(i)*H_il(i)
         Hil_3(i)=del_hi(i)*H_il(i)
         Hgl_1(i)=(1.d0-PsP(i))*H_gf(i)
      ENDDO 
!  
      DO i=1,ncell_fluid
         eg_1(i)=cell%eg_o(i)
         el_1(i)=cell%el_o(i)
         xn_1(i)=cell%quala_o(i)
      ENDDO
!
!.....Coefficients for convection and diffusion
!      
      DO i=1,ncell_fluid
!
         dtvol(i)=dt*volpr(i)
!
!........Vapor energy
!
         a(i,1)=(Hig_1(i)*(cell%dtsde(i)-cell%dtgde(i))+Hil_1(i)*cell%dtsde(i)+Hgl_1(i)*cell%dtgde(i))*dt
!
!........Liquid energy
!
         a(i,2)=(Hil_2(i)+Hgl_1(i))*cell%dtlde(i)*dt
!
!........Non-condensable gas
!
         xndt(i)=cell%quala_o(i)*dt
         gg(i)=Hig_3(i)*dtsg(i)+Hil_3(i)*dtsl(i)-gamma_wall(i)
!
         a(i,3)=-(Hig_3(i)*(cell%dtsdx(i)-cell%dtgdx(i))+Hil_3(i)*cell%dtsdx(i))*xndt(i)-gg(i)*dt
!
         a(i,1)=a(i,1)+(ar_gas(i)-parg(i)*cell%drhogde(i))
         a(i,2)=a(i,2)+(ar_liq(i)-parl(i)*cell%drholde(i))
         a(i,3)=a(i,3)+ar_gas(i)
!
         IF(cell%alphag_o(i).lt.1.0d-8) a(i,1)=1.0d30
         IF(cell%alphal_o(i).lt.1.0d-8) a(i,2)=1.0d30
         IF(cell%alphag_o(i).lt.1.0d-8) a(i,3)=1.0d30 
         !IF(cell%quala_o(i).lt.1.0d-8) a(i,3)=1.0d30
         IF(cell%quala_o(i).lt.1.0d-8) a(i,3)=DMAX1(a(i,3),1.0d0) !ST-pik
      ENDDO
!
      IF((imp_scalar_diff*imp_scalar_conv).eq.0) iter_scalar=1
! 
      DO iter=1,iter_scalar
!
!........Explicit convection
!
         IF(imp_scalar_conv.eq.0) THEN
            DO i=1,ncell_fluid
!
!..............Vapor energy, non-condensable gas, liquid energy convection
!               
               eg=-(arevg(i)-cell%eg_o(i)   *arvg(i))*dtvol(i)
               qu=-(arxvg(i)-cell%quala_o(i)*arvg(i))*dtvol(i)
               el=-(arevl(i)-cell%el_o(i)   *arvl(i))*dtvol(i)
!               
               eg_1(i)=eg_1(i)+eg/a(i,1)
               xn_1(i)=xn_1(i)+qu/a(i,3)
               xn_1(i)=MIN(MAX(0.d0,xn_1(i)),1.d0)
               el_1(i)=el_1(i)+el/a(i,2) 
            ENDDO
         ENDIF         
!
!........Explicit diffusion
!
         IF(imp_scalar_diff.eq.0)THEN
            DO i=1,ncell_fluid
               eg_1(i)=eg_1(i)+ediff_gas(i)*dtvol(i)/a(i,1)
               el_1(i)=el_1(i)+ediff_liq(i)*dtvol(i)/a(i,2)
               IF(ncg_diff.gt.0)THEN
                  eg_1(i)=eg_1(i)+mediff_gas(i)*dtvol(i)/a(i,1)
                  xn_1(i)=xn_1(i)+mdiff_gas(i)*dtvol(i)/a(i,3)
               ENDIF
            ENDDO
         ENDIF 
!
!........Implicit convection and diffusion
!
         IF((imp_scalar_diff+imp_scalar_conv).gt.0)THEN
            CALL implicit_scalar(eg_1,el_1,xn_1,eg_2,el_2,xn_2,a,iter)
            DO i=1,ncell_fluid
               xn_1(i)=MIN(MAX(0.d0,xn_1(i)),1.d0)
            ENDDO
         ENDIF         
!
!........Phase coupling
!
         IF(iter.eq.1)THEN
!
!...........Vapor energy,energy,Non-condensible gas
!
            DO i=1,ncell_fluid
               sm(i,1,1)=  a(i,1)
               sm(i,1,2)=-(Hil_1(i)+Hgl_1(i))*cell%dtlde(i)                                                     *dt
               sm(i,1,3)= (Hig_1(i)*(cell%dtsdx(i)-cell%dtgdx(i))+Hil_1(i)*cell%dtsdx(i)+Hgl_1(i)*cell%dtgdx(i))*dt-parg(i)*cell%drhogdx(i)
!               
               sm(i,2,1)=-(Hig_2(i)*(cell%dtsde(i)-cell%dtgde(i))+Hil_2(i)*cell%dtsde(i)+Hgl_1(i)*cell%dtgde(i))*dt
               sm(i,2,2)=  a(i,2)                                          
               sm(i,2,3)=-(Hig_2(i)*(cell%dtsdx(i)-cell%dtgdx(i))+Hil_2(i)*cell%dtsdx(i)+Hgl_1(i)*cell%dtgdx(i))*dt
!               
               sm(i,3,1)=-(Hig_3(i)*(cell%dtsde(i)-cell%dtgde(i))+Hil_3(i)*cell%dtsde(i))*xndt(i)
               sm(i,3,2)=  Hil_3(i)*cell%dtlde(i)                                        *xndt(i)
               sm(i,3,3)=  a(i,3)
            ENDDO
!         
            IF(p_work_face.eq.0) THEN
               DO i=1,ncell_fluid
                  rgp(i) =cell%p_o(i)/cell%rhog(i)*(arvg(i)  -cell%rhog(i)*avg(i))  *dtvol(i)
                  rlp(i) =cell%p_o(i)/cell%rhol(i)*(arvl(i)  -cell%rhol(i)*avl(i))  *dtvol(i)
               ENDDO
            ELSE
               DO i=1,ncell_fluid
                  rgp(i)=(arvg_p(i)-cell%rhog(i)*avg_p(i))/cell%rhog(i)*dtvol(i)
                  rlp(i)=(arvl_p(i)-cell%rhol(i)*avl_p(i))/cell%rhol(i)*dtvol(i)
               ENDDO
            ENDIF
!
            IF(limit_eng_src_opt.gt.0) CALL limit_energy_source(a,dtvol,hi_gas_w,h0_gas,hi_liq_w,h0_liq)
!
            DO i=1,ncell_fluid
               bm(i,1)=rgp(i)+qvol_gas(i)*dt+qporous_gas(i)*dtvol(i)+(hi_gas_w(i)-h0_gas(i))*gamma_wall(i)*dt-             &
                        (Hig_1(i)*dtsg(i)+Hil_1(i)*dtsl(i)+Hgl_1(i)*dtgl(i))*dt+parg(i)*cell%drhogdp(i)*pp(i)
               bm(i,2)=rlp(i)+qvol_liq(i)*dt+qporous_liq(i)*dtvol(i)-(hi_liq_w(i)-h0_liq(i))*gamma_wall(i)*dt+             &
                        (Hig_2(i)*dtsg(i)+Hil_2(i)*dtsl(i)+Hgl_1(i)*dtgl(i))*dt+parl(i)*cell%drholdp(i)*pp(i)
               bm(i,3)=gg(i)*xndt(i)
            ENDDO
!            
            IF(rad_model.ne.0) THEN
               DO i=1,ncell_fluid
                  bm(i,1)=bm(i,1)+rad_source(i)*dtvol(i)
               ENDDO
            ENDIF            
!
            DO i=1,ncell_fluid
               bm(i,1)=bm(i,1)+a(i,1)*(eg_1(i)-cell%eg_o(i))
               bm(i,2)=bm(i,2)+a(i,2)*(el_1(i)-cell%el_o(i))
               bm(i,3)=bm(i,3)+a(i,3)*(xn_1(i)-cell%quala_o(i))
            ENDDO   
         ELSE
            DO i=1,ncell_fluid
               bm(i,1)=a(i,1)*(eg_1(i)-eg_2(i))
               bm(i,2)=a(i,2)*(el_1(i)-el_2(i))
               bm(i,3)=a(i,3)*(xn_1(i)-xn_2(i))
            ENDDO   
         ENDIF
!
!.......EVVD
!         
         IF(l_subchannel)then
            !DO i=1,ncell_fluid     ! skip EVVD energy update (PSH edit)
            !   IF(npb(i).ne.0) CYCLE
            !   bm(i,1)= bm(i,1)+(tm_eng_g(i)+vd_eng_g(i))*dt            &
            !                   -(tm_mas_g(i)+vd_mas_g(i))*h0_gas(i)*dt
            !   bm(i,2)= bm(i,2)+(tm_eng_l(i)+vd_eng_l(i))*dt            &
            !                   -(tm_mas_l(i)+vd_mas_l(i))*h0_liq(i)*dt
            !ENDDO
            ! mixing_vane
            IF(l_mixing_vane)then
               DO i=1,ncell_fluid
                  IF(npb(i).ne.0) CYCLE
                  bm(i,2)= bm(i,2)+mixing_vane_l(3,i)*dtvol(i)           &
                                  -mixing_vane_l(2,i)*dtvol(i)*h0_liq(i)
               ENDDO
            ENDIF
         ENDIF                 
!
         IF(iter.eq.1) then
            CALL luinverse31(sm,ip,bm,ncell_fluid,npb)
         ELSE
            CALL solve31(sm,ip,bm,ncell_fluid,npb)
         ENDIF
!
         IF(iter.eq.1)THEN
            DO i=1,ncell_fluid
               eg_2(i)=cell%eg_o(i)+bm(i,1)
               el_2(i)=cell%el_o(i)+bm(i,2)
               xn_2(i)=cell%quala_o(i)+bm(i,3)
            ENDDO   
         ELSE
            DO i=1,ncell_fluid
               eg_2(i)=eg_2(i)+bm(i,1)
               el_2(i)=el_2(i)+bm(i,2)
               xn_2(i)=xn_2(i)+bm(i,3)
            ENDDO   
         ENDIF
         DO i=1,ncell_fluid
            xn_2(i)=MIN(MAX(0.d0,xn_2(i)),1.d0)
         ENDDO
!
         IF(relax_interface_dtemp.eq.1) THEN !somaflow.in
            IF(.not.repeat.and.iter.eq.2)THEN 
               DO i=1,ncell_fluid
                  eg_min=1.0d6 
                  el_min=1.0d4 
                  IF(eg_2(i).le.eg_min)THEN
                     deg=DMAX1((eg_2(i)),0.0d1)
                     deg=(deg-0.0d1)/eg_min
                     dtsg(i)=dtsg(i)*deg
                     dtsl(i)=dtsl(i)*deg
                     dtgl(i)=dtgl(i)*deg
                     repeat=.true.
                     CYCLE
                  ENDIF
                  IF(el_2(i).le.el_min)THEN 
                     del=DMAX1((el_2(i)),0.0d1) 
                     del=(del-0.0d1)/el_min
                     dtsg(i)=dtsg(i)*del
                     dtsl(i)=dtsl(i)*del
                     dtgl(i)=dtgl(i)*del
                     repeat=.true. 
                     CYCLE
                  ENDIF
               ENDDO
! 
               IF(np.gt.1) CALL allreducei_l1(repeat)
               IF(repeat) GOTO 100
!
            ENDIF
         ENDIF
!
!
      ENDDO !iter
!
      DO i=1,ncell_fluid
         IF(npb(i).eq.0)THEN
            cell%eg(i)=eg_2(i)
            cell%el(i)=el_2(i)
            cell%quala(i)=xn_2(i)
         ENDIF
      ENDDO
!
      IF(limit_iht_opt.gt.0)THEN
         DO i=1,ncell_fluid
            IF(cell%quala(i).le.1.0d-8) cell%quala(i)=0.0d0
         ENDDO
      ENDIF
!
      IF(wconden.ne.0.and.lwconden_alphal0)THEN
         DO nf_number=4,7
            istart=istart_nb1(1,nf_number)
            len   =istart_nb1(2,nf_number)
            DO nb=1,len
               i1=istart+nb
               i=icell_nb(i1)
               IF(icell_type(i).eq.1)THEN
                  cell%alphag(i)=1.d0
                  cell%alphal(i)=0.d0
                  cell%alphad(i)=0.d0
               ENDIF
            ENDDO
         ENDDO
      ENDIF      
!
!.....Save scalar changes
!
      DO i=1,ncell_fluid
         del_eg(i)=cell%eg(i)-cell%eg_o(i)
         del_el(i)=cell%el(i)-cell%el_o(i)
         del_ag(i)=cell%alphag(i)-cell%alphag_o(i)
         del_x(i)=cell%quala(i)-cell%quala_o(i)
      ENDDO
!       
      IF(prn_div_eng.gt.0.and.myrank.eq.0)THEN
         DO i=1,ncell_fluid
            IF(ABS(del_eg(i)/cell%eg_o(i)).gt.dsrc)THEN
               print *, 'gas', i, ABS(del_eg(i)/cell%eg_o(i)), cell%alphag_o(i)
            ENDIF
            IF(ABS(del_el(i)/cell%el_o(i)).gt.dsrc)THEN
               print *, 'liq', i, ABS(del_el(i)/cell%el_o(i)), cell%alphag_o(i)
            ENDIF
         ENDDO
      ENDIF
!
      CALL ncg_transport
!
!.....Update material properties
!
      CALL property_calc(1)
!
      repeat_smac=.false.
      IF(suspend_iht_opt.gt.0)THEN
         is=0
         DO i=1,ncell_fluid
            IF(err_stm(i).gt.0) is=1
        ENDDO
        IF(np.gt.1) CALL allreducei_max_i1(is)
        IF(is.gt.0)THEN
           repeat_smac=.true.
           RETURN
        ENDIF
      ENDIF
!
!.....Vapor generation rate
!
      CALL set_vapor_generation
!
      DEALLOCATE(a,sm,bm,ip)
!
      END SUBROUTINE scalar_update_smac3
!
!-------------------------------------------------------------------------
!
      SUBROUTINE limit_energy_source(a,dtvol,hi_gas_w,h0_gas,hi_liq_w,h0_liq)
!
      USE VOL_DATA      , ONLY: cell
      USE Zconst2       , ONLY: dt
      USE Zdel_scalar   , ONLY: dsrc
      USE Zqvol         , ONLY: qvol_gas,qvol_liq,gamma_wall,qporous_gas,qporous_liq
      USE Zzone         , ONLY: ncell_fluid
!
      IMPLICIT NONE
      INTEGER i
      REAL(8) eg,el,delta,dth
      REAL(8) a(ncell_fluid,3),dtvol(ncell_fluid),hi_gas_w(ncell_fluid),h0_gas(ncell_fluid), &
              hi_liq_w(ncell_fluid),h0_liq(ncell_fluid)
!
      DO i=1,ncell_fluid
!
         eg=dsrc*cell%eg_o(i)
         el=dsrc*cell%el_o(i)
!
         delta=qvol_gas(i)*dt/a(i,1)
         IF(delta.gt.eg) qvol_gas(i)=eg*a(i,1)/dt
         IF(delta.lt.-eg) qvol_gas(i)=-eg*a(i,1)/dt
         delta=qvol_liq(i)*dt/a(i,2)
         IF(delta.gt.el) qvol_liq(i)=el*a(i,2)/dt
         IF(delta.lt.-el) qvol_liq(i)=-el*a(i,2)/dt
!
         delta=qporous_gas(i)*dtvol(i)/a(i,1)
         IF(delta.gt.eg) qporous_gas(i)=eg*a(i,1)/dtvol(i)
         IF(delta.lt.-eg) qporous_gas(i)=-eg*a(i,1)/dtvol(i)
         delta=qporous_liq(i)*dtvol(i)/a(i,2)
         IF(delta.gt.el) qporous_liq(i)=el*a(i,2)/dtvol(i)
         IF(delta.lt.-el) qporous_liq(i)=-el*a(i,2)/dtvol(i)
!
         dth=(hi_gas_w(i)-h0_gas(i))*dt
         delta=dth*gamma_wall(i)/a(i,1)
         IF(delta.gt.eg) gamma_wall(i)=eg*a(i,1)/dth
         IF(delta.lt.-eg) gamma_wall(i)=-eg*a(i,1)/dth
         dth=(hi_liq_w(i)-h0_liq(i))*dt
         delta=dth*gamma_wall(i)/a(i,2)
         IF(delta.gt.el) gamma_wall(i)=el*a(i,2)/dth
         IF(delta.lt.-el) gamma_wall(i)=-el*a(i,2)/dth
!
      ENDDO    
!
      RETURN
      END SUBROUTINE limit_energy_source
