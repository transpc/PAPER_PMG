!
!!!!!!!
         ALLOCATE(cell%alphag(ncell_fp),cell%alphal(ncell_fp),cell%alphad(ncell_fp),cell%quala(ncell_fp))
         ALLOCATE(cell%rhog(ncell_fp),cell%rhol(ncell_fp),cell%rhod(ncell_fp),cell%rhom(ncell_fp))
         ALLOCATE(cell%rhomr(ncell_fp),cell%rhoa(ncell_fp))
!
         ALLOCATE(cell%eg(ncell_fp),cell%el(ncell_fp),cell%ed(ncell_fp),cell%p(ncell_fp),cell%pps(ncell_fp))
         ALLOCATE(cell%hg(ncell_fp),cell%hl(ncell_fp),cell%hgsat(ncell_fp),cell%hlsat(ncell_fp),cell%egsat(ncell_fp),cell%elsat(ncell_fp),cell%quals(ncell_fp))
         ALLOCATE(cell%tg(ncell_fp),cell%tl(ncell_fp),cell%td(ncell_fp),cell%ts(ncell_fp),cell%tst(ncell_fp))
         ALLOCATE(cell%drhogdp(ncell_fp),cell%drhogde(ncell_fp),cell%drhogdx(ncell_fp))
         ALLOCATE(cell%drholdp(ncell_fp),cell%drholde(ncell_fp))
         ALLOCATE(cell%dtgdp(ncell_fp),cell%dtgde(ncell_fp),cell%dtgdx(ncell_fp))
         ALLOCATE(cell%dtldp(ncell_fp),cell%dtlde(ncell_fp),cell%dtsdp(ncell_fp),cell%dtsde(ncell_fp),cell%dtsdx(ncell_fp))

         ALLOCATE(cell%lviscosg(ncell_fp),cell%lviscosl(ncell_fp),cell%lviscosd(ncell_fp))
         ALLOCATE(cell%tviscosg(ncell_fp),cell%tviscosl(ncell_fp),cell%tviscosd(ncell_fp))
         ALLOCATE(cell%eviscosg(ncell_fp),cell%eviscosl(ncell_fp),cell%eviscosd(ncell_fp))

         ALLOCATE(cell%vFgl(ncell_fp),cell%vFgd(ncell_fp))
         ALLOCATE(cell%entr(ncell_fp),cell%dentr(ncell_fp),cell%yeta(ncell_fp))

         ALLOCATE(cell%lcondg(ncell_fp),cell%lcondl(ncell_fp),cell%sigma(ncell_fp),cell%betag(ncell_fp),cell%betal(ncell_fp),cell%cpg(ncell_fp),cell%cpl(ncell_fp))
         ALLOCATE(cell%condg(ncell_fp),cell%condl(ncell_fp))

         ALLOCATE(cell%alphag_o(ncell_fp),cell%alphal_o(ncell_fp),cell%alphad_o(ncell_fp),cell%quala_o(ncell_fp))
         ALLOCATE(cell%eg_o(ncell_fp),cell%el_o(ncell_fp),cell%ed_o(ncell_fp),cell%p_o(ncell_fp))
         ALLOCATE(cell%tg_o(ncell_fp),cell%tl_o(ncell_fp),cell%td_o(ncell_fp),cell%ts_o(ncell_fp))

         ALLOCATE(cell%aint1(ncell_fp),cell%aint2(ncell_fp),cell%aint3(ncell_fp),cell%D1(ncell_fp),cell%D2(ncell_fp),cell%Ddepart(ncell_fp),cell%Dlift(ncell_fp))
         ALLOCATE(cell%limiter(ncell_fp))
         ALLOCATE(cell%vfwg(ncell_fp),cell%vfwl(ncell_fp))
         ALLOCATE(cell%vfwg_x(ncell_fp),cell%vfwg_y(ncell_fp),cell%vfwg_z(ncell_fp),cell%vfwl_x(ncell_fp),cell%vfwl_y(ncell_fp),cell%vfwl_z(ncell_fp)) !ibundle
         ALLOCATE(cell%rhog_o(ncell_fp))
!        ALLOCATE(cell%tcondg(ncell_fp),cell%tcondl !pik-tub-20080509-ins)

         ALLOCATE(cell%cboron(ncell_fp),cell%cboron_o(ncell_fp),cell%rhol_o(ncell_fp))      !pik-boron-2010-04-20-ins
         ALLOCATE(cell%twall(ncell_fp))
         ALLOCATE(cell%fwkl(ncell_fp),cell%fwkg(ncell_fp))
         ALLOCATE(cell%estm(ncell_fp),cell%estm_o(ncell_fp),cell%pps_o(ncell_fp))
         ALLOCATE(cell%film_thickness(ncell_fp),cell%film_shear(ncell_fp))
         ALLOCATE(cell%regime(ncell_fp),cell%idummyV(ncell_fp))
         ALLOCATE(cell%ced33(ncell_fp),cell%T_top(ncell_fp),cell%T_bot(ncell_fp))
         ALLOCATE(cell%mdiff(ncell_fp),cell%ha(ncell_fp))
!
!........RV model
!
         ALLOCATE(cell%vst(ncell_fp))
         ALLOCATE(cell%wf_dry(ncell_fp),cell%wf_VST(ncell_fp))
         ALLOCATE(cell%alpha_bs(ncell_fp),cell%alpha_de(ncell_fp),cell%alpha_sa(ncell_fp),cell%alpha_cd(ncell_fp),cell%alpha_gs(ncell_fp))
         ALLOCATE(cell%ia_bubbly(ncell_fp))
         ALLOCATE(cell%ia_slug_tb(ncell_fp),cell%ia_slug_sb(ncell_fp))         
         ALLOCATE(cell%ia_churn(ncell_fp))
         ALLOCATE(cell%ia_annular_drp(ncell_fp),cell%ia_annular_ann(ncell_fp))         
         ALLOCATE(cell%ia_mpr(ncell_fp))
         ALLOCATE(cell%ia_invann_ann(ncell_fp),cell%ia_invann_sb(ncell_fp))         
         ALLOCATE(cell%ia_invchn(ncell_fp))
         ALLOCATE(cell%ia_invslg_drp(ncell_fp),cell%ia_invslg_ann(ncell_fp))         
         ALLOCATE(cell%ia_mist(ncell_fp))
         ALLOCATE(cell%ia_mpo(ncell_fp))
         ALLOCATE(cell%ia_VST(ncell_fp))
         ALLOCATE(cell%ia_vst_st(ncell_fp),cell%ia_vst_sb(ncell_fp))

         ALLOCATE(cell%dbb(ncell_fp))
         ALLOCATE(cell%dsb(ncell_fp),cell%dtb(ncell_fp))
         ALLOCATE(cell%ddrp(ncell_fp))

         ALLOCATE(cell%length(ncell_fp))
         ALLOCATE(cell%fdir(ncell_fp))
         
         ALLOCATE(cell%c1(ncell_fp),cell%c0(ncell_fp))
         
         ALLOCATE(cell%ireflod(ncell_fp))
         
         ALLOCATE(cell%alphagf(ncell_fp),cell%alphalf(ncell_fp))