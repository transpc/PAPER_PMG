!        ALLOCATE(cell(ncell_fluid))
!!!!!!!
         ALLOCATE(cell%alphag(ncell_fluid),cell%alphal(ncell_fluid),cell%alphad(ncell_fluid),cell%quala(ncell_fluid))
         ALLOCATE(cell%rhog(ncell_fluid),cell%rhol(ncell_fluid),cell%rhod(ncell_fluid),cell%rhom(ncell_fluid))
         ALLOCATE(cell%rhomr(ncell_fluid),cell%rhoa(ncell_fluid))
!
         ALLOCATE(cell%eg(ncell_fluid),cell%el(ncell_fluid),cell%ed(ncell_fluid),cell%p(ncell_fluid),cell%pps(ncell_fluid))
         ALLOCATE(cell%hg(ncell_fluid),cell%hl(ncell_fluid),cell%hgsat(ncell_fluid),cell%hlsat(ncell_fluid),cell%egsat(ncell_fluid),cell%elsat(ncell_fluid),cell%quals(ncell_fluid))
         ALLOCATE(cell%tg(ncell_fluid),cell%tl(ncell_fluid),cell%td(ncell_fluid),cell%ts(ncell_fluid),cell%tst(ncell_fluid))
         ALLOCATE(cell%drhogdp(ncell_fluid),cell%drhogde(ncell_fluid),cell%drhogdx(ncell_fluid))
         ALLOCATE(cell%drholdp(ncell_fluid),cell%drholde(ncell_fluid))
         ALLOCATE(cell%dtgdp(ncell_fluid),cell%dtgde(ncell_fluid),cell%dtgdx(ncell_fluid))
         ALLOCATE(cell%dtldp(ncell_fluid),cell%dtlde(ncell_fluid),cell%dtsdp(ncell_fluid),cell%dtsde(ncell_fluid),cell%dtsdx(ncell_fluid))

         ALLOCATE(cell%lviscosg(ncell_fluid),cell%lviscosl(ncell_fluid),cell%lviscosd(ncell_fluid))
         ALLOCATE(cell%tviscosg(ncell_fluid),cell%tviscosl(ncell_fluid),cell%tviscosd(ncell_fluid))
         ALLOCATE(cell%eviscosg(ncell_fluid),cell%eviscosl(ncell_fluid),cell%eviscosd(ncell_fluid))

         ALLOCATE(cell%vFgl(ncell_fluid),cell%vFgd(ncell_fluid))
         ALLOCATE(cell%entr(ncell_fluid),cell%dentr(ncell_fluid),cell%yeta(ncell_fluid))

         ALLOCATE(cell%lcondg(ncell_fluid),cell%lcondl(ncell_fluid),cell%sigma(ncell_fluid),cell%betag(ncell_fluid),cell%betal(ncell_fluid),cell%cpg(ncell_fluid),cell%cpl(ncell_fluid))
         ALLOCATE(cell%condg(ncell_fluid),cell%condl(ncell_fluid))

         ALLOCATE(cell%alphag_o(ncell_fluid),cell%alphal_o(ncell_fluid),cell%alphad_o(ncell_fluid),cell%quala_o(ncell_fluid))
         ALLOCATE(cell%eg_o(ncell_fluid),cell%el_o(ncell_fluid),cell%ed_o(ncell_fluid),cell%p_o(ncell_fluid))
         ALLOCATE(cell%tg_o(ncell_fluid),cell%tl_o(ncell_fluid),cell%td_o(ncell_fluid),cell%ts_o(ncell_fluid))

         ALLOCATE(cell%aint1(ncell_fluid),cell%aint2(ncell_fluid),cell%aint3(ncell_fluid),cell%D1(ncell_fluid),cell%D2(ncell_fluid),cell%Ddepart(ncell_fluid),cell%Dlift(ncell_fluid))
         ALLOCATE(cell%limiter(ncell_fluid))
         ALLOCATE(cell%vfwg(ncell_fluid),cell%vfwl(ncell_fluid))
         ALLOCATE(cell%vfwg_x(ncell_fluid),cell%vfwg_y(ncell_fluid),cell%vfwg_z(ncell_fluid),cell%vfwl_x(ncell_fluid),cell%vfwl_y(ncell_fluid),cell%vfwl_z(ncell_fluid)) !ibundle
         ALLOCATE(cell%rhog_o(ncell_fluid))
!        ALLOCATE(cell%tcondg(ncell_fluid),cell%tcondl !pik-tub-20080509-ins)

         ALLOCATE(cell%cboron(ncell_fluid),cell%cboron_o(ncell_fluid),cell%rhol_o(ncell_fluid))      !pik-boron-2010-04-20-ins
         ALLOCATE(cell%twall(ncell_fluid))
         ALLOCATE(cell%fwkl(ncell_fluid),cell%fwkg(ncell_fluid))
         ALLOCATE(cell%estm(ncell_fluid),cell%estm_o(ncell_fluid),cell%pps_o(ncell_fluid))
         ALLOCATE(cell%film_thickness(ncell_fluid),cell%film_shear(ncell_fluid))
         ALLOCATE(cell%regime(ncell_fluid),cell%idummyV(ncell_fluid))
         ALLOCATE(cell%ced33(ncell_fluid),cell%T_top(ncell_fluid),cell%T_bot(ncell_fluid))
         ALLOCATE(cell%mdiff(ncell_fluid),cell%ha(ncell_fluid))
!
!........RV model
!
         ALLOCATE(cell%vst(ncell_fluid))
         ALLOCATE(cell%wf_dry(ncell_fluid),cell%wf_VST(ncell_fluid))
         ALLOCATE(cell%alpha_bs(ncell_fluid),cell%alpha_de(ncell_fluid),cell%alpha_sa(ncell_fluid),cell%alpha_cd(ncell_fluid),cell%alpha_gs(ncell_fluid))
         ALLOCATE(cell%ia_bubbly(ncell_fluid))
         ALLOCATE(cell%ia_slug_tb(ncell_fluid),cell%ia_slug_sb(ncell_fluid))         
         ALLOCATE(cell%ia_churn(ncell_fluid))
         ALLOCATE(cell%ia_annular_drp(ncell_fluid),cell%ia_annular_ann(ncell_fluid))         
         ALLOCATE(cell%ia_mpr(ncell_fluid))
         ALLOCATE(cell%ia_invann_ann(ncell_fluid),cell%ia_invann_sb(ncell_fluid))         
         ALLOCATE(cell%ia_invchn(ncell_fluid))
         ALLOCATE(cell%ia_invslg_drp(ncell_fluid),cell%ia_invslg_ann(ncell_fluid))         
         ALLOCATE(cell%ia_mist(ncell_fluid))
         ALLOCATE(cell%ia_mpo(ncell_fluid))
         ALLOCATE(cell%ia_VST(ncell_fluid))
         ALLOCATE(cell%ia_vst_st(ncell_fluid),cell%ia_vst_sb(ncell_fluid))

         ALLOCATE(cell%dbb(ncell_fluid))
         ALLOCATE(cell%dsb(ncell_fluid),cell%dtb(ncell_fluid))
         ALLOCATE(cell%ddrp(ncell_fluid))

         ALLOCATE(cell%length(ncell_fluid))
         ALLOCATE(cell%fdir(ncell_fluid))
         
         ALLOCATE(cell%ireflod(ncell_fluid))
         
         ALLOCATE(cell%alphagf(ncell_fluid),cell%alphalf(ncell_fluid))