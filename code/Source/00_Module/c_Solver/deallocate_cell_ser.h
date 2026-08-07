!        DEALLOCATE(cell)
!!!!!!!
         DEALLOCATE(cell%alphag,cell%alphal,cell%alphad,cell%quala)
         DEALLOCATE(cell%rhog,cell%rhol,cell%rhod,cell%rhom)
         DEALLOCATE(cell%rhomr,cell%rhoa)
!
         DEALLOCATE(cell%eg,cell%el,cell%ed,cell%p,cell%pps)
         DEALLOCATE(cell%hg,cell%hl,cell%hgsat,cell%hlsat,cell%egsat,cell%elsat,cell%quals)
         DEALLOCATE(cell%tg,cell%tl,cell%td,cell%ts,cell%tst)
         DEALLOCATE(cell%drhogdp,cell%drhogde,cell%drhogdx)
         DEALLOCATE(cell%drholdp,cell%drholde)
         DEALLOCATE(cell%dtgdp,cell%dtgde,cell%dtgdx)
         DEALLOCATE(cell%dtldp,cell%dtlde,cell%dtsdp,cell%dtsde,cell%dtsdx)

         DEALLOCATE(cell%lviscosg,cell%lviscosl,cell%lviscosd)
         DEALLOCATE(cell%tviscosg,cell%tviscosl,cell%tviscosd)
         DEALLOCATE(cell%eviscosg,cell%eviscosl,cell%eviscosd)

         DEALLOCATE(cell%vFgl,cell%vFgd)
         DEALLOCATE(cell%entr,cell%dentr,cell%yeta)

         DEALLOCATE(cell%lcondg,cell%lcondl,cell%sigma,cell%betag,cell%betal,cell%cpg,cell%cpl)
         DEALLOCATE(cell%condg,cell%condl)

         DEALLOCATE(cell%alphag_o,cell%alphal_o,cell%alphad_o,cell%quala_o)
         DEALLOCATE(cell%eg_o,cell%el_o,cell%ed_o,cell%p_o)
         DEALLOCATE(cell%tg_o,cell%tl_o,cell%td_o,cell%ts_o)

         DEALLOCATE(cell%aint1,cell%aint2,cell%aint3,cell%D1,cell%D2,cell%Ddepart,cell%Dlift)
         DEALLOCATE(cell%limiter)
         DEALLOCATE(cell%vfwg,cell%vfwl)
         DEALLOCATE(cell%vfwg_x,cell%vfwg_y,cell%vfwg_z,cell%vfwl_x,cell%vfwl_y,cell%vfwl_z) !ibundle
         DEALLOCATE(cell%rhog_o)
!        DEALLOCATE(cell%tcondg,cell%tcondl !pik-tub-20080509-ins)

         DEALLOCATE(cell%cboron,cell%cboron_o,cell%rhol_o)      !pik-boron-2010-04-20-ins
         DEALLOCATE(cell%twall)
         DEALLOCATE(cell%fwkl,cell%fwkg)
         DEALLOCATE(cell%estm,cell%estm_o,cell%pps_o)
         DEALLOCATE(cell%film_thickness,cell%film_shear)
         DEALLOCATE(cell%regime,cell%idummyV)
         DEALLOCATE(cell%mdiff,cell%ha)
!
!........RV model
!
         DEALLOCATE(cell%wf_dry,cell%wf_VST)
         DEALLOCATE(cell%alpha_bs,cell%alpha_de,cell%alpha_sa,cell%alpha_cd,cell%alpha_gs)
         DEALLOCATE(cell%ia_bubbly)
         DEALLOCATE(cell%ia_slug_tb,cell%ia_slug_sb)         
         DEALLOCATE(cell%ia_churn)
         DEALLOCATE(cell%ia_annular_drp,cell%ia_annular_ann)         
         DEALLOCATE(cell%ia_mpr)
         DEALLOCATE(cell%ia_invann_ann,cell%ia_invann_sb)         
         DEALLOCATE(cell%ia_invchn)
         DEALLOCATE(cell%ia_invslg_drp,cell%ia_invslg_ann)         
         DEALLOCATE(cell%ia_mist)
         DEALLOCATE(cell%ia_mpo)
         DEALLOCATE(cell%ia_VST)
         DEALLOCATE(cell%ia_vst_st,cell%ia_vst_sb)

         DEALLOCATE(cell%dbb)
         DEALLOCATE(cell%dsb,cell%dtb)
         DEALLOCATE(cell%ddrp)

         DEALLOCATE(cell%dh,cell%length)
         DEALLOCATE(cell%fdir)
          
         DEALLOCATE(cell%ireflod)
         
         DEALLOCATE(cell%alphagf,cell%alphalf)