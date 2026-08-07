!
      SUBROUTINE allocate_user_var
!
!     This routine allocates user variables.
!
      USE Zmpi            , ONLY: ncell_fp
      USE Zzone           , ONLY: ncell_fluid
      USE Zparam          , ONLY: ndim,nb_max
      USE Zconst1         , ONLY: cplmaster,vv_prob,wconden
      USE Zconst2         , ONLY: vl_origin,vg_origin,vd_origin                                       
      USE Zmodel          , ONLY: coef_diff,qconden,rad_source,molefr,rad_model,qrad
      USE Zndforce        , ONLY: relax_hik
      USE Zporous         , ONLY: fric_model_gas,fric_model_liq,fric_model_drp, & 
                                  tm_mas_l,tm_mas_g,tm_eng_l,tm_eng_g,vd_mas_l,vd_mas_g,vd_eng_l,vd_eng_g              
      USE Zrocom_specific , ONLY: tl_space_min,tl_space_avg,out_mdot,in_mdot,mdot_bc,                  &
                                   vtop_liq,vtop_drp,vtop_gas,etop_liq,etop_drp,etop_gas,               &
                                   ttop_liq,ttop_drp,ttop_gas,rhotop_liq,rhotop_drp,rhotop_gas,         &
                                   alphatop_liq,alphatop_drp,alphatop_gas,qualatop,closeloop,injection,rad_group
      USE Zporous         , ONLY: mixing_vane_l
      USE Zqvol           , ONLY: hil_o,hig_o
      USE Zporous         , ONLY: l_subchannel,l_mixing_vane
      USE Zrv_ncell       , ONLY: p3d_cupid,qvol_mas,ncell_fuel_rod,dnbr_cupid1,dnbr_cupid2
      USE Zrv_model       , ONLY: rv_ht_str,rv_model
      USE Zturb           ,ONLY: yplus,yplusg,utau,tauw,turb_ke,turb_keg
      USE MASTER4         , ONLY: nxy_th,nxyf,nz_th,npinx, &
                                  tfc,tfs,tcoo,dcoo,bcoo,  &
                                  p3d_th,vol_th,pin3d_th
!
      IMPLICIT NONE
!
      INTEGER n,n1
!
      n=ncell_fp
      n1=ncell_fluid

!
!.....Zrocom_specific
!
      IF(vv_prob.eq.'rocom' .or. vv_prob.eq.'rocom_mc') THEN
         ALLOCATE(tl_space_min(8),tl_space_avg(8),out_mdot(nb_max),in_mdot(nb_max),rad_group(n))
         ALLOCATE(mdot_bc(nb_max),vtop_liq(nb_max,ndim),vtop_drp(nb_max,ndim),  &
                   vtop_gas(nb_max,ndim),etop_liq(nb_max),etop_drp(nb_max),      &
                   etop_gas(nb_max),ttop_liq(nb_max),ttop_drp(nb_max),ttop_gas(nb_max),&
                   rhotop_liq(nb_max),rhotop_drp(nb_max),rhotop_gas(nb_max),&
                   alphatop_liq(nb_max),alphatop_drp(nb_max),alphatop_gas(nb_max),&
                   qualatop(nb_max))             
         closeloop=0
         injection=0
         tl_space_min(:)=0.0d0
         tl_space_avg(:)=0.0d0
         out_mdot(:)=0.0d0
         in_mdot(:)=0.0d0
         mdot_bc(:)=0.0d0
         vtop_liq(:,:)=0.0d0
         vtop_drp(:,:)=0.0d0
         vtop_gas(:,:)=0.0d0
         etop_liq(:)=0.0d0
         etop_drp(:)=0.0d0
         etop_gas(:)=0.0d0
         ttop_liq(:)=0.0d0
         ttop_drp(:)=0.0d0
         ttop_gas(:)=0.0d0
         rhotop_liq(:)=0.0d0
         rhotop_drp(:)=0.0d0
         rhotop_gas(:)=0.0d0
         alphatop_liq(:)=0.0d0
         alphatop_drp(:)=0.0d0
         alphatop_gas(:)=0.0d0
         qualatop(:)=0.0d0
         rad_group(:)=0.0d0
      ENDIF
!
!.....Zinlet_bc
!
      ALLOCATE(vl_origin(nb_max,ndim),vg_origin(nb_max,ndim),vd_origin(nb_max,ndim))                   
      vl_origin(:,:)=0.0d0
      vg_origin(:,:)=0.0d0 
      vd_origin(:,:)=0.0d0            
!
!.....Relaxation of hik
!      
      IF(relax_hik.gt.1.0d-10.or.rv_model.gt.0)THEN
         ALLOCATE(hil_o(n),hig_o(n))
         hil_o(:)=0.0d0
         hig_o(:)=0.0d0           
      ENDIF         
!
!.....Porous: Friction Coefficient
!    
      ALLOCATE(fric_model_gas(n1,ndim),fric_model_liq(n1,ndim),fric_model_drp(n1,ndim))
      fric_model_gas(:,:)=0.d0
      fric_model_liq(:,:)=0.d0
      fric_model_drp(:,:)=0.d0   
!
!.....Porous: EVVD
!
      IF(l_subchannel)then
         ALLOCATE(tm_mas_l(n),tm_mas_g(n),tm_eng_l(n),tm_eng_g(n))      
         ALLOCATE(vd_mas_l(n),vd_mas_g(n),vd_eng_l(n),vd_eng_g(n))  
         tm_mas_l(:)=0.d0
         tm_mas_g(:)=0.d0      
         tm_eng_l(:)=0.d0
         tm_eng_g(:)=0.d0
         vd_mas_l(:)=0.d0
         vd_mas_g(:)=0.d0      
         vd_eng_l(:)=0.d0
         vd_eng_g(:)=0.d0                  
      ENDIF   
      IF(l_mixing_vane)then
         ALLOCATE(mixing_vane_l(3,n))      
         mixing_vane_l=0.0d0
      ENDIF   
!
!.....Wall condensation
!
      IF(wconden.ne.0)THEN
         ALLOCATE(coef_diff(n)) 
         coef_diff(:)=0.0d0
      ENDIF
      ALLOCATE(qconden(n)) 
      qconden(:)=0.0d0   
!
!.....Radiation model
!
      IF(rad_model.ne.0)THEN
         ALLOCATE(rad_source(n1))
         rad_source(:)=0.0d0
      ENDIF
      ALLOCATE(qrad(n)) 
      qrad(:)=0.0d0 
!
!.....MASTER power output
!
      IF(rv_ht_str.eq.1)then
         ALLOCATE(p3d_cupid(ncell_fuel_rod))
         IF(cplmaster.gt.0) then
            p3d_cupid=0.0d0
         ELSE
            p3d_cupid=1.0d0
         ENDIF
      ENDIF
      IF(cplmaster.gt.0)then
         ALLOCATE(tfc(nxy_th,nz_th),tfs(nxy_th,nz_th))
         ALLOCATE(tcoo(nxy_th,nz_th),dcoo(nxy_th,nz_th),bcoo(nxy_th,nz_th))
         ALLOCATE(p3d_th(nxy_th,nz_th),vol_th(nxy_th,nz_th))
         ALLOCATE(pin3d_th(npinx,npinx,nxyf,nz_th))
         tfc=0.0d0
         tfs=0.0d0
         tcoo=0.0d0
         dcoo=0.0d0
         bcoo=0.0d0
         p3d_th=0.0d0
         vol_th=0.0d0
         pin3d_th=0.0d0
      ENDIF
      IF(cplmaster.gt.0 .or. l_subchannel)then
         ALLOCATE(qvol_mas(n))
         ALLOCATE(dnbr_cupid1(n),dnbr_cupid2(n))
         qvol_mas=0.0d0
         dnbr_cupid1=0.0d0
         dnbr_cupid2=0.0d0      
      ENDIF   
!
!.....SG
!
      IF(vv_prob.eq.'sgp_separator')THEN
         ALLOCATE(turb_ke(n),turb_keg(n),yplus(n),yplusg(n),utau(n),tauw(n))
         turb_ke(:)=0.0d0
         turb_keg(:)=0.0d0
         yplus(:)=0.0d0
         yplusg(:)=0.0d0
         utau(:)=0.0d0  
         tauw(:)=0.0d0
      ENDIF
!
!.....Molar fraction (for HYMERES)
!
      IF(vv_prob.eq.'h2p1_0'.or.vv_prob.eq.'h2p1_0x'.or. &
             vv_prob.eq.'h2p1_1'.or.vv_prob.eq.'h2p1_1x'.or. &
             vv_prob.eq.'h2p1_2'.or.vv_prob.eq.'h2p1_2x'.or. &
             vv_prob.eq.'h2p1_3'.or.vv_prob.eq.'h2p1_3x'.or. &
             vv_prob.eq.'h2p1_4'.or.vv_prob.eq.'h2p1_4x'.or. &
             vv_prob.eq.'VD_h2p1_0')THEN
         ALLOCATE(molefr(n))
         molefr(:)=0.0d0
      ENDIF
!      
      RETURN
      END SUBROUTINE allocate_user_var
