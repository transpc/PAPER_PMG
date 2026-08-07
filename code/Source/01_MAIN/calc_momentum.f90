!
      SUBROUTINE calc_momentum
!
!     This routine discretizes and solves the linearized momentum equations
!     explicitly (imp_mom_diff+imp_mom_conv=0 or iter_mom=1) or implicitly.
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: np
      USE Zparam       , ONLY: ndim
      USE Zvec_param   , ONLY: nf_non
      USE Zare         , ONLY: ar_gas,ar_liq,ar_drp
      USE Zbicg        , ONLY: relax_u
      USE Zconst1      , ONLY: mdiffoff,nwlf,ntdf,fric_face
      USE Zconst2      , ONLY: grav,gfactor,dt
      USE Zcoord3      , ONLY: volpr
      USE Zgradoption  , ONLY: grav_grad
      USE Zimplicit    , ONLY: imp_mom_conv,imp_mom_diff,ag_min_m,al_min_m,iter_mom
      USE Zm_src       , ONLY: src_gas,src_liq,src_drp
      USE Zmodel       , ONLY: drift_c0,drift_c1,cb_bubble
      USE Zporous      , ONLY: fric_model_gas,fric_model_liq,fric_model_drp
      USE Zpress       , ONLY: p,dpdx,dpdx_o
      USE Zpress_coeff , ONLY: coefp_g,coefp_l,coefp_d,coefm_g,coefm_l
      USE Zqvol        , ONLY: gamma,gamma_wall
      USE Zuserdefined , ONLY: udfl_mom_press_source,udfl_mom_loss,udfl_mom_film_shear
      USE Zporous      , ONLY: l_subchannel      
      USE Zvector      , ONLY: vl_n,vg_n,vd_n,vl_o,vg_o,vd_o,vg_t,vl_t,vd_t,vl_f_non,vg_f_non
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,ix,iter,nrhs
      REAL(8) :: bm1,bm2,bm3
      REAL(8) :: gamma_EV,gamma_CD,gamma_EV_w,gamma_CD_w
      REAL(8) :: yeta1
      REAL(8) :: a_g,a_l,a_d
      REAL(8) :: avm_gl_vg,avm_gl_vl,dvlg
      REAL(8) :: denom_g1,denom_g2,denom_g3         
      REAL(8) :: denom_l1,denom_l2,denom_l3         
      REAL(8) :: denom_d1,denom_d2,denom_d3         
!.....Local arrays
      INTEGER :: ip(ncell_fluid,2)
      REAL(8) :: am(ncell_fluid,3,3),bm(ncell_fluid,3,ndim+2)
      REAL(8) :: ag(ncell_fluid),al(ncell_fluid),ad(ncell_fluid)
      REAL(8) :: avm_gl(ncell_fluid),avm_gd(ncell_fluid)
!
!     Initialize variables
!
      IF(ndim.eq.2) THEN
         DO i=1,ncell_fluid
            fric_model_gas(i,1)=0.d0
            fric_model_gas(i,2)=0.d0
            fric_model_liq(i,1)=0.d0
            fric_model_liq(i,2)=0.d0
            fric_model_drp(i,1)=0.d0
            fric_model_drp(i,2)=0.d0
         ENDDO
      ELSE
         DO i=1,ncell_fluid
            fric_model_gas(i,1)=0.d0
            fric_model_gas(i,2)=0.d0
            fric_model_gas(i,3)=0.d0
            fric_model_liq(i,1)=0.d0
            fric_model_liq(i,2)=0.d0
            fric_model_liq(i,3)=0.d0
            fric_model_drp(i,1)=0.d0
            fric_model_drp(i,2)=0.d0
            fric_model_drp(i,3)=0.d0
         ENDDO
      ENDIF
!
!.....Set boundary conditions
!
      CALL set_inlet_flux
      CALL set_outlet_property
      CALL check_mass
!
      IF(ndim.eq.2) THEN
         DO i=1,ncell_fp
            vl_t(i,1)=vl_o(i,1)
            vl_t(i,2)=vl_o(i,2)
            vg_t(i,1)=vg_o(i,1)
            vg_t(i,2)=vg_o(i,2)
            vd_t(i,1)=vd_o(i,1)
            vd_t(i,2)=vd_o(i,2)
         ENDDO
      ELSE
         DO i=1,ncell_fp
            vl_t(i,1)=vl_o(i,1)
            vl_t(i,2)=vl_o(i,2)
            vl_t(i,3)=vl_o(i,3)
            vg_t(i,1)=vg_o(i,1)
            vg_t(i,2)=vg_o(i,2)
            vg_t(i,3)=vg_o(i,3)
            vd_t(i,1)=vd_o(i,1)
            vd_t(i,2)=vd_o(i,2)
            vd_t(i,3)=vd_o(i,3)
         ENDDO
      ENDIF
!
!.....Mixture, virtual mass and slip coefficient
!
      CALL virtual_mass(cell%alphag_o,cell%alphal_o,cell%alphad_o,cell%rhog,cell%rhol,cell%rhod, &
                        ar_gas,ar_liq,ar_drp,ag,al,ad,avm_gl,avm_gd,ncell_fluid)
!
!.....Calculate gradient vector components at cell center
!
      IF(grav_grad.eq.0)THEN
         CALL grad_press(p,dpdx,0)
         IF(np.gt.1) CALL communicate_2d(dpdx)
      ELSEIF(grav_grad.eq.1)THEN
         CALL grad_pressK1(p,dpdx,0)
         IF(np.gt.1) CALL communicate_2d(dpdx)
      ELSEIF(grav_grad.ge.2)THEN
         CALL grad_frink(p,dpdx,grav_grad)
         IF(np.gt.1) CALL communicate_2d(dpdx)
      ENDIF
!
      IF(np.gt.1) CALL communicate_1d(ar_liq,        &
                                      ar_gas,        &
                                      ar_drp,        &
                                      cell%eviscosl, &
                                      cell%eviscosg, &
                                      cell%eviscosd)
!
!.....Explicit diffusion
!
      IF(mdiffoff.eq.0.and.imp_mom_diff.eq.0) THEN
         CALL momentum_diffusion(src_liq,src_gas,src_drp)
      ELSE
         IF(ndim.eq.2) THEN
            DO i=1,ncell_fluid
               src_liq(i,1)=0.d0
               src_liq(i,2)=0.d0
               src_gas(i,1)=0.d0
               src_gas(i,2)=0.d0
               src_drp(i,1)=0.d0
               src_drp(i,2)=0.d0
            ENDDO
         ELSE
            DO i=1,ncell_fluid
               src_liq(i,1)=0.d0
               src_liq(i,2)=0.d0
               src_liq(i,3)=0.d0
               src_gas(i,1)=0.d0
               src_gas(i,2)=0.d0
               src_gas(i,3)=0.d0
               src_drp(i,1)=0.d0
               src_drp(i,2)=0.d0
               src_drp(i,3)=0.d0
            ENDDO
         ENDIF
      ENDIF
!
!.....Explicit convection
!
      IF(imp_mom_conv.eq.0) CALL momentum_convection(src_liq,src_gas,src_drp)
!
!.....Gravity source contribution
!
      IF(imp_mom_conv.eq.0 .or. (mdiffoff.eq.0.and.imp_mom_diff.eq.0)) THEN
         IF(ndim.eq.2) THEN
            DO i=1,ncell_fluid
               a_g=ar_gas(i)*gfactor(i) 
               a_l=ar_liq(i)*gfactor(i)
               a_d=ar_drp(i)*gfactor(i) 
               src_gas(i,1)=src_gas(i,1)*volpr(i)+a_g*grav(1)
               src_gas(i,2)=src_gas(i,2)*volpr(i)+a_g*grav(2)
               src_liq(i,1)=src_liq(i,1)*volpr(i)+a_l*grav(1)
               src_liq(i,2)=src_liq(i,2)*volpr(i)+a_l*grav(2)
               src_drp(i,1)=src_drp(i,1)*volpr(i)+a_d*grav(1)
               src_drp(i,2)=src_drp(i,2)*volpr(i)+a_d*grav(2)
            ENDDO
         ELSE
            DO i=1,ncell_fluid
               a_g=ar_gas(i)*gfactor(i) 
               a_l=ar_liq(i)*gfactor(i)
               a_d=ar_drp(i)*gfactor(i) 
               src_gas(i,1)=src_gas(i,1)*volpr(i)+a_g*grav(1)
               src_gas(i,2)=src_gas(i,2)*volpr(i)+a_g*grav(2)
               src_gas(i,3)=src_gas(i,3)*volpr(i)+a_g*grav(3)
               src_liq(i,1)=src_liq(i,1)*volpr(i)+a_l*grav(1)
               src_liq(i,2)=src_liq(i,2)*volpr(i)+a_l*grav(2)
               src_liq(i,3)=src_liq(i,3)*volpr(i)+a_l*grav(3)
               src_drp(i,1)=src_drp(i,1)*volpr(i)+a_d*grav(1)
               src_drp(i,2)=src_drp(i,2)*volpr(i)+a_d*grav(2)
               src_drp(i,3)=src_drp(i,3)*volpr(i)+a_d*grav(3)
            ENDDO
         ENDIF
      ELSE
         IF(ndim.eq.2) THEN
            DO i=1,ncell_fluid
               a_g=ar_gas(i)*gfactor(i) 
               a_l=ar_liq(i)*gfactor(i)
               a_d=ar_drp(i)*gfactor(i) 
               src_gas(i,1)=a_g*grav(1)
               src_gas(i,2)=a_g*grav(2)
               src_liq(i,1)=a_l*grav(1)
               src_liq(i,2)=a_l*grav(2)
               src_drp(i,1)=a_d*grav(1)
               src_drp(i,2)=a_d*grav(2)
            ENDDO
         ELSE
            DO i=1,ncell_fluid
               a_g=ar_gas(i)*gfactor(i) 
               a_l=ar_liq(i)*gfactor(i)
               a_d=ar_drp(i)*gfactor(i) 
               src_gas(i,1)=a_g*grav(1)
               src_gas(i,2)=a_g*grav(2)
               src_gas(i,3)=a_g*grav(3)
               src_liq(i,1)=a_l*grav(1)
               src_liq(i,2)=a_l*grav(2)
               src_liq(i,3)=a_l*grav(3)
               src_drp(i,1)=a_d*grav(1)
               src_drp(i,2)=a_d*grav(2)
               src_drp(i,3)=a_d*grav(3)
            ENDDO
         ENDIF
      ENDIF
!
!.....Explicit wall shear force
!
      IF(udfl_mom_film_shear)THEN
         CALL udfn_mom_film_shear
      ENDIF 
!
!.....Pressure source contribution
!
      IF(udfl_mom_press_source)THEN
         CALL udfn_mom_press_source
      ELSE      
         IF(ndim.eq.2) THEN
            DO i=1,ncell_fluid
               src_gas(i,1)=src_gas(i,1)-cell%alphag(i)*dpdx(i,1)
               src_gas(i,2)=src_gas(i,2)-cell%alphag(i)*dpdx(i,2)
               src_liq(i,1)=src_liq(i,1)-cell%alphal(i)*dpdx(i,1)
               src_liq(i,2)=src_liq(i,2)-cell%alphal(i)*dpdx(i,2)
               src_drp(i,1)=src_drp(i,1)-cell%alphad(i)*dpdx(i,1)
               src_drp(i,2)=src_drp(i,2)-cell%alphad(i)*dpdx(i,2)
            ENDDO
         ELSE
            DO i=1,ncell_fluid
               src_gas(i,1)=src_gas(i,1)-cell%alphag(i)*dpdx(i,1)
               src_gas(i,2)=src_gas(i,2)-cell%alphag(i)*dpdx(i,2)
               src_gas(i,3)=src_gas(i,3)-cell%alphag(i)*dpdx(i,3)
               src_liq(i,1)=src_liq(i,1)-cell%alphal(i)*dpdx(i,1)
               src_liq(i,2)=src_liq(i,2)-cell%alphal(i)*dpdx(i,2)
               src_liq(i,3)=src_liq(i,3)-cell%alphal(i)*dpdx(i,3)
               src_drp(i,1)=src_drp(i,1)-cell%alphad(i)*dpdx(i,1)
               src_drp(i,2)=src_drp(i,2)-cell%alphad(i)*dpdx(i,2)
               src_drp(i,3)=src_drp(i,3)-cell%alphad(i)*dpdx(i,3)
            ENDDO
         ENDIF
      ENDIF
!
!.....Force loss by form loss and turbulence mixing
!
      IF(udfl_mom_loss) CALL udfn_mom_source
!
!.....Anisotropic Friction Model - subchannel model (explicit)
!
      IF(l_subchannel) CALL udfn_subchannel_mom_source      
!
!.....Non-drag forces
!
      CALL int_non_drag_coeff
      IF(nwlf .ne.-1.d0) CALL int_non_drag_lub
      IF(ntdf .ne.-1.d0) CALL int_non_drag_turb(Cb_bubble)
!
!.....vg_t,vl_t,vd_t
!    
      IF(ndim.eq.2) THEN
         DO i=1,ncell_fluid
            IF(cell%alphag_o(i).gt.ag_min_m) THEN
               bm1=src_gas(i,1)
               bm2=src_gas(i,2)
               denom_g1=ar_gas(i)-fric_model_gas(i,1)*dt
               denom_g2=ar_gas(i)-fric_model_gas(i,2)*dt
               IF(ABS(denom_g1).lt.1.d-8) denom_g1=SIGN(1.d0,denom_g1)*MAX(1.d-8,ABS(denom_g1))
               IF(ABS(denom_g2).lt.1.d-8) denom_g2=SIGN(1.d0,denom_g2)*MAX(1.d-8,ABS(denom_g2))
               vg_t(i,1)=(ar_gas(i)*vg_o(i,1)+bm1*dt)/denom_g1
               vg_t(i,2)=(ar_gas(i)*vg_o(i,2)+bm2*dt)/denom_g2
            ENDIF
            IF(cell%alphal_o(i).gt.al_min_m) THEN
               bm1=src_liq(i,1)
               bm2=src_liq(i,2)
               denom_l1=ar_liq(i)-fric_model_liq(i,1)*dt
               denom_l2=ar_liq(i)-fric_model_liq(i,2)*dt
               IF(ABS(denom_l1).lt.1.d-8) denom_l1=SIGN(1.d0,denom_l1)*MAX(1.d-8,ABS(denom_l1))
               IF(ABS(denom_l2).lt.1.d-8) denom_l2=SIGN(1.d0,denom_l2)*MAX(1.d-8,ABS(denom_l2))
               vl_t(i,1)=(ar_liq(i)*vl_o(i,1)+bm1*dt)/denom_l1
               vl_t(i,2)=(ar_liq(i)*vl_o(i,2)+bm2*dt)/denom_l2
            ENDIF
            IF(cell%alphad_o(i).gt.al_min_m) THEN
               bm1=src_drp(i,1)
               bm2=src_drp(i,2)
               denom_d1=ar_drp(i)-fric_model_drp(i,1)*dt
               denom_d2=ar_drp(i)-fric_model_drp(i,2)*dt
               IF(ABS(denom_d1).lt.1.d-8) denom_d1=SIGN(1.d0,denom_d1)*MAX(1.d-8,ABS(denom_d1))
               IF(ABS(denom_d2).lt.1.d-8) denom_d2=SIGN(1.d0,denom_d2)*MAX(1.d-8,ABS(denom_d2))
               vd_t(i,1)=(ar_drp(i)*vd_o(i,1)+bm1*dt)/denom_d1
               vd_t(i,2)=(ar_drp(i)*vd_o(i,2)+bm2*dt)/denom_d2
            ENDIF
!
         ENDDO
      ELSE
         DO i=1,ncell_fluid
            IF(cell%alphag_o(i).gt.ag_min_m) THEN
               bm1=src_gas(i,1)
               bm2=src_gas(i,2)
               bm3=src_gas(i,3)
               denom_g1=ar_gas(i)-fric_model_gas(i,1)*dt
               denom_g2=ar_gas(i)-fric_model_gas(i,2)*dt
               denom_g3=ar_gas(i)-fric_model_gas(i,3)*dt
               IF(ABS(denom_g1).lt.1.d-8) denom_g1=SIGN(1.d0,denom_g1)*MAX(1.d-8,ABS(denom_g1))
               IF(ABS(denom_g2).lt.1.d-8) denom_g2=SIGN(1.d0,denom_g2)*MAX(1.d-8,ABS(denom_g2))
               IF(ABS(denom_g3).lt.1.d-8) denom_g3=SIGN(1.d0,denom_g3)*MAX(1.d-8,ABS(denom_g3))
               vg_t(i,1)=(ar_gas(i)*vg_o(i,1)+bm1*dt)/denom_g1
               vg_t(i,2)=(ar_gas(i)*vg_o(i,2)+bm2*dt)/denom_g2
               vg_t(i,3)=(ar_gas(i)*vg_o(i,3)+bm3*dt)/denom_g3
            ENDIF
            IF(cell%alphal_o(i).gt.al_min_m) THEN
               bm1=src_liq(i,1)
               bm2=src_liq(i,2)
               bm3=src_liq(i,3)
               denom_l1=ar_liq(i)-fric_model_liq(i,1)*dt
               denom_l2=ar_liq(i)-fric_model_liq(i,2)*dt
               denom_l3=ar_liq(i)-fric_model_liq(i,3)*dt
               IF(ABS(denom_l1).lt.1.d-8) denom_l1=SIGN(1.d0,denom_l1)*MAX(1.d-8,ABS(denom_l1))
               IF(ABS(denom_l2).lt.1.d-8) denom_l2=SIGN(1.d0,denom_l2)*MAX(1.d-8,ABS(denom_l2))
               IF(ABS(denom_l3).lt.1.d-8) denom_l3=SIGN(1.d0,denom_l3)*MAX(1.d-8,ABS(denom_l3))
               vl_t(i,1)=(ar_liq(i)*vl_o(i,1)+bm1*dt)/denom_l1
               vl_t(i,2)=(ar_liq(i)*vl_o(i,2)+bm2*dt)/denom_l2
               vl_t(i,3)=(ar_liq(i)*vl_o(i,3)+bm3*dt)/denom_l3
            ENDIF
            IF(cell%alphad_o(i).gt.al_min_m) THEN
               bm1=src_drp(i,1)
               bm2=src_drp(i,2)
               bm3=src_drp(i,3)
               denom_d1=ar_drp(i)-fric_model_drp(i,1)*dt
               denom_d2=ar_drp(i)-fric_model_drp(i,2)*dt
               denom_d3=ar_drp(i)-fric_model_drp(i,3)*dt
               IF(ABS(denom_d1).lt.1.d-8) denom_d1=SIGN(1.d0,denom_d1)*MAX(1.d-8,ABS(denom_d1))
               IF(ABS(denom_d2).lt.1.d-8) denom_d2=SIGN(1.d0,denom_d2)*MAX(1.d-8,ABS(denom_d2))
               IF(ABS(denom_d3).lt.1.d-8) denom_d3=SIGN(1.d0,denom_d3)*MAX(1.d-8,ABS(denom_d3))
               vd_t(i,1)=(ar_drp(i)*vd_o(i,1)+bm1*dt)/denom_d1
               vd_t(i,2)=(ar_drp(i)*vd_o(i,2)+bm2*dt)/denom_d2
               vd_t(i,3)=(ar_drp(i)*vd_o(i,3)+bm3*dt)/denom_d3
            ENDIF
!
         ENDDO
      ENDIF
!
!.....1st order Euler time advancing; do 1st block if iter_mom=1, otherwise, do 1st and 2nd block
!
      IF((imp_mom_diff+imp_mom_conv).eq.0) iter_mom=1
!
      nrhs=ndim+2
!     
      DO iter=1,iter_mom
!
         IF(iter.eq.1)THEN
!DIR$ SIMD
            DO i=1,ncell_fluid
               gamma_EV= MAX(gamma(i),0.d0)
               gamma_CD=-MIN(gamma(i),0.d0)
               gamma_EV_w= MAX(gamma_wall(i),0.d0)
               gamma_CD_w=-MIN(gamma_wall(i),0.d0)               
               yeta1=1.d0-cell%yeta(i)
!
               am(i,1,1)=     ag(i)+(cell%vFgl(i)*drift_c1(i)+cell%vFgd(i) +             gamma_EV+gamma_EV_w   +cell%vfwg(i))*dt
               am(i,1,2)=-avm_gl(i)-(cell%vFgl(i)*drift_c0(i)              +       yeta1*gamma_EV+gamma_EV_w                )*dt
               am(i,1,3)=-avm_gd(i)-(                         cell%vFgd(i) +cell%yeta(i)*gamma_EV                           )*dt
!
               am(i,2,1)=-avm_gl(i)-(cell%vFgl(i)*drift_c0(i)              +       yeta1*gamma_CD+gamma_CD_w                )*dt
               am(i,2,2)=     al(i)+(cell%vFgl(i)*drift_c1(i)+cell%dentr(i)+       yeta1*gamma_CD+gamma_CD_w   +cell%vfwl(i))*dt
               am(i,2,3)=                                    -cell%dentr(i)                                                  *dt
!
               am(i,3,1)=-avm_gd(i)-(cell%vFgd(i)                          +cell%yeta(i)*gamma_CD)*dt
               am(i,3,2)=                                    -cell%entr(i)                        *dt
               am(i,3,3)=     ad(i)+(cell%vFgd(i)            +cell%entr(i) +cell%yeta(i)*gamma_CD)*dt
!
               DO ix=1,ndim
                  bm(i,1,ix)=ar_gas(i)*vg_t(i,ix)-(avm_gl(i)*(vl_o(i,ix)-vg_o(i,ix))+avm_gd(i)*(vd_o(i,ix)-vg_o(i,ix)))
                  bm(i,2,ix)=ar_liq(i)*vl_t(i,ix)-(avm_gl(i)*(vg_o(i,ix)-vl_o(i,ix)))
                  bm(i,3,ix)=ar_drp(i)*vd_t(i,ix)-(avm_gd(i)*(vg_o(i,ix)-vd_o(i,ix)))
               ENDDO
               bm(i,1,ndim+1)=cell%alphag(i)
               bm(i,2,ndim+1)=cell%alphal(i)
               bm(i,3,ndim+1)=cell%alphad(i)
!
               bm(i,1,ndim+2)=dt
               bm(i,2,ndim+2)=dt
               bm(i,3,ndim+2)=dt
!
            ENDDO
            CALL luinverse3m(am,ip,bm,ncell_fluid,nrhs)
!
            IF(ndim.eq.2) THEN
               DO i=1,ncell_fluid
                  vg_t(i,1)=bm(i,1,1)
                  vg_t(i,2)=bm(i,1,2)
                  vl_t(i,1)=bm(i,2,1)
                  vl_t(i,2)=bm(i,2,2)
                  vd_t(i,1)=bm(i,3,1)
                  vd_t(i,2)=bm(i,3,2)
                  coefp_g(i)=bm(i,1,ndim+1)
                  coefp_l(i)=bm(i,2,ndim+1)
                  coefp_d(i)=bm(i,3,ndim+1)
!
                  coefm_g(i)=bm(i,1,ndim+2)
                  coefm_l(i)=bm(i,2,ndim+2)
!
               ENDDO
            ELSE
               DO i=1,ncell_fluid
                  vg_t(i,1)=bm(i,1,1)
                  vg_t(i,2)=bm(i,1,2)
                  vg_t(i,3)=bm(i,1,3)
                  vl_t(i,1)=bm(i,2,1)
                  vl_t(i,2)=bm(i,2,2)
                  vl_t(i,3)=bm(i,2,3)
                  vd_t(i,1)=bm(i,3,1)
                  vd_t(i,2)=bm(i,3,2)
                  vd_t(i,3)=bm(i,3,3)
                  coefp_g(i)=bm(i,1,ndim+1)
                  coefp_l(i)=bm(i,2,ndim+1)
                  coefp_d(i)=bm(i,3,ndim+1)
!
                  coefm_g(i)=bm(i,1,ndim+2)
                  coefm_l(i)=bm(i,2,ndim+2)
!
               ENDDO
            ENDIF
!            
  100       FORMAT(i3,3x,100(e14.7,1x))    
            CLOSE(18)                     
!
         ELSE
!
            DO i=1,ncell_fluid
               gamma_EV= MAX(gamma(i),0.d0)
               gamma_CD=-MIN(gamma(i),0.d0)
!
               DO ix=1,ndim
                  avm_gl_vg=avm_gl(i)*vg_t(i,ix)
                  avm_gl_vl=avm_gl(i)*vl_t(i,ix)
                  dvlg=cell%vFgl(i)*(vl_t(i,ix)-vg_t(i,ix))*dt
                  bm(i,1,ix)=ar_gas(i)*vg_n(i,ix)+avm_gl_vg-avm_gl_vl-dvlg+(gamma_EV+gamma_wall(i)+cell%vfwg(i))*vg_t(i,ix)*dt
                  bm(i,2,ix)=ar_liq(i)*vl_n(i,ix)-avm_gl_vg+avm_gl_vl+dvlg+(gamma_CD              +cell%vfwl(i))*vl_t(i,ix)*dt
                  bm(i,3,ix)=0.d0
               ENDDO
            ENDDO
!               
            CALL solve3m(am,ip,bm,ncell_fluid,ndim)
!
            IF(ndim.eq.2) THEN
               DO i=1,ncell_fluid
                  vg_t(i,1)=bm(i,1,1)
                  vg_t(i,2)=bm(i,1,2)
                  vl_t(i,1)=bm(i,2,1)
                  vl_t(i,2)=bm(i,2,2)
               ENDDO
            ELSE
               DO i=1,ncell_fluid
                  vg_t(i,1)=bm(i,1,1)
                  vg_t(i,2)=bm(i,1,2)
                  vg_t(i,3)=bm(i,1,3)
                  vl_t(i,1)=bm(i,2,1)
                  vl_t(i,2)=bm(i,2,2)
                  vl_t(i,3)=bm(i,2,3)
               ENDDO
            ENDIF 
!
         ENDIF 
!
!........Implicit momentum convection and diffsuion
!
         IF((imp_mom_diff+imp_mom_conv).gt.0)THEN
            CALL implicit_momentum(iter)
         ELSE
            IF(ndim.eq.2) THEN
               DO i=1,ncell_fp
                  vl_n(i,1)=vl_t(i,1)
                  vl_n(i,2)=vl_t(i,2)
                  vg_n(i,1)=vg_t(i,1)
                  vg_n(i,2)=vg_t(i,2)
                  vd_n(i,1)=vd_t(i,1)
                  vd_n(i,2)=vd_t(i,2)
               ENDDO
            ELSE
               DO i=1,ncell_fp
                  vl_n(i,1)=vl_t(i,1)
                  vl_n(i,2)=vl_t(i,2)
                  vl_n(i,3)=vl_t(i,3)
                  vg_n(i,1)=vg_t(i,1)
                  vg_n(i,2)=vg_t(i,2)
                  vg_n(i,3)=vg_t(i,3)
                  vd_n(i,1)=vd_t(i,1)
                  vd_n(i,2)=vd_t(i,2)
                  vd_n(i,3)=vd_t(i,3)
               ENDDO
            ENDIF
         ENDIF
!
      ENDDO 
!
!.....Body forces defined at cell face
!
      IF(ndim.eq.2) THEN
         DO i=1,nf_non
            vl_f_non(i,1)=0.d0
            vl_f_non(i,2)=0.d0
            vg_f_non(i,1)=0.d0
            vg_f_non(i,2)=0.d0
         ENDDO
      ELSE
         DO i=1,nf_non
            vl_f_non(i,1)=0.d0
            vl_f_non(i,2)=0.d0
            vl_f_non(i,3)=0.d0
            vg_f_non(i,1)=0.d0
            vg_f_non(i,2)=0.d0
            vg_f_non(i,3)=0.d0
         ENDDO
      ENDIF
!
      IF(fric_face.gt.0) CALL f_fric_face
!
      CALL f_nd_face
!
      IF(relax_u.gt.0.d0) THEN
         IF(ndim.eq.2) THEN
            DO i=1,ncell_fp
               vl_n(i,1)=vl_o(i,1)+(1.d0-relax_u)*(vl_n(i,1)-vl_o(i,1))
               vl_n(i,2)=vl_o(i,2)+(1.d0-relax_u)*(vl_n(i,2)-vl_o(i,2))
               vg_n(i,1)=vg_o(i,1)+(1.d0-relax_u)*(vg_n(i,1)-vg_o(i,1))
               vg_n(i,2)=vg_o(i,2)+(1.d0-relax_u)*(vg_n(i,2)-vg_o(i,2))
            ENDDO
         ELSE
            DO i=1,ncell_fp
               vl_n(i,1)=vl_o(i,1)+(1.d0-relax_u)*(vl_n(i,1)-vl_o(i,1))
               vl_n(i,2)=vl_o(i,2)+(1.d0-relax_u)*(vl_n(i,2)-vl_o(i,2))
               vl_n(i,3)=vl_o(i,3)+(1.d0-relax_u)*(vl_n(i,3)-vl_o(i,3))
               vg_n(i,1)=vg_o(i,1)+(1.d0-relax_u)*(vg_n(i,1)-vg_o(i,1))
               vg_n(i,2)=vg_o(i,2)+(1.d0-relax_u)*(vg_n(i,2)-vg_o(i,2))
               vg_n(i,3)=vg_o(i,3)+(1.d0-relax_u)*(vg_n(i,3)-vg_o(i,3))
            ENDDO
         ENDIF
      ENDIF
!
!.....fluxBC model: choke model, mcp model, valve_model
!      
      dpdx_o=dpdx       
!         
      END SUBROUTINE calc_momentum
