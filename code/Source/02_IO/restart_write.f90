!
      SUBROUTINE restart_write(nout)
!
!     This routine writes the restart file
!
      USE VOL_DATA                 
      USE Wall_DATA    , ONLY: face
      USE SOLID_DATA   , ONLY: solid
      USE Zparam       , ONLY: ndim
      USE Zb_condition , ONLY: rhob_gas,eb_gas
      USE Zbc_index    , ONLY: nvin
      USE Zconst1      , ONLY: lrestart_overwrite,iat,iheatpart,iturb,nlift,nwlf,ntdf
      USE Zconst2      , ONLY: dt
      USE Zio_unit     , ONLY: unit_restart
      USE Zcore        , ONLY: myrank
      USE Zdel_scalar  , ONLY: del_eg,del_el,del_x,del_ag,del_ad,del_rhog,del_rhol
      USE Zface        , ONLY: qqcell,qecell,qclcell,qcgcell
      USE Ziat         , ONLY: iat_nucl
      USE Zndforce     , ONLY: cwlf,clift,ctd
      USE Zpress       , ONLY: p
      USE Zqvol        , ONLY: gamma,h_il,h_ig,h_gf,qvol_ice_solid,qvol_liq  
      USE Ztimecon     , ONLY: time,itim
      USE Zturb        , ONLY: turb_ke_o,turb_dp_o,turb_ke,turb_dp,turb_keg_o,turb_dpg_o,turb_keg,turb_dpg      
      USE Zvector      , ONLY: vl_n,vd_n,vg_n,vl_o,vd_o,vg_o
      USE Zzone        , ONLY: ncell_fluid,ncell_cond,ncell_cond_all
      USE Zuserdefined , ONLY: user_iary,user_rary
!
      USE Zrv_model    , ONLY: rv_model                              
!
      USE Zvec_major    , ONLY: flux_l_nf,flux_g_nf,flux_d_nf,rad_ir_nbcon
      USE Zvec_param    , ONLY: nf_flux,nf_tot_nbcon
!
      USE Zbc_index    , ONLY: nvin,npin,vin_norm
      USE Zb_condition , ONLY: e_gas_nd,e_liq_nd,                       &
                               eb_gas,eb_liq,                           &
                               rho_gas_nd,rho_liq_nd,                   &
                               rhob_gas,rhob_liq,                       &
                               t_gas_nd,t_liq_nd,                       &
                               tb_gas,tb_liq,                           &
                               qualab,quala_nd,pbnd,p_fb,               &
                               vb_liq,vb_gas,vb_drp,                    &
                               vin_liq,vin_gas,vin_Drp
      USE Zncg         , ONLY: tao,cvao_nvin,uao_nvin,dcva_nvin,ra_nvin,  &
                                cvao_npin,uao_npin,dcva_npin,ra_npin,qn_nvin,qn_npin
      USE Zmars        , ONLY: n_marsbc                          
      USE Zqvol        , ONLY: qrv_gas,qrv_liq,qrv_gamma   
      USE unitManager  , ONLY: createUnit
!
      IMPLICIT NONE
!
      INCLUDE '../10_LinkToMARS/c3com.h'
!
!     input
      INTEGER :: nout
!     local variables
      INTEGER :: i,ix,i1,j
      INTEGER :: na,na_c
!      
      CHARACTER*16 myrank_chr
      CHARACTER(50) f_restart      
!      
      LOGICAL, SAVE::INITIAL_r
!
      DATA INITIAL_r /.TRUE./
!
!     return
      na=ncell_fluid
      na_c=ncell_cond 
!
      cell%ed(:)=cell%el(:)      
      cell%ed_o(:)=cell%el_o(:)        
!
!.....Restart for RV
!
      IF(rv_model.eq.1)THEN
        CALL restart_write_rv(nout)      
      ENDIF
!
!.....Write restart file (open separate files by myrank number )
!
      IF(INITIAL_r)THEN
         WRITE(myrank_chr,*)myrank
         myrank_chr=adjustl(myrank_chr)
         f_restart='restart'//trim(myrank_chr)//'.dat'
         ! OPEN(700+myrank,file=f_restart,form='unformatted')
         unit_restart=createUnit("restart write")
         unit_restart=700 !+myrank
         OPEN(unit_restart,file=f_restart,form='unformatted')
         INITIAL_r=.FALSE.         
      ELSE
         IF(lrestart_overwrite)THEN
            CLOSE(unit_restart)
            WRITE(myrank_chr,*)myrank
            myrank_chr=adjustl(myrank_chr)
            f_restart='restart'//trim(myrank_chr)//'.dat'
            OPEN(unit_restart,file=f_restart,form='unformatted')
            INITIAL_r=.FALSE.   
         ENDIF                    
      ENDIF

!
!.....key parameter
!
         WRITE(unit_restart)nout,itim,time,dt
!
!.....typed-cell data
!
         DO i=1,na
            WRITE(unit_restart)cell%alphag(i),cell%alphal(i),cell%alphad(i),cell%quala(i)
            WRITE(unit_restart)cell%rhog(i),cell%rhol(i),cell%rhod(i),cell%rhom(i)
            WRITE(unit_restart)cell%eg(i),cell%el(i),cell%ed(i),cell%p(i),cell%pps(i)
            WRITE(unit_restart)cell%hg(i),cell%hl(i),cell%hgsat(i),cell%hlsat(i),cell%egsat(i),cell%elsat(i),cell%quals(i),cell%ha(i)
            WRITE(unit_restart)cell%tg(i),cell%tl(i),cell%td(i),cell%ts(i),cell%tst(i)
            WRITE(unit_restart)cell%drhogdp(i),cell%drhogde(i),cell%drhogdx(i)
            WRITE(unit_restart)cell%drholdp(i),cell%drholde(i)
            WRITE(unit_restart)cell%dtgdp(i),cell%dtgde(i),cell%dtgdx(i)
            WRITE(unit_restart)cell%dtldp(i),cell%dtlde(i),cell%dtsdp(i),cell%dtsde(i),cell%dtsdx(i)
            WRITE(unit_restart)cell%eviscosg(i),cell%eviscosl(i),cell%eviscosd(i)
            WRITE(unit_restart)cell%lviscosg(i),cell%lviscosl(i),cell%lviscosd(i)
            WRITE(unit_restart)cell%vFgl(i),cell%vFgd(i)
            WRITE(unit_restart)cell%entr(i),cell%dentr(i),cell%yeta(i)
            WRITE(unit_restart)cell%condg(i),cell%condl(i),cell%sigma(i),cell%betag(i),cell%betal(i),cell%cpg(i),cell%cpl(i)
            WRITE(unit_restart)cell%alphag_o(i),cell%alphal_o(i),cell%alphad_o(i),cell%quala_o(i)
            WRITE(unit_restart)cell%eg_o(i),cell%el_o(i),cell%ed_o(i),cell%p_o(i)
            WRITE(unit_restart)cell%tg_o(i),cell%tl_o(i),cell%td_o(i)
            WRITE(unit_restart)cell%tviscosg(i),cell%tviscosl(i),cell%tviscosd(i)
            WRITE(unit_restart)cell%lcondg(i),cell%lcondl(i)
            WRITE(unit_restart)cell%aint1(i),cell%aint2(i),cell%aint3(i),cell%D1(i),cell%D2(i),cell%Ddepart(i),cell%Dlift(i)
            WRITE(unit_restart)cell%regime(i)
            WRITE(unit_restart)cell%vfwg(i),cell%vfwl(i)
            WRITE(unit_restart)cell%rhog_o(i)
            WRITE(unit_restart)cell%cboron(i)
            WRITE(unit_restart)cell%estm(i),cell%estm_o(i),cell%pps_o(i)
         ENDDO
!
!.....TYPEd-face data
!
         DO i=1,na
            WRITE(unit_restart)face%twall_partition(i)
            WRITE(unit_restart)face%wall_fluxl_diff(i),face%wall_fluxg_diff(i),face%wall_fluxd_diff(i)
            WRITE(unit_restart)face%ddepartw(i),face%ratio_evap(i)
         ENDDO
!
!.....TYPEd-solid data
!
      IF(na_c.gt.0)THEN
         DO i=1,na_c
            WRITE(unit_restart)solid%tsol(i)
            WRITE(unit_restart)solid%tsol_o(i)
            WRITE(unit_restart)solid%tsol_max(i),solid%tpellet_surf(i)
            WRITE(unit_restart)qvol_ice_solid(i)
         ENDDO
      ENDIF
!
!.....normal cell data                                                                            
!
         DO i=1,na
            WRITE(unit_restart)p(i)
            WRITE(unit_restart)gamma(i)          
            WRITE(unit_restart)H_il(i),H_ig(i),H_gf(i)
            WRITE(unit_restart)del_eg(i),del_el(i),del_x(i)
            WRITE(unit_restart)del_ag(i),del_ad(i),del_rhog(i),del_rhol(i)
         ENDDO
      IF(iturb.ge.0)THEN
         DO i=1,na
            WRITE(unit_restart)turb_ke_o(i),turb_dp_o(i),turb_ke(i),turb_dp(i)   
         ENDDO
      ENDIF
      IF(iturb.gt.0)THEN
         DO i=1,na
            WRITE(unit_restart)turb_keg_o(i),turb_dpg_o(i),turb_keg(i),turb_dpg(i)    
         ENDDO
      ENDIF      
      IF(nlift.ne.-1.0d0)THEN
         DO i=1,na
            WRITE(unit_restart)Clift(i)
         ENDDO         
      ENDIF
      IF(ntdf.ne.-1.0d0)THEN
         DO i=1,na
            WRITE(unit_restart)Ctd(i)
         ENDDO         
      ENDIF
      IF(nwlf.ne.-1.0d0)THEN
         DO i=1,na
            WRITE(unit_restart)Cwlf(i)
         ENDDO         
      ENDIF
      IF(iat.gt.0)THEN
         DO i=1,na
            WRITE(unit_restart)iat_nucl(i)
         ENDDO
      ENDIF
      IF(iheatpart.gt.0)THEN
         DO i=1,na
            WRITE(unit_restart)qqcell(i),qclcell(i),qcgcell(i)
         ENDDO              
      ENDIF            
      IF(iheatpart.gt.0.or.ncell_cond_all.gt.0)THEN
         DO i=1,na
            WRITE(unit_restart)qecell(i)
         ENDDO 
      ENDIF       
!
!.....i,j array data                                                                              
!
      WRITE(unit_restart)nf_flux
      DO i1=1,nf_flux
         WRITE(unit_restart)flux_l_nf(i1),flux_g_nf(i1),flux_d_nf(i1)
      ENDDO
!
      WRITE(unit_restart)nf_tot_nbcon
      DO i1=1,nf_tot_nbcon
         WRITE(unit_restart)rad_ir_nbcon(i1)
      ENDDO
!
!.....i,ix array data                                                                             
!
         DO i=1,na
            DO ix=1,ndim
               WRITE(unit_restart)vl_n(i,ix),vd_n(i,ix),vg_n(i,ix),vl_o(i,ix),vd_o(i,ix),vg_o(i,ix)
            ENDDO
         ENDDO
!
!.....flow & pressure boundary
!
      WRITE(unit_restart)tao
      WRITE(unit_restart)nvin
      DO i=1,nvin
         WRITE(unit_restart)vin_norm(i),&
                          vb_gas(i,1),vb_gas(i,2),vb_gas(i,ndim),vin_gas(i),&
                          vb_liq(i,1),vb_liq(i,2),vb_liq(i,ndim),vin_liq(i),&
                          vb_drp(i,1),vb_drp(i,2),vb_drp(i,ndim),vin_drp(i),&
                          p_fb(i),tb_liq(i),tb_gas(i),qualab(i),eb_liq(i),eb_gas(i),rhob_liq(i),rhob_gas(i),&
                          cvao_nvin(i),uao_nvin(i),dcva_nvin(i),ra_nvin(i),                                 &
                          qn_nvin(i,1),qn_nvin(i,2),qn_nvin(i,3),qn_nvin(i,4),qn_nvin(i,5),qn_nvin(i,6),qn_nvin(i,7),qn_nvin(i,8)
      ENDDO
!
      WRITE(unit_restart)npin
      DO i=1,npin
         WRITE(unit_restart)pbnd(i),t_liq_nd(i),t_gas_nd(i),quala_nd(i),e_liq_nd(i), e_gas_nd(i),rho_liq_nd(i),rho_gas_nd(i), &
                          cvao_npin(i),uao_npin(i),dcva_npin(i),ra_npin(i),                                                 &
                          qn_npin(i,1),qn_npin(i,2),qn_npin(i,3),qn_npin(i,4),qn_npin(i,5),qn_npin(i,6),qn_npin(i,7),qn_npin(i,8)
      ENDDO
!      
!.....user array
!
      DO i=1,100
         WRITE(unit_restart)user_iary(i),user_rary(i)
      ENDDO
!
!....MARS interface
!
      IF(n_marsbc.gt.0)THEN
        DO i=1,n_marsbc
            WRITE(unit_restart)c3vg(1,i),c3vl(1,i),c3delp(1,i),c3alphf(1,i),c3betaf(1,i),c3alphg(1,i),c3betag(1,i),c3xi(1,i)  
        ENDDO
        DO i=1,n_marsbc
           DO j=1,n_marsbc
               WRITE(unit_restart)c3yeta(1,i,j)
           ENDDO 
        ENDDO         
      ENDIF 
!      
!.....rv core power                                                                            
!
         DO i=1,na
            WRITE(unit_restart)qrv_liq(i),qrv_gas(i),qrv_gamma(i)
         ENDDO
!      
!.....core power
!    
        DO i=1,na
            WRITE(unit_restart)qvol_liq(i)
        ENDDO   
!           
      RETURN
      END SUBROUTINE restart_write
