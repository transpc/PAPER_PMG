!
      SUBROUTINE restart_read(nout)
!
!     This routine reads the restart file
!
      USE VOL_DATA                 
      USE Wall_DATA    , ONLY: face
      USE SOLID_DATA   , ONLY: solid
      USE Zzone        , ONLY: ncell_fluid,ncell_cond,ncell_cond_all
      USE Zparam       , ONLY: ndim
      USE Zcore        , ONLY: myrank
      USE Zconst1      , ONLY: restart,lrestart_changed_nbcon,iat,iheatpart,iturb,nlift,nwlf,ntdf
      USE Zio_unit     , ONLY: unit_restart,unit_saveout,unit_log
      USE Zconst2      , ONLY: dt
      USE Zbc_index    , ONLY: nvin
      USE Zb_condition , ONLY: rhob_gas,eb_gas
      USE Zdel_scalar  , ONLY: del_eg,del_el,del_x,del_ag,del_ad,del_rhog,del_rhol
      USE Zface        , ONLY: qqcell,qecell,qclcell,qcgcell      
      USE Ziat         , ONLY: iat_nucl
      USE Zndforce     , ONLY: cwlf,clift,ctd
      USE Zpress       , ONLY: p
      USE Zqvol        , ONLY: gamma,h_il,h_ig,h_gf,qvol_ice_solid,qvol_liq  
      USE Ztimecon     , ONLY: time,itim,itim_last,nbline
      USE Zturb        , ONLY: turb_ke_o,turb_dp_o,turb_ke,turb_dp,turb_keg_o,turb_dpg_o,turb_keg,turb_dpg      
      USE Zvector      , ONLY: vl_n,vd_n,vg_n,vl_o,vd_o,vg_o
      USE Zuserdefined , ONLY: user_iary,user_rary
      USE Zrv_model    , ONLY: rv_model      
!
      USE Zvec_param   , ONLY: nf_flux,nf_tot_nbcon
      USE Zvec_major   , ONLY: flux_l_nf,flux_g_nf,flux_d_nf,rad_ir_nbcon
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
      use unitManager  , only: createUnit                          
!   
      IMPLICIT NONE
!
      INCLUDE '../10_LinkToMARS/c3com.h'
!
      INTEGER nout,i,ix,i1,restart_start,nf_flux_trans,j,loop,nf_tot_nbcon_trans
      INTEGER na,na_c
      INTEGER restart_step
      REAL(8) tmp1,tmp2,tmp3
!
      CHARACTER*16 myrank_chr
      CHARACTER(50) f_restart  
!
      DATA restart_step/1/
!      
      IF(lrestart_changed_nbcon)CALL user_def_inp(-1) 
!
!.....update time before reading restart file
!
      IF(myrank.eq.0)WRITE(*,"(11x,a)")'Read saveout.dat...'
      OPEN(unit_saveout,file="saveout.dat",status='old',iostat=restart_step)
      IF(restart_step.ne.0)then
         WRITE(*,*)'          No saveout.dat at myrank= ',myrank
         WRITE(unit_log,*)'          No saveout.dat file at myrank= ',myrank
         STOP
      ENDIF     
      DO i=1,nbline
      READ(unit_saveout,*)dt,time,nout,restart_start      
         IF(restart.eq.2) THEN 
            IF(restart_start.eq.itim_last)EXIT
         ELSE
            IF(restart_start.eq.itim)EXIT
         ENDIF  
      ENDDO  
      CLOSE(unit_saveout)
!
!.....Restart for RV
!
      IF(rv_model.eq.1)THEN
          CALL restart_read_rv
      ENDIF      
!
      na=ncell_fluid
      na_c=ncell_cond       
!
!.....Read restart file (open separate files by myrank number )
!
      IF(myrank.eq.0)WRITE(*,"(11x,a)")'Read restart?.dat...'
      WRITE(myrank_chr,*)myrank
      myrank_chr=adjustl(myrank_chr)
      f_restart='restart'//trim(myrank_chr)//'.dat'
!      write(*,*) f_restart
      unit_restart=createUnit("restart read")
      unit_restart=700 !+myrank
      OPEN(unit_restart,file=f_restart,status='old',form='unformatted',iostat=restart_step)
      ! OPEN(700+myrank,file=f_restart,status='old',form='unformatted',iostat=restart_step)
      IF(restart_step.ne.0)then
         WRITE(*,*)'          No restart data file at myrank= ',myrank
         WRITE(unit_log,*)'          No restart data file at myrank= ',myrank
         STOP
      ENDIF   
!
      loop=0
      DO
         loop=loop+1
!
!........key parameter
!
         READ(unit_restart,iostat=restart_step)nout,restart_start,time,dt
         ! READ(700+myrank,iostat=restart_step)nout,restart_start,time,dt
         IF(restart_step.ne.0)then
            WRITE(*,"(11x,a,1i3,1i10,1i5,1i5)")'Check restart option,step,nout,set!!! ',restart,itim,nout,loop
            IF(myrank.eq.0) WRITE(unit_log,*)'          Check restart option & step &nout !!! ',restart,itim,nout
            STOP
         ENDIF
!
!........TYPEd-cell data
!
         DO i=1,na
            READ(unit_restart)cell%alphag(i),cell%alphal(i),cell%alphad(i),cell%quala(i)      
            READ(unit_restart)cell%rhog(i),cell%rhol(i),cell%rhod(i),cell%rhom(i)
!
            cell%rhomr(i)=1.0d0/cell%rhom(i)
!            
            READ(unit_restart)cell%eg(i),cell%el(i),cell%ed(i),cell%p(i),cell%pps(i)
            READ(unit_restart)cell%hg(i),cell%hl(i),cell%hgsat(i),cell%hlsat(i),cell%egsat(i),cell%elsat(i),cell%quals(i), cell%ha(i)
            READ(unit_restart)cell%tg(i),cell%tl(i),cell%td(i),cell%ts(i),cell%tst(i)
            READ(unit_restart)cell%drhogdp(i),cell%drhogde(i),cell%drhogdx(i)
            READ(unit_restart)cell%drholdp(i),cell%drholde(i)
            READ(unit_restart)cell%dtgdp(i),cell%dtgde(i),cell%dtgdx(i)
            READ(unit_restart)cell%dtldp(i),cell%dtlde(i),cell%dtsdp(i),cell%dtsde(i),cell%dtsdx(i)
            READ(unit_restart)cell%eviscosg(i),cell%eviscosl(i),cell%eviscosd(i)
            READ(unit_restart)cell%lviscosg(i),cell%lviscosl(i),cell%lviscosd(i)
            READ(unit_restart)cell%vFgl(i),cell%vFgd(i)
            READ(unit_restart)cell%entr(i),cell%dentr(i),cell%yeta(i)
            READ(unit_restart)cell%condg(i),cell%condl(i),cell%sigma(i),cell%betag(i),cell%betal(i),cell%cpg(i),cell%cpl(i)
            READ(unit_restart)cell%alphag_o(i),cell%alphal_o(i),cell%alphad_o(i),cell%quala_o(i)         
            READ(unit_restart)cell%eg_o(i),cell%el_o(i),cell%ed_o(i),cell%p_o(i)
            READ(unit_restart)cell%tg_o(i),cell%tl_o(i),cell%td_o(i)
            READ(unit_restart)cell%tviscosg(i),cell%tviscosl(i),cell%tviscosd(i)
            READ(unit_restart)cell%lcondg(i),cell%lcondl(i)
            READ(unit_restart)cell%aint1(i),cell%aint2(i),cell%aint3(i),cell%D1(i),cell%D2(i),cell%Ddepart(i),cell%Dlift(i)
            READ(unit_restart)cell%regime(i)
            READ(unit_restart)cell%vfwg(i),cell%vfwl(i)
            READ(unit_restart)cell%rhog_o(i)
            READ(unit_restart)cell%cboron(i)
            READ(unit_restart)cell%estm(i),cell%estm_o(i),cell%pps_o(i)
         ENDDO
!
!........TYPEd-face data
!
         DO i=1,na
            READ(unit_restart)face%twall_partition(i)
            READ(unit_restart)face%wall_fluxl_diff(i),face%wall_fluxg_diff(i),face%wall_fluxd_diff(i)
            READ(unit_restart)face%ddepartw(i),face%ratio_evap(i)
         ENDDO
!
!........TYPEd-solid data
!
         IF(na_c.gt.0)THEN 
            DO i=1,na_c
               READ(unit_restart)solid%tsol(i)
               READ(unit_restart)solid%tsol_o(i)
               READ(unit_restart)solid%tsol_max(i),solid%tpellet_surf(i)
               READ(unit_restart)qvol_ice_solid(i)
            ENDDO
         ENDIF
!
!........normal cell data                                                                            
!
         DO i=1,na
            READ(unit_restart)p(i)
            READ(unit_restart)gamma(i)
            READ(unit_restart)H_il(i),H_ig(i),H_gf(i)
            READ(unit_restart)del_eg(i),del_el(i),del_x(i)
            READ(unit_restart)del_ag(i),del_ad(i),del_rhog(i),del_rhol(i)
         ENDDO
      IF(iturb.ge.0)THEN
         DO i=1,na
            READ(unit_restart)turb_ke_o(i),turb_dp_o(i),turb_ke(i),turb_dp(i) 
         ENDDO
      ENDIF
      IF(iturb.gt.0)THEN
         DO i=1,na
            READ(unit_restart)turb_keg_o(i),turb_dpg_o(i),turb_keg(i),turb_dpg(i)    
         ENDDO
      ENDIF      
      IF(nlift.ne.-1.0d0)THEN
         DO i=1,na
            READ(unit_restart)Clift(i)
         ENDDO         
      ENDIF
      IF(ntdf.ne.-1.0d0)THEN
         DO i=1,na
            READ(unit_restart)Ctd(i)
         ENDDO         
      ENDIF
      IF(nwlf.ne.-1.0d0)THEN
         DO i=1,na
            READ(unit_restart)Cwlf(i)
         ENDDO         
      ENDIF        
      IF(iat.gt.0)THEN
         DO i=1,na
            READ(unit_restart)iat_nucl(i)
         ENDDO
      ENDIF
      IF(iheatpart.gt.0)THEN
         DO i=1,na
            READ(unit_restart)qqcell(i),qclcell(i),qcgcell(i)
         ENDDO              
      ENDIF            
      IF(iheatpart.gt.0.or.ncell_cond_all.gt.0)THEN
         DO i=1,na
            READ(unit_restart)qecell(i)
         ENDDO 
      ENDIF           
!
!........i,j array data                                                                              
!
         READ(unit_restart)nf_flux_trans
         DO i1=1,nf_flux_trans
            READ(unit_restart)tmp1,tmp2,tmp3
            IF(i1.le.nf_flux)THEN
               flux_l_nf(i1)=tmp1
               flux_g_nf(i1)=tmp2
               flux_d_nf(i1)=tmp3
            ENDIF
         ENDDO
!
         READ(unit_restart)nf_tot_nbcon_trans
         DO i1=1,nf_tot_nbcon_trans
            READ(unit_restart)tmp1
            IF(i1.le.nf_tot_nbcon)THEN
               rad_ir_nbcon(i1)=tmp1
            ENDIF
         ENDDO
!
!........i,ix array data                                                                             
!
         DO i=1,na
            DO ix=1,ndim
               READ(unit_restart)vl_n(i,ix),vd_n(i,ix),vg_n(i,ix),vl_o(i,ix),vd_o(i,ix),vg_o(i,ix)
            ENDDO
         ENDDO
!
!........flow & pressure boundary
!
         READ(unit_restart)tao
         READ(unit_restart)nvin
         DO i=1,nvin
            READ(unit_restart)vin_norm(i),&
                            vb_gas(i,1),vb_gas(i,2),vb_gas(i,ndim),vin_gas(i),&
                            vb_liq(i,1),vb_liq(i,2),vb_liq(i,ndim),vin_liq(i),&
                            vb_drp(i,1),vb_drp(i,2),vb_drp(i,ndim),vin_drp(i),&
                            p_fb(i),tb_liq(i),tb_gas(i),qualab(i),eb_liq(i),eb_gas(i),rhob_liq(i),rhob_gas(i),&
                            cvao_nvin(i),uao_nvin(i),dcva_nvin(i),ra_nvin(i),                                 &
                          qn_nvin(i,1),qn_nvin(i,2),qn_nvin(i,3),qn_nvin(i,4),qn_nvin(i,5),qn_nvin(i,6),qn_nvin(i,7),qn_nvin(i,8)
         ENDDO
!      
         READ(unit_restart)npin
         DO i=1,npin
            READ(unit_restart)pbnd(i),t_liq_nd(i),t_gas_nd(i),quala_nd(i),e_liq_nd(i), e_gas_nd(i),rho_liq_nd(i),rho_gas_nd(i), &
                            cvao_npin(i),uao_npin(i),dcva_npin(i),ra_npin(i),                                                 &
                          qn_npin(i,1),qn_npin(i,2),qn_npin(i,3),qn_npin(i,4),qn_npin(i,5),qn_npin(i,6),qn_npin(i,7),qn_npin(i,8)
         ENDDO
!      
!.....user array
!
         DO i=1,100
            READ(unit_restart)user_iary(i),user_rary(i)
         ENDDO
!
!....MARS interface
!
      IF(n_marsbc.gt.0)THEN
        DO i=1,n_marsbc
            READ(unit_restart)c3vg(1,i),c3vl(1,i),c3delp(1,i),c3alphf(1,i),c3betaf(1,i),c3alphg(1,i),c3betag(1,i),c3xi(1,i)  
        ENDDO
        DO i=1,n_marsbc
           DO j=1,n_marsbc
               READ(unit_restart)c3yeta(1,i,j)
           ENDDO 
        ENDDO         
      ENDIF 
!
!.....rv core power                                                                            
!
         DO i=1,na
            READ(unit_restart)qrv_liq(i),qrv_gas(i),qrv_gamma(i)
         ENDDO
!
!.....core power
!    
        DO i=1,na
            READ(unit_restart)qvol_liq(i)
        ENDDO 
!-------------------------------------------------------------------------------                 
!
!........Stop reading at designated restarting point                                                                             
!         
         IF(restart.eq.2) THEN
            IF(restart_start.eq.itim_last)EXIT
         ELSE
            IF(restart_start.eq.itim)EXIT
         ENDIF
      ENDDO
      CLOSE(unit_restart)
!
!.....Print screen
!
      IF(restart_start.eq.itim)THEN
         IF(myrank.eq.0)THEN
            WRITE(*,*)'          Restart from itim=',restart_start
            WRITE(unit_log,*)'          Restart from itim=',restart_start
         ENDIF
      ELSE
         IF(myrank.eq.0)THEN      
             WRITE(*,*)'          Restart from the final restart number, itim=',restart_start
             WRITE(unit_log,*)'          Restart from the final restart number, itim=',restart_start
         ENDIF    
         itim=restart_start
      ENDIF
!
!.....Save the first step of restarting point
!
    1 CONTINUE
      IF(myrank.eq.0)OPEN(unit_saveout,file="saveout.dat")  
      ! IF(myrank.eq.0)OPEN(96,file="saveout.dat")  
       CALL restart_write(nout)
      IF(myrank.eq.0)WRITE(unit_saveout,*)dt,time,nout,itim      
!            
      IF(myrank.eq.0)WRITE(*,*)'          End of restart read.'
      IF(myrank.eq.0)WRITE(unit_log,*)'          End of restart read.'
!      
      RETURN
      END SUBROUTINE restart_read
