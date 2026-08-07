!
      SUBROUTINE shift_solutions
!
!     This routine shift solutions fo time-marching.
!
      USE VOL_DATA   , ONLY: cell
      USE SOLID_DATA , ONLY: solid
      USE Zmpi       , ONLY: ncell_fp
      USE Zzone      , ONLY: ncell_fluid,ncell_cond
      USE Zparam     , ONLY: ndim
      USE Zbicg      , ONLY: pbcgind,pbcgsig,pbcgind_max,pbcgsig_max
      USE Ziat       , ONLY: ia,ia_old
      USE Zncg       , ONLY: qn_cell,qn_cell_o
      USE Zndforce   , ONLY: vfgl_o,relax_cd,relax_hik
      USE Zvector    , ONLY: vg_o,vl_o,vd_o,vg_n,vl_n,vd_n,ug_o,ul_o,vrel_o
      USE Zqvol      , ONLY: H_ig,H_il,hil_o,hig_o
!
      USE Zvec_param    , ONLY: nf_flux
      USE Zvec_major    , ONLY: flux_l_nf,flux_g_nf,flux_d_nf,       &
                                flux_l_nf_o,flux_g_nf_o,flux_d_nf_o
!      
      IMPLICIT NONE
!
      INTEGER i,ix,i1
      REAL(8) vr
!      
      cell%alphal_o(:)=cell%alphal(:)
      cell%alphag_o(:)=cell%alphag(:)
      cell%alphad_o(:)=cell%alphad(:)
      cell%eg_o(:)    =cell%eg(:)
      cell%ed_o(:)    =cell%ed(:)
      cell%el_o(:)    =cell%el(:)
      cell%tl_o(:)    =cell%tl(:)
      cell%tg_o(:)    =cell%tg(:)
      cell%quala_o(:) =cell%quala(:)
      cell%estm_o(:)  =cell%estm(:)
      cell%p_o(:)     =cell%p(:)
      cell%pps_o(:)   =cell%pps(:)
      cell%ts_o(:)    =cell%ts(:)
      cell%rhol_o(:)  =cell%rhol(:)
      cell%rhog_o(:)  =cell%rhog(:)
      qn_cell_o(:,:)  =qn_cell(:,:)
      IF(ndim.eq.2) THEN
         DO i=1,ncell_fp
            vl_o(i,1)=vl_n(i,1)
            vl_o(i,2)=vl_n(i,2)
            vg_o(i,1)=vg_n(i,1)
            vg_o(i,2)=vg_n(i,2)
            vd_o(i,1)=vd_n(i,1)
            vd_o(i,2)=vd_n(i,2)
         ENDDO
      ELSE
         DO i=1,ncell_fp
            vl_o(i,1)=vl_n(i,1)
            vl_o(i,2)=vl_n(i,2)
            vl_o(i,3)=vl_n(i,3)
            vg_o(i,1)=vg_n(i,1)
            vg_o(i,2)=vg_n(i,2)
            vg_o(i,3)=vg_n(i,3)
            vd_o(i,1)=vd_n(i,1)
            vd_o(i,2)=vd_n(i,2)
            vd_o(i,3)=vd_n(i,3)
         ENDDO
      ENDIF
!
      IF(ncell_cond.gt.0)solid%tsol_o(:)=solid%tsol(:)
!
      ia_old(:)=ia(:)
      cell%cboron_o(:)=cell%cboron(:)
!
      DO i1=1,nf_flux
         flux_l_nf_o(i1)=flux_l_nf(i1)
         flux_g_nf_o(i1)=flux_g_nf(i1)
         flux_d_nf_o(i1)=flux_d_nf(i1)
      ENDDO
!
      IF(relax_cd.ne.0.0d0)vfgl_o(:)=cell%vfgl(:)
!
!.....Relative velocity for K-3, calc_momentum, heat_partition
!
      ul_o(:)=0.0d0
      ug_o(:)=0.0d0
      vrel_o(:)=0.0d0
      DO i=1,ncell_fluid
         DO ix=1,ndim
            ul_o(i)=ul_o(i)+vl_o(i,ix)*vl_o(i,ix)
            ug_o(i)=ug_o(i)+vg_o(i,ix)*vg_o(i,ix)
            vr=vl_o(i,ix)-vg_o(i,ix)
            vrel_o(i)=vrel_o(i)+vr*vr
         ENDDO
         ul_o(i)=DSQRT(ul_o(i))
         ug_o(i)=DSQRT(ug_o(i))
         vrel_o(i)=DSQRT(vrel_o(i))
      ENDDO
!
!.....Reset pressure iteration indicator 
!
      pbcgind=0
      pbcgind_max=0
      pbcgsig=0
      pbcgsig_max=0
!
!.....Heat transfer coefficients
!    
      IF(relax_hik.gt.1.0d-10)THEN
         hil_o(:)=H_il(:)
         hig_o(:)=H_ig(:)
      ENDIF
!
      END SUBROUTINE shift_solutions
