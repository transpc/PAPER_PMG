!
      SUBROUTINE set_outlet_property
!
!     This routine calculate the physical properties of the outlet pressure boundary cells
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: ndim      
      USE Zvec_param   , ONLY: nf_nonk,nf_non
      USE Znum_cell    , ONLY: istart_nf, &
                               right_nb_k, &
                               nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Zbc_index    , ONLY: npb
      USE Zb_condition , ONLY: alpha_gas_nd,alpha_liq_nd,alpha_drp_nd,e_gas_nd,e_liq_nd,e_drp_nd,quala_nd
      USE Zuserdefined , ONLY: udfl_outlet_property
      USE Zbc_index    , ONLY: i_horizontal_outlet
      USE Zvector      , ONLY: vg_o,vl_o,vd_o
      USE Zvec_geo     , ONLY: djir_non
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,k
      INTEGER :: ii,kk,ix
      INTEGER :: nv,nf_number,istart0,istart,len,i0,i1
      REAL(8) :: wsum,wsumr,el,eg,al,ag,ed,ad,quala,ha
      REAL(8) :: vl,vg,vd
      REAL(8) :: rdji
!.....Local arrays
!      REAL(8),DIMENSION(ncell_fluid) :: vel,veg,val,vag,vad,vquala,vha,vwsum
!      REAL(8),DIMENSION(ncell_fluid,ndim) :: vvl,vvg,vvd
!.....Local vector arrays
!      REAL(8),DIMENSION(nf_nonk+nf_non) :: el_nf,eg_nf,al_nf,ag_nf,ad_nf,quala_nf,rdji_nf,ha_nf
!      REAL(8),DIMENSION(nf_nonk+nf_non,ndim) :: vlo_nf,vgo_nf,vdo_nf
!
! sang test
      REAL(8),DIMENSION(:), ALLOCATABLE::vel,veg,val,vag,vad,vquala,vha,vwsum,       &
              el_nf,eg_nf,al_nf,ag_nf,ad_nf,quala_nf,rdji_nf,ha_nf   
      REAL(8),DIMENSION(:,:), ALLOCATABLE:: vvl,vvg,vvd,                             &  
               vlo_nf,vgo_nf,vdo_nf
! 
      ALLOCATE(vel(ncell_fluid),veg(ncell_fluid),val(ncell_fluid),vag(ncell_fluid),    &
              vad(ncell_fluid),vquala(ncell_fluid),vha(ncell_fluid),vwsum(ncell_fluid))
              
      ALLOCATE(vvl(ncell_fluid,ndim),vvg(ncell_fluid,ndim),vvd(ncell_fluid,ndim))
      ALLOCATE(el_nf(nf_nonk+nf_non),eg_nf(nf_nonk+nf_non),al_nf(nf_nonk+nf_non),           &
               ag_nf(nf_nonk+nf_non),ad_nf(nf_nonk+nf_non),quala_nf(nf_nonk+nf_non),        &
               rdji_nf(nf_nonk+nf_non),ha_nf(nf_nonk+nf_non))     
      ALLOCATE(vlo_nf(nf_nonk+nf_non,ndim),vgo_nf(nf_nonk+nf_non,ndim),vdo_nf(nf_nonk+nf_non,ndim))
      
          
!.....Build summation info for non
!
      nf_number_nb=0
      nf_number_id(-1)=-1
      nf_number_id(0)=0
      istart_nfs(-1)=0
      istart_nfs(0)=istart_nfs(-1)+nf_nonk
      lens         =istart_nfs(0)+nf_non
!
      nv=0
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)
         IF(npb(ii).gt.0 .and. npb(kk).eq.0) THEN
            rdji=djir_non(i)
            el_nf(i0)   =rdji*cell%el(kk)
            eg_nf(i0)   =rdji*cell%eg(kk)
            al_nf(i0)   =rdji*cell%alphal(kk)
            ag_nf(i0)   =rdji*cell%alphag(kk)
            ad_nf(i0)   =rdji*cell%alphad(kk)
            quala_nf(i0)=rdji*cell%quala(kk)
            ha_nf(i0)   =rdji*cell%ha(kk)
            rdji_nf(i0) =rdji
            DO ix=1,ndim
               vlo_nf(i0,ix)=rdji*vl_o(kk,ix)
               vgo_nf(i0,ix)=rdji*vg_o(kk,ix)
               vdo_nf(i0,ix)=rdji*vd_o(kk,ix)
            ENDDO
         ELSE
            el_nf(i0)   =0.d0
            eg_nf(i0)   =0.d0
            al_nf(i0)   =0.d0
            ag_nf(i0)   =0.d0
            ad_nf(i0)   =0.d0
            quala_nf(i0)=0.d0
            ha_nf(i0)   =0.d0
            rdji_nf(i0) =0.d0
            DO ix=1,ndim
               vlo_nf(i0,ix)=0.d0
               vgo_nf(i0,ix)=0.d0
               vdo_nf(i0,ix)=0.d0
            ENDDO
         ENDIF
      ENDDO
!
      nv=-1
      nf_number=nf_number_id(nv)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         k=right_nb_k(i)
         ii=right_non(k)
         kk=left_nf(k)
         IF(npb(ii).gt.0 .and. npb(kk).eq.0) THEN
            rdji=djir_non(k)
            el_nf(i)   =rdji*cell%el(kk)
            eg_nf(i)   =rdji*cell%eg(kk)
            al_nf(i)   =rdji*cell%alphal(kk)
            ag_nf(i)   =rdji*cell%alphag(kk)
            ad_nf(i)   =rdji*cell%alphad(kk)
            quala_nf(i)=rdji*cell%quala(kk)
            ha_nf(i)   =rdji*cell%ha(kk)
            rdji_nf(i) =rdji
            DO ix=1,ndim
               vlo_nf(i,ix)=rdji*vl_o(kk,ix)
               vgo_nf(i,ix)=rdji*vg_o(kk,ix)
               vdo_nf(i,ix)=rdji*vd_o(kk,ix)
            ENDDO
         ELSE
            el_nf(i)   =0.d0
            eg_nf(i)   =0.d0
            al_nf(i)   =0.d0
            ag_nf(i)   =0.d0
            ad_nf(i)   =0.d0
            quala_nf(i)=0.d0
            ha_nf(i)   =0.d0
            rdji_nf(i) =0.d0
            DO ix=1,ndim
               vlo_nf(i,ix)=0.d0
               vgo_nf(i,ix)=0.d0
               vdo_nf(i,ix)=0.d0
            ENDDO
         ENDIF
      ENDDO
!
      CALL sum_nf(0,0,             &
                  el_nf   ,vel,    &
                  eg_nf   ,veg,    &
                  al_nf   ,val,    &
                  ag_nf   ,vag,    &
                  ad_nf   ,vad,    &
                  quala_nf,vquala, &
                  ha_nf   ,vha,    &
                  rdji_nf ,vwsum)
!
      CALL sum_nf_ndim(0,0,ncell_fluid, &
                       vlo_nf,vvl,      &
                       vgo_nf,vvg,      &
                       vdo_nf,vvd)
!
!..... Neumann momentum B.C
!
      DO ix=1,ndim
         DO i=1,ncell_fluid
            IF(npb(i).eq.0) cycle
            wsum=vwsum(i)
            vl=vvl(i,ix)
            vg=vvg(i,ix)
            vd=vvd(i,ix)
            IF(wsum.ne.0.0d0)THEN
               wsumr=1.d0/wsum
               vl_o(i,ix)=vl*wsumr
               vg_o(i,ix)=vg*wsumr
               vd_o(i,ix)=vd*wsumr
            ELSE
               vl_o(i,ix)=vl
               vg_o(i,ix)=vg
               vd_o(i,ix)=vd
            ENDIF
         ENDDO
      ENDDO
!
      DO i=1,ncell_fluid
         IF(npb(i).eq.0) cycle
         wsum=vwsum(i)
         el=vel(i)
         eg=veg(i)
         ed=vel(i)
         al=val(i)
         ag=vag(i)
         ad=vad(i)
         quala=vquala(i)
         ha=vha(i)
!
         IF(wsum.ne.0.0d0)THEN
            wsumr=1.d0/wsum
         ELSE
            wsumr=1.d0
         ENDIF
!
!........alpha_gas_nd(1) > 1.0d1 : Press. B.C. + Neumann scalar B.C.
!
         k=npb(i)
         IF(alpha_gas_nd(k).gt.1.d1.and.wsum.ne.0.d0)THEN
            cell%el(i)    =el*wsumr
            cell%eg(i)    =eg*wsumr
            cell%ed(i)    =ed*wsumr
            cell%alphal(i)=al*wsumr
            cell%alphag(i)=ag*wsumr
            cell%alphad(i)=ad*wsumr
            cell%quala(i) =quala*wsumr
            cell%ha(i)    =ha*wsumr
         ELSEIF(alpha_gas_nd(k).gt.1.d1.and.wsum.eq.0.d0)THEN
            cell%el(i)    =e_liq_nd(k)
            cell%eg(i)    =e_gas_nd(k)
            cell%ed(i)    =e_drp_nd(k)
            cell%alphal(i)=cell%alphal(i)
            cell%alphag(i)=cell%alphag(i)
            cell%alphad(i)=cell%alphad(i)
            cell%quala(i) =quala_nd(k)
         ELSE
            cell%el(i)    =e_liq_nd(k)
            cell%eg(i)    =e_gas_nd(k)
            cell%ed(i)    =e_drp_nd(k)
            cell%alphal(i)=alpha_liq_nd(k)
            cell%alphag(i)=alpha_gas_nd(k)
            cell%alphad(i)=alpha_drp_nd(k)
            cell%quala(i) =quala_nd(k)
         ENDIF
      ENDDO
!
!.....Must loop on ncell_fp as cell%.._o never communicated
!
      DO i=1,ncell_fp
         cell%el_o(i)    =cell%el(i)
         cell%eg_o(i)    =cell%eg(i)
         cell%ed_o(i)    =cell%ed(i)
         cell%alphal_o(i)=cell%alphal(i)
         cell%alphag_o(i)=cell%alphag(i)
         cell%alphad_o(i)=cell%alphad(i)
         cell%quala_o(i) =cell%quala(i)
      ENDDO
!
!.....UDFs for outlet boundary control
!.....udfl_outlet_property=problem-dependent outlet condition (cavitation,plume,rocom)
!.....i_horizontal_outlet==problem-dependent outlet condition (horizontal_flow)
!
      IF(udfl_outlet_property) CALL udfn_outlet_property
      IF(i_horizontal_outlet.ge.1) CALL udfn_horizontal_outlet
!
      
! sang test
      DEALLOCATE(vel,veg,val,vag,vad,vquala,vha,vwsum,                 &
              el_nf,eg_nf,al_nf,ag_nf,ad_nf,quala_nf,rdji_nf,ha_nf)   
      DEALLOCATE(vvl,vvg,vvd,                                          &  
               vlo_nf,vgo_nf,vdo_nf)
!
! 
      END SUBROUTINE set_outlet_property
