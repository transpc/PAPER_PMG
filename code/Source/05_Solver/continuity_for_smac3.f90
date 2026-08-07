!
      SUBROUTINE continuity_for_smac3
!
!     This routine solves non-conservative form of continuity equation.
!     This routine is activated when smac=3.
!
      USE Zinterface
      USE VOL_DATA    , ONLY: cell               
      USE Zzone       , ONLY: ncell_fluid
      USE Zconst2     , ONLY: dt
      USE Zvec_param  , ONLY: nf_nonk,nf_non,nf_mcc,nf_inl,nf_out,nf_fluxk
      USE Znum_cell   , ONLY: right_nb_k,istart_nf,                     &
                              nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zvec_index  , ONLY: left_nf,right_non
      USE Zare        , ONLY: ar_liq,ar_gas
      USE Zcoord3     , ONLY: volpr
      USE Zdecoupled  , ONLY: al_min_c,ag_min_c
      USE Zdel_scalar , ONLY: del_el,del_eg,del_x
      USE Zpress      , ONLY: pp
      USE Zqvol       , ONLY: gamma,gamma_wall
      USE Ztimecon    , ONLY: iso_thermal
      USE Zbc_index   , ONLY: npb
      USE Zvec_major  , ONLY: flux_l_nf,flux_g_nf,    &
                              liq_conv_nf,vap_conv_nf
      !OPR1000 rod-scale (EVVD)                              
      USE Zporous     , ONLY: tm_mas_l,tm_mas_g, &
                              vd_mas_l,vd_mas_g, &
                              mixing_vane_l      
      USE Zporous , ONLY: l_subchannel,l_mixing_vane
                                                        
!
      IMPLICIT NONE
!
!.....Local variable
      INTEGER :: i,k
      INTEGER :: ii,kk
      INTEGER :: nv,nf_number,len,istart0,istart,i0,i1
      REAL(8) :: rhog,rhol,vt
      REAL(8) :: dr
      REAL(8) :: ag,al
!.....Local arrays
      REAL(8) :: aal(ncell_fluid),aag(ncell_fluid)
      REAL(8) :: rhog_v(ncell_fluid),rhol_v(ncell_fluid)
      !OPR1000 rod-scale (EVVD)                              
      REAL(8) :: evvd(ncell_fluid),evvdg(ncell_fluid),evvdl(ncell_fluid)
!.....Local vector arrays
      REAL(8),DIMENSION(nf_fluxk) :: al_nf,ag_nf
!
!.....Build summation info for non,mcc,inl,out
!
      nf_number_nb=3
      nf_number_id(-1)=-1
      nf_number_id(0)=0
      nf_number_id(1)=1
      nf_number_id(2)=2
      nf_number_id(3)=3
      istart_nfs(-1)=0
      istart_nfs(0)=istart_nfs(-1)+nf_nonk
      istart_nfs(1)=istart_nfs(0) +nf_non
      istart_nfs(2)=istart_nfs(1) +nf_mcc
      istart_nfs(3)=istart_nfs(2) +nf_inl
      lens         =istart_nfs(3) +nf_out
!
!.....Cells non
!
      nv=0
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len =istart_nf(2,nf_number)
      DO i=1,len
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)
         IF(npb(ii).eq.0) THEN
            al_nf(i0)=-liq_conv_nf(i1)*flux_l_nf(i1)
            ag_nf(i0)=-vap_conv_nf(i1)*flux_g_nf(i1)
         ELSE
            al_nf(i0)=0.d0
            ag_nf(i0)=0.d0
         ENDIF
      ENDDO
!
      nv=-1
      nf_number=nf_number_id(nv)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         k=right_nb_k(i)
         ii=right_non(k)
         IF(npb(ii).eq.0) THEN
            al_nf(i)=liq_conv_nf(k)*flux_l_nf(k)
            ag_nf(i)=vap_conv_nf(k)*flux_g_nf(k)
         ELSE
            al_nf(i)=0.d0
            ag_nf(i)=0.d0
         ENDIF
      ENDDO
!
!.....The rest
!
      DO nv=1,3
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(nv)
         istart=istart_nf(1,nf_number)
         len =istart_nf(2,nf_number)
         DO i=1,len
            i0=istart0+i
            i1=istart+i
            ii=left_nf(i1)
            IF(npb(ii).eq.0) THEN
               al_nf(i0)=-liq_conv_nf(i1)*flux_l_nf(i1)
               ag_nf(i0)=-vap_conv_nf(i1)*flux_g_nf(i1)
            ELSE
               al_nf(i0)=0.d0
               ag_nf(i0)=0.d0
            ENDIF
         ENDDO
      ENDDO
!
      CALL sum_nf(0,0,       &
                  al_nf,aal, &
                  ag_nf,aag)
!
      DO i=1,ncell_fluid
         IF(npb(i).ne.0) cycle
         vt=dt*volpr(i)
         aal(i)=aal(i)*vt
         aag(i)=aag(i)*vt
      ENDDO
!
      DO i=1,ncell_fluid
         IF(npb(i).ne.0) CYCLE
         rhog_v(i)=cell%rhog(i)+cell%drhogdp(i)*pp(i) 
         rhol_v(i)=cell%rhol(i)+cell%drholdp(i)*pp(i) 
      ENDDO
      IF(iso_thermal.eq.0)THEN
         DO i=1,ncell_fluid
            IF(npb(i).ne.0) CYCLE
            rhog_v(i)=rhog_v(i)+cell%drhogde(i)*del_eg(i)+cell%drhogdx(i)*del_x(i)
            rhol_v(i)=rhol_v(i)+cell%drholde(i)*del_el(i)
         ENDDO
      ENDIF
!
      IF(l_subchannel)then
         DO i=1,ncell_fluid
            IF(npb(i).ne.0) CYCLE
            rhog=rhog_v(i)   
            rhol=rhol_v(i)   
!
!...........Convection
!
            al=aal(i)
            ag=aag(i)
!
!...........Calculate liquid fraction
!
!...........EVVD    
            evvdg(i)=tm_mas_g(i)+vd_mas_g(i)
            evvdl(i)=tm_mas_l(i)+vd_mas_l(i)
            IF(l_mixing_vane)then
               evvdl(i)=evvdl(i)+mixing_vane_l(2,i)*volpr(i)
            ENDIF
            evvd(i)=evvdg(i)+evvdl(i)

            IF(cell%alphag_o(i).lt.ag_min_c)THEN
               ag=(ag+(gamma(i)+gamma_wall(i)+evvdg(i))*dt)/cell%rhog(i) !EVVD            
               cell%alphag(i)=cell%alphag_o(i)+ag
               cell%alphag(i)=DMIN1(DMAX1(0.0d0,cell%alphag(i)),1.0d0)
               cell%alphal(i)=1.d0-cell%alphag(i)
            ELSEIF(cell%alphal_o(i).lt.al_min_c)THEN
               al=(al-(gamma(i)+gamma_wall(i))*dt+evvdl(i)*dt)/cell%rhol(i) !EVVD           
               cell%alphal(i)=cell%alphal_o(i)+al
               cell%alphal(i)=DMIN1(DMAX1(0.0d0,cell%alphal(i)),1.0d0)
               cell%alphag(i)=1.d0-cell%alphal(i)
            ELSE
               dr=DMAX1(1.d-5,rhol-rhog)
!              al(i)=(ar_liq(i)+ar_gas(i)+al(i)+ag(i)-rhog)/(rhol-rhog)
               al=(ar_liq(i)+ar_gas(i)+al+ag-rhog)/dr
               al=al+evvd(i)*dt/dr    ! EVVD
               cell%alphal(i)=al
               cell%alphal(i)=DMIN1(DMAX1(0.0d0,cell%alphal(i)),1.0d0)
               cell%alphag(i)=1.d0-cell%alphal(i)
            ENDIF
         ENDDO
      ELSE
         DO i=1,ncell_fluid
            IF(npb(i).ne.0) CYCLE
            rhog=rhog_v(i)   
            rhol=rhol_v(i)   
!
!...........Convection
!
            al=aal(i)
            ag=aag(i)
!
!...........Calculate liquid fraction
!
            IF(cell%alphag_o(i).lt.ag_min_c)THEN
               ag=(ag+(gamma(i)+gamma_wall(i))*dt)/cell%rhog(i)
               cell%alphag(i)=cell%alphag_o(i)+ag
               cell%alphag(i)=DMIN1(DMAX1(0.0d0,cell%alphag(i)),1.0d0)
               cell%alphal(i)=1.d0-cell%alphag(i)
            ELSEIF(cell%alphal_o(i).lt.al_min_c)THEN
               al=(al-(gamma(i)+gamma_wall(i))*dt)/cell%rhol(i)
               cell%alphal(i)=cell%alphal_o(i)+al
               cell%alphal(i)=DMIN1(DMAX1(0.0d0,cell%alphal(i)),1.0d0)
               cell%alphag(i)=1.d0-cell%alphal(i)
            ELSE
               dr=DMAX1(1.d-5,rhol-rhog)
!              al(i)=(ar_liq(i)+ar_gas(i)+al(i)+ag(i)-rhog)/(rhol-rhog)
               al=(ar_liq(i)+ar_gas(i)+al+ag-rhog)/dr
               cell%alphal(i)=al
               cell%alphal(i)=DMIN1(DMAX1(0.0d0,cell%alphal(i)),1.0d0)
               cell%alphag(i)=1.d0-cell%alphal(i)
            ENDIF
         ENDDO
      ENDIF
!
!........Assume no droplet fraction 
!
      DO i=1,ncell_fluid
         cell%alphad(i)=0.0d0
      ENDDO
!
!........Isothermal option to update temperatures
!
      IF(iso_thermal.eq.0)THEN
         DO i=1,ncell_fluid
            IF(npb(i).ne.0) CYCLE
            cell%ts(i)=cell%ts_o(i)+cell%dtsdp(i)*pp(i)
            cell%tg(i)=cell%tg_o(i)+cell%dtgdp(i)*pp(i)
            cell%tl(i)=cell%tl_o(i)+cell%dtldp(i)*pp(i)
         ENDDO
      ENDIF
!
      END SUBROUTINE continuity_for_smac3
