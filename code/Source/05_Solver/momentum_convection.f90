!
      SUBROUTINE momentum_convection(cnvc_liq,cnvc_gas,cnvc_drp)
!
!     This routine calculates convective fluxes through the cell face between cells i and j
!
      USE Zinterface
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: ndim
      USE Zvec_param   , ONLY: nf_nonk,nf_non,nf_mcc,nf_inl,nf_out,nf_flux,nf_fluxk
      USE Znum_cell    , ONLY: istart_nf,istart_nbcon_nf,                      &
                               nf_number_nb,lens,                              &
                               right_nb_k,istart_nfs,nf_number_id
      USE Zvec_index   , ONLY: left_nf,right_non,nbcon_nf
      USE Z2nd_order   , ONLY: mom_conv_2nd
      USE Zare         , ONLY: ar_gas,ar_liq,ar_drp
      USE c3com_cupid  , ONLY: i3invtbl,mcdirect,c3dpv,mcgdirect
      USE Zb_condition , ONLY: alphab_gas,alphab_liq,alphab_drp,     &
                               rhob_gas,rhob_liq,rhob_drp , &
                               vb_liq,vb_gas,vb_drp, &
                               vin_liq,vin_gas,vin_drp
      USE Zbc_index    , ONLY: vin_norm
      USE Zuserdefined , ONLY: vel_bc_profile_inl
      USE Zvector      , ONLY: vg_o,vl_o,vd_o
      USE Zvec_major   , ONLY: flux_l_nf,flux_g_nf,flux_d_nf
      USE Zvec_geo     , ONLY: xn_nf
      USE Zrv_choke    
      USE Zmcp
      USE Zcore   
      USE Zrv_model    , ONLY: rv_mcp,rv_choke,rv_valve
!
      IMPLICIT NONE
!      
      INCLUDE '../10_LinkToMARS/c3com.h'
!      
!.....Output
      REAL(8),DIMENSION(ncell_fluid,ndim) :: cnvc_liq,cnvc_gas,cnvc_drp
!.....Local variables
      INTEGER :: i,ii,kk,k,ix,idx
      INTEGER :: nv,nf_number,istart0,istart,istart2,len,i0,i1,i2
      REAL(8) :: a_l,a_g,a_d
      REAL(8) :: b_l,b_g,b_d
      REAL(8) :: vl1,vl2,vl3
      REAL(8) :: vg1,vg2,vg3
      REAL(8) :: vd1,vd2,vd3
      REAL(8) :: f_profile
!.....Local arrays
      REAL(8),DIMENSION(nf_flux) :: mfluxl_cup_nf,mfluxg_cup_nf,mfluxd_cup_nf
      REAL(8),DIMENSION(nf_flux,ndim) :: fluxl_c_nf,fluxg_c_nf,fluxd_c_nf
      REAL(8),DIMENSION(nf_fluxk,ndim) :: cnvc_l_nf,cnvc_g_nf,cnvc_d_nf
!
!.....Cells non
!
      IF(ndim.eq.2) THEN 
         nf_number=0
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len  
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)           
            a_l=min(flux_l_nf(i1),0.d0)
            b_l=max(flux_l_nf(i1),0.d0)
            a_g=min(flux_g_nf(i1),0.d0)
            b_g=max(flux_g_nf(i1),0.d0)
            a_d=min(flux_d_nf(i1),0.d0)
            b_d=max(flux_d_nf(i1),0.d0)
            mfluxl_cup_nf(i1)=a_l*ar_liq(kk)+b_l*ar_liq(ii)
            mfluxg_cup_nf(i1)=a_g*ar_gas(kk)+b_g*ar_gas(ii)
            mfluxd_cup_nf(i1)=a_d*ar_drp(kk)+b_d*ar_drp(ii)
            fluxl_c_nf(i1,1)=(a_l*ar_liq(kk))*vl_o(kk,1)+(b_l*ar_liq(ii))*vl_o(ii,1)
            fluxg_c_nf(i1,1)=(a_g*ar_gas(kk))*vg_o(kk,1)+(b_g*ar_gas(ii))*vg_o(ii,1)
            fluxd_c_nf(i1,1)=(a_d*ar_drp(kk))*vd_o(kk,1)+(b_d*ar_drp(ii))*vd_o(ii,1)
            fluxl_c_nf(i1,2)=(a_l*ar_liq(kk))*vl_o(kk,2)+(b_l*ar_liq(ii))*vl_o(ii,2)
            fluxg_c_nf(i1,2)=(a_g*ar_gas(kk))*vg_o(kk,2)+(b_g*ar_gas(ii))*vg_o(ii,2)
            fluxd_c_nf(i1,2)=(a_d*ar_drp(kk))*vd_o(kk,2)+(b_d*ar_drp(ii))*vd_o(ii,2)
         ENDDO
!         
!........fluxBC: choke model, mcp model, valve model
!          
         IF(rv_valve.eq.1.or.rv_choke.eq.1.or.rv_mcp.eq.1) CALL fluxBC_convection_ice(fluxl_c_nf,fluxg_c_nf,fluxd_c_nf,mfluxl_cup_nf,mfluxg_cup_nf,mfluxd_cup_nf)
!
!.....Cells mcc
!
         nf_number=1
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len  
            i1=istart+i
            ii=left_nf(i1)
            idx=i3invtbl(i)
            a_l=min(flux_l_nf(i1),0.d0)
            b_l=max(flux_l_nf(i1),0.d0)
            a_g=min(flux_g_nf(i1),0.d0)
            b_g=max(flux_g_nf(i1),0.d0)
            a_d=min(flux_d_nf(i1),0.d0)
            b_d=max(flux_d_nf(i1),0.d0)
            IF(mcdirect(idx).lt.0)THEN
               mfluxl_cup_nf(i1)=a_l*c3dpv(idx,5)+b_l*ar_liq(ii)
               mfluxd_cup_nf(i1)=a_d*c3dpv(idx,4)+b_d*ar_drp(ii)
               vl1=c3vl(1,idx)*xn_nf(i1,1)
               vl2=c3vl(1,idx)*xn_nf(i1,2)
               vd1=c3vl(1,idx)*xn_nf(i1,1)
               vd2=c3vl(1,idx)*xn_nf(i1,2)
               fluxl_c_nf(i1,1)=(a_l*c3dpv(idx,5))*vl1+(b_l*ar_liq(ii))*vl_o(ii,1)
               fluxl_c_nf(i1,2)=(a_l*c3dpv(idx,5))*vl2+(b_l*ar_liq(ii))*vl_o(ii,2)
               fluxd_c_nf(i1,1)=(a_d*c3dpv(idx,4))*vd1+(b_d*ar_drp(ii))*vd_o(ii,1)
               fluxd_c_nf(i1,2)=(a_d*c3dpv(idx,4))*vd2+(b_d*ar_drp(ii))*vd_o(ii,2)
            ELSE
               mfluxl_cup_nf(i1)=flux_l_nf(i1)*ar_liq(ii)
               mfluxd_cup_nf(i1)=flux_d_nf(i1)*ar_drp(ii)
               fluxl_c_nf(i1,1)=(flux_l_nf(i1)*ar_liq(ii))*vl_o(ii,1)
               fluxl_c_nf(i1,2)=(flux_l_nf(i1)*ar_liq(ii))*vl_o(ii,2)
               fluxd_c_nf(i1,1)=(flux_d_nf(i1)*ar_drp(ii))*vd_o(ii,1)
               fluxd_c_nf(i1,2)=(flux_d_nf(i1)*ar_drp(ii))*vd_o(ii,2)
            ENDIF
            IF(mcgdirect(idx).lt.0)THEN
               mfluxg_cup_nf(i1)=a_g*c3dpv(idx,3)+b_g*ar_gas(ii)
               vg1=c3vg(1,idx)*xn_nf(i1,1)
               vg2=c3vg(1,idx)*xn_nf(i1,2)
               fluxg_c_nf(i1,1)=(a_g*c3dpv(idx,3))*vg1+(b_g*ar_gas(ii))*vg_o(ii,1)
               fluxg_c_nf(i1,2)=(a_g*c3dpv(idx,3))*vg2+(b_g*ar_gas(ii))*vg_o(ii,2)
            ELSE
               mfluxg_cup_nf(i1)=flux_g_nf(i1)*ar_gas(ii)
               fluxg_c_nf(i1,1)=(flux_g_nf(i1)*ar_gas(ii))*vg_o(ii,1)
               fluxg_c_nf(i1,2)=(flux_g_nf(i1)*ar_gas(ii))*vg_o(ii,2)
            ENDIF
         ENDDO
!
!.....Cells inl
!
         nf_number=2
         istart=istart_nf(1,nf_number)
         istart2=istart_nbcon_nf(nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len  
            i1=istart+i
            i2=istart2+i
            ii=left_nf(i1)
            k=nbcon_nf(i2)
            a_l=min(flux_l_nf(i1),0.d0)
            b_l=max(flux_l_nf(i1),0.d0)
            a_g=min(flux_g_nf(i1),0.d0)
            b_g=max(flux_g_nf(i1),0.d0)
            a_d=min(flux_d_nf(i1),0.d0)
            b_d=max(flux_d_nf(i1),0.d0)
            mfluxl_cup_nf(i1)=a_l*(alphab_liq(k)*rhob_liq(k))+b_l*ar_liq(ii)
            mfluxg_cup_nf(i1)=a_g*(alphab_gas(k)*rhob_gas(k))+b_g*ar_gas(ii)
            mfluxd_cup_nf(i1)=a_d*(alphab_drp(k)*rhob_drp(k))+b_d*ar_drp(ii)
            f_profile=vel_bc_profile_inl(i)
            IF(vin_norm(k).eq.0)THEN
               vl1=vb_liq(k,1)*f_profile
               vl2=vb_liq(k,2)*f_profile
               vg1=vb_gas(k,1)*f_profile
               vg2=vb_gas(k,2)*f_profile
               vd1=vb_drp(k,1)*f_profile
               vd2=vb_drp(k,2)*f_profile
            ELSE
               vl1=vin_liq(k)*xn_nf(i1,1)*f_profile
               vl2=vin_liq(k)*xn_nf(i1,2)*f_profile
               vg1=vin_gas(k)*xn_nf(i1,1)*f_profile
               vg2=vin_gas(k)*xn_nf(i1,2)*f_profile
               vd1=vin_drp(k)*xn_nf(i1,1)*f_profile
               vd2=vin_drp(k)*xn_nf(i1,2)*f_profile
            ENDIF
            fluxl_c_nf(i1,1)=(a_l*(alphab_liq(k)*rhob_liq(k)))*vl1+(b_l*ar_liq(ii))*vl_o(ii,1)
            fluxl_c_nf(i1,2)=(a_l*(alphab_liq(k)*rhob_liq(k)))*vl2+(b_l*ar_liq(ii))*vl_o(ii,2)
            fluxg_c_nf(i1,1)=(a_g*(alphab_gas(k)*rhob_gas(k)))*vg1+(b_g*ar_gas(ii))*vg_o(ii,1)
            fluxg_c_nf(i1,2)=(a_g*(alphab_gas(k)*rhob_gas(k)))*vg2+(b_g*ar_gas(ii))*vg_o(ii,2)
            fluxd_c_nf(i1,1)=(a_d*(alphab_drp(k)*rhob_drp(k)))*vd1+(b_d*ar_drp(ii))*vd_o(ii,1)
            fluxd_c_nf(i1,2)=(a_d*(alphab_drp(k)*rhob_drp(k)))*vd2+(b_d*ar_drp(ii))*vd_o(ii,2)
         ENDDO
!
!.....Cells out
!
         nf_number=3
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len  
            i1=istart+i
            ii=left_nf(i1)
            mfluxl_cup_nf(i1)=flux_l_nf(i1)*ar_liq(ii)
            mfluxg_cup_nf(i1)=flux_g_nf(i1)*ar_gas(ii)
            mfluxd_cup_nf(i1)=flux_d_nf(i1)*ar_drp(ii)
            fluxl_c_nf(i1,1)=(flux_l_nf(i1)*ar_liq(ii))*vl_o(ii,1)
            fluxg_c_nf(i1,1)=(flux_g_nf(i1)*ar_gas(ii))*vg_o(ii,1)
            fluxd_c_nf(i1,1)=(flux_d_nf(i1)*ar_drp(ii))*vd_o(ii,1)
            fluxl_c_nf(i1,2)=(flux_l_nf(i1)*ar_liq(ii))*vl_o(ii,2)
            fluxg_c_nf(i1,2)=(flux_g_nf(i1)*ar_gas(ii))*vg_o(ii,2)
            fluxd_c_nf(i1,2)=(flux_d_nf(i1)*ar_drp(ii))*vd_o(ii,2)
         ENDDO
      ELSE
!
!.....Cells non
!
         nf_number=0
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len  
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)           
            a_l=min(flux_l_nf(i1),0.d0)
            b_l=max(flux_l_nf(i1),0.d0)
            a_g=min(flux_g_nf(i1),0.d0)
            b_g=max(flux_g_nf(i1),0.d0)
            a_d=min(flux_d_nf(i1),0.d0)
            b_d=max(flux_d_nf(i1),0.d0)
            mfluxl_cup_nf(i1)=a_l*ar_liq(kk)+b_l*ar_liq(ii)
            mfluxg_cup_nf(i1)=a_g*ar_gas(kk)+b_g*ar_gas(ii)
            mfluxd_cup_nf(i1)=a_d*ar_drp(kk)+b_d*ar_drp(ii)
            fluxl_c_nf(i1,1)=(a_l*ar_liq(kk))*vl_o(kk,1)+(b_l*ar_liq(ii))*vl_o(ii,1)
            fluxg_c_nf(i1,1)=(a_g*ar_gas(kk))*vg_o(kk,1)+(b_g*ar_gas(ii))*vg_o(ii,1)
            fluxd_c_nf(i1,1)=(a_d*ar_drp(kk))*vd_o(kk,1)+(b_d*ar_drp(ii))*vd_o(ii,1)
            fluxl_c_nf(i1,2)=(a_l*ar_liq(kk))*vl_o(kk,2)+(b_l*ar_liq(ii))*vl_o(ii,2)
            fluxg_c_nf(i1,2)=(a_g*ar_gas(kk))*vg_o(kk,2)+(b_g*ar_gas(ii))*vg_o(ii,2)
            fluxd_c_nf(i1,2)=(a_d*ar_drp(kk))*vd_o(kk,2)+(b_d*ar_drp(ii))*vd_o(ii,2)
            fluxl_c_nf(i1,3)=(a_l*ar_liq(kk))*vl_o(kk,3)+(b_l*ar_liq(ii))*vl_o(ii,3)
            fluxg_c_nf(i1,3)=(a_g*ar_gas(kk))*vg_o(kk,3)+(b_g*ar_gas(ii))*vg_o(ii,3)
            fluxd_c_nf(i1,3)=(a_d*ar_drp(kk))*vd_o(kk,3)+(b_d*ar_drp(ii))*vd_o(ii,3)
         ENDDO
!         
!........fluxBC: choke model, mcp model, valve model
!          
         IF(rv_valve.eq.1.or.rv_choke.eq.1.or.rv_mcp.eq.1) CALL fluxBC_convection_ice(fluxl_c_nf,fluxg_c_nf,fluxd_c_nf,mfluxl_cup_nf,mfluxg_cup_nf,mfluxd_cup_nf)            
!
!.....Cells mcc
!
         nf_number=1
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len  
            i1=istart+i
            ii=left_nf(i1)
            idx=i3invtbl(i)
            a_l=min(flux_l_nf(i1),0.d0)
            b_l=max(flux_l_nf(i1),0.d0)
            a_g=min(flux_g_nf(i1),0.d0)
            b_g=max(flux_g_nf(i1),0.d0)
            a_d=min(flux_d_nf(i1),0.d0)
            b_d=max(flux_d_nf(i1),0.d0)
            IF(mcdirect(idx).lt.0)THEN
               mfluxl_cup_nf(i1)=a_l*c3dpv(idx,5)+b_l*ar_liq(ii)
               mfluxd_cup_nf(i1)=a_d*c3dpv(idx,4)+b_d*ar_drp(ii)
               vl1=c3vl(1,idx)*xn_nf(i1,1)
               vl2=c3vl(1,idx)*xn_nf(i1,2)
               vl3=c3vl(1,idx)*xn_nf(i1,3)
               vd1=c3vl(1,idx)*xn_nf(i1,1)
               vd2=c3vl(1,idx)*xn_nf(i1,2)
               vd3=c3vl(1,idx)*xn_nf(i1,3)
               fluxl_c_nf(i1,1)=(a_l*c3dpv(idx,5))*vl1+(b_l*ar_liq(ii))*vl_o(ii,1)
               fluxl_c_nf(i1,2)=(a_l*c3dpv(idx,5))*vl2+(b_l*ar_liq(ii))*vl_o(ii,2)
               fluxl_c_nf(i1,3)=(a_l*c3dpv(idx,5))*vl3+(b_l*ar_liq(ii))*vl_o(ii,3)
               fluxd_c_nf(i1,1)=(a_d*c3dpv(idx,4))*vd1+(b_d*ar_drp(ii))*vd_o(ii,1)
               fluxd_c_nf(i1,2)=(a_d*c3dpv(idx,4))*vd2+(b_d*ar_drp(ii))*vd_o(ii,2)
               fluxd_c_nf(i1,3)=(a_d*c3dpv(idx,4))*vd3+(b_d*ar_drp(ii))*vd_o(ii,3)
            ELSE
               mfluxl_cup_nf(i1)=flux_l_nf(i1)*ar_liq(ii)
               mfluxd_cup_nf(i1)=flux_d_nf(i1)*ar_drp(ii)
               fluxl_c_nf(i1,1)=(flux_l_nf(i1)*ar_liq(ii))*vl_o(ii,1)
               fluxl_c_nf(i1,2)=(flux_l_nf(i1)*ar_liq(ii))*vl_o(ii,2)
               fluxl_c_nf(i1,3)=(flux_l_nf(i1)*ar_liq(ii))*vl_o(ii,3)
               fluxd_c_nf(i1,1)=(flux_d_nf(i1)*ar_drp(ii))*vd_o(ii,1)
               fluxd_c_nf(i1,2)=(flux_d_nf(i1)*ar_drp(ii))*vd_o(ii,2)
               fluxd_c_nf(i1,3)=(flux_d_nf(i1)*ar_drp(ii))*vd_o(ii,3)
            ENDIF
            IF(mcgdirect(idx).lt.0)THEN
               mfluxg_cup_nf(i1)=a_g*c3dpv(idx,3)+b_g*ar_gas(ii)
               vg1=c3vg(1,idx)*xn_nf(i1,1)
               vg2=c3vg(1,idx)*xn_nf(i1,2)
               vg3=c3vg(1,idx)*xn_nf(i1,3)
               fluxg_c_nf(i1,1)=(a_g*c3dpv(idx,3))*vg1+(b_g*ar_gas(ii))*vg_o(ii,1)
               fluxg_c_nf(i1,2)=(a_g*c3dpv(idx,3))*vg2+(b_g*ar_gas(ii))*vg_o(ii,2)
               fluxg_c_nf(i1,3)=(a_g*c3dpv(idx,3))*vg3+(b_g*ar_gas(ii))*vg_o(ii,3)
            ELSE
               mfluxg_cup_nf(i1)=flux_g_nf(i1)*ar_gas(ii)
               fluxg_c_nf(i1,1)=(flux_g_nf(i1)*ar_gas(ii))*vg_o(ii,1)
               fluxg_c_nf(i1,2)=(flux_g_nf(i1)*ar_gas(ii))*vg_o(ii,2)
               fluxg_c_nf(i1,3)=(flux_g_nf(i1)*ar_gas(ii))*vg_o(ii,3)
            ENDIF
         ENDDO
!
!.....Cells inl
!
         nf_number=2
         istart=istart_nf(1,nf_number)
         istart2=istart_nbcon_nf(nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len  
            i1=istart+i
            i2=istart2+i
            ii=left_nf(i1)
            k=nbcon_nf(i2)
            a_l=min(flux_l_nf(i1),0.d0)
            b_l=max(flux_l_nf(i1),0.d0)
            a_g=min(flux_g_nf(i1),0.d0)
            b_g=max(flux_g_nf(i1),0.d0)
            a_d=min(flux_d_nf(i1),0.d0)
            b_d=max(flux_d_nf(i1),0.d0)
            mfluxl_cup_nf(i1)=a_l*(alphab_liq(k)*rhob_liq(k))+b_l*ar_liq(ii)
            mfluxg_cup_nf(i1)=a_g*(alphab_gas(k)*rhob_gas(k))+b_g*ar_gas(ii)
            mfluxd_cup_nf(i1)=a_d*(alphab_drp(k)*rhob_drp(k))+b_d*ar_drp(ii)
            f_profile=vel_bc_profile_inl(i)
            IF(vin_norm(k).eq.0)THEN
               vl1=vb_liq(k,1)*f_profile
               vl2=vb_liq(k,2)*f_profile
               vl3=vb_liq(k,3)*f_profile
               vg1=vb_gas(k,1)*f_profile
               vg2=vb_gas(k,2)*f_profile
               vg3=vb_gas(k,3)*f_profile
               vd1=vb_drp(k,1)*f_profile
               vd2=vb_drp(k,2)*f_profile
               vd3=vb_drp(k,3)*f_profile
            ELSE
               vl1=vin_liq(k)*xn_nf(i1,1)*f_profile
               vl2=vin_liq(k)*xn_nf(i1,2)*f_profile
               vl3=vin_liq(k)*xn_nf(i1,3)*f_profile
               vg1=vin_gas(k)*xn_nf(i1,1)*f_profile
               vg2=vin_gas(k)*xn_nf(i1,2)*f_profile
               vg3=vin_gas(k)*xn_nf(i1,3)*f_profile
               vd1=vin_drp(k)*xn_nf(i1,1)*f_profile
               vd2=vin_drp(k)*xn_nf(i1,2)*f_profile
               vd3=vin_drp(k)*xn_nf(i1,3)*f_profile
            ENDIF
            fluxl_c_nf(i1,1)=(a_l*(alphab_liq(k)*rhob_liq(k)))*vl1+(b_l*ar_liq(ii))*vl_o(ii,1)
            fluxl_c_nf(i1,2)=(a_l*(alphab_liq(k)*rhob_liq(k)))*vl2+(b_l*ar_liq(ii))*vl_o(ii,2)
            fluxl_c_nf(i1,3)=(a_l*(alphab_liq(k)*rhob_liq(k)))*vl3+(b_l*ar_liq(ii))*vl_o(ii,3)
            fluxg_c_nf(i1,1)=(a_g*(alphab_gas(k)*rhob_gas(k)))*vg1+(b_g*ar_gas(ii))*vg_o(ii,1)
            fluxg_c_nf(i1,2)=(a_g*(alphab_gas(k)*rhob_gas(k)))*vg2+(b_g*ar_gas(ii))*vg_o(ii,2)
            fluxg_c_nf(i1,3)=(a_g*(alphab_gas(k)*rhob_gas(k)))*vg3+(b_g*ar_gas(ii))*vg_o(ii,3)
            fluxd_c_nf(i1,1)=(a_d*(alphab_drp(k)*rhob_drp(k)))*vd1+(b_d*ar_drp(ii))*vd_o(ii,1)
            fluxd_c_nf(i1,2)=(a_d*(alphab_drp(k)*rhob_drp(k)))*vd2+(b_d*ar_drp(ii))*vd_o(ii,2)
            fluxd_c_nf(i1,3)=(a_l*(alphab_drp(k)*rhob_drp(k)))*vd3+(b_l*ar_drp(ii))*vd_o(ii,3)
         ENDDO
!
!.....Cells out
!
         nf_number=3
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len  
            i1=istart+i
            ii=left_nf(i1)
            mfluxl_cup_nf(i1)=flux_l_nf(i1)*ar_liq(ii)
            mfluxg_cup_nf(i1)=flux_g_nf(i1)*ar_gas(ii)
            mfluxd_cup_nf(i1)=flux_d_nf(i1)*ar_drp(ii)
            fluxl_c_nf(i1,1)=(flux_l_nf(i1)*ar_liq(ii))*vl_o(ii,1)
            fluxg_c_nf(i1,1)=(flux_g_nf(i1)*ar_gas(ii))*vg_o(ii,1)
            fluxd_c_nf(i1,1)=(flux_d_nf(i1)*ar_drp(ii))*vd_o(ii,1)
            fluxl_c_nf(i1,2)=(flux_l_nf(i1)*ar_liq(ii))*vl_o(ii,2)
            fluxg_c_nf(i1,2)=(flux_g_nf(i1)*ar_gas(ii))*vg_o(ii,2)
            fluxd_c_nf(i1,2)=(flux_d_nf(i1)*ar_drp(ii))*vd_o(ii,2)
            fluxl_c_nf(i1,3)=(flux_l_nf(i1)*ar_liq(ii))*vl_o(ii,3)
            fluxg_c_nf(i1,3)=(flux_g_nf(i1)*ar_gas(ii))*vg_o(ii,3)
            fluxd_c_nf(i1,3)=(flux_d_nf(i1)*ar_drp(ii))*vd_o(ii,3)
         ENDDO
      ENDIF
!      
      IF(mom_conv_2nd.gt.0) CALL mom_2nd_conv(fluxl_c_nf,fluxg_c_nf,fluxd_c_nf)
!
!...........Define momentum convection: cnvc_liq,cnvc_gas,cnvc_drp
!
!
!.....Build summation info for all non,mcc,inl,out
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
      nv=-1
      nf_number=nf_number_id(nv)
      len   =istart_nf(2,nf_number)
      DO ix=1,ndim
         DO i=1,len
            k=right_nb_k(i)
            ii=right_non(k)
            cnvc_l_nf(i,ix)=fluxl_c_nf(k,ix)-vl_o(ii,ix)*mfluxl_cup_nf(k)
            cnvc_g_nf(i,ix)=fluxg_c_nf(k,ix)-vg_o(ii,ix)*mfluxg_cup_nf(k)
            cnvc_d_nf(i,ix)=fluxd_c_nf(k,ix)-vd_o(ii,ix)*mfluxd_cup_nf(k)
         ENDDO
      ENDDO
!
      DO nv=0,3
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(nv)
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO ix=1,ndim
            DO i=1,len
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               cnvc_l_nf(i0,ix)=-fluxl_c_nf(i1,ix)+vl_o(ii,ix)*mfluxl_cup_nf(i1)
               cnvc_g_nf(i0,ix)=-fluxg_c_nf(i1,ix)+vg_o(ii,ix)*mfluxg_cup_nf(i1)
               cnvc_d_nf(i0,ix)=-fluxd_c_nf(i1,ix)+vd_o(ii,ix)*mfluxd_cup_nf(i1)
            ENDDO
         ENDDO
      ENDDO
!
      CALL sum_nf_ndim(1,0,ncell_fluid,    &
                       cnvc_l_nf,cnvc_liq, &
                       cnvc_g_nf,cnvc_gas, &
                       cnvc_d_nf,cnvc_drp)
!      
      END SUBROUTINE momentum_convection
