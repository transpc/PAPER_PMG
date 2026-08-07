!
      SUBROUTINE set_inlet_flux
!
!     This routine calculate the inlet and outlet mass flow
!     to check the system mass error
!
      USE Zparam       , ONLY: ndim
      USE Zb_condition , ONLY: vb_liq,vb_drp,vb_gas,vin_liq,vin_drp,vin_gas
      USE c3com_cupid  , ONLY: i3invtbl
      USE Zbc_index    , ONLY: vin_norm
!
      USE Znum_cell    , ONLY: istart_nf,istart_nbcon_nf
      USE Zvec_index   , ONLY: left_nf,nbcon_nf
      USE Zvec_major   , ONLY: flux_l_nf,flux_g_nf,flux_d_nf
      USE Zvec_geo     , ONLY: xn_nf,sv_nf,svp_nf
!
      IMPLICIT NONE
!      
      INCLUDE '../10_LinkToMARS/c3com.h' 
!
      INTEGER :: i,k,idx
      INTEGER :: ii
      INTEGER :: nf_number,istart,isize,istart2,i1,i2
!      
      REAL(8) vli,vgi,vdi
!      
      IF(ndim.eq.2)THEN
!
         nf_number=2
         istart=istart_nf(1,nf_number)
         istart2=istart_nbcon_nf(nf_number)
         isize =istart_nf(2,nf_number)
         DO i=1,isize
            i1=istart+i
            i2=istart2+i
            ii=left_nf(i1)
            k=nbcon_nf(i2)
            IF(vin_norm(k).eq.0)THEN
               vli=vb_liq(k,1)*svp_nf(i1,1)+vb_liq(k,2)*svp_nf(i1,2)
               vgi=vb_gas(k,1)*svp_nf(i1,1)+vb_gas(k,2)*svp_nf(i1,2)
               vdi=vb_drp(k,1)*svp_nf(i1,1)+vb_drp(k,2)*svp_nf(i1,2)
            ELSE
               vli=vin_liq(k)*xn_nf(i1,1)*svp_nf(i1,1)+vin_liq(k)*xn_nf(i1,2)*svp_nf(i1,2)
               vgi=vin_gas(k)*xn_nf(i1,1)*svp_nf(i1,1)+vin_gas(k)*xn_nf(i1,2)*svp_nf(i1,2)
               vdi=vin_drp(k)*xn_nf(i1,1)*svp_nf(i1,1)+vin_drp(k)*xn_nf(i1,2)*svp_nf(i1,2)
            ENDIF 
            flux_l_nf(i1)=vli
            flux_g_nf(i1)=vgi
            flux_d_nf(i1)=vdi
         ENDDO
!MARS interface
         nf_number=1
         istart=istart_nf(1,nf_number)
         isize =istart_nf(2,nf_number)
         DO i=1,isize
            i1=istart+i
            ii=left_nf(i1)
            idx=i3invtbl(i)
            vli=0.0d0
            vgi=0.0d0
            vdi=0.0d0
            IF(sv_nf(i1,1).ne.0.0d0)THEN
               vli=c3vl(1,idx)*xn_nf(i1,1)*svp_nf(i1,1)
               vgi=c3vg(1,idx)*xn_nf(i1,1)*svp_nf(i1,1)
               vdi=c3vl(1,idx)*xn_nf(i1,1)*svp_nf(i1,1)
            ENDIF
            IF(sv_nf(i1,2).ne.0.0d0)THEN
               vli=vli+c3vl(1,idx)*xn_nf(i1,2)*svp_nf(i1,2)
               vgi=vgi+c3vg(1,idx)*xn_nf(i1,2)*svp_nf(i1,2)
               vdi=vdi+c3vl(1,idx)*xn_nf(i1,2)*svp_nf(i1,2)
            ENDIF
            flux_l_nf(i1)=vli
            flux_g_nf(i1)=vgi
            flux_d_nf(i1)=vdi
         ENDDO
      ELSE
!
         nf_number=2
         istart=istart_nf(1,nf_number)
         istart2=istart_nbcon_nf(nf_number)
         isize =istart_nf(2,nf_number)
         DO i=1,isize
            i1=istart+i
            i2=istart2+i
            ii=left_nf(i1)
            k=nbcon_nf(i2)
            IF(vin_norm(k).eq.0)THEN
               vli=vb_liq(k,1)*svp_nf(i1,1)+vb_liq(k,2)*svp_nf(i1,2)+vb_liq(k,3)*svp_nf(i1,3)
               vgi=vb_gas(k,1)*svp_nf(i1,1)+vb_gas(k,2)*svp_nf(i1,2)+vb_gas(k,3)*svp_nf(i1,3)
               vdi=vb_drp(k,1)*svp_nf(i1,1)+vb_drp(k,2)*svp_nf(i1,2)+vb_drp(k,3)*svp_nf(i1,3)
            ELSE
               vli=vin_liq(k)*xn_nf(i1,1)*svp_nf(i1,1)+vin_liq(k)*xn_nf(i1,2)*svp_nf(i1,2)+vin_liq(k)*xn_nf(i1,3)*svp_nf(i1,3)
               vgi=vin_gas(k)*xn_nf(i1,1)*svp_nf(i1,1)+vin_gas(k)*xn_nf(i1,2)*svp_nf(i1,2)+vin_gas(k)*xn_nf(i1,3)*svp_nf(i1,3)
               vdi=vin_drp(k)*xn_nf(i1,1)*svp_nf(i1,1)+vin_drp(k)*xn_nf(i1,2)*svp_nf(i1,2)+vin_drp(k)*xn_nf(i1,3)*svp_nf(i1,3)
            ENDIF 
            flux_l_nf(i1)=vli
            flux_g_nf(i1)=vgi
            flux_d_nf(i1)=vdi
         ENDDO
!
!MARS interface
         nf_number=1
         istart=istart_nf(1,nf_number)
         isize =istart_nf(2,nf_number)
         DO i=1,isize
            i1=istart+i
            ii=left_nf(i1)
            idx=i3invtbl(i)
            vli=0.0d0
            vgi=0.0d0
            vdi=0.0d0
            IF(sv_nf(i1,1).ne.0.0d0)THEN
               vli=c3vl(1,idx)*xn_nf(i1,1)*svp_nf(i1,1)
               vgi=c3vg(1,idx)*xn_nf(i1,1)*svp_nf(i1,1)
               vdi=c3vl(1,idx)*xn_nf(i1,1)*svp_nf(i1,1)
            ENDIF
            IF(sv_nf(i1,2).ne.0.0d0)THEN
               vli=vli+c3vl(1,idx)*xn_nf(i1,2)*svp_nf(i1,2)
               vgi=vgi+c3vg(1,idx)*xn_nf(i1,2)*svp_nf(i1,2)
               vdi=vdi+c3vl(1,idx)*xn_nf(i1,2)*svp_nf(i1,2)
            ENDIF
            IF(sv_nf(i1,3).ne.0.0d0)THEN
               vli=vli+c3vl(1,idx)*xn_nf(i1,3)*svp_nf(i1,3)
               vgi=vgi+c3vg(1,idx)*xn_nf(i1,3)*svp_nf(i1,3)
               vdi=vdi+c3vl(1,idx)*xn_nf(i1,3)*svp_nf(i1,3)
            ENDIF
            flux_l_nf(i1)=vli
            flux_g_nf(i1)=vgi
            flux_d_nf(i1)=vdi
         ENDDO
      ENDIF
!
      RETURN
      END SUBROUTINE set_inlet_flux
