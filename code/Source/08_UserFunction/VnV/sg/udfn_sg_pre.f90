!
      SUBROUTINE udfn_sg_pre
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell 
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: np
      USE Zvec_param   , ONLY: nf_nonk,nf_non
      USE Znum_cell    , ONLY: istart_nf,istart_nbcon_nf,                    &
                               nf_number_nb,lens,                           &
                               right_nb_k,istart_nfs,nf_number_id,istart_nf
      USE Zvec_index   , ONLY: left_nf,right_non,nbcon_nf
      USE Ztimecon     , ONLY: time
      USE Zb_condition , ONLY: rhob_liq,vb_liq,vb_gas,alphab_liq,                     &
                               alphab_gas,qualab,rhob_gas,eb_gas,tb_gas,eb_liq,tb_liq
      USE Zbc_index    , ONLY: npb,vin_norm
      USE Zvec_major   , ONLY: flux_l_nf,flux_g_nf
      USE Zvec_geo     , ONLY: sa_nf
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,k,ii,kk
      INTEGER :: nv,nf_number,istart,istart0,istart2,len,i0,i1,i2
      LOGICAL, SAVE :: initial=.true.
      REAL(8) :: flux_liq,flux_gas,mflux_liq,mflux_gas,eflux_liq,eflux_gas,tflux_liq,tflux_gas
      REAL(8) :: vlb,vgb,rholb,rhogb,elb,egb,tlb,tgb
      REAL(8) :: fw_flow
      REAL(8),SAVE :: sb1,sb2,sb3,sb4
      REAL(8) :: tmp(8)
!.....Local arrays 
      REAL(8),DIMENSION(ncell_fluid) :: flux_liqv,flux_gasv,mflux_liqv,mflux_gasv,  &
                                        eflux_liqv,eflux_gasv,tflux_liqv,tflux_gasv 
!.....Local vector arrays 
      REAL(8),DIMENSION(nf_nonk+nf_non) :: flux_liqa_nf,mflux_liqa_nf,mflux_gasa_nf,eflux_liqa_nf, &
                                           flux_gasa_nf,eflux_gasa_nf,tflux_liqa_nf,tflux_gasa_nf
!
      IF(initial)THEN
!
!........Cells inl
!
         nf_number=2
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         istart2=istart_nbcon_nf(nf_number)
         sb1=0.d0
         sb2=0.d0
         sb3=0.d0
         sb4=0.d0
         DO i=1,len  
            i1=istart+i
            i2=istart2+i
            ii=left_nf(i1)
            IF    (nbcon_nf(i2).eq.1) THEN
               sb1=sb1+sa_nf(i1)
            ELSEIF(nbcon_nf(i2).eq.2) THEN
               sb2=sb2+sa_nf(i1)
            ELSEIF(nbcon_nf(i2).eq.3) THEN
               sb3=sb3+sa_nf(i1)
            ELSEIF(nbcon_nf(i2).eq.4) THEN
               sb4=sb4+sa_nf(i1)
            ENDIF
         ENDDO
!
         IF(np.gt.1)THEN
            tmp(1)=sb1
            tmp(2)=sb2
            tmp(3)=sb3
            tmp(4)=sb4
            CALL allreducei_r(tmp,4)
            sb1=tmp(1)
            sb2=tmp(2)
            sb3=tmp(3)
            sb4=tmp(4)
         ENDIF
!
         fw_flow=1130.6d0   ! should be an input
         vb_liq(1,1)=0.9d0*fw_flow/rhob_liq(1)/sb1
         vb_liq(2,1)=0.1d0*fw_flow/rhob_liq(2)/sb2
!
         initial=.false.
      ENDIF    
!
!.....Build summation info for non
!
      nf_number_nb=0
      nf_number_id(-1)=-1
      nf_number_id(0)=0
      istart_nfs(-1)=0
      istart_nfs(0)=istart_nfs(-1)+nf_nonk
      lens         =istart_nfs(0) +nf_non
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
         IF((npb(ii).eq.2 .and. npb(kk).eq.0) .and. &
             flux_l_nf(i1).lt.0.d0) THEN
             flux_liqa_nf(i0) =flux_l_nf(i1)*cell%alphal(kk)
             mflux_liqa_nf(i0)=flux_l_nf(i1)*cell%rhol(kk)*cell%alphal(kk)
             mflux_gasa_nf(i0)=flux_g_nf(i1)*cell%rhog(kk)*cell%alphag(kk)
             eflux_liqa_nf(i0)=flux_l_nf(i1)*cell%el(kk)*cell%alphal(kk)
         ELSE
             flux_liqa_nf(i0) =0.d0
             mflux_liqa_nf(i0)=0.d0
             mflux_gasa_nf(i0)=0.d0
             eflux_liqa_nf(i0)=0.d0
         ENDIF 
         IF((npb(ii).eq.2 .and. npb(kk).eq.0) .and. &
             flux_g_nf(i1).lt.0.d0) THEN
             flux_gasa_nf(i0) =flux_g_nf(i1)*cell%alphag(kk)
             eflux_gasa_nf(i0)=flux_g_nf(i1)*cell%eg(kk)*cell%alphag(kk)
             tflux_liqa_nf(i0)=flux_l_nf(i1)*cell%tl(kk)*cell%alphal(kk)
             tflux_gasa_nf(i0)=flux_g_nf(i1)*cell%tg(kk)*cell%alphag(kk)
         ELSE
             flux_gasa_nf(i0) =0.d0
             eflux_gasa_nf(i0)=0.d0
             tflux_liqa_nf(i0)=0.d0
             tflux_gasa_nf(i0)=0.d0
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
         IF((npb(ii).eq.2 .and. npb(kk).eq.0) .and. &
             flux_l_nf(k).gt.0.d0) THEN
             flux_liqa_nf(i) =-flux_l_nf(k)*cell%alphal(kk)
             mflux_liqa_nf(i)=-flux_l_nf(k)*cell%rhol(kk)*cell%alphal(kk)
             mflux_gasa_nf(i)=-flux_g_nf(k)*cell%rhog(kk)*cell%alphag(kk)
             eflux_liqa_nf(i)=-flux_l_nf(k)*cell%el(kk)*cell%alphal(kk)
         ELSE
             flux_liqa_nf(i) =0.d0
             mflux_liqa_nf(i)=0.d0
             mflux_gasa_nf(i)=0.d0
             eflux_liqa_nf(i)=0.d0
         ENDIF 
         IF((npb(ii).eq.2 .and. npb(kk).eq.0) .and. &
             flux_g_nf(k).gt.0.d0) THEN
             flux_gasa_nf(i) =-flux_g_nf(k)*cell%alphag(kk)
             eflux_gasa_nf(i)=-flux_g_nf(k)*cell%eg(kk)*cell%alphag(kk)
             tflux_liqa_nf(i)=-flux_l_nf(k)*cell%tl(kk)*cell%alphal(kk)
             tflux_gasa_nf(i)=-flux_g_nf(k)*cell%tg(kk)*cell%alphag(kk)
         ELSE
             flux_gasa_nf(i) =0.d0
             eflux_gasa_nf(i)=0.d0
             tflux_liqa_nf(i)=0.d0
             tflux_gasa_nf(i)=0.d0
         ENDIF 
      ENDDO
!
      CALL sum_nf(0,0,                      &
                  flux_liqa_nf ,flux_liqv,  &
                  mflux_liqa_nf,mflux_liqv, & 
                  mflux_gasa_nf,mflux_gasv, &
                  eflux_liqa_nf,eflux_liqv, &
                  flux_gasa_nf ,flux_gasv,  &
                  eflux_gasa_nf,eflux_gasv, &
                  tflux_liqa_nf,tflux_liqv, &
                  tflux_gasa_nf,tflux_gasv)
!
      flux_liq =0.d0
      flux_gas =0.d0
      mflux_liq=0.d0
      mflux_gas=0.d0
      eflux_liq=0.d0
      eflux_gas=0.d0
      tflux_liq=0.d0
      tflux_gas=0.d0
      DO i=1,ncell_fluid
         IF(npb(i).eq.2)THEN
            flux_liq =flux_liq +flux_liqv(i)
            mflux_liq=mflux_liq+mflux_liqv(i)
            mflux_gas=mflux_gas+mflux_gasv(i)
            eflux_liq=eflux_liq+eflux_liqv(i)
!
            flux_gas =flux_gas +flux_gasv(i)
            eflux_gas=eflux_gas+eflux_gasv(i)
            tflux_liq=tflux_liq+tflux_liqv(i)
            tflux_gas=tflux_gas+tflux_gasv(i)
         ENDIF
      ENDDO
!
      IF(np.gt.1)THEN
         tmp(1)=flux_liq
         tmp(2)=mflux_liq
         tmp(3)=mflux_gas
         tmp(4)=eflux_liq
         tmp(5)=flux_gas 
         tmp(6)=eflux_gas
         tmp(7)=tflux_liq
         tmp(8)=tflux_gas
         CALL allreducei_r(tmp,8)
         flux_liq =tmp(1)
         mflux_liq=tmp(2)
         mflux_gas=tmp(3)
         eflux_liq=tmp(4)
         flux_gas =tmp(5)
         eflux_gas=tmp(6)
         tflux_liq=tmp(7)
         tflux_gas=tmp(8)
      ENDIF
!
!     IF(.not.initial)THEN
!        IF(time.gt.10.)THEN
      IF(time.gt.2.)THEN
!        IF(time.ge.0.0)THEN
!
         IF(flux_liq.lt.0.d0)THEN
            vlb=-flux_liq/sb3
            rholb=mflux_liq/flux_liq
            elb=eflux_liq/flux_liq
            tlb=tflux_liq/flux_liq
         ELSE
            vlb  =0.d0
            rholb=0.d0
            elb  =0.d0
            tlb  =0.d0
         ENDIF
         IF(flux_gas.lt.0.d0)THEN
            vgb=-flux_gas/sb4
            rhogb=mflux_gas/flux_gas
            egb=eflux_gas/flux_gas
            tgb=tflux_gas/flux_gas
         ELSE
            vgb  =0.d0
            rhogb=0.d0
            egb  =0.d0
            tgb  =0.d0
         ENDIF
!
if(1)then
         vin_norm(3)=0
         alphab_gas(3)=0.d0
         alphab_liq(3)=1.d0 
         qualab(3)=0.d0
         rhob_liq(3)=rholb
         eb_liq(3)=elb
         tb_liq(3)=tlb
         vb_liq(3,3)=-vlb
endif
!
if(1)then
         vin_norm(4)=0
         alphab_gas(4)=1.d0
         alphab_liq(4)=0.d0
         qualab(3)=0.d0
         rhob_gas(4)=rhogb
         eb_gas(4)=egb
         tb_gas(4)=tgb
         vb_gas(4,3)=vgb
endif
      ENDIF
!
      END SUBROUTINE udfn_sg_pre
