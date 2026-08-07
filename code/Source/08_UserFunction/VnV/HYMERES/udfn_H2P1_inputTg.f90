!
      SUBROUTINE udfn_H2P1_inputTg
!
      USE Zb_condition , ONLY: p_fb,tb_liq,tb_gas,qualab,eb_liq,eb_gas,   &
                                rhob_liq,rhob_gas,vb_gas,turb_kegb,turb_dpgb,lvisb_gas
      USE Zconst1      , ONLY: vv_prob,iturb,iVisRatio,vis_ratio
      USE Zconst2      , ONLY: grav
      USE Zcore        , ONLY: np
      USE Zncg         , ONLY: tao,cvao_nvin,uao_nvin,dcva_nvin,ra_nvin,qn_nvin
      USE Znum_cell    , ONLY: istart_nf,istart_nbcon_nf
      USE Zparam       , ONLY: ndim,cmu
      USE Ztimecon     , ONLY: time
      USE Zmodel       , ONLY: h2p1_tp1,h2p1_tp2
      USE Zvec_geo     , ONLY: saa_nf
      USE Zvec_index   , ONLY: nbcon_nf
!
      IMPLICIT NONE
!
      INTEGER i,k,ix
      INTEGER nf_number,istart,isize,i1,istart2,i2
      INTEGER n_injectionT
      INTEGER,SAVE:: n_injectionT1
      REAL(8) FACT,vb_gas_size
      REAL(8),SAVE:: area_in
      REAL(8),SAVE,ALLOCATABLE:: injectionT(:,:),inputTp(:,:)
      REAL(8) tmp1,tmp2 
!
      LOGICAL, SAVE::initialT
      DATA initialT/.TRUE./
!
      IF(initialT)THEN
         initialT=.FALSE.
         IF(vv_prob.eq.'h2p1_0'.or.vv_prob.eq.'h2p1_0x'.or.vv_prob.eq.'VD_h2p1_0')THEN       !2000s
            n_injectionT=1044
         ELSEIF(vv_prob.eq.'h2p1_1'.or.vv_prob.eq.'h2p1_1x')THEN   !2400s
            n_injectionT=1206
         ELSEIF(vv_prob.eq.'h2p1_2'.or.vv_prob.eq.'h2p1_2x')THEN   !2400s
            n_injectionT=1214
         ELSEIF(vv_prob.eq.'h2p1_3'.or.vv_prob.eq.'h2p1_3x')THEN   !4500s
            n_injectionT=2296
         ELSEIF(vv_prob.eq.'h2p1_4'.or.vv_prob.eq.'h2p1_4x')THEN   !5000s
            n_injectionT=2502
         ENDIF
         n_injectionT1=n_injectionT-1
         ALLOCATE(injectionT(2,n_injectionT))
!
         IF(vv_prob.eq.'VD_h2p1_0')THEN
            OPEN(870,file='injectionT.dat',status='old',form='unformatted')
            DO i=1,n_injectionT
               READ(870) injectionT(1,i),injectionT(2,i)
            ENDDO
         ELSE
            OPEN(870,file='injectionT.dat',status='old')
            READ(870,*) injectionT(:,:)
         ENDIF
         CLOSE(870)
!
         IF(vv_prob.eq.'h2p1_1'.or.vv_prob.eq.'h2p1_1x'.or. &
            vv_prob.eq.'h2p1_2'.or.vv_prob.eq.'h2p1_2x'.or. &
            vv_prob.eq.'h2p1_3'.or.vv_prob.eq.'h2p1_3x'.or. &
            vv_prob.eq.'h2p1_4'.or.vv_prob.eq.'h2p1_4x')THEN
            ALLOCATE(inputTp(3,n_injectionT))
            OPEN(876,file='inputTp.dat',status='old')
            READ(876,*) inputTp(:,:)
            CLOSE(876)
         ENDIF
!
         IF(vv_prob.eq.'h2p1_0'.or.vv_prob.eq.'h2p1_1'.or.vv_prob.eq.'h2p1_2'.or. &
            vv_prob.eq.'h2p1_3'.or.vv_prob.eq.'h2p1_4'.or.vv_prob.eq.'VD_h2p1_0')THEN
            area_in=0.0d0
            nf_number=2
            istart=istart_nf(1,nf_number)
            istart2=istart_nbcon_nf(nf_number)
            isize =istart_nf(2,nf_number)
            DO i=1,isize
               i1=istart+i
               i2=istart2+i
               k=nbcon_nf(i2)
               IF(k.eq.1) area_in=area_in+saa_nf(i1)
            ENDDO
            IF(np.gt.1) CALL allreducei_r1(area_in)
         ENDIF
      ENDIF
!
      h2p1_tp1=0.0d0
      h2p1_tp2=0.0d0
      IF(time.le.injectionT(1,1))THEN
         tb_gas(1)=injectionT(2,1)+273.16d0
         CALL convert_temp2erg(p_fb(1),tb_liq(1),tb_gas(1),qualab(1),eb_liq(1),eb_gas(1),rhob_liq(1),rhob_gas(1),tmp1,tmp2, &
                         tao,cvao_nvin(1),uao_nvin(1),dcva_nvin(1),ra_nvin(1),qn_nvin(1,:))
         IF(vv_prob.eq.'h2p1_1'.or.vv_prob.eq.'h2p1_1x'.or. &
            vv_prob.eq.'h2p1_2'.or.vv_prob.eq.'h2p1_2x'.or. &
            vv_prob.eq.'h2p1_3'.or.vv_prob.eq.'h2p1_3x'.or. &
            vv_prob.eq.'h2p1_4'.or.vv_prob.eq.'h2p1_4x')THEN
            h2p1_tp1=inputTp(2,1)+273.16d0
            h2p1_tp2=inputTp(3,1)+273.16d0
         ENDIF
      ELSE
         DO i=1,n_injectionT1
            IF(time.gt.injectionT(1,i).and.time.le.injectionT(1,i+1))THEN
               FACT=(time-injectionT(1,i))/(injectionT(1,i+1)-injectionT(1,i))
               tb_gas(1)=injectionT(2,i)+FACT*(injectionT(2,i+1)-injectionT(2,i))
               tb_gas(1)=tb_gas(1)+273.16d0
               CALL convert_temp2erg(p_fb(1),tb_liq(1),tb_gas(1),qualab(1),eb_liq(1),eb_gas(1),rhob_liq(1),rhob_gas(1),tmp1,tmp2, &
                               tao,cvao_nvin(1),uao_nvin(1),dcva_nvin(1),ra_nvin(1),qn_nvin(1,:))
               IF(vv_prob.eq.'h2p1_1'.or.vv_prob.eq.'h2p1_1x'.or. &
                  vv_prob.eq.'h2p1_2'.or.vv_prob.eq.'h2p1_2x'.or. &
                  vv_prob.eq.'h2p1_3'.or.vv_prob.eq.'h2p1_3x'.or. &
                  vv_prob.eq.'h2p1_4'.or.vv_prob.eq.'h2p1_4x')THEN
                  h2p1_tp1=inputTp(2,i)+FACT*(inputTp(2,i+1)-inputTp(2,i))
                  h2p1_tp1=h2p1_tp1+273.16d0
                  h2p1_tp2=inputTp(3,i)+FACT*(inputTp(3,i+1)-inputTp(3,i))
                  h2p1_tp2=h2p1_tp2+273.16d0
               ENDIF
            ENDIF
         ENDDO
      ENDIF
!
      IF(vv_prob.eq.'h2p1_0x'.or.vv_prob.eq.'h2p1_1x'.or.vv_prob.eq.'h2p1_2x'.or. &
         vv_prob.eq.'h2p1_3x'.or.vv_prob.eq.'h2p1_4x')THEN
!         udfl_vel_bc_profile=.TRUE.
!         udfl_ke_bc_profile=.TRUE.
!         udfl_dp_bc_profile=.FALSE.
      ELSEIF(vv_prob.eq.'h2p1_0'.or.vv_prob.eq.'h2p1_1'.or.vv_prob.eq.'h2p1_2'.or. &
             vv_prob.eq.'h2p1_3'.or.vv_prob.eq.'h2p1_4'.or.vv_prob.eq.'VD_h2p1_0')THEN
         DO ix=1,ndim
            vb_gas(1,ix)=-(grav(ix)/9.81d0)*30.0d0/(rhob_gas(1)*area_in*1.0d3)
         ENDDO
!
         vb_gas_size=0.d0
         DO ix=1,ndim
            vb_gas_size=vb_gas_size+(vb_gas(1,ix)**2)
         ENDDO
         turb_kegb(1)=1.5d0 * 0.05d0**2.d0*vb_gas_size 
!
         IF(iVisRatio.eq.0)THEN       !Turbulence Intensity and Length Scale
            IF(iturb.eq.1)THEN
               turb_dpgb(1)=cmu**(-0.25)*turb_kegb(1)**0.5d0/(0.038d0 * 0.01d0)
            ELSEIF(iturb.ge.2)THEN
               turb_dpgb(1)=cmu**0.75*turb_kegb(1)**1.5d0/(0.07d0 *0.01d0)
            ENDIF
         ELSEIF(iVisRatio.eq.1)THEN
            IF(iturb.eq.1)THEN
               turb_dpgb(1)=rhob_gas(1)*turb_kegb(1)/(lvisb_gas(1)*vis_ratio)
            ELSEIF(iturb.ge.2.and.iturb.le.4)THEN
               turb_dpgb(1)=rhob_gas(1)*cmu*turb_kegb(1)**2/(lvisb_gas(1)*vis_ratio)
            ENDIF
         ENDIF
      ENDIF
!
      RETURN
      END SUBROUTINE udfn_H2P1_inputTg
