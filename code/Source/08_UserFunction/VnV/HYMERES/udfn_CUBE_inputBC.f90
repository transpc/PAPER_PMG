!
      SUBROUTINE udfn_CUBE_inputBC
!
      USE Zb_condition , ONLY: p_fb,tb_liq,tb_gas,qualab,eb_liq,eb_gas,   &
                                rhob_liq,rhob_gas,vb_gas,turb_kegb,turb_dpgb,lvisb_gas
      USE Zconst1      , ONLY: vv_prob,iturb,iVisRatio,vis_ratio
      USE Zcore        , ONLY: np
      USE Zncg         , ONLY: tao,cvao_nvin,uao_nvin,dcva_nvin,ra_nvin,qn_nvin
      USE Znum_cell    , ONLY: istart_nf,istart_nbcon_nf
      USE Zparam       , ONLY: ndim,cmu
      USE Ztimecon     , ONLY: time
      USE Zvec_geo     , ONLY: saa_nf
      USE Zvec_index   , ONLY: nbcon_nf
!
      IMPLICIT NONE
!
      INTEGER i,k,ix
      INTEGER nf_number,istart,isize,i1,istart2,i2
      INTEGER n_inputBC
      INTEGER,SAVE:: n_inputBC1
      REAL(8) FACT,vb_gas_size
      REAL(8),SAVE:: area_in
      REAL(8),SAVE,ALLOCATABLE:: inputBC(:,:)
      REAL(8) tmp1,tmp2,inputMFR
      REAL(8) qn_nvin0(8)
!
      LOGICAL, SAVE::initialBC
      DATA initialBC/.TRUE./
!
      IF(initialBC)THEN
         initialBC=.FALSE.
         IF(vv_prob.eq.'ST2-CT-01')THEN
            n_inputBC=8154
         ELSEIF(vv_prob.eq.'ST2-CT-02')THEN
            n_inputBC=8270
         ELSEIF(vv_prob.eq.'ST2-CT-03')THEN
            n_inputBC=11898
         ENDIF
         n_inputBC1=n_inputBC-1
         ALLOCATE(inputBC(4,n_inputBC))
!
         OPEN(870,file='CUBE_BC.dat',status='old')
         READ(870,*) inputBC(:,:)
         CLOSE(870)
!
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
         IF(np.gt.1)THEN
            CALL allreducei_r1(area_in)
         ENDIF
      ENDIF
!
      IF(time.le.inputBC(1,1))THEN
         inputMFR=inputBC(2,1)
         p_fb(1)=inputBC(3,1)*1.0d6
         tb_gas(1)=inputBC(4,1)+273.15d0
         tb_liq(1)=tb_gas(1)
         qn_nvin0(:)=qn_nvin(1,:)
         CALL convert_temp2erg(p_fb(1),tb_liq(1),tb_gas(1),qualab(1),eb_liq(1),eb_gas(1),rhob_liq(1),rhob_gas(1),tmp1,tmp2, &
                         tao,cvao_nvin(1),uao_nvin(1),dcva_nvin(1),ra_nvin(1),qn_nvin0)
      ELSE
         DO i=1,n_inputBC1
            IF(time.gt.inputBC(1,i).and.time.le.inputBC(1,i+1))THEN
               FACT=(time-inputBC(1,i))/(inputBC(1,i+1)-inputBC(1,i))
               inputMFR=inputBC(2,i)+FACT*(inputBC(2,i+1)-inputBC(2,i))
               p_fb(1)=inputBC(3,i)+FACT*(inputBC(3,i+1)-inputBC(3,i))
               p_fb(1)=p_fb(1)*1.0d6
               tb_gas(1)=inputBC(4,i)+FACT*(inputBC(4,i+1)-inputBC(4,i))
               tb_gas(1)=tb_gas(1)+273.15d0
               tb_liq(1)=tb_gas(1)
               qn_nvin0(:)=qn_nvin(1,:)
               CALL convert_temp2erg(p_fb(1),tb_liq(1),tb_gas(1),qualab(1),eb_liq(1),eb_gas(1),rhob_liq(1),rhob_gas(1),tmp1,tmp2, &
                               tao,cvao_nvin(1),uao_nvin(1),dcva_nvin(1),ra_nvin(1),qn_nvin0)
            ENDIF
         ENDDO
      ENDIF
!
      vb_gas(1,1)=0.0d0
      vb_gas(1,2)=-inputMFR/(rhob_gas(1)*area_in)
      vb_gas(1,3)=0.0d0
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
!
      RETURN
      END SUBROUTINE udfn_CUBE_inputBC
