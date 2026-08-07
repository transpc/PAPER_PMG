!
      SUBROUTINE check_scalar(s_flg,repeat)
!
!     This routine check amount of scalar difference between old time step
!     and new time step, and back to old time when the difference is big.
!
      USE Zinterface
      USE VOL_DATA   , ONLY: cell 
      USE Zzone      , ONLY: ncell_fluid
      USE Zcore      , ONLY: np,myrank
      USE Zvec_param , ONLY: nf_nonk,nf_non,nf_mcc,nf_inl,nf_out,nf_fluxk
      USE Znum_cell  , ONLY: right_nb_k,istart_nf, &
                             nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zvec_index , ONLY: left_nf,right_non
      USE Zconst2    , ONLY: dt,dt_old,i_repeat,dtr
      USE Ztimecon   , ONLY: dt_opt,cfl_ratio,cfl_ratio_max,time
      USE Zbc_index  , ONLY: npb
      USE Zcoord3    , ONLY: volr
      USE Zpress     , ONLY: pp
      USE Zare       , ONLY: ar_gas,ar_liq,ar_drp
      USE Zvec_major , ONLY: flux_l_nf,flux_g_nf,flux_d_nf,     &
                             liq_conv_nf,vap_conv_nf,drp_conv_nf
      USE Zio_unit   , ONLY: unit_log
      USE Zcheck_scalar, ONLY:eps_rho,eps_p  
!
      IMPLICIT NONE 
!
!.....Input 
      LOGICAL :: s_flg(3),repeat
!.....Local variables
      INTEGER :: i,k
      INTEGER :: ii,kk
      INTEGER :: nv,nf_number,istart0,istart,len,i0,i1
      LOGICAL :: l_rho,l_dp
      REAL(8) :: drho,dtvol,drho_max,delp
!      REAL(8) :: eps_rho=8.d-3
!      REAL(8) :: eps_p=30.d3
      REAL(8) :: rhos
!.....Local arrays
      REAL(8) :: rho(ncell_fluid)
!.....Local vector arrays
      REAL(8) :: rho_nf(nf_fluxk)
!
      drho_max=0.0d0
      l_rho=.FALSE.
      l_dp=.FALSE.
!      
      IF(dt_opt.eq.2.or.dt_opt.eq.3)THEN
         repeat=s_flg(1).or.s_flg(2).or.s_flg(3)
      ELSEIF(dt_opt.eq.4.or.dt_opt.eq.5)THEN
         repeat=s_flg(1).or.s_flg(2)
      ELSEIF(dt_opt.eq.6.or.dt_opt.eq.7)THEN
         repeat=s_flg(3)
      ENDIF
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
         IF(npb(ii).eq.0 ) THEN
            rho_nf(i0)= vap_conv_nf(i1)*flux_g_nf(i1)  &
                       +liq_conv_nf(i1)*flux_l_nf(i1)  &
                       +drp_conv_nf(i1)*flux_d_nf(i1)
         ELSE
            rho_nf(i0)=0.d0
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
            rho_nf(i)=-vap_conv_nf(k)*flux_g_nf(k)  &
                      -liq_conv_nf(k)*flux_l_nf(k)  &
                      -drp_conv_nf(k)*flux_d_nf(k)
         ELSE
            rho_nf(i)=0.d0
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
               rho_nf(i0) = vap_conv_nf(i1)*flux_g_nf(i1)  &
                           +liq_conv_nf(i1)*flux_l_nf(i1)  &
                           +drp_conv_nf(i1)*flux_d_nf(i1)
            ELSE
               rho_nf(i0)=0.d0
            ENDIF
         ENDDO
      ENDDO
!
      CALL sum_nf(0,0,        &
                  rho_nf,rho)
!
!DIR$ NOVECTOR
      DO i=1,ncell_fluid
         IF(npb(i).gt.0) cycle
         dtvol=dt*volr(i)
         rhos=(ar_gas(i)+ar_liq(i)+ar_drp(i))-rho(i)*dtvol
         drho=dabs(rhos-cell%rhom(i))
         drho_max=MAX(drho_max,drho/cell%rhom(i))
         delp=ABS(pp(i))
         IF(delp .gt. eps_p) l_dp=.TRUE.
      ENDDO
!
      l_rho=drho_max.gt.eps_rho
!
      IF(dt_opt.eq.2.or.dt_opt.eq.5.or.dt_opt.eq.6)THEN
         repeat=repeat.or.l_rho.or.l_dp
      ELSEIF(dt_opt.eq.8)THEN
         repeat=l_rho.or.l_dp
      ELSEIF(dt_opt.eq.9)THEN
         repeat=l_dp
      ENDIF
!
      IF(np.gt.1)THEN
         CALL allreducei_l(s_flg,3)
         CALL allreducei_l1(l_rho)
         CALL allreducei_l1(l_dp)
         CALL allreducei_l1(repeat)
      ENDIF
!
      IF(repeat)THEN
         cfl_ratio=0.5d0*cfl_ratio
         i_repeat=0
      ELSE
         i_repeat=i_repeat+1
         IF(i_repeat.gt.20)THEN
            cfl_ratio=2.d0*cfl_ratio
            cfl_ratio=DMIN1(cfl_ratio,cfl_ratio_max)
         ENDIF
      ENDIF
!
      IF(repeat)THEN
         CALL scalar_reset
         time=time-dt
         dt=dt_old*0.5d0
         dtr=1.0d0/dt
         CALL set_dt
         time=time+dt
         IF(myrank.eq.0) WRITE(*,100) time,dt,s_flg(1),s_flg(2),s_flg(3),l_rho,l_dp
         IF(myrank.eq.0) WRITE(unit_log,100) time,dt,s_flg(1),s_flg(2),s_flg(3),l_rho,l_dp
      ENDIF
!
      IF(dt.lt.1.0d-15)THEN
         IF(myrank.eq.0)WRITE(*,*)'dt is less than 1.0d-15!'
         STOP
      ENDIF
!
  100 FORMAT(e11.4,'dt is reduced to ',e11.4,' : l_el,l_eg,l_vol,l_rho,l_dp = ',5(l1,1x))
!
      END SUBROUTINE check_scalar
!
      SUBROUTINE check_iteration
      USE Zbicg    , ONLY: pbcgsig
      USE Zcore    , ONLY: myrank
      USE Zconst2  , ONLY: dt,dt_old,dtr
      USE Ztimecon , ONLY: time   
!      
      IMPLICIT NONE
!      
!      
      pbcgsig=1
      CALL scalar_reset         
      time=time-dt
      dt=dt_old*0.5d0
      dtr=1.0d0/dt
      CALL set_dt
      time=time+dt
      IF(myrank.eq.0) WRITE(*,"(a,2e20.10)")'          reduce dt:',dt,dt_old
!      
      END SUBROUTINE check_iteration
