!
      SUBROUTINE set_dt
!
!.....Set a current time step given a CFL number
!
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: np
      USE Zparam       , ONLY: ndim      
!!!   USE Zb_condition , ONLY: vb_gas
      USE Zbicg        , ONLY: pbcgsig,pbcgsig_max      
      USE Zconst2      , ONLY: dt,dt_old,dtr
      USE Ztimecon     , ONLY: ctrl_opt,cfl_ratio,cfl_ratio_max,dt_opt,dt_max,dxmin
      USE Zvector      , ONLY: vg_n,vl_n,vd_n
!
      IMPLICIT NONE
!
      INTEGER i
!
      REAL(8) vmax_liq,vmax_gas,vmax_drp,vel_liq,vel_drp,vel_gas,vmax
      REAL(8) dxmin_vmax 
!
!.....Calculate maximum fluid velocity
!
      vmax_liq=0.0d0
      vmax_drp=0.0d0
      vmax_gas=0.0d0
      IF(ndim.eq.2) THEN
         DO i=1,ncell_fluid
            vel_liq=vl_n(i,1)**2+vl_n(i,2)**2
            vel_drp=vd_n(i,1)**2+vd_n(i,2)**2
            vel_gas=vg_n(i,1)**2+vg_n(i,2)**2
            vmax_liq=max(vmax_liq,vel_liq)
            vmax_drp=max(vmax_drp,vel_drp)
            vmax_gas=max(vmax_gas,vel_gas)
         ENDDO
      ELSE
         DO i=1,ncell_fluid
            vel_liq=vl_n(i,1)**2+vl_n(i,2)**2+vl_n(i,3)**2
            vel_drp=vd_n(i,1)**2+vd_n(i,2)**2+vd_n(i,3)**2
            vel_gas=vg_n(i,1)**2+vg_n(i,2)**2+vg_n(i,3)**2
            vmax_liq=max(vmax_liq,vel_liq)
            vmax_drp=max(vmax_drp,vel_drp)
            vmax_gas=max(vmax_gas,vel_gas)
         ENDDO
      ENDIF
!      
      vmax=max(vmax_liq,vmax_drp,vmax_gas)
      vmax=dsqrt(vmax)
!!!      DO ix=1,ndim
!!!         vmax=dmax1(vmax,dabs(vb_gas(ix,1)))
!!!      ENDDO
!      
      IF(np.gt.1) CALL allreducei_max_r1(vmax)
!
      IF(vmax.gt.1.0d-3)THEN
         dxmin_vmax=dxmin/vmax
      ELSE
         dxmin_vmax=0.001d0
      ENDIF      
!
!.....Set dt depending on dt_opt
!
      IF(ctrl_opt.eq.1) cfl_ratio=cfl_ratio_max
      IF(dt_opt.eq.0) dt_opt=1
      IF(dt_opt.eq.1.or.dt_opt.eq.3.or.dt_opt.eq.4) THEN
         dt_old=dt
         dt=min(cfl_ratio*dxmin_vmax,1.2d0*dt_old)
         dt=min(dt,dt_max)
         dtr=1.d0/dt
!         
      ELSEIF(dt_opt.eq.2) THEN
         dt_old=dt
         dt=min(cfl_ratio*dxmin_vmax,1.2d0*dt_old)
         dt=min(dt,dt_max)
         dtr=1.d0/dt
         IF(dt.lt.1.d-12) STOP '### minimum dt !'
      ELSE
         dt_opt = 1
      ENDIF
!
      IF(pbcgsig.gt.0.and.pbcgsig.le.pbcgsig_max)THEN
         pbcgsig=pbcgsig+1
         dt=min(dt,dt_old)
         dtr=1.0d0/dt
      ELSE
         pbcgsig=0
      ENDIF
!
      END SUBROUTINE set_dt
