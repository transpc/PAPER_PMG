!
      SUBROUTINE udfn_vel_bc_ramp
!
!     This routine defines velocity boundary condition
!
      USE Zparam          , ONLY: ndim
      USE Zconst2         , ONLY: stime_vup,stime_vflat,vl_origin,vg_origin,vd_origin
      USE Zb_condition    , ONLY: vb_liq,vb_gas,vb_drp,vin_liq,vin_gas,vin_drp
      USE Zbc_index       , ONLY: nvin,vin_norm      
      USE Ztimecon        , ONLY: time
!
      IMPLICIT NONE
!
      INTEGER i,ix
!
      LOGICAL, SAVE::initial
!
      DATA INITIAL /.TRUE./ 
!      
      IF(initial)THEN
         initial=.FALSE.
         DO ix=1,ndim
            DO i=1,nvin
               vl_origin(i,ix)=vb_liq(i,ix)
               vg_origin(i,ix)=vb_gas(i,ix)
               vd_origin(i,ix)=vb_drp(i,ix)
            ENDDO
         ENDDO     
      ENDIF
!
!.....Linearly increase the heat flux in time
!
      IF(time.ge.stime_vup.and.time.lt.stime_vflat)THEN
         DO ix=1,ndim
            DO i=1,nvin
               vb_liq(i,ix)=vl_origin(i,ix)*(time-stime_vup)/(stime_vflat-stime_vup)
               vb_gas(i,ix)=vg_origin(i,ix)*(time-stime_vup)/(stime_vflat-stime_vup)
               vb_drp(i,ix)=vd_origin(i,ix)*(time-stime_vup)/(stime_vflat-stime_vup)
            ENDDO
         ENDDO       
      ELSEIF(time.lt.stime_vup)THEN
         DO i=1,nvin
            DO ix=1,ndim
               vb_liq(i,ix)=0.d0
               vb_gas(i,ix)=0.d0
               vb_drp(i,ix)=0.d0
            ENDDO
         ENDDO       
      ELSE
         DO ix=1,ndim
            DO i=1,nvin
               vb_liq(i,ix)=vl_origin(i,ix)
               vb_gas(i,ix)=vg_origin(i,ix)
               vb_drp(i,ix)=vd_origin(i,ix)
            ENDDO
         ENDDO       
      ENDIF
      DO i=1, nvin
         IF(vin_norm(i).gt.0)THEN
            vin_liq(i)=-vb_liq(i,1)
            vin_gas(i)=-vb_gas(i,1)
            vin_drp(i)=-vb_drp(i,1)
         ENDIF      
      ENDDO
!
      RETURN
      END SUBROUTINE udfn_vel_bc_ramp 
