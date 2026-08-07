!
      SUBROUTINE turb_ke_vis_gas
!
!     This routine calculates liquid and gas turbulent viscosity 
!     using the solutions of k-e equations
!
      USE VOL_DATA              
      USE Zparam    , ONLY: cmu,ke_small,cappa,clog
      USE Zbc_index , ONLY: icell_type
      USE Zconst1   , ONLY: iturb,lowreynolds,turb_phase
      USE Zface     , ONLY: Kepsilon,Kepsilon_real,gas_only
      USE Zturb     , ONLY: turb_keg,turb_dpg,tauwg,veltg,yplusg, &
                             walln,utaug,wvis_liq,wvis_gas,ustar_keg,w_real_keg,cmug_real
      USE Zzone     , ONLY: ncell_fluid
!
      IMPLICIT NONE
!
      INTEGER i
!      
      REAL(8) fmu
      REAL(8) tv_max,tv_mult_max
      REAL(8) pi_Real,a_s
      REAL(8) Ret,a_square
!
      DATA tv_mult_max/2.0d3/
!
!.....Calculation Cmu for Realizable k-e
!
      IF(iturb.eq.Kepsilon_real)THEN
         DO i=1,ncell_fluid
            pi_Real=1.d0/3.d0*DACOS(DMAX1(DMIN1(DSQRT(6.d0)*w_real_keg(i),1.0d0),-1.0d0))
            a_s=DSQRT(6.d0)*DCOS(pi_real)
            cmug_real(i)=1.0d0/(4.04d0+a_s*ustar_keg(i)*turb_keg(i)/turb_dpg(i))          
         ENDDO
      ENDIF        
!      
!.....Wall sheer stress and viscosity
!
      DO i=1,ncell_fluid
!          
         IF(icell_type(i).eq.1)THEN
!
!...........Single phase wall sheer and viscosity
!
            IF(iturb.ne.Kepsilon_real)THEN      !std k-e & RNG k-e model
               tauwg(i)=cell%rhog(i)*cmu**0.25d0*DSQRT(turb_keg(i))*cappa*veltg(i)/DLOG(clog*yplusg(i))
               wvis_gas(i)=cell%rhog(i)*cmu**0.25d0*DSQRT(turb_keg(i))*cappa/DLOG(clog*yplusg(i))*walln(i)
            ELSE                      !Realizable k-e model
               tauwg(i)=cell%rhog(i)*cmug_real(i)**0.25d0*DSQRT(turb_keg(i))*cappa*veltg(i)/DLOG(clog*yplusg(i))
               wvis_gas(i)=cell%rhog(i)*cmug_real(i)**0.25d0*DSQRT(turb_keg(i))*cappa/DLOG(clog*yplusg(i))*walln(i)         
            ENDIF
            utaug(i)=DSQRT(DABS(tauwg(i))/cell%rhog(i))
            wvis_gas(i)=DMAX1(cell%lviscosg(i),wvis_gas(i))
            IF(yplusg(i).le.11.63d0) wvis_gas(i)=cell%lviscosg(i)
!
!...........LIquid viscosity
!
            IF(turb_phase.eq.gas_only)THEN
               wvis_liq(i)=wvis_gas(i)*cell%rhol(i)/cell%rhog(i)
               wvis_liq(i)=DMAX1(cell%lviscosl(i),wvis_liq(i))
            ENDIF
!
         ENDIF       
!         
      ENDDO
!
!.....Effective viscosity and conductivity
!
      DO i=1,ncell_fluid
!
!........Turbulence viscosity
!
         IF(turb_dpg(i).gt.ke_small) THEN
            IF(iturb.ne.Kepsilon_real)THEN       !std k-e & RNG k-e model
               cell%tviscosg(i)=cell%rhog(i)*cmu*turb_keg(i)*turb_keg(i)/turb_dpg(i)
            ELSE                       !Realizable k-e model
               cell%tviscosg(i)=cell%rhog(i)*cmug_real(i)*turb_keg(i)*turb_keg(i)/turb_dpg(i)
            ENDIF            
            IF(lowreynolds.eq.1) THEN
               fmu=1.d0-DEXP(-0.0115d0*yplusg(i))    
               cell%tviscosg(i)=fmu*cell%tviscosg(i) 
            ELSEIF(lowreynolds.eq.2)THEN
               Ret=turb_keg(i)*turb_keg(i)*cell%rhog(i)/cell%lviscosg(i)/turb_dpg(i)
               a_square=(1.0d0+Ret/50.0d0)*(1.0d0+Ret/50.0d0)
               fmu=DEXP(-3.4d0/a_square)
               cell%tviscosg(i)=fmu*cell%tviscosg(i)
            ENDIF             
         ELSE
            cell%tviscosg(i)=0.0d0
         ENDIF
!
!........Control the maximum turbulent viscosity by tv_max
!
         tv_max=cell%lviscosg(i)*tv_mult_max
         wvis_gas(i)=DMIN1(wvis_gas(i),tv_max)
         cell%tviscosg(i)=DMIN1(cell%tviscosg(i),tv_max)
!         
!........Effective viscosity
!         
         cell%eviscosg(i)=cell%lviscosg(i)+cell%tviscosg(i)
         IF(icell_type(i).eq.1.and.lowreynolds.eq.0) THEN
            IF(yplusg(i).le.11.63d0) THEN
               cell%tviscosg(i)=0.d0
            ELSE
               cell%tviscosg(i)=wvis_gas(i)-cell%lviscosg(i)
            ENDIF   
            cell%eviscosg(i)=wvis_gas(i)
         ENDIF 
!
!........Liquid viscosity
!
         IF(turb_phase.eq.gas_only) THEN
            cell%tviscosl(i)=cell%tviscosg(i)*cell%rhol(i)/cell%rhog(i)
            cell%eviscosl(i)=cell%eviscosg(i)*cell%rhol(i)/cell%rhog(i)
         ENDIF
!
      ENDDO
!
      RETURN
      END SUBROUTINE turb_ke_vis_gas
