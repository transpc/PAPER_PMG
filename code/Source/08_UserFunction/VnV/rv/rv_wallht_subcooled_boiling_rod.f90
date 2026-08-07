!
      SUBROUTINE subcooled_boiling_rod(i,mul_o,qe,p,rhog,rhol,alphal)
!
!     This routine calculates wall heat transfer coefficient and heat
!     flux for single phase convection
!     
      USE VOL_DATA
      USE Zconst1      , ONLY: restart
      USE Zconst2      , ONLY: dt
      USE Zwall_HTC    , ONLY: qflux_l,mflux_tota,dia_rod,cp,cond,sat_hfp,h_liq,hfg_p     
     !rod-scale
      USE Zparam       , ONLY: pi
      USE Zrv_ncell    , ONLY: channel_cell_hts2d,cupid_cell_hts2d         
      USE Zrv_hts_2d   , ONLY: nr_2d,ri_2d
      USE Zcoord2      , ONLY: cell_leng
      USE Zcoord3      , ONLY: volp
      USE Zwall_HTC    , ONLY: gamma_wall_rod
      USE Zio_unit     , ONLY: unit_log

!
      IMPLICIT NONE
!
      INTEGER i
!
      REAL(8) mul_o,ht_area
      REAL(8) f_pump,qe
      REAL(8) St,Pe,Nu,h_min,h_crit,Mul
      REAL(8) pbr,fprg_min,crenhr,p5r,fpreps,gmmllp,gmmlep    !SRL model
!
      REAL(8) tau_t
      ! Additional connectivity - jrlee
      INTEGER k,m
      REAL(8) p
      REAL(8) rhog,rhol
      REAL(8) alphal

      DATA tau_t/0.1d0/      !time constant for Mul relaxation
!
      k=channel_cell_hts2d(i)
      m=cupid_cell_hts2d(i)

      Pe=mflux_tota*dia_rod*cp/cond
      Nu=qflux_l*dia_rod/cond
      h_min=DMIN1(h_liq,sat_hfp)
      IF(Pe.gt.7.0d4)THEN
         St=Nu/Pe
         IF(.true.)THEN       ! SRL model
            pbr=p/6.894757d3            !psia=6.894757d3
            fprg_min=1.0782d0/(1.015d0+DEXP((pbr-145.75d0)/28.d0))         
            h_crit=sat_hfp-St*cp/(0.0055d0-0.0009d0*fprg_min)
         ELSE            ! Default Saha model
            h_crit=sat_hfp-St*cp*153.85d0    !1/0.0065=153.85            
         ENDIF
      ELSE
         h_crit=sat_hfp-Nu*cp*0.0021978      !1/445=0.0021978
      ENDIF

      IF(h_min.gt.h_crit)THEN
         IF(.true.)THEN       ! SRL model
            pbr=p/6.894757d3            !psia=6.894757d3
!            fprg_min=1.0782d0/(1.015d0+DEXP((pbr-145.75d0)/28.d0)) 
            fprg_min=1.0782d0/(1.015d0+DEXP((pbr-140.75d0)/28.d0)) !yjm modification 145.75 -> 140.75
            crenhr=(h_min-h_crit)/(sat_hfp-h_crit)
!            p5r=0.0022d0+crenhr*(0.11262d0+crenhr*(-0.53224d0+crenhr*(8.68227d0+crenhr*(-11.29044d0+crenhr*(4.253448d0)))))
            p5r=0.0022d0+crenhr*(0.11262d0+crenhr*(-0.59224d0+crenhr*(8.68227d0+crenhr*(-11.29044d0+crenhr*(4.253448d0))))) ! yjm mod 0.53224 -> 0.59224
            p5r=DMIN1(1.0d0,p5r)
            fpreps=1.0d0/(0.97d0+38.0d0*DEXP(-(pbr+60.0d0)/42.0d0))
            fpreps=DMIN1(1.0d0,fpreps)
            gmmllp=crenhr+fprg_min*(p5r-crenhr)
            gmmlep=(1.0d0+rhol*(sat_hfp-h_min)*fpreps/(rhog*hfg_p))
            IF(DABS(gmmlep).gt.1.0d-10)THEN
               gmmlep=gmmllp/gmmlep                                              
            ELSE
               gmmlep=1.d10
            ENDIF
            Mul=gmmlep
         ELSE            ! Default Saha model
            Mul=(h_min-h_crit)/(sat_hfp-h_crit)
            f_pump=rhol*(sat_hfp-h_min)/(rhog*hfg_p)
            Mul=Mul/(1.0d0+f_pump)
         ENDIF
         IF(restart.ne.0)Mul_o=Mul
!
!........Relaxation of Mul
!    
         Mul=(Mul_o+dt*Mul/tau_t)/(1.0d0+dt/tau_t)
!
!........Calculate vapor generation rate (kg/s)
!        
!........Origin
         !gamma_wall(i)=gamma_wall(i)+qflux_l*ht_area*volpr(i)/DMAX1(1.0d4,hfg_p)*Mul
!........rod-scale         
         ht_area=2.d0*pi*ri_2d(nr_2d)*cell_leng(m,3)*0.25d0
         gamma_wall_rod(i)=qflux_l*ht_area/DMAX1(1.0d4,hfg_p)*Mul
!
!........Limit vapor generate rate. The rate cannot exceed '90% of of liquid in the cell'
!
!........Origin
         !gamma_wall(i)=DMIN1(gamma_wall(i),0.9d0*cell%alphal(i)*cell%rhol(i)/dt)
!........rod-scale
         gamma_wall_rod(i)=DMIN1(gamma_wall_rod(i),0.9d0*alphal*rhol/dt*volp(m)*0.25d0)
!         
         qe=qflux_l*Mul
      ELSE
         qe=0.0d0
         !gamma_wall(i)=0.0d0
         Mul=Mul_o
      ENDIF
      Mul_o=Mul
!
     
      IF(ISNAN(qe).or.ISNAN(Mul))THEN
         WRITE(*,*)'subcooled_boiling:qe,Mul=',qe,Mul
         WRITE(unit_log,*)'subcooled_boiling:qe,MUl=',qe,Mul
      ENDIF
!
      RETURN
    END SUBROUTINE subcooled_boiling_rod
!
!======================================================================
!======================================================================
!
      SUBROUTINE subcooled_boiling_rod_sahazuber(i,mul_o,qe,p,rhog,rhol,alphal)
!
!     This routine calculates wall heat transfer coefficient and heat
!     flux for single phase convection
!     
      USE VOL_DATA
      USE Zconst1      , ONLY: restart
      USE Zconst2      , ONLY: dt
      USE Zwall_HTC    , ONLY: qflux_l,mflux_tota,dia_rod,cp,cond,sat_hfp,h_liq,hfg_p     
     !rod-scale
      USE Zparam       , ONLY: pi
      USE Zrv_ncell    , ONLY: channel_cell_hts2d,cupid_cell_hts2d         
      USE Zrv_hts_2d   , ONLY: nr_2d,ri_2d
      USE Zcoord2      , ONLY: cell_leng
      USE Zcoord3      , ONLY: volp
      USE Zwall_HTC    , ONLY: gamma_wall_rod
      USE Zio_unit     , ONLY: unit_log

!
      IMPLICIT NONE
!
      INTEGER i
!
      REAL(8) mul_o,ht_area
      REAL(8) f_pump,qe
      REAL(8) St,Pe,Nu,h_min,h_crit,Mul
      REAL(8) pbr,fprg_min,crenhr,p5r,fpreps,gmmllp,gmmlep    !SRL model
!
      REAL(8) tau_t
      ! Additional connectivity - jrlee
      INTEGER k,m
      REAL(8) p
      REAL(8) rhog,rhol
      REAL(8) alphal

      DATA tau_t/0.1d0/      !time constant for Mul relaxation
!
      k=channel_cell_hts2d(i)
      m=cupid_cell_hts2d(i)

      Pe=mflux_tota*dia_rod*cp/cond
      Nu=qflux_l*dia_rod/cond
      h_min=DMIN1(h_liq,sat_hfp)

      IF(Pe.gt.7.0d4)THEN
         St=Nu/Pe
         IF(.true.)THEN       ! SRL model
            pbr=p/6.894757d3            !psia=6.894757d3
            fprg_min=1.0782d0/(1.015d0+DEXP((pbr-145.75d0)/28.d0))         
            h_crit=sat_hfp-St*cp/(0.0055d0-0.0009d0*fprg_min)
         ELSE            ! Default Saha model
            h_crit=sat_hfp-St*cp*153.85d0    !1/0.0065=153.85            
         ENDIF
      ELSE
         h_crit=sat_hfp-Nu*cp*0.0021978      !1/445=0.0021978
      ENDIF

      IF(h_min.gt.h_crit)THEN
         IF(.true.)THEN       ! SRL model
            pbr=p/6.894757d3            !psia=6.894757d3
!            fprg_min=1.0782d0/(1.015d0+DEXP((pbr-145.75d0)/28.d0)) 
            fprg_min=1.0782d0/(1.015d0+DEXP((pbr-140.75d0)/28.d0)) !yjm modification 145.75 -> 140.75
            crenhr=(h_min-h_crit)/(sat_hfp-h_crit)
!            p5r=0.0022d0+crenhr*(0.11262d0+crenhr*(-0.53224d0+crenhr*(8.68227d0+crenhr*(-11.29044d0+crenhr*(4.253448d0)))))
            p5r=0.0022d0+crenhr*(0.11262d0+crenhr*(-0.59224d0+crenhr*(8.68227d0+crenhr*(-11.29044d0+crenhr*(4.253448d0))))) ! yjm mod 0.53224 -> 0.59224
            p5r=DMIN1(1.0d0,p5r)
            fpreps=1.0d0/(0.97d0+38.0d0*DEXP(-(pbr+60.0d0)/42.0d0))
            fpreps=DMIN1(1.0d0,fpreps)
            gmmllp=crenhr+fprg_min*(p5r-crenhr)
            gmmlep=(1.0d0+rhol*(sat_hfp-h_min)*fpreps/(rhog*hfg_p))
            IF(DABS(gmmlep).gt.1.0d-10)THEN
               gmmlep=gmmllp/gmmlep                                              
            ELSE
               gmmlep=1.d10
            ENDIF
            Mul=gmmlep
         ELSE            ! Default Saha model
            Mul=(h_min-h_crit)/(sat_hfp-h_crit)
            f_pump=rhol*(sat_hfp-h_min)/(rhog*hfg_p)
            Mul=Mul/(1.0d0+f_pump)
         ENDIF
         IF(restart.ne.0)Mul_o=Mul
!
!........Relaxation of Mul
!    
         Mul=(Mul_o+dt*Mul/tau_t)/(1.0d0+dt/tau_t)
!
!........Calculate vapor generation rate (kg/s)
!        
!........Origin
         !gamma_wall(i)=gamma_wall(i)+qflux_l*ht_area*volpr(i)/DMAX1(1.0d4,hfg_p)*Mul
!........rod-scale         
         ht_area=2.d0*pi*ri_2d(nr_2d)*cell_leng(m,3)*0.25d0
         gamma_wall_rod(i)=qflux_l*ht_area/DMAX1(1.0d4,hfg_p)*Mul
!
!........Limit vapor generate rate. The rate cannot exceed '90% of of liquid in the cell'
!
!........Origin
         !gamma_wall(i)=DMIN1(gamma_wall(i),0.9d0*cell%alphal(i)*cell%rhol(i)/dt)
!........rod-scale
         gamma_wall_rod(i)=DMIN1(gamma_wall_rod(i),0.9d0*alphal*rhol/dt*volp(m)*0.25d0)
!         
         qe=qflux_l*Mul
      ELSE
         qe=0.0d0
         !gamma_wall(i)=0.0d0
         Mul=Mul_o
      ENDIF
      Mul_o=Mul
!
     
      IF(ISNAN(qe).or.ISNAN(Mul))THEN
         WRITE(*,*)'subcooled_boiling:qe,Mul=',qe,Mul
         WRITE(unit_log,*)'subcooled_boiling:qe,MUl=',qe,Mul
      ENDIF
!
      RETURN
      END SUBROUTINE subcooled_boiling_rod_sahazuber
!
