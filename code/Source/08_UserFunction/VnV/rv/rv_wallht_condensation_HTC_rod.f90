!
      SUBROUTINE condensation_HTC_rod(i,mode,ts,quala,hg,hlsat,qc, &
                                      p,alphag,alphal,rhog,rhol,tg,tl,vg,vl,condg,condl,visg,visl,cpg,cpl,betag,betal)
!
!     This routine calculates wall heat transfer coefficient for condensation near wall
!
      USE STM_TBL_cupid  , ONLY: pcrit
      USE Zwall_HTC      , ONLY: HTC_T,HTC_tl,HTC_tst,HTC_tg,HTC_tgp,HTC_tgt,qflux_t,qflux_l,qflux_g, &
                                tw,qual_eq,mflux_tota,dia_rod
     !rod-scale
      USE Zparam         , ONLY:pi
      USE Zrv_ncell      , ONLY: cupid_cell_hts2d         
      USE Zrv_hts_2d     , ONLY: nr_2d,ri_2d
      USE Zcoord2        , ONLY: cell_leng
      USE Zwall_HTC      , ONLY: gamma_wall_rod
      USE Zio_unit       , ONLY: unit_log
!
      IMPLICIT NONE
!
!.....Input
      INTEGER i,mode
      REAL(8) ts,quala,hg,hlsat
      ! arguments for sinlge_HTC_rod
      REAL(8) p
      REAL(8) alphag,rhog,condg,visg,cpg,betag,tg,vg
      REAL(8) alphal,rhol,condl,visl,cpl,betal,tl,vl
!
!.....Output
      REAL(8) qc      
!
!.....local      
      INTEGER m
      REAL(8) dt_pps,dt_liq
      REAL(8) HTC_ditt,HTC_cond
      REAL(8) qual_max,Re,Pr
      REAL(8) hdb,hf,z,ftr
      REAL(8) tdiff
      REAL(8) ht_area
!
!.....CUPID local cell number,m
!      
      m=cupid_cell_hts2d(i)

      dt_pps=tw-ts
      dt_liq=tw-tl
!      
      IF(alphag.le.0.3d0.or.dt_liq.gt.0.0d0)THEN
         mode=2
         CALL single_phase_HTC_rod(i,mode,alphag,alphal,rhog,rhol,tg,tl,vg,vl,condg,condl,visg,visl,cpg,cpl,betag,betal)
         HTC_ditt=HTC_tl
         qflux_g=0.0d0
      ELSE
         HTC_ditt=0.0d0
         qflux_t=0.0d0
      ENDIF
!
!.....Calculate condensation HTC when NC gas is absent
!
      IF(quala.lt.0.001d0)THEN
         qual_max=DMAX1(1.0d-9,DMIN1(1.0d0,qual_eq))
         Re=DMIN1(4000.0d0,mflux_tota)*dia_rod/visl
         Pr=visl*cpl/condl
         hdb=0.023d0*condl*Re**0.8d0*Pr**0.4d0/dia_rod
         hf=hdb*(1.0d0-qual_max)**0.8d0
         z=(p/pcrit)**0.04d0*(1.0d0/qual_max-1.0d0)**0.8d0
         ftr=1.0d0
         IF(z.ne.0.0d0)ftr=1.0d0+3.8d0/z**0.95d0
         HTC_cond=hf*ftr
         !HTC_cond=DMAX1(HTC_cond,hf*ftr)    !If horizontal model is calculated above 
      ELSE
         WRITE(*,*) '        Condensation model with NC gas is not implemented yet!!', i
         PAUSE 
      ENDIF
!
!.....Calculate HTC and qflux
!
      IF(dt_liq.lt.0.0d0)THEN
         qflux_t=HTC_cond*dt_pps
         qflux_l=HTC_cond*dt_liq
         qflux_l=DMAX1(qflux_t,qflux_l)
         qflux_g=qflux_t-qflux_l
         qflux_g=DMIN1(qflux_g,0.0d0)
         HTC_tgp=qflux_g/dt_pps
         HTC_tl=HTC_cond
      ELSE
         HTC_tl=HTC_ditt
         qflux_l=HTC_tl*dt_liq
         HTC_tgp=HTC_cond
         qflux_g=HTC_tgp*dt_pps
      ENDIF      
!
!.....Interpolete HTC&qflux if ag < 0.3
!     
      IF(alphag.lt.0.3d0)THEN
!         HTC_tl=HTC_tl+(HTC_tl-HTC_ditt)*(cell%alphal(i)-0.1d0)/0.2d0
         HTC_tl=HTC_ditt+(HTC_tl-HTC_ditt)*(alphal-0.1d0)/0.2d0   !lsj
!         HTC_tgp=HTC_tgp+HTC_tgp*(cell%alphal(i)-0.1d0)/0.2d0
         HTC_tgp=HTC_tgp*(alphal-0.1d0)/0.2d0                     !lsj
         qflux_l=HTC_tl*dt_liq
         qflux_g=HTC_tgp*dt_pps
      ENDIF
!
!.....Save HTC and heat flux
!     
      mode=11
      qflux_t=qflux_l+qflux_g
      HTC_t=HTC_tl+HTC_tgp
      HTC_tg=0.0d0
      HTC_tgt=0.0d0
      HTC_tst=0.0d0
     !origin
     !gamma_wall(i)=gamma_wall(i)+qflux_g*ht_area*volpr(i)/DMAX1(1.0d4,cell%hg(i)-cell%hlsat(i))
!.....Rod-scale
      ht_area=2.d0*pi*ri_2d(nr_2d)*cell_leng(m,3)*0.25d0
      gamma_wall_rod(i)=gamma_wall_rod(i)+qflux_g*ht_area/DMAX1(1.0d4,hg-hlsat)
!
      IF(ISNAN(qflux_l).or.ISNAN(qflux_g))THEN
         WRITE(*,*)'condensation_HTC:',qflux_l,qflux_g
         WRITE(unit_log,*)'condensation_HTC:',qflux_l,qflux_g
      ENDIF
!       
      qc=DABS(qflux_g)
!      
      tdiff=(tw-tl)
      IF(DABS(tdiff).gt.1.d-1)THEN
         HTC_tl=qflux_l/tdiff
      ELSE
         qflux_l=0.0d0
         HTC_tl=0.0d0
      ENDIF   
      HTC_tst=0.0d0
      HTC_tg=0.0d0
      IF(DABS(dt_pps).gt.1.d-1)THEN
          HTC_tgp=qflux_g/dt_pps
      ELSE
          HTC_tgp=0.0d0
          qflux_g=0.0d0
      ENDIF    
!         
      END SUBROUTINE condensation_HTC_rod
            
