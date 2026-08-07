!
      SUBROUTINE condensation_HTC(i,mode,ht_area,qc)
!
!     This routine calculates wall heat transfer coefficient for condensation near wall
!
      USE VOL_DATA      , ONLY: cell
      USE STM_TBL_cupid , ONLY: pcrit
      USE Zcoord3       , ONLY: volpr
      USE Zqvol         , ONLY: gamma_wall
      USE Zwall_HTC     , ONLY: HTC_T,HTC_tl,HTC_tst,HTC_tg,HTC_tgp,HTC_tgt,qflux_t,qflux_l,qflux_g,          &
                                tw,rho,cond,cp,viscos,tl,beta,qual_eq,mflux_tota,mflux_liqa, dia_rod,HTC_cond
      USE Zmodel        , ONLY: qconden
      USE Zio_unit      , ONLY: unit_log
!
      IMPLICIT NONE
!      
      INTEGER i,mode
!
      REAL(8) ht_area
      REAL(8) dt_pps,dt_liq,qual_max,Re,Pr,hdb,hf,z,ftr
      REAL(8) HTC_ditt,HTC_Shah,qc
      REAL(8) gravity,Re_film,film_thick,fr_x
      
      DATA gravity/9.81d0/      
!      
      dt_pps=tw-cell%ts(i)
      dt_liq=tw-cell%tl(i)
!
      rho=cell%rhol(i) 
      tl=cell%tl(i) 
      cond=cell%lcondl(i) 
      viscos=cell%lviscosl(i) 
      cp=cell%cpl(i) 
      beta=cell%betal(i)       
      IF(cell%alphag(i).le.0.3d0.or.dt_liq.gt.0.0d0)THEN  !to smooth transition from condensation to other modes  
         mode=2
         CALL single_phase_HTC(i,mode)                    ! Calculation of HTC_tl(i), qflux_l
         HTC_ditt=HTC_tl
         qflux_g=0.0d0
      ELSE
         HTC_ditt=0.0d0
         qflux_t=0.0d0
      ENDIF
!
!.....Laminar film condensation
!.....1.Vertical pipe or plate
!     Base it on film thickness to make it a local form instead of the average value used in MOD3 up to 3.1.1.1              
!     Local Form of Nusselt for Vertical Surfaces.                         
!      -  calculate film thickness from film Re.                         
!      -  gravity should be replaced into the direction cosine for the gravity vector to consider the angle of grid
! 
      Re_film=DABS(mflux_liqa)*dia_rod/cell%lviscosl(i)
      film_thick=0.9086d0*(Re_film*(cell%lviscosl(i)/cell%rhol(i))**2/gravity)**0.333333d0
      film_thick=DMAX1(film_thick,1.0d-5)
      HTC_cond=cell%lcondl(i)/film_thick  !Nusselt HTC.                                                                   
      HTC_cond=DMAX1(HTC_cond,4.36d0*cell%condl(i)/dia_rod) !do not allow the HTC to be less than that for laminar flow.      
!
!.....2.Horizontal stratified condensation --- !!!cyj: temporary Inactivation because CUPID cannot automatically determine the angle of each cell
!   
      !HTC_cond=0.296d0*((cell%rhol(i)*DMAX1(cell%rhol(i)-cell%rhog(i)),0.0d0)*gravity*hfg*cell%conl(i)**3) &
      !         /(dia_rod*cell%lviscosl(i)*DMAX1(-dt_pps,1.0d0))**0.25d0
!
!.....turbulent condensation heat transfer correlation
!
      IF(cell%quala(i).lt.0.001d0)THEN
         qual_max=DMAX1(1.0d-9,DMIN1(1.0d0,qual_eq))
         Re=DMIN1(4000.0d0,mflux_tota)*dia_rod/viscos
         Pr=viscos*cp/cond
         hdb=0.023d0*cond*Re**0.8d0*Pr**0.4d0/dia_rod
         hf=hdb*(1.0d0-qual_max)**0.8d0
         z=(cell%p(i)/pcrit)**0.04d0*(1.0d0/qual_max-1.0d0)**0.8d0
         ftr=1.0d0
         IF(z.ne.0.0d0)ftr=1.0d0+3.8d0/z**0.95d0
         HTC_Shah=hf*ftr
         HTC_cond=DMAX1(HTC_cond,HTC_Shah) 
      ENDIF
!
!.....Interoplation according to the range of Quala
!      
      IF(cell%quala(i).lt.0.0001d0)THEN
         !do nothing
      ELSEIF(cell%quala(i).gt.0.001d0)THEN
         CALL condensation_ncg(i)
         qflux_g=0.0d0
         HTC_cond=HTC_cond*qflux_t/HTC_cond/dt_pps
         HTC_tl=HTC_cond
      ELSE ! 0.0001<cell%quala(i)<0.001
         HTC_Shah=HTC_cond
         CALL condensation_ncg(i)
         fr_x=(0.001d0-cell%quala(i))*1111.1111d0       ! 1111.1111=1/(0.001-0.0001)       
         HTC_cond=HTC_Shah*fr_x+HTC_cond*(1.0d0-fr_x)
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
         qflux_g=HTC_cond*dt_pps
      ENDIF      
!
!.....Interpolete HTC&qflux if ag < 0.3
!     
      IF(cell%alphag(i).lt.0.3d0)THEN
         HTC_tl=HTC_ditt+(HTC_tl-HTC_ditt)*(cell%alphag(i)-0.1d0)*5.0d0    ! 5.0 = 1/(0.3-0.1)
         HTC_tgp=HTC_tgp*(cell%alphag(i)-0.1d0)*5.0d0
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
      gamma_wall(i)=gamma_wall(i)+qflux_g*ht_area*volpr(i)/DMAX1(1.0d4,cell%hg(i)-cell%hlsat(i))
      gamma_wall(i)=DMAX1(gamma_wall(i),-0.9d0*cell%alphag(i)*cell%rhog(i))
      qconden(i)=qflux_t
!
      IF(ISNAN(qflux_l).or.ISNAN(qflux_g))THEN
         WRITE(*,*)'condensation_HTC:',qflux_l,qflux_g
         WRITE(unit_log,*)'condensation_HTC:',qflux_l,qflux_g
      ENDIF
!       
      qc=DABS(qflux_g)      ! qflux_g+qc=0, mass/energy are already considered in gamma_wall
!        
      RETURN
      END SUBROUTINE condensation_HTC
            
