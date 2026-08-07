!
      SUBROUTINE nucl_boiling_HTC_rod(i,mode,p,alphag,alphal,rhog,rhol,tg,tl,vg,vl,condg,condl,visg,visl,cpg,cpl,betag,betal)
!
!     This routine calculates wall heat transfer coefficient and heat
!     flux for single phase convection
!     
      USE STM_TBL_cupid  , ONLY: pcrit
      USE Zwall_HTC    , ONLY: qflux_t,qflux_l,qflux_g,HTC_t,HTC_tl,HTC_tst,HTC_tg,HTC_tgp, &
                               f_direc,dia_rod,mflux_liqa,mflux_gasa, &
                                hfg_p,rho,cond,viscos,cp,beta,sigma, &
                                qual_eq,tw,tsat_t,dt_sat
      USE Zio_unit     , ONLY: unit_log
!
      IMPLICIT NONE
!
      INTEGER i,mode
!
      REAL(8) tl_min
!     REAL(8) gravity
      REAL(8) h_mac,h_mic,q_mac,q_mic
      REAL(8) Re,Re_tp,f,s_f,x_tt,p_x,dpdt,p_mic,alpha
!
!.....Additional connectivity - jrlee
      REAL(8) p
      REAL(8) alphag,rhog,condg,visg,cpg,betag,tg,vg
      REAL(8) alphal,rhol,condl,visl,cpl,betal,tl,vl
      REAL(8) hinter,alo,ht1,ht2,alp,diff
!
!.....Interpolation function
!
      hinter(alo,ht1,ht2,alp,diff)=ht1+(ht2-ht1)*(alp-alo)*diff       
!      
      LOGICAL err
!      
!     DATA gravity/9.81d0/
!
!.....Set fluid property
!
      rho=rhol
      rhog=rhog
      cond=condl
      viscos=visl
      cp=cpl
      beta=betal                 
!
!.....Initial values
!
      tl_min=DMIN1(tsat_t,tl)
      mflux_liqa=DMAX1(1.0d0,alphal*rho*0.06d0,mflux_liqa) !!cyj: Set minimum mass flux
      Re=mflux_liqa*dia_rod/viscos
      f=1.d0
!
!.....calculate 'f' factor
!      
      IF(tl_min.ge.tsat_t-5.0d0)THEN      
         IF(DABS(mflux_liqa).gt.1.0d-3)THEN
            x_tt=(mflux_gasa/mflux_liqa)**0.9d0 *(rhol/rhog)**0.5d0*(visg/visl)**0.1d0
         ELSE
            x_tt=100.d0
         ENDIF 
         IF(x_tt.gt.0.1d0)f=2.35d0*(x_tt+0.213d0)**0.736d0
!
!........Interpolate f_factor from t_subcooled to t_saturation
!         
          f=f-0.2d0*(f-1.0d0)*DMAX1(0.0d0,DMIN1(5.0d0,tsat_t-tl_min))
      ENDIF
!
!.....Interpolate f_factor when dt_sat<1.0      
!
      IF(dt_sat.le.1.0d0)f=hinter(0.d0, 1.d0, f, dt_sat, 1.d0)
!
!.....Calculate 'h_mac' via single phase HTC (Max of Dittus-Boelter, Churchil-Chu, McAdams)
!     
      mode=2
      CALL single_phase_HTC_rod(i,mode,alphag,alphal,rhog,rhol,tg,tl,vg,vl,condg,condl,visg,visl,cpg,cpl,betag,betal)

      h_mac=HTC_t*f
!
!.....calculate 's' factor
!        
      Re_tp=DMIN1(70.d0,Re*f**1.25d0 *1.0d-4)
      IF(Re_tp.lt.32.5d0)THEN
         s_f=1.0d0/(1.0d0+0.12d0*Re_tp**1.14d0)
      ELSEIF(Re_tp.lt.70.0d0)THEN
         s_f=1.0d0/(1.0d0+0.42d0*Re_tp**0.78d0)
      ELSE
         s_f=0.0797d0
      ENDIF
!
!.....calculate 'h_mic' 
!
      CALL psatpd_cupid(tw,p_x,dpdt,1,err)
      IF(err)p_x=pcrit
      p_mic=p_x-p
      IF(p_mic.gt.0.and.dt_sat.ge.0)THEN
         p_mic=p_mic**0.75d0
         h_mic=0.00122d0*cond**0.79d0*cp**0.45d0*rho**0.49d0*(dt_sat/hfg_p/rhog)**0.24d0*p_mic*s_f/DSQRT(sigma)/viscos**0.29d0
      ENDIF
!
!.....IF equilibrium quality is greater than 0.95, consider heat transfer to gas phase
!      
      IF(qual_eq.ge.0.95d0)THEN
         mode=9
         rho=rhog
         !tl=tg
         cond=condg
         viscos=visg
         cp=cpg
         beta=betag
         CALL single_phase_HTC_rod(i,mode,alphag,alphal,rhog,rhol,tg,tl,vg,vl,condg,condl,visg,visl,cpg,cpl,betag,betal)
         alpha=DMIN1(qual_eq,0.999d0)     
         h_mac=hinter(0.95d0,h_mac,0.d0,alpha,25.d0)   
         h_mic=hinter(0.95d0,h_mic,0.d0,alpha,25.d0)  
         HTC_tg=hinter(0.95d0,0.d0,HTC_tg,alpha,25.d0)  
         qflux_g=HTC_tg*(tw-tl)
      ENDIF
!
!.....Save HTC and heat flux
!
      HTC_t=h_mac+h_mic+HTC_tg
      qflux_t=HTC_t*tw-h_mic*tsat_t-h_mac*tl-HTC_tg*tg
      qflux_l=qflux_t-qflux_g
!
!.....Compensate 'Negatice heat flux' (when liquid if superheat)      
!
      IF(qflux_t.lt.0.d0.and.qflux_g.eq.0)THEN
         q_mic=h_mic*(tw-tsat_t)
         q_mac=-(q_mic+qflux_g)
         h_mac=q_mac/DMAX1(1.d-9,(tw-tl))
         qflux_t=0.0d0
         qflux_l=q_mac+q_mic
         HTC_t=h_mac+h_mic+HTC_tg
      ENDIF
      HTC_tl=h_mac
      HTC_tst=h_mic
      HTC_tgp=0.0d0      
!
      IF(ISNAN(qflux_l).or.ISNAN(qflux_g))THEN
         WRITE(*,*)'nucl_boiling_HTC:',qflux_l,qflux_g
         WRITE(unit_log,*)'nucl_boiling_HTC:',qflux_l,qflux_g
      ENDIF
!       
      qflux_l=HTC_tl*(tw-tl)+HTC_tst*dt_sat    
!         
      RETURN     
      END SUBROUTINE nucl_boiling_HTC_rod
!
