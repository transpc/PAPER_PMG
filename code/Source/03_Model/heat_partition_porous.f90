!
      SUBROUTINE heat_partition_porous
!
!     This routine defines heat partitioning model and calculate heat flux & 
!     using given solid wall temperature in porous cells.  
!
      USE VOL_DATA                 
      USE SOLID_DATA      , ONLY: solid
      USE Zparam          , ONLY: pi,sqrt_pi
      USE Zconst1         , ONLY: wconden
      USE Zcoord3         , ONLY: volp,aporous
      USE Zface           , ONLY: qecell
      USE Zflowregime     , ONLY: alphag_bc,alphag_cm
      USE Zheat_partition , ONLY: kfactor      
      USE Znum_cell       , ONLY: n_fluid
      USE Zqvol           , ONLY: qporous_gas,qporous_liq,qporous_gamma,gamma_wall
      USE Zvector         , ONLY: ug_o,ul_o
      USE Zzone           , ONLY: ncell_cond,nmaterial_c
!
      IMPLICIT NONE
!
      INTEGER i,ii
!      
      REAL(8) ndensity,bfreq,d_depart,deltarho,deltahlg
      REAL(8) deltaTs,deltaTl,deltaTg,deltaTsol
      REAL(8) twait,hconvl,hconvg,a_single,a_two
      REAL(8) qq,qe,qcl,qcg,weight
      REAL(8) st
      REAL(8) twall_porous
      REAL(8) theta
!   
      IF(wconden.lt.0)RETURN
!DIR$ SIMD
      DO i=1,ncell_cond
         IF(nmaterial_c(i).lt.0)THEN
      !
            ii=n_fluid(i)
      !
      !.....Assign solid wall temperature in a porous cell
      !
            twall_porous=solid%tsol(i)
      !
            deltaTs=DMAX1(0.0d0,twall_porous-cell%ts(ii))
            deltaTl=(twall_porous-cell%tl(ii))
            deltaTg=(twall_porous-cell%tg(ii))
            deltarho=DMAX1(0.0d0,cell%rhol(ii)-cell%rhog(ii))
            deltahlg=DMAX1(0.0d0,cell%hg(ii)-cell%hl(ii))
            deltaTsol=DMAX1(0.0d0,solid%tsol(i)-twall_porous)
      !
      !.....Nucleation site density (CFX4)
      !  
            ndensity=(185d0*deltaTs)**1.805d0
            theta=38.0d0
      !
      !.....Departure diameter model (Fritz)	      
      !
            d_depart=0.0208d0*theta*sqrt(cell%sigma(ii)/9.806d0/deltarho)
      !
      !.....Bubble departure frequency and buble waiting time
      ! 	    
            bfreq=DSQRT(4.0d0 * 9.806d0*deltarho/(3.0d0*d_depart*cell%rhol(ii)))
            twait=0.8d0/bfreq
      !
      !.....Bubble influence factor (area)
      !            
      !      kfactor=4.0d0
      !
      !.....Area fraction for two-phase and nucleate site density
      !            
            a_two=DMAX1(0.0d0,DMIN1(1.0d0,ndensity*pi*d_depart**2/4.0d0*kfactor))
            IF(a_two.ge.1.0d0)THEN
               ndensity=1.0d0/(pi*d_depart**2/4.0d0*kfactor)
            ENDIF
      !
      !.....Area fraction for single-phase, and convective heat transfer coefficients
      !           
            a_single=DMAX1(0.0d0,DMIN1(1.0d0,1.d0-a_two))
            st=0.0045d0
            hconvl=DMAX1(10.0d0,st*cell%rhol(ii)*cell%cpl(ii)*ul_o(ii))
            hconvg=st*cell%rhog(ii)*cell%cpg(ii)*ug_o(ii)
      !
      !.....Calculate partitioned heat
      !
            qq=2.0d0/sqrt_pi*DSQRT(twait*cell%lcondl(ii)*cell%rhol(ii)*cell%cpl(ii))*bfreq*a_two*deltaTl
            qe=ndensity*bfreq*pi/6.0d0*d_depart**3*cell%rhog(ii)*deltahlg
            qcl=hconvl*a_single*deltaTl
            qcg=hconvg * 1.0d0*deltaTg
      !
            weight=DMAX1(0.0d0,DMIN1(1.0d0,(alphag_cm-cell%alphag(ii))/(alphag_cm-alphag_bc)))
            qq=weight*qq
            qe=weight*qe
            qcl=weight*qcl
            qcg=(1.0d0-weight)*qcg
      !
      !.....Assign partitioned heats to equation variables
      !
            qporous_liq(ii)=qporous_liq(ii)+(qq+qcl)*aporous(ii)
            qporous_gas(ii)=qporous_gas(ii)+qcg*aporous(ii)
            qporous_gamma(ii)=qporous_gamma(ii)+qe*aporous(ii)
            gamma_wall(ii)=gamma_wall(ii)+qe/deltahlg*aporous(ii)/volp(ii) 
            qecell(ii)=qe     
      !   
            solid%hconv_rod_g(i)=hconvg
            solid%hconv_rod_l(i)=hconvl
!
         ENDIF
      ENDDO
!
      RETURN
      END SUBROUTINE heat_partition_porous
