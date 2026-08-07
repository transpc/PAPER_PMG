!
      SUBROUTINE dpdx_htfs(i, voidf, fwf, fwg, rhof, rhog, visf, visg, diameter, iflag, dpdx)
!
!     This routine calculates pressure drop, H.F.F.S correlation, Claxton et al. (1972).
!
      USE VOL_DATA                 
      USE Zvector      , ONLY: ul_o,ug_o
      USE Zrv_wall_fric, ONLY: rough_wall
      USE Zconst2      , ONLY: hydraulicd      
!      
!
      IMPLICIT NONE
!
      INTEGER i
      INTEGER iflag !pik
!      
      REAL(8) dpdx, dp_l, dp_v
      REAL(8) voidg, voidf
      REAL(8) c,mf,t1,f1,gamma
      REAL(8) rhof,v_liq,visf,visg,diameter,roughness,fwg,rhog,fwf,epsilon,v_gas !pik   
!
!   -------------- variables -------------------
!    voidf : liquid fraction
!    rhof : liquid density
!    v_liq : absolute value of liquid velocity 
!    visf : liquid viscosity
!    voidg : gas fraction
!    rhog : gas density
!    v_gas : absolute value of gas velocity
!    visg : gas viscosity
!    diameter : hydraulic diameter
!    roughness : wall roughness (input of SPACE)
!-------------------------------------------------
!
      voidf=cell%alphal(i)
      rhof=cell%rhol(i)
      v_liq=ul_o(i)
      visf=cell%lviscosl(i)
      voidg=cell%alphag(i)
      v_gas=ug_o(i)
      visg=cell%lviscosg(i)
      diameter=hydraulicd(i)
      roughness=rough_wall
!  
      epsilon = 1.0e-9
!
!   calculate pressure drop caused by wall friction
!
      voidg=1.0e0-voidf
	  dp_v=fwg*rhog*(voidg*v_gas)*(voidg*v_gas)
	  dp_l=fwf*rhof*(voidf*v_liq)*(voidf*v_liq)
!
	  mf=voidg*rhog*v_gas+voidf*rhof*v_liq
	  f1=dmax1(EPSILON,28.0-0.3*dsqrt(dabs(mf)))
!	
	  IF(iflag==0)THEN
		  gamma=rhog/rhof*(visf/visg)**0.2
      END IF
!    
	  IF(iflag==1)THEN
		  gamma=rhof/rhog*(visg/visf)**0.2
      END IF
!
	  t1=dexp(-(dlog10(gamma)+2.5)**2/dmax1(2.4-1.0e-4*mf,EPSILON))
	  c=dmax1(2.0,-2.0+f1*t1)
!	
	  dpdx=2.0*(dp_l + c*dsqrt(dp_l*dp_v) + dp_v)/diameter     
!
      RETURN
      END SUBROUTINE dpdx_htfs
