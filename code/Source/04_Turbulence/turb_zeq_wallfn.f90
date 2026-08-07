!
      SUBROUTINE Wall_const
!
!     This routine assigns the constant coefficients
!
      USE wallfunction , ONLY: one,kappa,elog,ctrans,utau_max,utau_tol,  &
                               uplus_max,uplus_tol,sigt,rsigt,cdes,rcap, &
                               emkb,c16,c124,small20
!
      IMPLICIT NONE
!
      REAL(8) b_wfn
      REAL(8) fctrans
!
      one=1.d0
      kappa=0.41d0
      b_wfn=4.9d0
      elog=DEXP(kappa*b_wfn)
      ctrans=fctrans(kappa,elog)
!
      utau_max=10.0d0
      utau_tol=1.0d-5
      uplus_max=50.d0
      uplus_tol=1.d-5
!
!.....Turbulence Prandtl number
!
      sigt=1.d0
      cdes=0.1d0
!
!.....1/kappa
!
      rcap=one/kappa
      rsigt=one/sigt
!
      emkb=DEXP(-kappa*b_wfn)
      c16=1.d0/16.d0
      c124=1.d0/24.d0
      small20=1d-20
!      
      RETURN
      END SUBROUTINE Wall_const
!      
      FUNCTION fctrans(kappa,elog)
      USE Zio_unit        , ONLY: unit_log
!
!     This routine calculates
!
      IMPLICIT NONE
!      
      INTEGER:: i
!      
      REAL(8):: fctrans,kappa,elog,yplus_tol,yplus0,yplus
!
      yplus_tol=1.0d-5
      yplus0=11.0d0
      i=0
      DO
         i=i+1
         yplus=yplus0-fwall_layer(yplus0)/dfwall_layer(yplus0)
         IF(DABS(yplus-yplus0)<yplus_tol) EXIT
         yplus0=yplus
!         
         IF(i>100) THEN
            WRITE(*,*)"Iteration number for ctrans exceeds",i,', and yplus is',yplus  
            WRITE(unit_log,*)"Iteration number for ctrans exceeds",i,', and yplus is',yplus   
            EXIT
         ENDIF
!         
      ENDDO
      fctrans=yplus
!
      CONTAINS
!
         FUNCTION fwall_layer(yplus)
!      
            REAL(8):: fwall_layer,yplus
!         
            fwall_layer=yplus-1.0d0/kappa*dlog(elog*yplus)
         END FUNCTION fwall_layer
!
         FUNCTION dfwall_layer(yplus)
!      
            REAL(8):: dfwall_layer,yplus
!         
            dfwall_layer=1.0d0-1.0d0/(kappa*yplus)
         END FUNCTION dfwall_layer
!
      END FUNCTION fctrans
!
      SUBROUTINE Wall_Zeq_i(icell,inbcon,vgn_i,vln_i,lviscosl_mi,tviscosl_mi,   &
                             lcondl_mi,cpl_mi,rhol_mi,lviscosg_mi,tviscosg_mi,   &
                             lcondg_mi,cpg_mi,rhog_mi,csv,csa,dji_n,wviscosl,    &
                             wcondl,wviscosg,wcondg,utaul,yplusl)
!
!     This routine calculates turbulence viscosity at the wall
!     TWM_NO       = 0, NO wall function
!     TWM_Peric0   = 2, wall function for zero-eq model
!     TWM_Utau     = 3, wall function based on the utau
!
      USE wallfunction , ONLY: rcap,ctrans,elog,one,sigt,rsigt, &
                               small20 
!
      IMPLICIT NONE
!
      INTEGER:: icell,inbcon
      INTEGER:: wall_model
!      
      INTEGER, PARAMETER:: TWM_NO=0,TWM_Peric0=2,TWM_Utau=3      
!
      REAL(8):: vgn_i(3),vln_i(3)
      REAL(8):: lviscosl_mi,tviscosl_mi,lcondl_mi,cpl_mi,rhol_mi
      REAL(8):: lviscosg_mi,tviscosg_mi,lcondg_mi,cpg_mi,rhog_mi
      REAL(8):: csv(3),csa,dji_n
      REAL(8):: wviscosl,wcondl,wviscosg,wcondg
      REAL(8):: sx,sy,sz,sa,dn,vnorm,ut,vt,wt,sigl,rsigl,sig_rat
      REAL(8):: wall_fn,jayatillaka
      REAL(8):: uc_bcell=0.d0,vc_bcell=0.d0,wc_bcell=0.d0
      REAL(8):: veltl,veltg,utaul,utaug
      REAL(8):: uplusl,uplusg,tplusl,tplusg,yplusl,yplusg
      REAL(8):: yplusl_min,yplusl_max,yplusg_min,yplusg_max
!
      CALL Wall_const
!
!.....Select type of wall function
!
      wall_model=TWM_Peric0
!
      SELECT CASE(wall_model)
         CASE(TWM_NO)
            wviscosl=0.d0
            wcondl=0.d0
            wviscosg=0.d0
            wcondg=0.d0
         CASE(TWM_Peric0)
!
!...........For liquid
!
            yplusl_min=10000.0d0
            yplusl_max=0.0d0
            yplusl=rcap*(lviscosl_mi+tviscosl_mi)/lviscosl_mi
!
            IF(yplusl<=ctrans)THEN
               uplusl=yplusl
               wviscosl=lviscosl_mi
               wcondl=lcondl_mi
            ELSE
               uplusl=rcap*DLOG(elog*yplusl)
               wviscosl=lviscosl_mi*DMAX1(one,yplusl/uplusl)
               sigl=cpl_mi*lviscosl_mi/lcondl_mi
               rsigl=1.0d0/sigl
               sig_rat=sigl*rsigt
               tplusl=(sigt*rsigl)*(uplusl+jayatillaka(sig_rat))
               wcondl=(cpl_mi*rsigl*lviscosl_mi)*DMAX1(one,yplusl/tplusl)
            ENDIF
!
            yplusl_min=DMIN1(yplusl,yplusl_min)
            yplusl_max=DMAX1(yplusl,yplusl_max)
!               
!...........For gas
!
            yplusg_min=10000.0d0
            yplusg_max=0.0d0
            yplusg=rcap*(lviscosg_mi+tviscosg_mi)/lviscosg_mi
!
            IF(yplusg<=ctrans)THEN
               uplusg=yplusg
               wviscosg=lviscosg_mi
               wcondg=lcondg_mi
            ELSE
               uplusg=rcap*DLOG(elog*yplusg)
               wviscosg=lviscosg_mi*DMAX1(one,yplusg/uplusg)
               sigl=cpg_mi*lviscosg_mi/lcondg_mi
               rsigl=1.0d0/sigl
               sig_rat=sigl*rsigt
               tplusg=(sigt*rsigl)*(uplusg+jayatillaka(sig_rat))
               wcondg=(cpg_mi*rsigl*lviscosg_mi)*DMAX1(one,yplusg/tplusg)
            ENDIF
!
            yplusg_min=DMIN1(yplusg,yplusg_min)
            yplusg_max=DMAX1(yplusg,yplusg_max)
!
         CASE(TWM_Utau)
            sx=csv(1)
            sy=csv(2)
            sz=csv(3)
            sa=csa
            dn=dji_n
!         
!...........For liquid
!
            yplusl_min=10000.0d0
            yplusl_max=0.0d0
!            
!...........Calculate liquid tangential velocity at icell
!
            vnorm=sx*vln_i(1)+sy*vln_i(2)+sz*vln_i(3)
            ut=vln_i(1)-(sx*vnorm+uc_bcell)
            vt=vln_i(2)-(sy*vnorm+vc_bcell)
            wt=vln_i(3)-(sz*vnorm+wc_bcell)
            veltl=DSQRT(ut*ut+vt*vt+wt*wt)
            veltl=DMAX1(veltl,small20)
!            
!...........Calculate liquid shear velocity using the wall function
!
            utaul=wall_fn(lviscosl_mi,rhol_mi,dn,veltl,icell,inbcon)
            yplusl=rhol_mi*utaul*dn/lviscosl_mi
!            
!...........Liquid wall viscosity
!
            IF(yplusl<=ctrans) THEN
               wviscosl=lviscosl_mi
               wcondl=lcondl_mi
            ELSE
               uplusl=veltl/utaul
               wviscosl=lviscosl_mi*DMAX1(one,yplusl/uplusl)
               sigl=cpl_mi*lviscosl_mi/lcondl_mi
               rsigl=1.0d0/sigl
               sig_rat=sigl*rsigt
               tplusl=(sigt*rsigl)*(uplusl+jayatillaka(sig_rat))
               wcondl=(cpl_mi*rsigl*lviscosl_mi)*DMAX1(one,yplusl/tplusl)
            ENDIF
            yplusl_min=DMIN1(yplusl,yplusl_min)
            yplusl_max=DMAX1(yplusl,yplusl_max)
!
!...........For gas
!
            IF(0) THEN
               yplusg_min=10000.0d0
               yplusg_max=0.0d0
!               
!..............Calculate gas tangential velocity at icell
!
               vnorm=sx*vln_i(1)+sy*vln_i(2)+sz*vln_i(3)
               ut=vgn_i(1)-(sx*vnorm+uc_bcell)
                vt=vgn_i(2)-(sy*vnorm+vc_bcell)
               wt=vgn_i(3)-(sz*vnorm+wc_bcell)
               veltg=DSQRT(ut*ut+vt*vt+wt*wt)
               veltg=DMAX1(veltg,small20)
!               
!..............Calculate gas shear velocity using the wall function
!
               utaug=wall_fn(lviscosg_mi,rhog_mi,dn,veltg,icell,inbcon)
               yplusg=rhog_mi*utaug*dn/lviscosg_mi
!               
!..............Calculate gas wall viscosity
!
               IF(yplusg<=ctrans) THEN
                  wviscosg=lviscosg_mi
                  wcondg=lcondg_mi
               ELSE
                  uplusg=veltg/utaug
                  wviscosg=lviscosg_mi*DMAX1(one,yplusg/uplusg)
                  sigl=cpg_mi*lviscosg_mi/lcondg_mi
                  rsigl=1.0d0/sigl
                  sig_rat=sigl*rsigt
                  tplusg=(sigt*rsigl)*(uplusg+jayatillaka(sig_rat))
                  wcondg=(cpg_mi*rsigl*lviscosl_mi)*DMAX1(one,yplusg/tplusg)
               ENDIF
!            
               yplusg_min=dmin1(yplusg,yplusg_min)
               yplusg_max=dmax1(yplusg,yplusg_max)
            ENDIF
      END SELECT

      RETURN
      END
!
      FUNCTION jayatillaka(sig_rat)
!
!     This function is Jayatillaka's function for turbulence wall function
!     heat transfer coefficient
!
      IMPLICIT NONE
!
      REAL(8):: jayatillaka,sig_rat

      jayatillaka=0.9d0*(sig_rat-1.0d0)*sig_rat**(-0.25d0)

      END FUNCTION jayatillaka
!
      FUNCTION wall_fn(visi,rhoi,delyi,veli,icell,inbcon)
!
!     This function is wall function model to find shear velocity utau 
!     using Newtom-Raphson method
!
      USE wallfunction , ONLY: small20,utau_max,utau_tol,rcap
      USE Zio_unit     , ONLY: unit_log
!
      IMPLICIT NONE
!      
      INTEGER:: iter,inbcon,icell
!      
      REAL(8):: wall_fn,visi,rhoi,delyi,veli,utau,utau0,re_local
!
      re_local=rhoi*veli*delyi/visi
!
      IF(re_local<=5.0d0) THEN
         utau=veli/DSQRT(re_local+small20)
      ELSE
         utau0=0.045d0*veli
         iter=0
         DO
            iter=iter+1
            IF(inbcon.ge.-1) THEN
               utau=utau0-fn_wall(utau0)/dfn_wall(utau0)
            ELSE
               utau=utau0-fn_wall_yun(utau0)/dfn_wall_yun(utau0)
            ENDIF
            IF(DABS(utau-utau0)<utau_tol) EXIT
            utau0=DMIN1(utau,utau_max)
            IF(iter.gt.1000) THEN
               WRITE(*,*)'          Iteration number for utau exceeds, max=',iter
               WRITE(unit_log,*)'          Iteration number for utau exceeds, max=',iter
               EXIT
            ENDIF
         ENDDO

      ENDIF
!
      wall_fn=utau
!
      CONTAINS
!
         FUNCTION fn_wall(utau)
!
!        This function calculates wall function value 
!        with m-iterative step utau
!      
         USE wallfunction , ONLY: elog
!         
         IMPLICIT NONE
!      
         REAL(8):: fn_wall,utau
!      
         fn_wall=veli/utau-rcap*DLOG(elog*rhoi*delyi/visi*utau)
!      
         END FUNCTION fn_wall
!
         FUNCTION dfn_wall(utau)
!
!        This function calculates differential wall function value 
!        with m-iterative step utau
!         
         REAL(8):: dfn_wall,utau
!      
         dfn_wall=-veli/(utau*utau)-rcap/utau
!      
         END FUNCTION dfn_wall
!
         FUNCTION fn_wall_yun(utau)
!
!        This function calculates modified wall function by Yun's correction
!        B.J.Yun, et.al., Prediction of a subcooled boiling flow with advanced
!        two-phase flow models
!
         USE VOL_DATA                 
         USE Wall_DATA    , ONLY: face
!
         IMPLICIT NONE
!
         REAL(8):: fn_wall_yun,utau,kr_plus,y_plus,E,calc_E_yun,Dbubble
         REAL(8):: kr
!      
         REAL(8), PARAMETER:: kappa_loc=0.4d0
         REAL(8), PARAMETER:: yplus_cr=11.23d0
         REAL(8), PARAMETER:: eta=1.0d0
         REAL(8), PARAMETER:: ksi=0.174d0
!
         Dbubble=cell%D1(icell)
         kr=eta*Dbubble*face%ratio_evap(icell)**ksi
         kr_plus=rhoi*kr*utau/visi
         E = calc_E_yun(kappa_loc,kr_plus)
         y_plus=rhoi*delyi/visi*utau
         fn_wall_yun=veli/utau-DMIN1(y_plus,1.0/kappa_loc*DLOG(E*y_plus))
!      
         IF(y_plus.le.yplus_cr) THEN
            fn_wall_yun=veli/utau-1.0/kappa_loc*DLOG(E*yplus_cr)*y_plus/yplus_cr
         ENDIF    
!      
         END FUNCTION fn_wall_yun
!
         FUNCTION dfn_wall_yun(utau)
!
!        This function calculates modified differential wall function value
!        by Yun's correction
!      
         REAL(8):: dfn_wall_yun,utau
         REAL(8), PARAMETER:: kappa=0.4d0
!      
         dfn_wall_yun=-veli/(utau*utau)-1.0d0/kappa/utau
!      
         END FUNCTION dfn_wall_yun      
!      
      END FUNCTION wall_fn
!
      FUNCTION calc_E_yun(kappa,kr_plus)
!
!     This function correct E of wall fuction with kr_plus
!
      IMPLICIT NONE
!
	  REAL(8):: kr_plus
      REAL(8):: calc_E_yun,rough_B_yun,B,kappa
 !     
      IF(kr_plus.le.4.0d0) THEN
         calc_E_yun = 9.0d0
	  ELSEIF(kr_plus.ge.4.0d0.and.kr_plus.le.70.0d0) THEN
         B=rough_B_yun(kr_plus)
		 calc_E_yun = DEXP(kappa*B)/kr_plus
      ELSE
         calc_E_yun = DEXP(kappa*8.5)/kr_plus
	  ENDIF
!	  
      RETURN
	  END FUNCTION calc_E_yun
!
      FUNCTION rough_B_yun(kr_plus)
!
!     This function calculates wall roughness
!      
      IMPLICIT NONE
!      
	  REAL(8):: rough_B_yun,kr_plus,x
!	  
	  REAL(8), PARAMETER:: A  = 6.1063d0
	  REAL(8), PARAMETER:: B1 =-1.8692d0
	  REAL(8), PARAMETER:: B2 =28.54543d0
	  REAL(8), PARAMETER:: B3 =-44.08548d0
	  REAL(8), PARAMETER:: B4 =28.54657d0
	  REAL(8), PARAMETER:: B5 =-8.7074d0
	  REAL(8), PARAMETER:: B6 =1.03894d0
!
      x=DLOG10(kr_plus)   
      rough_B_yun = A+B1*x+B2*x**2.0d0+B3*x**3.0d0+B4*x**4.0d0 &
                + B5*x**5.0d0 + B6*x**6.0d0
!
      RETURN
	  END FUNCTION rough_B_yun
