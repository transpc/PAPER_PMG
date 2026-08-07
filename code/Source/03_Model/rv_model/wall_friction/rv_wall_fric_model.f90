!
      SUBROUTINE rv_wall_fric_model
!
!     This routine calculates total pressure drop by wall friction and partition it.
!
      USE VOL_DATA                 
      USE Zvector      , ONLY: ul_o,ug_o
      USE Zparam       , ONLY: pi !pik 
      USE Zzone        , ONLY: ncell_fluid
      USE Zrv_wall_fric, ONLY: rough_wall
      USE Zconst2      , ONLY: hydraulicd      
!
      IMPLICIT NONE
!
      INTEGER i
      INTEGER hsf,vsf             
!      
      REAL(8) fwf,fwg,pd_htfs
      REAL(8) dpdx_nor, dpdx_fb
      REAL(8) a_liq,r_liq,v_liq,vis_liq,diameter,roughness !pik
      REAL(8) r_gas,v_gas,vis_gas     
      REAL(8) rey_gas,rey_liq,pfnrgj,dpdx_tot
      REAL(8) frac_liq, frac_gas
      REAL(8) epsilon, tolerance      
      REAL(8) a_vw,hst_angle,a_lw !pik
      REAL(8) ss_liq,ss_gas
      REAL(8) wall_liq,wall_gas !pik    
!
      DO i=1,ncell_fluid
!     -------------- variables -------------------
!     
!        a_liq : liquid fraction
!        r_liq : liquid density
!        v_liq : absolute value of liquid velocity 
!        vis_liq : liquid viscosity
!        r_gas : gas density
!        v_gas : absolute value of gas velocity
!        vis_gas : gas viscosity
!        diameter : hydraulic diameter
!        roughness : wall roughness (input of SPACE)
!     
!     -------------------------------------------------
         roughness=rough_wall
         epsilon = 1.0e-9
         tolerance = 1.0e-4
         hsf=1
         vsf=0
         hst_angle=pi/6.0d0    
!               
         a_liq=cell%alphal(i)  
         r_liq=cell%rhol(i)
         v_liq=ul_o(i)
         vis_liq=cell%lviscosl(i)
         r_gas=cell%rhog(i)
         v_gas=ug_o(i)
         vis_gas=cell%lviscosg(i)
         diameter=hydraulicd(i)
!         diameter=0.1d0
!
         a_lw=a_liq
         a_vw=1.0-a_liq
!      
         pfnrgj = cell%wf_dry(i)
!       
         pd_htfs=0.0     ! pressure drop, H.F.F.S correlation, Claxton et al. (1972)
!    
!        Calculate pressure drop
!    
!        define reynolds number (alpha * rho * |v| * hydraulic dia / mu)
!    
!         rey_gas=DMAX1(50.d0,DMAX1(1.0e-8, 1.0e0-a_liq)*r_gas*DMAX1(1.0e-5,v_gas)*diameter/vis_gas)
!         rey_liq=DMAX1(50.d0,DMAX1(1.0e-8, a_liq)*r_liq*DMAX1(1.0e-5,v_liq)*diameter/vis_liq)
         rey_gas=DMAX1(1.0d-3,DMAX1(1.0e-8, 1.0e0-a_liq)*r_gas*DMAX1(1.0e-5,v_gas)*diameter/vis_gas)
         rey_liq=DMAX1(1.0d-3,DMAX1(1.0e-8, a_liq)*r_liq*DMAX1(1.0e-5,v_liq)*diameter/vis_liq)         
!       
!        calculate fanning factor
!    
         CALL fanning_factor(rey_liq,roughness,diameter,fwf)
         CALL fanning_factor(rey_gas,roughness,diameter,fwg)
!    
	     dpdx_nor=0.0
	     dpdx_fb=0.0
!    
	     IF(pfnrgj<1.0e0)THEN
!    
!           pd_htfs=dpdx_htfs(0)	// 0 means normal flows
!    
            CALL dpdx_htfs(i,a_liq,fwf,fwg,r_liq,r_gas,vis_liq,vis_gas,diameter,0,pd_htfs)	 ! 0 means normal flows  
!    
!           no annular flow for cross flow
!    
            dpdx_nor=pd_htfs
         ENDIF
!    
     	 IF(pfnrgj>0.0e0)THEN	
	   	    CALL dpdx_htfs(i,a_liq,fwf,fwg,r_liq,r_gas,vis_liq,vis_gas,diameter,1,pd_htfs)	! 1 means film boiling flows
!    
           dpdx_fb=pd_htfs
	     ENDIF
!	     
         dpdx_tot=dpdx_nor**(1.0e0-pfnrgj)*dpdx_fb**(pfnrgj)
!      --------------------------------------------------------------------------------------------------      
!    
!        partition wall drage
!    
!        call wall_drag_partition
!    
!        calculate liquid/gas fraction on the wall for horizontally stratified flow
!    
!         IF (hsf) THEN
!            a_vw=dmax1(1.0e-5, hst_angle/PI)
!            a_lw=dmax1(1.0e-5, 1.0-a_vw)
!!    
!!        calculate liquid/gas fraction on the wall for vertically stratified flow
!!    
!         ELSE IF (vsf) THEN
!            a_lw=a_liq
!            a_vw=1.0-a_liq
!!    
!!        calculate liquid/gas fraction on the wall for annular flow
!!    
!         ELSE
!            aliq_min=0.01
!            fwet=sqrt(dmax1(1.0e-9,a_liq)/aliq_min)
!            a_lw=fwet
!            a_vw=1.0-fwet
!         END IF

!         CALL rv_ihtc_vst(i,vsf)     
         !IF (vsf) THEN
!            a_lw=a_liq
!            a_vw=1.0-a_liq
         !ELSE
         !   aliq_min=0.01
         !   fwet=DSQRT(DMAX1(1.0e-9,a_liq)/aliq_min)
         !   a_lw=fwet
         !   a_vw=1.0-fwet
         !END IF
!    
!        calculate fanning factor
!    
!	     CALL fanning_factor(rey_liq,roughness,diameter,fwf)
!	     CALL fanning_factor(rey_gas,roughness,diameter,fwg)
!        
!        calculate liquid/gas fraction
!        
!        ss_liq = a_lw*fwf_hst*r_liq*v_liq*v_liq
!        ss_vap = a_vw*fwg_hst*r_gas*v_vap*v_vap
!        frac_liq_hst=ss_liq/max(ss_liq+ss_vap,EPSILON)		
!        frac_vap_hst=(1.0-frac_liq_hst)*a_vap/max(EPSILON,a_vap+a_dis)
!        
!        ss_liq = a_lw*fwf_vst*r_liq*v_liq*v_liq
!        ss_vap = a_vw*fwg_vst*r_gas*v_vap*v_vap
!        frac_liq_vst=ss_liq/max(ss_liq+ss_vap,EPSILON)		
!        frac_vap_vst=(1.0-frac_liq_vst)*a_vap/max(EPSILON,a_vap+a_dis)
!        
         ss_liq = a_lw*fwf*r_liq*v_liq*v_liq
         ss_gas = a_vw*fwg*r_gas*v_gas*v_gas
         frac_liq = ss_liq/DMAX1(ss_liq+ss_gas,EPSILON)		
         frac_gas = 1.0-frac_liq
!       
!        calculate wall friction
!       
         wall_liq=frac_liq*dpdx_tot/DMAX1(tolerance,v_liq) ! force/velocity
	   	 wall_gas=frac_gas*dpdx_tot/DMAX1(tolerance,v_gas) ! force/velocity
!           
!        wall_liq=frac_liq*dpdx_tot/max(1.0e-8,a_liq)/max(tolerance,v_liq)	!phase intensive form
!        wall_vap=frac_vap*dpdx_tot/max(1.0e-8,a_vap)/max(tolerance,v_vap)	!phase intensive form
!    
         cell%vfwl(i)=wall_liq		
         cell%vfwg(i)=wall_gas
!	   	
      ENDDO
!
!      f_cal=diameter*(cell%p(1)-cell%p(ncell_fluid))/(2.d0*cell%rhol(ncell_fluid)*v_liq*v_liq*30.d0)
!      WRITE(123,100) time,v_liq,rey_liq,rey_gas,fwf,fwg,f_cal
!100   FORMAT(10(e14.7,1x))
!
      RETURN
      END SUBROUTINE rv_wall_fric_model
