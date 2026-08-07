!
      SUBROUTINE int_htc_full_topology(ag)
!
!     This routine calculates the heat transfer coefficient at the interface 
!      when iHTC=1 and udfl_calc_HTC_int_i=.true..
!
      USE VOL_DATA                
      USE Zparam      , ONLY: ndim      
      USE Zdaint_ag   , ONLY: daint1_ag,daint2_ag,daint1_ag_bc,daint2_ag_cm,  &
                               D1_01,D1_bc,D1_09,                             &
                               aint_bc,aint_cm
      USE Zdhda       , ONLY: dHldag,dHgdag,dHfgdag
      USE Zflowregime , ONLY: alphag_bc,alphag_cm,gamma_1,gamma_2
      USE Zmodel      , ONLY: H_il_min,H_ig_min
      USE Zmpi        , ONLY: ncell_fp
      USE Zqvol       , ONLY: H_ig,H_il,H_gf
      USE Ztimecon    , ONLY: alpha_min
      USE Zvector     , ONLY: vrel_o,vl_o,vg_o,ul_o,ug_o
      USE Zzone       , ONLY: ncell_fluid
      USE Zvoid       , ONLY: dagdx,gamma_void
!
      IMPLICIT NONE
!
      INTEGER i
!      
      REAL(8) ag(ncell_fp)
      REAL(8) Nu,Re,Pr,st,hig1,hig2,hig3,hil1,hil2,hil3,delv
      REAL(8) H_il_i,H_il_1,H_il_2,H_il_3,H_ig_1,H_ig_2,H_ig_3,H_ig_i
      REAL(8) dHldag_1,dHldag_2,dHldag_3,dHgdag_1,dHgdag_2,dHgdag_3
      REAL(8) hil_01,hil_bc,hil_09
      REAL(8) hig_01,hig_09
      REAL(8) utangent(ndim)  
      REAL(8) weight    
!
!DIR$ SIMD
      DO i=1,ncell_fluid
         delv=vrel_o(i)
!
!........For bubbly flow
!
         Re=cell%rhol(i)*delv*cell%D1(i)/cell%lviscosl(i)
         Pr=cell%lviscosl(i)*cell%cpl(i)/cell%lcondl(i)
         Nu=2.0d0+0.6d0*Re**0.5d0*Pr**0.3d0
         hil1=cell%lcondl(i)*Nu
         hil_01=hil1/D1_01(i)
         hil_bc=hil1/D1_bc(i)
         hil_09=hil1/D1_09(i)
         hil1=hil1/cell%D1(i)
         hig1=1000.0d0
!
!........For mist flow
!
         Re=cell%rhog(i)*delv*cell%D2(i)/cell%lviscosg(i)
         Pr=cell%lviscosg(i)*cell%cpg(i)/cell%lcondg(i)
         Nu=2.0d0+0.6d0*Re**0.5d0*Pr**0.3d0
         hil2=cell%lcondl(i)*Nu/cell%D2(i)
         hig2=cell%lcondg(i)*Nu/cell%D2(i)
!         hig2=1000.0d0
!
!........For sharp interface
!         
         utangent(:)=(vl_o(i,:)-vg_o(i,:))-DOT_PRODUCT(vl_o(i,:)-vg_o(i,:),dagdx(i,:))  &
                    /DMAX1(alpha_min,DSQRT(DOT_PRODUCT(dagdx(i,:),dagdx(i,:))))             &
                    *dagdx(i,:)/DMAX1(alpha_min,DSQRT(DOT_PRODUCT(dagdx(i,:),dagdx(i,:))))
         st=0.0045d0*((cell%rhog(i)*ug_o(i)*cell%lviscosl(i))/(cell%rhol(i)*DMAX1(ul_o(i),1.d-8)*cell%lviscosg(i)))**(1.0d0/3.0d0)
         hil3=st*cell%rhol(i)*cell%cpl(i)*DSQRT(DOT_PRODUCT(utangent(:),utangent(:)))
         hig3=st*cell%rhog(i)*cell%cpg(i)*DSQRT(DOT_PRODUCT(utangent(:),utangent(:)))
!
!........Limit value
!
         hig1=DMAX1(hig1,1.0d0)
         hig2=DMAX1(hig2,1.0d0)
         !hig3=DMAX1(hig1,1.0d0) ! for sharp interface: case1
         hig3=DMAX1(hig3,1.0d0) ! for sharp interface: case2
         
         hil1=DMAX1(hil1,1.0d0)
         hil2=DMAX1(hil2,1.0d0)
         !hil3=DMAX1(hil1,1.0d0) ! for sharp interface: case1
         hil3=DMAX1(hil3,1.0d0) ! for sharp interface: case2
!
!........Bubbly flow regime
!
         IF(cell%regime(i).eq.11)THEN
            H_il(i)=hil1*cell%aint1(i)
            H_ig(i)=hig1*cell%aint1(i)
            dHldag(i)=hil1*daint1_ag(i)
            dHgdag(i)=hig1*daint1_ag(i)
!
!........Mist flow regime
!
         ELSEIF(cell%regime(i).eq.13)THEN
            H_il(i)=hil2*cell%aint2(i)
            H_ig(i)=hig2*cell%aint2(i)
            dHldag(i)=hil2*daint2_ag(i)
            dHgdag(i)=hig2*daint2_ag(i)
!
!........Bubbly-Mist transition regime
!           
         ELSEIF(cell%regime(i).eq.12)THEN
            H_il_1=hil_bc*aint_bc(i)
            H_il_2=hil2*aint_cm(i)
            H_ig_1=hig1*aint_bc(i)
            H_ig_2=hig2*aint_cm(i)
            dHldag_1=hil_bc*daint1_ag_bc(i)
            dHldag_2=hil2*daint2_ag_cm(i)
            dHgdag_1=hig1*daint1_ag_bc(i)
            dHgdag_2=hig2*daint2_ag_cm(i)
!
!...........Cubic interpolation for the trasition region
!
            CALL cub_interp2(ag(i),alphag_bc,alphag_cm,H_il_1,H_il_2,dHldag_1,dHldag_2,H_il(i),dHldag(i))
            CALL cub_interp2(ag(i),alphag_bc,alphag_cm,H_ig_1,H_ig_2,dHgdag_1,dHgdag_2,H_ig(i),dHgdag(i))
!
!........Sharp interface regime
!        
         ELSEIF(cell%regime(i).eq.3)THEN
            H_il(i)=hil3*cell%aint3(i)
            H_ig(i)=hig3*cell%aint3(i)
            dHldag(i)=0.d0 
            dHgdag(i)=0.d0
!
!........Bubbly-interface transition regime
!             
         ELSEIF(cell%regime(i).eq.21)THEN
            H_il_1=hil1*cell%aint1(i)
            H_il_3=DMAX1(H_il_min,hil3*cell%aint3(i))
            H_ig_1=hig1*cell%aint1(i)
            H_ig_3=DMAX1(H_ig_min,hig3*cell%aint3(i))
            dHldag_1=hil1*daint1_ag(i)
            dHldag_3=0.d0
            dHgdag_1=hig1*daint1_ag(i)
            dHgdag_3=0.d0
!            
            weight=(gamma_2-gamma_void(i))/(gamma_2-gamma_1)
            H_il_i=weight*H_il_1+(1.0d0-weight)*H_il_3
            H_ig_i=weight*H_ig_1+(1.0d0-weight)*H_ig_3
            dHldag(i)=weight*dHldag_1+(1.0d0-weight)*dHldag_3
            dHgdag(i)=weight*dHgdag_1+(1.0d0-weight)*dHgdag_3  
            H_il(i)=H_il_i            
            H_ig(i)=H_ig_i
!
!........Mist-interface transition regime
!  
         ELSEIF(cell%regime(i).eq.23)THEN
            H_il_2=hil2*cell%aint2(i)
            H_il_3=DMAX1(H_il_min,hil3*cell%aint3(i))
            H_ig_2=hig2*cell%aint2(i)
            H_ig_3=DMAX1(H_ig_min,hig3*cell%aint3(i))
            dHldag_2=hil2*daint2_ag(i)
            dHldag_3=0.d0
            dHgdag_2=hig2*daint2_ag(i)
            dHgdag_3=0.d0   
            weight=(gamma_2-gamma_void(i))/(gamma_2-gamma_1)
            H_il_i=weight*H_il_2+(1.0d0-weight)*H_il_3
            H_ig_i=weight*H_ig_2+(1.0d0-weight)*H_ig_3
            dHldag(i)=weight*dHldag_2+(1.0d0-weight)*dHldag_3
            dHgdag(i)=weight*dHgdag_2+(1.0d0-weight)*dHgdag_3  
            H_il(i)=H_il_i            
            H_ig(i)=H_ig_i
!
!........Churn-interface transition regime
!           
         ELSEIF(cell%regime(i).eq.22)THEN
            H_il_1=hil_bc*aint_bc(i)
            H_il_2=hil2*aint_cm(i)
            H_ig_1=hig1*aint_bc(i)
            H_ig_2=hig2*aint_cm(i)
            dHldag_1=hil_bc*daint1_ag_bc(i)
            dHldag_2=hil2*daint2_ag_cm(i)
            dHgdag_1=hig1*daint1_ag_bc(i)
            dHgdag_2=hig2*daint2_ag_cm(i)
!
!...........Cubic interpolation for the trasition region
!
            CALL cub_interp2(ag(i),alphag_bc,alphag_cm,H_il_1,H_il_2,dHldag_1,dHldag_2,H_il_1,dHldag_1)
            CALL cub_interp2(ag(i),alphag_bc,alphag_cm,H_ig_1,H_ig_2,dHgdag_1,dHgdag_2,H_ig_1,dHgdag_1)         
!         
            H_il_2=DMAX1(H_il_min,hil3*cell%aint3(i))
            H_ig_2=DMAX1(H_ig_min,hig3*cell%aint3(i))
            dHldag_2=0.d0
            dHgdag_2=0.d0
 !           
            weight=(gamma_2-gamma_void(i))/(gamma_2-gamma_1)
            H_il_i=weight*H_il_1+(1.0d0-weight)*H_il_2
            H_ig_i=weight*H_ig_1+(1.0d0-weight)*H_ig_2
            dHldag(i)=weight*dHldag_1+(1.0d0-weight)*dHldag_2
            dHgdag(i)=weight*dHgdag_1+(1.0d0-weight)*dHgdag_2  
            H_il(i)=H_il_i            
            H_ig(i)=H_ig_i
         ENDIF
!
         hil_01=H_il(i)
         hig_01=H_ig(i)
         hil_09=H_il(i)
         hig_09=H_ig(i)
!
!........Linear interpolation of Hik for temperature when alphag<0.1 or alphag>0.9
!
         CALL int_htc_interp_temp2(i,ag(i),hil_01,hig_01,hil_09,hig_09)
!
!........Relaxation of Hik
!
!!       H_il(i)=(1.0d0-relax_hik)*H_il(i)+relax_hik*hil_o(i)
!!       H_ig(i)=(1.0d0-relax_hik)*H_ig(i)+relax_hik*hig_o(i)
!!       hil_o(i)=H_il(i)
!!       hig_o(i)=H_ig(i)
!         
         H_gf(i)=H_ig(i)
         dHfgdag(i)=dHgdag(i)
!
      ENDDO
!
      RETURN
      END SUBROUTINE int_htc_full_topology
!
!DEC$ ATTRIBUTES INLINE :: int_htc_interp_temp2
      SUBROUTINE int_htc_interp_temp2(i,ag,hil_01,hig_01,hil_09,hig_09)
!
!     Linear interpolation of Hik for temperature when alphag < 0.1 or  alphag > 0.9
!
      USE VOL_DATA                
      USE Zdhda       , ONLY: dHldag,dHgdag,dHldtl,dHgdtg,dHfgdtg
      USE Zmodel      , ONLY: H_il_min,H_ig_min,dtl,dtg
      USE Ztimecon    , ONLY: alpha_min
      USE Zqvol       , ONLY: H_ig,H_il
!
      IMPLICIT NONE
!      
      INTEGER i
!
      REAL(8) ag
      REAL(8) Hila,Higa
      REAL(8) hil_01,hig_01,hil_09,hig_09
!
      dHldtl(i)=0.0d0
      dHgdtg(i)=0.0d0
!
!.....Subcooled liquid only
!
      IF(ag.le.2.0d0*alpha_min)THEN
         IF(cell%tl(i).lt.cell%ts(i)) H_il(i)=0.0d0
      ENDIF
!
!.....Superheated steam only
!
      IF(ag.ge.(1.0d0-2.0d0*alpha_min))THEN
         IF(cell%tg(i).gt.cell%ts(i)) H_ig(i)=0.0d0
      ENDIF
!
      IF(ag.le.0.1d0)THEN
!
!........Linear interpolation for liquid temperature when alphag < 0.1
!
         IF(cell%tl(i).lt.cell%ts(i))THEN
            Hila=0.0d0
         ELSEIF(cell%tl(i).lt.cell%ts(i)+dtl)THEN
            dHldtl(i)=H_il_min/dtl
            Hila=dHldtl(i)*(cell%tl(i)-cell%ts(i))
         ELSE
            Hila=H_il_min
         ENDIF
         dHldag(i)=10.d0*(hil_01-Hila)
         H_il(i)=Hila+dHldag(i)*ag
         Higa=H_ig_min
         dHgdag(i)=10.d0*(hig_01-Higa)
         H_ig(i)=Higa+dHgdag(i)*ag
!
      ELSEIF(ag.ge.0.9d0)THEN
!
!........Linear interpolation for gas temperature when alphag > 0.9
!
         IF(cell%tg(i).ge.cell%ts(i))THEN
            Higa=0.0d0
         ELSEIF(cell%tg(i).gt.cell%ts(i)-dtg)THEN
            dHgdtg(i)=H_ig_min/dtg
            Higa=dHgdtg(i)*(cell%ts(i)-cell%tg(i))
         ELSE
            Higa=H_ig_min
         ENDIF
!  
         dHgdag(i)=10.d0*(Higa-hig_09)
         H_ig(i)=Higa+dHgdag(i)*(ag-1.0d0)
         Hila=H_il_min
         dHldag(i)=10.d0*(Hila-hil_09)
         H_il(i)=Hila+dHldag(i)*(ag-1.0d0)
!
      ENDIF
      dHfgdtg(i)=dHgdtg(i)
!
      RETURN
      END SUBROUTINE int_htc_interp_temp2
!
!DEC$ ATTRIBUTES INLINE :: cub_interp2
      SUBROUTINE cub_interp2(x,x1,x2,y1,y2,dydx1,dydx2,c,dcdx)
!
!     Interploation using cubic polynomial
!
      IMPLICIT NONE
!      
      REAL(8) c,dcdx,x,x1,x2,y1,y2,dydx1,dydx2
      REAL(8) a1,a2,a3,a4,dx,dy,ddy,dxx,dxxx,xx1,xx,xxx
!
      dx=x2-x1
      dy=y2-y1
      ddy=dydx2-dydx1
      dxx=dx*dx
      dxxx=dxx*dx
      xx1=x1*x1
      xx=x*x
      xxx=xx*x
!
      a1=ddy/dxx+2.0d0*(dydx1/dxx-dy/dxxx)
      a2=(ddy/dx-3.0d0*(x2+x1)*a1)/2.0d0
      a3=dydx1-3.0d0*a1*xx1-2.0d0*a2*x1
      a4=y1-a1*xx1*x1-a2*xx1-a3*x1
!
      c=a1*xxx+a2*xx+a3*x+a4
      dcdx=3.0d0*xx*a1+2.0d0*x*a2+a3
!
      RETURN
      END SUBROUTINE cub_interp2
!



