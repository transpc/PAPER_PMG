!
      SUBROUTINE int_htc_simple_topology(ag)
!
!     This routine calculates the heat transfer coefficient at the interface 
!     when iHTC=1 and and udfl_calc_HTC_int_i=.false..
!
      USE VOL_DATA                  
      USE Zdaint_ag     , ONLY: daint1_ag,daint2_ag,daint1_ag_bc,daint2_ag_cm,  &
                                 D1_01,D1_bc,D1_09,aint_01b,aint_09b,aint_09d,   &
                                 aint_bc,aint_cm
      USE Zdhda         , ONLY: dHldag,dHgdag,dHfgdag
      USE Zflowregime   , ONLY: alphag_bc,alphag_cm
      USE Zmpi          , ONLY: ncell_fp
      USE Zqvol         , ONLY: H_ig,H_il,H_gf
      USE Zscalar_coeff , ONLY: l_min_hik
      USE Zvector       , ONLY: vrel_o
      USE Zzone         , ONLY: ncell_fluid
!
      IMPLICIT NONE
!
      INTEGER i
!      
      REAL(8) ag(ncell_fp)
      REAL(8) Nu,Re,Pr,hig1,hig2,hil1,hil2,delv
      REAL(8) H_il_1,H_il_2,H_ig_1,H_ig_2
      REAL(8) dHldag_1,dHldag_2,dHgdag_1,dHgdag_2
      REAL(8) hil_01,hil_bc,hil_09
      REAL(8) hig_01,hig_09
!
!DIR$ SIMD
      DO i=1,ncell_fluid
!
         delv=vrel_o(i)
!
!........For bubble
!
         Re=cell%rhol(i)*delv*cell%D1(i)/cell%lviscosl(i)
         Pr=cell%lviscosl(i)*cell%cpl(i)/cell%lcondl(i)
         Nu=2.0d0+0.6d0*Re**0.5d0*Pr**0.3d0
!
         hil1=cell%lcondl(i)*Nu
         hil_01=hil1/D1_01(i)
         hil_bc=hil1/D1_bc(i)
         hil_09=hil1/D1_09(i)
         hil1=hil1/cell%D1(i)
         hig1=1000.0d0
!
!........For drop
!
         Re=cell%rhog(i)*delv*cell%D2(i)/cell%lviscosg(i)
         Pr=cell%lviscosg(i)*cell%cpg(i)/cell%lcondg(i)
         Nu=2.0d0+0.6d0*Re**0.5d0*Pr**0.3d0
         hil2=cell%lcondg(i)*Nu/cell%D2(i)
         hig2=cell%lcondg(i)*Nu/cell%D2(i)
         hig2=1000.0d0
!
!........Limit value
!
         hig1=DMAX1(hig1,1.0d0)
         hig2=DMAX1(hig2,1.0d0)
         hil1=DMAX1(hil1,1.0d0)
         hil2=DMAX1(hil2,1.0d0)
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
!
!...........Cubic interpolation for the trasition region
!
            H_il_1=hil_bc*aint_bc(i)
            H_il_2=hil2*aint_cm(i)
            H_ig_1=hig1*aint_bc(i)
            H_ig_2=hig2*aint_cm(i)
            dHldag_1=hil_bc*daint1_ag_bc(i)
            dHldag_2=hil2*daint2_ag_cm(i)
            dHgdag_1=hig1*daint1_ag_bc(i)
            dHgdag_2=hig2*daint2_ag_cm(i)
            CALL cub_interp1(ag(i),alphag_bc,alphag_cm,H_il_1,H_il_2,dHldag_1,dHldag_2,H_il(i),dHldag(i))
            CALL cub_interp1(ag(i),alphag_bc,alphag_cm,H_ig_1,H_ig_2,dHgdag_1,dHgdag_2,H_ig(i),dHgdag(i))
         ENDIF
!
         hil_01=hil_01*aint_01b(i)
         hig_01=hig1*aint_01b(i)
         IF(cell%regime(i).eq.13)THEN
            hil_09=hil2*aint_09d(i)
            hig_09=hig2*aint_09d(i)
         ELSE
            hil_09=hil_09*aint_09b(i)
            hig_09=hig1*aint_09b(i)
         ENDIF
!
!........Linear interpolation of Hik for temperature when alphag<0.1 or alphag>0.9
!
         CALL int_htc_interp_temp1(i,ag(i),hil_01,hig_01,hil_09,hig_09)
!
!........Relaxation of Hik
!
!!!         H_il(i)=(1.0d0-relax_hik)*H_il(i)+relax_hik*hil_o(i)
!!!         H_ig(i)=(1.0d0-relax_hik)*H_ig(i)+relax_hik*hig_o(i)
!!!         hil_o(i)=H_il(i)
!!!         hig_o(i)=H_ig(i)
!         
         H_gf(i)=H_ig(i)
         dHfgdag(i)=dHgdag(i)
      ENDDO
!
      IF(l_min_hik.gt.0) CALL set_min_hik
!
      RETURN
      END SUBROUTINE int_htc_simple_topology
!
!
!DEC$ ATTRIBUTES NOINLINE :: set_min_hik 
      SUBROUTINE set_min_hik
!
!.....This routine sets the minimum value of Hik not to allow negative cell phase mass due to phase chage
!
      USE VOL_DATA                
      USE Zconst2     , ONLY: dt
      USE Zqvol       , ONLY: H_ig,H_il
      USE Zzone       , ONLY: ncell_fluid
!
      IMPLICIT NONE
!
      INTEGER i
!
      REAL(8) PsP,hi_liq,hi_gas,gg,gm,dalphag,agn,fr
!
      DO i=1,ncell_fluid
         PsP=cell%pps(i)/cell%p(i)
         gg=-(H_ig(i)*PsP*(cell%ts(i)-cell%tg(i))+H_il(i)*(cell%ts(i)-cell%tl(i)))
         IF(gg.ge.0.d0)THEN
            hi_gas=cell%hgsat(i)
            hi_liq=cell%hl(i)
         ELSE
            hi_gas=cell%hg(i)
            hi_liq=cell%hlsat(i)
         ENDIF
         gm=gg/(hi_gas-hi_liq)
         dalphag=gm*dt/cell%rhog(i)
         agn=cell%alphag(i)+dalphag
         IF(agn.lt.0.0d0)THEN
            fr=0.5d0*(dalphag-agn)/dalphag
            IF(cell%tg(i).lt.cell%ts(i)) H_ig(i)=fr*H_ig(i)
            IF(cell%tl(i).lt.cell%ts(i)) H_il(i)=fr*H_il(i)
         ENDIF
         IF(agn.gt.1.0d0)THEN
            fr=0.5d0*(dalphag-agn-1.0d0)/dalphag
            IF(cell%tg(i).gt.cell%ts(i)) H_ig(i)=fr*H_ig(i)
            IF(cell%tl(i).gt.cell%ts(i)) H_il(i)=fr*H_il(i)            
         ENDIF
      ENDDO
!
      RETURN
      END SUBROUTINE set_min_hik
!
!DEC$ ATTRIBUTES INLINE :: cub_interp1
      SUBROUTINE cub_interp1(x,x1,x2,y1,y2,dydx1,dydx2,c,dcdx)
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
      END SUBROUTINE cub_interp1
!
!DEC$ ATTRIBUTES INLINE :: int_htc_interp_temp1
      SUBROUTINE int_htc_interp_temp1(i,ag,hil_01,hig_01,hil_09,hig_09)
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
      END SUBROUTINE int_htc_interp_temp1
