!----------------------------------------------------------------------------------
!
      SUBROUTINE udfn_sg_heat_partition
!
      USE VOL_DATA   ,ONLY: cell
      USE Zzone      ,ONLY: ncell_fluid
      USE Zparam     ,ONLY: ndim
      USE Zcoord3    ,ONLY: volp
      USE Zqvol      ,ONLY: qporous_liq
      USE Zvector    ,ONLY: vl_o,ul_o
      USE Zsg        ,ONLY: mult_cell,mult_3d_cell1,mult_3d_cell2,f1_mult,f2_mult,  &
                            q_sd,t_wall,hyd_d,tl_sg,htcl_sec,htcb_sec,ts_sg,ht_area,do_tube,     &
                            igr,j1d,relax_hb,tune_hb,iboil_model
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER g,m,k,i,m1,m2
      LOGICAL EPRI78,LELE81,CONST
      REAL(8),PARAMETER::X_DB=0.75d0
      REAL(8)::pr,RE,gl
      REAL(8)::ps_sg,vp_sg
      REAL(8)::al_sg,qu_sg,vsl_sg,cpl_sg,cdl_sg,rl_sg,vl_sg(ndim),ul_sg,ql_sg
      REAL(8)::hsv,hsh,hsl,q_flux,hb,cb,htmp
!
      CONST=.false.
      EPRI78=.false.
      LELE81=.false.
!
      IF(iboil_model.eq.1)THEN
         EPRI78=.true.
      ELSEIF(iboil_model.eq.2)THEN
         LELE81=.true.
      ELSE
         CONST=.true.
      ENDIF
!
      DO i=1,ncell_fluid
         g=igr(i)
         IF(g.eq.0) CYCLE
         m=j1d(i)
!
         k=mult_cell(g,m)
         IF(k.eq.0)THEN
!
            ts_sg(g,m)=cell%ts(i)
            ps_sg=cell%p(i)
            vp_sg=volp(i)
!
            tl_sg(g,m)=cell%tl(i)
            al_sg=cell%alphal(i)
            qu_sg=cell%quals(i)
            vsl_sg=cell%lviscosl(i)
            cpl_sg=cell%cpl(i)
            cdl_sg=cell%lcondl(i)
            rl_sg=cell%rhol(i)
            ul_sg=ul_o(i)
            vl_sg(:)=vl_o(i,:)
            ql_sg=qporous_liq(i)
!
         ELSE
            m1=mult_3d_cell1(k)
            m2=mult_3d_cell2(k)
!
            ts_sg(g,m)=f1_mult(k)*cell%ts(m1)+f2_mult(k)*cell%ts(m2)
            ps_sg=f1_mult(k)*cell%p(m1)+f2_mult(k)*cell%p(m2)
            vp_sg=volp(m1)+volp(m2)
!
            tl_sg(g,m)=f1_mult(k)*cell%tl(m1)+f2_mult(k)*cell%tl(m2)
            al_sg=f1_mult(k)*cell%alphal(m1)+f2_mult(k)*cell%alphal(m2)
            qu_sg=f1_mult(k)*cell%quals(m1)+f2_mult(k)*cell%quals(m2)
            vsl_sg=f1_mult(k)*cell%lviscosl(m1)+f2_mult(k)*cell%lviscosl(m2)
            cpl_sg=f1_mult(k)*cell%cpl(m1)+f2_mult(k)*cell%cpl(m2)
            cdl_sg=f1_mult(k)*cell%lcondl(m1)+f2_mult(k)*cell%lcondl(m2)
            rl_sg=f1_mult(k)*cell%rhol(m1)+f2_mult(k)*cell%rhol(m2)
            ul_sg=f1_mult(k)*ul_o(m1)+f2_mult(k)*ul_o(m2)
            vl_sg(:)=f1_mult(k)*vl_o(m1,:)+f2_mult(k)*vl_o(m2,:)
            ql_sg=qporous_liq(m1)+qporous_liq(m2)
!
         ENDIF
!
!........Single phase liquid heat transfer coefficient (ATHOS3: EPRI78 Option)
!
         IF(ps_sg.le.5.5158d6)THEN
            cb=2.6399d0*(1.554d0+ps_sg*(1.2425d-7+5.6798d-14*ps_sg))
         ELSE
            cb=2.6399d0*(-4.557d0+ps_sg*(2.273d-6-1.2584d-13*ps_sg))
         ENDIF
!
         pr=vsl_sg*cpl_sg/cdl_sg
         gl=rl_sg*al_sg*ul_sg
         re=gl*hyd_d(1)/vsl_sg
         hsv=0.03d0*re**0.8d0*pr**0.333d0*cdl_sg/hyd_d(1)
         hsh=0.29d0*re**0.6d0*pr**0.333d0*cdl_sg/do_tube
         hsl=DSQRT(hsv*hsv+hsh*hsh)
!
!........Biling heat transfer coefficient (ATHOS3: EPRI78 Option)
!
         IF(EPRI78)THEN
            q_flux=ql_sg/ht_area(g,m)
            hb=cb*q_flux**0.6667d0
            IF(t_wall(g,m).le.tl_sg(g,m))THEN
               htcl_sec(i)=0.0d0
               htcb_sec(i)=0.0d0
            ELSE
               htmp=hsl*(t_wall(g,m)-tl_sg(g,m))/(t_wall(g,m)-ts_sg(g,m))
               IF(hb.gt.htmp)THEN
                  htcb_sec(i)=hb
                  htcl_sec(i)=0.0d0
               ELSE
                  htcb_sec(i)=0.0d0
                  htcl_sec(i)=hsl
               ENDIF
            ENDIF
         ENDIF
!
!........Biling heat transfer coefficient (ATHOS3: LELE81 Option)
!
         IF(LELE81)THEN
            IF(t_wall(g,m).ge.tl_sg(g,m))THEN
               htcl_sec(i)=hsl
            ELSE
               htcl_sec(i)=0.0d0
            ENDIF
            IF(t_wall(g,m).ge.ts_sg(g,m))THEN
               hb=1792.653d0*DEXP(ps_sg/4342590.0d0)*tune_hb                                   ! "tune_hb" has to be tunned to get a design recirculation flow ratio
               htcb_sec(i)=relax_hb*htcb_sec(i)+(1.0d0-relax_hb)*hb*(t_wall(g,m)-ts_sg(g,m))   ! Relaxation in time
            ELSE
               htcb_sec(i)=0.0d0
            ENDIF
         ENDIF
!
!........Spline at high void fraction
!
         IF(qu_sg.gt.X_DB)THEN
            htcl_sec(i)=htcl_sec(i)*(1.0d0-qu_sg)/(1.0d0-X_DB)
            htcb_sec(i)=htcb_sec(i)*(1.0d0-qu_sg)/(1.0d0-X_DB)
         ENDIF
!
!........Constant heat transfer coefficient for a test
!
         IF(CONST)THEN
            htcb_sec(i)=0.0d0
            htcl_sec(i)=al_sg*35000.d0
         ENDIF
!
         q_sd(g,m)=(htcl_sec(i)*(t_wall(g,m)-tl_sg(g,m))+htcb_sec(i)*(t_wall(g,m)-ts_sg(g,m)))*ht_area(g,m)
!
         IF(k.eq.0)THEN
            qporous_liq(i)=q_sd(g,m)
         ELSE
            qporous_liq(m1)=f1_mult(k)*q_sd(g,m)
            qporous_liq(m2)=f2_mult(k)*q_sd(g,m)
         ENDIF
!
      ENDDO
!
      END SUBROUTINE udfn_sg_heat_partition
!----------------------------------------------------------------------------------
