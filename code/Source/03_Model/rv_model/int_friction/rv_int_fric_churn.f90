      SUBROUTINE bubblyslug_annular_for_churn(flag)
!      
      USE VOL_DATA       , ONLY: cell
      USE Zconst2        , ONLY: hydraulicd
      USE Zparam         , ONLY: pi,ndim
      USE Zvector        , ONLY: vg_o,vl_o
      USE STM_TBL_cupid  , ONLY: pcrit
      USE Zrv_int_friction , ONLY: red_min,reb_min,dh,drift,slug,annular,alphag,alphag,       &
                                   drift_model,alphade,alphal,epri,bestion,wf_bestion,pres,   &
                                   wf_epri,sideward,fore,upward,asali,rhog,rhol,mug,mul,      &
                                   a1,b1,b2,c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,d2,             &
                                   gravity,vgj,ci,sigma,vr_drift,vrd,cp,L,reg,rel,jg,         &
                                   re,k0,jl,rhogoverl,r,d1,fdebug,rhol_g,vgj0,                &
                                   vfgl_drift,mdrag_drift,bubbly,vr,xa,dcon,dht,dbub,mum,     &
                                   reb,cd,ia_bubbly,vfgl_bubbly,iam,alphags,alphabs,alphab,   &
                                   ia_slug_tb,lcell,vfgl_slug_tb,vfgl_slug_sb,vfgl_slug,      &
                                   alphasa,ia_slug_sb,drop_max,ia_annular_drp,ia_annular_ann, &
                                   vfgl_annular,dgas,vg,alphalf,f1,f2,f3,fi,cann,             &
                                   vfgl_annular_ann,alphald,ddrp,dh,red,vfgl_annular_drp,     &
                                   i
      USE Zrv_model      , ONLY: rv_fric_i      
      USE Zwall_HTC      , ONLY: reflod,inline_bundle
!
      IMPLICIT NONE
!
      LOGICAL, SAVE :: mars_coding_only      
!     
      INTEGER ix 
      INTEGER flag 
      INTEGER drift_store,slug_store,annular_store
      REAL(8) alphag_store
      REAL(8) alphag_save,alphal_save
      REAL(8) c0_epri,c1_epri,ci_epri,vfgl_epri,mdrag_epri
      REAL(8) c0_bestion,c1_bestion,ci_bestion,vfgl_bestion,mdrag_bestion
      REAL(8) gamma_Star,alphaef,alphaad,G_star,Ce,sigma_star,v_crit,lambda,F11,gamma,alphaad_star !annular-alphalf
      REAL(8) vfg_star,vfg_star2,vfg2,vfg      !annular-vr2
!     REAL(8) alpha_bub !bubbly
      REAL(8) F9 !slug
      REAL(8) c5_under,c6_under,c7_under,c8_under  !LSJ      
      REAL(8) velo
      REAL(8) alpan,alpanf,alphac,alphad,alphef
      REAL(8) velo_g,velo_l,rhofg,vcritl,rvcrtj,centr,scrach,fluxm1,vfgbb,voidb,rdstar,alpdrp,rdiam,voidx,surfa,rey !,xxa
      REAL(8) fic,void,voidgs,sqvoid,sannu,fgf1,fgf2,fgf
      REAL(8) velg0,velf0,alpg,alpf,reynog,reynof,reyno,arg,rhorat,xk0,xr,xln,xld,xl,rfg,denom,jsubgs,jccfl,jsubfs,fract
      REAL(8) alphagi
!
!.....Choose MARS manual model/MARS code model
!
      IF(rv_fric_i.eq.1) THEN
         mars_coding_only=.FALSE.         !MARS manual model
      ELSEIF(rv_fric_i.eq.2)THEN
         mars_coding_only=.TRUE.          !MARS code model
      ENDIF 
!
      red_min=1.0d0
      reb_min=1.0d0
!
      dh=hydraulicd(i)       
!      
!.....store the original value
!     
      IF(flag.eq.1)THEN
         drift_store=drift
         slug_store=slug
         annular_store=annular
         alphag_store=alphag
         slug=1
         annular=1
         drift=1
!          
!........select drift flux model when drift_model is 1.
!
         IF(drift_model.eq.1)THEN
            slug=0
         ELSE
            drift=0
         ENDIF
      ENDIF
!
!.....set alphag to left-end value
!      
      IF(flag.eq.1)THEN
          alphag=alphade 
          alphal=1.0d0-alphag
      ENDIF    
!  
!.....drift flux model
!
      IF(drift.eq.1)THEN
         epri=0
         bestion=0
!         
!         IF(reflod)THEN
         IF(reflod.and.inline_bundle.eq.1)THEN
            CALL rv_int_fric_weight_interval(pres,10.0d5,20.0d5,wf_bestion,wf_epri)
         ELSE            
         wf_epri=1.0d0
         wf_bestion=0.0d0
         ENDIF
!         
         IF(wf_epri.gt.0.0d0.and.wf_epri.le.1.0d0)THEN
            epri=1
         ELSEIF(wf_bestion.gt.0.0d0.and.wf_bestion.le.1.0d0)THEN
            bestion=1 
         ENDIF    
         IF(upward.eq.1.or.sideward.eq.1)THEN
            fore=1
            asali=0
         ELSE
            fore=0
            asali=1
         ENDIF
         c0_epri=0.0d0
         c1_epri=0.0d0
         ci_epri=0.0d0
         vfgl_epri=0.0d0
         mdrag_epri=0.0d0
         c0_bestion=0.0d0
         c1_bestion=0.0d0
         ci_bestion=0.0d0
         vfgl_bestion=0.0d0
         mdrag_bestion=0.0d0
         alphag_save=alphag
         alphal_save=alphal
         !next alphag=DMAX1(1.0d-2,alphag)
         !next alphal=DMAX1(1.0d-2,1.0d0-alphag)
         alphag=DMAX1(2.0d-2,alphag)
         alphal=DMAX1(2.0d-2,1.0d0-alphag)
      ELSE
         epri=0
         bestion=0
      ENDIF          
!
!.....epri drift flux model for bubbly,slug
!
      IF(epri.eq.1)THEN
         IF(mars_coding_only) THEN
!
            velg0=DABS(vg_o(i,ndim))
            velf0=DABS(vl_o(i,ndim))
!
            alpg=DMAX1(1.d-2,alphag)         
            alpf=1.d0-alpg
!         
            reynog=rhog*alpg*velg0*dh/mug 
            reynof=rhol*alpf*velf0*dh/mul
            IF(reynog.gt.reynof.or.reynog.lt.0.0d0)then 
               reyno=reynog 
            ELSE 
               reyno=reynof 
            ENDIF 
!            arg=dmax1(-170.0d0,dmin1(170.0d0,-reyno/6.0d4))    !alternative coefficient but not effective
            arg=DMAX1(-85.0d0,DMIN1(85.0d0,-reyno/60000.0d0))
            a1=1.0d0/(1.0d0+dexp(arg)) 
            b1=dmin1(0.8d0,a1) 
            rhorat=rhog/rhol
!         
!...........distribution coefficient
!
!.....pcrit is constant read in stread.f90 why change?
!           pcrit=22.4d6
            xk0=b1+(1.0d0-b1)*rhorat**0.25 
            xr=(1.0d0+1.57d0*rhorat)/(1.0d0-b1) 
            c1=4.0d0*pcrit*pcrit/(cell%p(i)*(pcrit-cell%p(i))) 
            c1=dabs(c1)
            c1=DMAX1(1.d-5,c1) !Added by LSJ                
!
!            IF(c1*alpg.lt.170.0d0)then 
            IF(c1.lt.170.0d0)then 
               xln=1.0d0-dexp(-c1*alpg) 
               xln=dmax1(1.d-5,xln)
            ELSE 
               xln=1.0d0 
            ENDIF 
            IF(c1.lt.170.d0)then 
               xld=1.0d0-dexp(-c1) 
               xld=dmax1(1.d-5,xld)
            ELSE 
               xld=1.0d0 
            ENDIF 
            xl=xln/xld 
            c0=xl/(xk0+(1.0d0-xk0)*alpg**xr)
            c0=dmax1(0.d0,dmin1(1.d0,c0))     !This is MARS manual only, but makes big difference.
!
!...........drift velocity
!
            c7=(0.09144d0/dh)**0.6 
            IF(c7.ge.1.0d0)then 
               c4=1.0d0 
            ELSE 
               c7_under=1.d0-c7
!               IF(c7_under.lt.1.d-5.and.c7_under.ge.0.d0) c7_under=1.d-5
!               IF(c7_under.gt.-1.d-5.and.c7_under.lt.0.d0) c7_under=-1.d-5         
               IF(DABS(c7_under).lt.1.d-5) c7_under=1.d-5
               c8=c7/c7_under
!               
               c8_under=1.0d0-DEXP(-c8)
!               IF(c8_under.lt.1.d-5.and.c8_under.ge.0.d0) c8_under=1.d-5
!               IF(c8_under.gt.-1.d-5.and.c8_under.lt.0.d0) c8_under=-1.d-5            
               IF(DABS(c8_under).lt.1.d-5) c8_under=1.d-5
               c4=1.0d0/c8_under               
            ENDIF 
            c5=dsqrt(150.0d0*rhorat) 
            IF(c5.ge.1.0d0)then 
               c2=1.0d0 
            ELSE 
               c5_under=1.d0-c5
!               IF(c5_under.lt.1.d-5.and.c5_under.ge.0.d0) c5_under=1.d-5
!               IF(c5_under.gt.-1.d-5.and.c5_under.lt.0.d0) c5_under=-1.d-5
               IF(DABS(c5_under).lt.1.d-5) c5_under=1.d-5 
               c6=c5/c5_under               
!               
               c6_under=1.0d0-DEXP(-c6)
!               IF(c6_under.lt.1.d-5.and.c6_under.ge.0.d0) c6_under=1.d-5
!               IF(c6_under.gt.-1.d-5.and.c6_under.lt.0.d0) c6_under=-1.d-5
               IF(DABS(c6_under).lt.1.d-5) c6_under=1.d-5
               c2=1.0d0/c6_under
            
            ENDIF 
            IF(reynog.ge.0.0d0)then 
               c9=(1.d0-alpg)**b1 
            ELSE 
               c9=dmin1(0.7d0,(1.d0-alpg)**0.65d0) 
            ENDIF 
!
! actually, velocity can be negative accroding to MARS manual, but CUPID takes jg, jl as positive only. 
! thus, the following up-/down-flow should be checked later.
!
!           upflow
!            IF(velf0.ge.0.0d0.and.velg0.ge.0.0d0)then 
!            IF(jg.ge.0.d0.and.jl.ge.0.d0)THEN         
            IF(vg_o(i,ndim).ge.0.d0.and.vl_o(i,ndim).ge.0.d0)THEN
               c3=dmax1(0.5d0,2.0d0*dexp(-dabs(reynof)/6.0d4)) 
!
!           downflow
!            ELSEIF(velf0.le.0.0d0.and.velg0.le.0.0d0)then 
            ELSEIF(vg_o(i,ndim).le.0.d0.and.vl_o(i,ndim).le.0.d0)THEN
!               c10=2.0d0*(dexp((dabs(reynof)/3.5d5)**0.40d0))-1.75d0*(          &
!                   dabs(reynof))**0.03*dexp(-(dabs(reynof)/5.0d4)*              &
!                   (0.0381d0/dh)**2.0d0)+(0.0381d0/dh)**0.25*(dabs(reynof))**0.001
               c10=2.0d0*(dexp((dabs(reynof))**0.40d0/3.5d5))-1.75d0*(          &
                   dabs(reynof))**0.03*dexp(-(dabs(reynof)/5.0d4)*              &
                   (0.0381d0/dh)**2.0d0)+(0.0381d0/dh)**0.25*(dabs(reynof))**0.001
               b2=1.0d0/(1.0d0+0.05d0*dabs(reynof)/3.5d5)**0.4
               c3=2.0d0*(c10/2.0d0)**b2 
!
!           countercurrent flow
            ELSE 
!               c10=2.0d0*(dexp((dabs(reynof)/3.5d5)**0.40d0))-1.75d0*(          &
!                   dabs(reynof))**0.03*dexp(-(dabs(reynof)/5.0d4)*              &
!                   (0.0381d0/dh)**2.0d0)+(0.0381d0/dh)**0.25d0*(dabs(reynof))**0.001
               c10=2.0d0*(dexp((dabs(reynof)/3.5d5)**0.40d0))-1.75d0*(          &
                   dabs(reynof))**0.03*dexp(-(dabs(reynof)/5.0d4)*              &
                   (0.0381d0/dh)**2.0d0)+(0.0381d0/dh)**0.25d0*(dabs(reynof))**0.001 
               b2=1.0d0/(1.0d0+0.05d0*dabs(reynof)/3.5d5)**0.4
               c3=2.0d0*(c10/2.0d0)**b2 
!
               rfg=dmax1(1.0d-05,(rhol-rhog)) 
               denom=dsqrt(gravity*rfg*dh) 
               jsubgs=alpg*velg0*dsqrt(rhog)/denom 
               jsubgs=dmax1(0.0d0,dmin1(0.55d0,jsubgs)) 
               jccfl=(0.775d0-dsqrt(jsubgs))**2 
               jsubfs=alpf*dabs(velf0)*dsqrt(rhol)/denom 
               fract=dmin1(1.0d0,(jsubfs/jccfl))
!
               c3=c3*fract+2.0d0*(1.0d0-fract) 
            ENDIF        
!
!...........all flows
!
            rfg=dmax1((rhol-rhog),1.0d-05) 
            vgj=1.41d0*(rfg*sigma*gravity/(rhol*rhol))**0.25*c2*c3*c4*c9
            vgj=dmax1(vgj,1.0d-05) 
!
!            ci=alpg*alpf**3*rhorat*gravity/(DABS(vgj)*vgj)   !LSJ modification
            ci=alpg*alpf**3*rfg*gravity/(DABS(vgj)*vgj)   !LSJ modification
!            c1=(1.0d0-c0*alphag)/DMAX1(1.d-9,alpf)
            c1=(1.0d0-c0*alpg)/DMAX1(1.d-9,alpf)
            vr_drift(:)=c1*vg_o(i,:)-c0*vl_o(i,:)
            vrd=vr_drift(1)*vr_drift(1)+vr_drift(2)*vr_drift(2)
            IF(ndim.eq.3)vrd=vrd+vr_drift(ndim)*vr_drift(ndim)
            c0_epri=c0
            c1_epri=c1
            ci_epri=ci
            vfgl_epri=ci*vrd
            mdrag_epri=ci*vrd*vrd
      
         ELSE
      
!
!...........cp, L
!
            cp=DABS(4.0d0*pcrit**2.0d0/(cell%p(i)*DABS(pcrit-cell%p(i))))
            cp=DMAX1(1.d-5,cp) !Added by LSJ         
!         
            IF(cp.lt.170.0d0)THEN
               L=(1.0d0-DEXP(-cp*alphag))/(1.0d0-DEXP(-cp))
            ELSE
               L=1.0d0
            ENDIF
!
!...........Re
!
            reg=rhog*jg*dh/mug
            rel=rhol*jl*dh/mul
!         
            IF(reg.gt.rel .or. reg.lt.0.0d0)THEN
               re=reg
            ELSE
               re=rel
            ENDIF
!
!...........k0, r & b1
!
!            a1=1.0d0+DEXP(DMAX1(-85.0d0,DMIN1(-85.0d0,-re/60000.0d0)))
            a1=1.0d0+DEXP(DMAX1(-85.0d0,DMIN1(85.0d0,-re/60000.0d0)))    !LSJ modification
            a1=1.0d0/a1
            b1=DMIN1(0.8d0,a1)
            k0=b1+(1.d0-b1)*(rhogoverl)**(1.0d0/4.0d0)
            r=(1.0d0+1.57d0*(rhogoverl)/(1.0d0-b1))
!
!...........c0
!
!            IF(reg.ge.0.0d0)THEN
!               c0=L/(k0+(1.0d0-k0)*alphag**r)
!            ELSE
!               c0=DMAX1(L/(k0+(1.0d0-k0)*alphag**r),vgj0*(1.0d0-alphag)**0.2d0/(DABS(jg)+DABS(jl)))
!            ENDIF
!            c0=DMAX1(0.0d0,DMIN1(1.0d0,c0))
!
!...........c9
!
            IF(reg.ge.0.0d0)THEN
               c9=(1.0d0-alphag)**b1
            ELSE
               c9=DMIN1(0.7d0,(1.0d0-alphag)**0.65d0)
            ENDIF
!
!...........c2,c4
!
            c5=(150.0d0*(rhogoverl))**(1.0d0/2.0d0)
            c5_under=1.d0-c5
            IF(DABS(1.d0-c5).lt.1.d-5) c5_under=1.d-5
            c6=c5/c5_under
!         
            d2=0.09144d0
            c7=(d2/dh)**0.6d0
            c7_under=1.d0-c7
            IF(DABS(1.d0-c7).lt.1.d-5) c7_under=1.d-5
            c8=c7/c7_under
!         
            IF(c5.ge.1.0d0)THEN
               c2=1.0d0
            ELSE
               c6_under=1.0d0-DEXP(-c6)
               IF(DABS(1.0d0-DEXP(-c6)).lt.1.d-5) c6_under=1.d-5
               c2=1.0d0/c6_under
            ENDIF
!   
            IF(c7.ge.1.0d0)THEN
               c4=1.0d0
            ELSE
               c8_under=1.0d0-DEXP(-c8)
               IF(DABS(1.0d0-DEXP(-c8)).lt.1.d-5) c8_under=1.d-5
               c4=1.0d0/c8_under
            ENDIF
!
!...........c3
!
            IF(fore.eq.1)THEN !upward
               c3=DMAX1(0.5d0,2.0d0*DEXP(-DABS(rel)/60000.0d0))
            ELSE              !downward & countercurrent
               d1=0.0381d0
               b2=1.0d0/(1.0d0+0.05d0*DABS(rel/350000.0d0))**0.4d0
!               c10=2.0d0*DEXP(DABS(rel/350000.0d0)**0.40d0) &
!                  -1.75d0*DABS(rel)**0.03d0*(-DABS(rel)/50000.0d0*(d1/dh)**2.0d0)*DABS(rel)**0.001d0 &
!                  +(d1/dh)**0.25d0*DABS(rel)**0.001d0
               c10=2.0d0*DEXP(DABS(rel)**0.4/350000.0d0) &
                  -1.75d0*DABS(rel)**0.03*DEXP(-DABS(rel)/50000.0d0*d1*d1/(dh*dh))  &
                  +(d1/dh)**0.25*DABS(rel)**0.001                
               c3=2.0d0*(c10/2.0d0)**b2
            ENDIF
!
!...........vgj
!
!            vgj=DMAX1(1.d-5,rhol_g)*sigma*gravity/rhol**2.0d0
!            fdebug=vgj
!!            c3=DMIN1(c3,1.0d0) !next,critical
!            vgj=1.41d0*fdebug**(1.0d0/4.0d0)*c2*c3*c4*c9
!!            vgj=DMAX1(vgj,1.0d-1) !next,critical
            fdebug=DMAX1(1.d-5,rhol_g)*sigma*gravity/(rhol*rhol)
            vgj0=1.41d0*fdebug**0.25*c2*c3*c4
            vgj=vgj0*c9
            vgj=DMAX1(vgj,1.0d-1) !next,critical
!
!...........c0
!
            IF(reg.ge.0.0d0)THEN
               c0=L/(k0+(1.0d0-k0)*alphag**r)
            ELSE
               c0=DMAX1(L/(k0+(1.0d0-k0)*alphag**r),vgj0*(1.0d0-alphag)**0.2d0/(DABS(jg)+DABS(jl)))
            ENDIF
            c0=DMAX1(0.0d0,DMIN1(1.0d0,c0))           
!
!...........vfgl, vflg
!
!            phij=pi/2.0d0 !next
!            ci=alphag*alphal**3.0d0*rholoverg*gravity*DSIN(phij)/(DABS(vgj)*vgj)
            ci=alphag*alphal**3*rhol_g*gravity/(DABS(vgj)*vgj)   !LSJ modification
            c1=(1.0d0-c0*alphag)/DMAX1(1.d-9,alphal)
            vr_drift(:)=c1*vg_o(i,:)-c0*vl_o(i,:)
            vrd=vr_drift(1)*vr_drift(1)+vr_drift(2)*vr_drift(2)
            IF(ndim.eq.3)vrd=vrd+vr_drift(ndim)*vr_drift(ndim)
            c0_epri=c0
            c1_epri=c1
            ci_epri=ci
            vfgl_epri=ci*vrd
            mdrag_epri=ci*vrd*vrd
!
         ENDIF         
!         
      ENDIF 
!
!.....bestion drift flux model for bubbly,slug
!
      IF(bestion.eq.1)THEN 
!         c0=1.0d0
!!         vgj=0.124d0*(DABS(gravity*rhol_g*dh/rhog))**0.50d0
!         vgj=0.188d0*(DABS(gravity*rhol_g*dh/rhog))**0.5       !LSJ modification
!!         phij=pi/2.0d0 !next
!!         ci=alphag*alphal**3.0d0*rholoverg*gravity*DSIN(phij)/(DABS(vgj)*vgj)
!         ci=alphag*alphal**3*rhol_g*gravity/(DABS(vgj)*vgj)    !LSJ modification
!         c1=(1.0d0-c0*alphag)/alphal
!
         c0=1.2d0  ! yjm MARS source c0bes
         ci=65.0d0*DMAX1(alphag,1d-15)*DMAX1(alphal,1.d-15)**3.0d0*cell%rhog(i)/dh
         IF(alphag.gt.0.0d0) c0=DMIN1(c0,1.0d0/alphag)
         c1=(1.0d0-c0*alphag)/(1.0d0-alphag)         
!         
         vr_drift(:)=c1*vg_o(i,:)-c0*vl_o(i,:)
         vrd=vr_drift(1)*vr_drift(1)+vr_drift(2)*vr_drift(2)
         IF(ndim.eq.3)vrd=vrd+vr_drift(ndim)*vr_drift(ndim)
         c0_bestion=c0
         c1_bestion=c1
         ci_bestion=ci
         vfgl_bestion=ci*vrd
         mdrag_bestion=ci*vrd*vrd         
      ENDIF
!   
      IF(drift.eq.1)THEN
          c0=c0_epri*wf_epri+c0_bestion*wf_bestion
          c1=c1_epri*wf_epri+c1_bestion*wf_bestion
          ci=ci_epri*wf_epri+ci_bestion*wf_bestion
          vfgl_drift=vfgl_epri*wf_epri+vfgl_bestion*wf_bestion
          mdrag_drift=mdrag_epri*wf_epri+mdrag_bestion*wf_bestion
         alphag=alphag_save
         alphal=alphal_save          
      ENDIF
!
!.....drag coefficient model for bubbly
!
      IF(bubbly.eq.1)THEN
         vfg=vr
         vfg2=DMAX1(vfg*vfg,xa(1)/rhol*DMIN1(dcon(1)*alphag**0.333333333d0,dht))
         dbub=xa(1)/(rhol*vfg2)
         dbub=DMAX1(dbub,84.0d-6)
!                    
         mum=mul/DMAX1(1.d-9,1.0d0-alphag)
!         reb=(1.0d0-alphag)*rhol*dbub*vr/mum
         reb=rhol*dbub*vr/mum  ! yjm mod
         reb=DMAX1(reb,1.d-10)
         cd=24.0d0/reb*(1.0d0+0.1d0*reb**0.75d0)
         ia_bubbly=3.6d0*DMAX1(alphag,1.0d-5)/dbub         
         IF(iam.eq.1)ia_bubbly=cell%ia_bubbly(i)
         IF(iam.eq.0)cell%ia_bubbly(i)=ia_bubbly
!         
         vfgl_bubbly=1.0d0/8.0d0*rhol*cd*ia_bubbly*vr
      ENDIF      
!
!.....drag coefficient model for slug
!
      IF(slug.eq.1)THEN
         alphags=DMIN1(1.0d0,DMAX1(0.0d0,(alphag-alphabs)/(alphasa-alphabs)))
         F9=DEXP(-8.0d0*alphags)
         alphags=alphabs*F9
         alphab=(alphag-alphags)/(1.0d0-alphags)
         alphab=DMIN1(1.0d0,DMAX1(0.0d0,alphab))
         cd=10.9d0*alphab**0.5d0*(1.0d0-alphab)**3.0d0
         ia_slug_tb=alphab/lcell
         IF(iam.eq.1)ia_slug_tb=cell%ia_slug_tb(i)
         IF(iam.eq.0)cell%ia_slug_tb(i)=ia_slug_tb
         vfgl_slug_tb=0.5d0*rhol*cd*ia_slug_tb*vr
!         vfg=vr
! small bubble
         vfg=vr*F9*F9
!         
!         vfg2=DMAX1(vfg*vfg,xa(1)/rhol*DMIN1(dcon(1)*alphag**0.333333333d0,dht))
         vfg2=DMAX1(vfg*vfg,xa(1)/rhol*DMIN1(dcon(1)*alphags**0.333333333d0,dh))
!         
         dbub=xa(1)/(rhol*vfg2)
         dbub=DMAX1(dbub,84.0d-6)
!                    
!         mum=mul/(1.0d0-alphag) !next
         mum=mul
!         
         reb=(1.0d0-alphags)*rhol*dbub*vr/mum
         reb=DMAX1(1.0d0,reb) !next
         cd=24.0d0/reb*(1.0d0+0.1d0*reb**0.75d0)
!         ia_slug_sb=3.6d0*alphags/dbub*(1.0d0-alphab)
         ia_slug_sb=3.6d0*alphags/dbub*(1.0d0-alphab)*f9
         
         IF(iam.eq.1)ia_slug_sb=cell%ia_slug_sb(i)
         IF(iam.eq.0)cell%ia_slug_sb(i)=ia_slug_sb
         vfgl_slug_sb=1.0d0/8.0d0*rhol*cd*ia_slug_sb*vr
!   
         vfgl_slug=vfgl_slug_tb+vfgl_slug_sb
      ENDIF
!
!.....set alphag to right-end value
!      
      IF(flag.eq.1)THEN
         alphag=alphasa
         alphal=1.0d0-alphag
      ENDIF         
!
!.....drag coefficient model for annular
!
      IF(annular.eq.1)THEN
         IF(mars_coding_only) THEN      
            alpan=alphag
            alpanf=1.d0-alpan
            alphac=cell%alpha_sa(i)
            alphad=1.0d-4             
            alphef=DMAX1(2.0d0*alphad,DMIN1(2.0d-3*rhog/rhol,2.0d-4))
!          
            velo_g=DABS(vg_o(i,ndim))
            velo_l=DABS(vl_o(i,ndim))
!               
            rhofg=rhol_g
            vcritl=5.66d0*DSQRT(DSQRT(DMAX1(sigma,1.0d-7)*rhofg)/rhog)
            rvcrtj=alphag*velo_g/vcritl 
            centr=7.5d-5
            scrach=centr*rvcrtj*rvcrtj*rvcrtj*rvcrtj*rvcrtj*rvcrtj  
            fluxm1=1.0d-4*DSQRT(DSQRT(DABS(velo_l*alpanf*rhol*dh/mul)))          
            IF(scrach.le.200.d0) THEN
               vfgbb=DEXP(-scrach)*DMAX1(0.0d0,(1.0d0-fluxm1))
            ELSE
               vfgbb=0.d0
            ENDIF 
!              
!            IF(alphag.le.alphac.or.alpanf.ge.alphef) THEN
            IF(alphag.gt.alphac.and.alpanf.lt.alphef) THEN
                scrach=(alpanf-1.0d-7)/(alphef-1.0d-7) 
                alphad=alphad*scrach+1.0d-5*(1.0d0-scrach) 
                vfgbb=vfgbb*scrach 
            ENDIF             
!            
            voidb=DMIN1(0.999999999d0,DMAX1(0.d0,alpanf*vfgbb)) 
            rdstar=DSQRT(DSQRT(voidb)) 
!            alpdrp=DMAX1((alpanf-voidb)/DMAX1(1.d-5,(1.0d0-voidb)),alphad)
            alpdrp=DMAX1((alpanf-voidb)/(1.0d0-voidb),alphad)
            vfg=(velo_g-velo_l)*(1.0d0-vfgbb)
            IF(alphal.lt.1.0d-6) vfg=vfg*alphal*1.0d6
!
!           BubbleDropDrag: fic            
            vfg2=dmax1(vfg*vfg,xa(2)/(rhog*dmin1(dcon(2)*alpdrp**0.333333333d0,dh)))
!
!            rdiam=dmin1(rhog*vfg2/xa(2),1.0d0/84.0d-6)
            rdiam=DMAX1(1.0d0/dh,1.0d0/drop_max,DMIN1(rhog*vfg2/xa(2),1.0d0/84.0d-6))
!            
            voidx=1.0d0-alpdrp
            surfa=10.8d0*alpdrp*rdiam
            rey=xa(2)*voidx**3/(mug*dsqrt(vfg2*voidx)) 
            fic=rhog*surfa*(1.0d0+0.1d0*rey**0.75d0)/rey
!
!           sannu, void
            void=1.d0-voidb    
            voidgs=dsqrt(void) 
            sqvoid=dsqrt(alpan) 
            reg=dmax1(0.01d0,dabs(velo_g-velo_l))
            reg=rhog*reg*sqvoid*dh/mug
            sannu=6.1192775d0*sqrt(rdstar)*voidgs/dh
!
!           fg: Cd
            fgf1=64.0d0/reg 
            fgf2=0.02d0*(1.d0+150.0d0*(1.0d0-sqvoid)) 
            IF(reg.le.500.0d0)then 
               fgf=fgf1 
            ELSEIF(reg.lt.1500.0d0)then 
               fgf=fgf1*(1500.0d0-reg)*0.001d0+fgf2*(reg-500.0d0)*.001d0
            ELSE 
               fgf=fgf2 
            ENDIF                       
!            
!            fic=0.125d0*rhog*sannu*fgf*fic*void
            fic=0.125d0*rhog*sannu*fgf+fic*void  ! yjm
!
            ia_annular_drp=alpdrp
            ia_annular_ann=alpdrp
            IF(iam.eq.1) ia_annular_drp=cell%ia_annular_drp(i)
            IF(iam.eq.0) cell%ia_annular_drp(i)=ia_annular_drp
            IF(iam.eq.1)ia_annular_ann=cell%ia_annular_ann(i)
            IF(iam.eq.0)cell%ia_annular_ann(i)=ia_annular_ann            
!         
            vfgl_annular=fic*vr
!         
         ELSE  
!                  
!...........alphalf 
!   
            gamma_star=1.0d0
            alphaad=1.d-4
            alphaef=DMAX1(2.0d0*alphaad,DMIN1(2.0d-3*rhog/rhol,2.0d-4))
            gamma=(alphal-alphaad)/(alphaef-alphaad)
            IF(alphag.gt.alphasa.and.alphal.lt.alphaef)THEN
               gamma_star=gamma
            ELSE
               gamma_star=1.0d0
            ENDIF
            dgas=alphag**0.50d0*dh !equivalent wetted diameter
         
!            Rel=alphal*rhol*vr*dgas/mul !dgas?
            DO ix=1,ndim
               velo=vl_o(i,ix)*vl_o(i,ix) 
            ENDDO   
            velo=dsqrt(velo)
            Rel=alphal*rhol*velo*dh/mul
            G_star=1.d-4*Rel**0.25d0
!         
            Ce=7.5d0 !vertical
!
            sigma_star=DMAX1(sigma,1.d-7)
!            v_crit=3.2d0*(sigma_star*gravity*rhol_g)**0.25d0   
            v_crit=3.2d0*(sigma_star*gravity*rhol_g)**0.25/rhog**0.5  !LSJ modification           
            lambda=alphag*vg/v_crit
!         
            F11=gamma_star*DMAX1(0.0d0,(1-G_star))*DEXP(-Ce*1.0d-5*lambda**6.0d0)
            alphalf=DMAX1(0.0d0,alphal*F11)
!            
!............fi  
!    
!            reg=alphag*rhog*vr*dgas/mug
            reg=rhog*vr*dgas/mug   !LSJ modification
            reg=DMAX1(1.0d0,reg) !next
            f1=64.0d0/reg
            f3=0.02d0*(1.0d0+150.0d0*(1.0d0-(1.0d0-alphalf)**0.5d0))
            f2=(1500.d0-reg)/1000.0d0*f1+(reg-500.0d0)/1000.0d0*f3
            IF(reg.le.500.0d0)THEN
               fi=f1
            ELSEIF(reg.gt.500.0d0.and.reg.lt.1500.d0)THEN
               fi=f2
            ELSE
               fi=f3
            ENDIF
!            
!...........cann,ia_annular_ann         
!
!            cann=(30.0d0*alphalf)**1.8d0 !next
            cann=(30.0d0*alphalf)**0.125    !LSJ modification
            ia_annular_ann=4.0d0*cann/dh*(1.0d0-alphalf)**0.5d0
            IF(iam.eq.1)ia_annular_ann=cell%ia_annular_ann(i)
            IF(iam.eq.0)cell%ia_annular_ann(i)=ia_annular_ann
            vfgl_annular_ann=1.0d0/8.0d0*rhog*fi*ia_annular_ann*vr
!
!...........drp
!
!...........alphald
!
            alphaad_star=alphaad
            IF(alphag.gt.alphasa.and.alphal.lt.alphaef)alphaad_star=alphaad*gamma+(1.0d0-gamma)*1.d-5
            alphald=DMAX1((alphal-alphalf)/(1.0d0-alphalf),alphaad_star)
!            
!...........vfg2         
!   
            vfg=vr
            vfg_star=vfg*(1.0d0-F11)
!            IF(alphag.gt.alphasa.and.alphal.lt.alphaef)vfg_star=vfg*(1.0d0-F11*gamma_star)
            IF(alphag.gt.alphasa.and.alphal.lt.alphaef)vfg_star=vfg*(1.0d0-F11*gamma)   !LSJ modification   
            vfg_star2=vfg_star
            IF(alphal.lt.1.d-6)vfg_star2=vfg_star*alphal*1.d6
!
            vfg2=DMAX1(vfg*vfg,xa(2)/rhog*DMIN1(dcon(2)*alphal**0.333333333d0,dht))
            ddrp=xa(2)/(rhog*vfg2)
            ddrp=DMAX1(ddrp,84.0d-6)
            ddrp=DMIN1(ddrp,dh)
!           
                      
            alphagi=DMAX1(0.001d0,(1.0d0-alphal)) !pik_rv_debug           
            mum=mug/alphagi**0.25d0
!            red=(1.0d0-alphal)**2.5d0*rhog*ddrp*vr/mum
            red=alphagi**2.5d0*rhog*ddrp*dsqrt(vfg2)/mug    !LSJ modification
            red=DMAX1(red_min,red)
            cd=24.0d0/red*(1.0d0+0.1d0*red**0.75d0)
!
            ia_annular_drp=3.6d0*alphald/ddrp*(1.0d0-alphalf)
            IF(iam.eq.1)ia_annular_drp=cell%ia_annular_drp(i)
            IF(iam.eq.0)cell%ia_annular_drp(i)=ia_annular_drp
            vfgl_annular_drp=1.0d0/8.0d0*rhog*cd*ia_annular_drp*vr
!         
            vfgl_annular=vfgl_annular_ann+vfgl_annular_drp
         ENDIF           
!         
      ENDIF 
!
!.....restore the original value
!
      IF(flag.eq.1)THEN
         alphag=alphag_store
         alphal=1.0d0-alphag
         slug=slug_store
         annular=annular_store
         drift=drift_store         
      ENDIF
!      
      RETURN
      ENDSUBROUTINE bubblyslug_annular_for_churn     
