      SUBROUTINE invann_invslg_for_invchn(flag)
!
      USE VOL_DATA       , ONLY: cell
      USE Zparam         , ONLY: pi
      USE STM_TBL_cupid  , ONLY: pcrit
      USE Zrv_int_friction , ONLY: invann,invslg,invchn,alphag,alphal,alphabs,alphab,       &
                                   inertens,rhol_g,gravity,sigma,dh,dstar,dstar1,dstar2,    &
                                   delta,deltastar,cd,iam,ia_invann_ann,vfgl_invann_ann,    &
                                   rhog,vr,xa,rhol,dcon,dht,dbub,mum,mul,reb,reb_min,       &
                                   ia_invann_sb,vfgl_invann,vfgl_invann_sb,alphacd,alphasa, &
                                   alphadrp,alphamax,alphamin,alphab_1,ia_invslg_ann,mug,   &
                                   vfgl_invslg_ann,ddrp,ddrp,drop_max,red,red_min,          &
                                   ia_invslg_drp,vfgl_invslg,vfgl_invslg_drp,               &
                                   i
!
      IMPLICIT NONE      
!
      INTEGER flag
      INTEGER invann_store,invslg_store
      REAL(8) alphag_store
      REAL(8) F18,F17,F16,alpha_ian,alpha_bub !invann
      REAL(8) F21,alphastar,vfg !invslg
      REAL(8) vfg2,preduc,drmin 
!
!.....store the original value
!
      IF(flag.eq.1)THEN
          invann_store=invann
          invslg_store=invslg
          alphag_store=alphag
          invann=1
          invslg=1
      ENDIF      
!
!.....set alphag to left-end value
!
      IF(flag.eq.1)THEN
         alphag=alphabs
         alphal=1.0d0-alphag
      ENDIF         
!
!.....drag coefficient model for inverted annular
!
      IF(invann.eq.1)THEN
!
         F18=DMIN1(alphag/0.05d0,0.999999d0)
         F17=DEXP(-8.0d0*(alphabs-alphag)/alphabs)*F18
         alpha_ian=alphag
         IF(invchn.eq.1)alpha_ian=alphabs !next
         alphab=F17*alpha_ian
         alphab=DMIN1(0.99999d0,DMAX1(0.0d0,alphab)) !next, 0.99999d0
         alpha_bub=DMAX1((alpha_ian-alphab)/(1.0d0-alphab),1.d-7) !
         
         inertens=rhol_g*gravity/sigma
         inertens=inertens**0.5d0
         dstar=DMAX1(1.0d0/30.0d0,dh*inertens)
         dstar1=9.07d0/dstar
         dstar2=1.63d0+4.74d0/dstar
         delta=dh/2.0d0*(1.0d0-(1.0d0-alphab)**0.5d0)
         deltastar=DMAX1(delta*inertens,1.0d-8)
         cd=4.0d0*(0.005d0+0.2754d0*10.0d0**dstar1*deltastar*dstar2)
         ia_invann_ann=4.0d0/dh*(1.0d0-alphab)**0.5d0
         IF(iam.eq.1)ia_invann_ann=cell%ia_invann_ann(i)
         IF(iam.eq.0)cell%ia_invann_ann(i)=ia_invann_ann
         vfgl_invann_ann=1.0d0/8.0d0*rhog*cd*ia_invann_ann*vr
         vfg=vr
!!
!         dbub=DMAX1(dbub_min,5.0d0*sigma/(rhol*vr*vr)) !We=10.0
!         dbub=DMIN1(dbub_max,dbub)
!!
         vfg2=DMAX1(vfg*vfg,xa(1)/rhol*DMIN1(dcon(1)*alphag**0.333333333d0,dht))
         dbub=xa(1)/(rhol*vfg2)
         dbub=DMAX1(dbub,84.0d-6)      
!
!         mum=mul/(1.0d0-alphag)
         mum=mul/(1.0d0-alpha_bub)
!         
         reb=rhol*dbub*vr/mum
         reb=DMAX1(reb_min,reb)
         cd=24.0d0/reb*(1.0d0+0.1d0*reb**0.75d0)
         F16=1-F17
!         
!         ia_invann_sb=3.6d0*alpha_bub/dbub*(1.0d0-alphab)*F16
         ia_invann_sb=3.6d0*alpha_bub/dbub*(1.0d0-alphab)
!         
         IF(iam.eq.1)ia_invann_sb=cell%ia_invann_sb(i)
         IF(iam.eq.0)cell%ia_invann_sb(i)=ia_invann_sb        
         vfgl_invann_sb=1.0d0/8.0d0*rhol*cd*ia_invann_sb*vr         
!
         vfgl_invann=vfgl_invann_ann+vfgl_invann_sb         
!
      ENDIF
!
!.....set alphag to right-end value
!     
      IF(flag.eq.1)THEN
         alphag=alphacd
         alphal=1.0d0-alphag
      ENDIF    
!
!.....drag coefficient model for inverted slug
!
      IF(invslg.eq.1)THEN !next
!
         alphastar=(alphasa-alphag)/(alphasa-alphabs)
         F21=DEXP(-alphastar)
         alphadrp=(1.0d0-alphasa)*F21
         alphadrp=DMIN1(alphamax,(DMAX1(alphamin,alphadrp)))
!         alphadrp=DMAX1(1.0d-4,alphal)         
!         
         alphab=(alphal-alphadrp)/(1.0d0-alphadrp)
         alphab=DMIN1(1.0d0,DMAX1(0.0d0,alphab))
         alphab_1=1.0d0-alphab
         cd=10.9d0*alphab**0.5d0*alphab_1*alphab_1*alphab_1
         ia_invslg_ann=alphab
         IF(iam.eq.1)ia_invslg_ann=cell%ia_invslg_ann(i)
         IF(iam.eq.0)cell%ia_invslg_ann(i)=ia_invslg_ann
!
!         vfgl_invslg_ann=1.0d0/8.0d0*rhog*cd*ia_invslg_ann*vr !next-too small
         vfgl_invslg_ann=1.0d0/2.0d0*rhog*cd*ia_invslg_ann*vr
!         
         vfg=vr
!!         
!         We=6.d0
!         mum=mug
!         !next vfg=vr*F21*F21
!         ddrp=We*sigma/(rhog*vfg*vfg) 
!         pstar=pres/22.4d6
!         IF(pstar.lt.0.025d0)THEN
!            ddrp=DMAX1(0.0025d0,ddrp)
!         ELSEIF(pstar.gt.0.25d0)THEN
!            ddrp=DMAX1(0.0002d0,ddrp)
!         ENDIF
!         ddrp=DMIN1(ddrp_max,ddrp)
!!
         mum=mug
         vfg=DMAX1(DABS(vfg),0.001d0)
         vfg2=vfg*vfg !next
         ddrp=xa(3)/(rhog*vfg2)
         preduc=cell%p(i)/pcrit
         IF(preduc.lt.0.25d0)THEN
            IF(preduc.lt.0.025d0)THEN
               drmin=dcon(2)
            ELSE
               drmin=dcon(2)+(dcon(3)-dcon(2))*4.44444d0*(preduc-0.025d0)
            ENDIF
         ELSE
            drmin=dcon(3)
         ENDIF
!         ddrp=DMAX1(ddrp,drmin)
!         ddrp=DMIN1(dh,dcon(2),ddrp)
         ddrp=DMIN1(dh,ddrp)         ! yjm  ori: on
         ddrp=DMAX1(DMIN1(ddrp,drop_max),8.6d-5)
         cell%ddrp(i)=ddrp
!                    
!         red=rhog*ddrp*vr/mum 
         red=rhog*ddrp*vr*(1.0d0-alphadrp)**2.5d0/mum 
         red=DMAX1(red_min,red)
!         cd=DMIN1(24.0d0/red*(1.0d0+0.1d0*red**0.75d0),0.05625d0*8.0d0) 
         cd=DMAX1(24.0d0/red*(1.0d0+0.1d0*red**0.75d0),0.05625d0*8.0d0) 
         ia_invslg_drp=3.6d0*alphadrp/ddrp*(1.0d0-alphab)
         IF(iam.eq.1)ia_invslg_drp=cell%ia_invslg_drp(i)
         IF(iam.eq.0)cell%ia_invslg_drp(i)=ia_invslg_drp
         vfgl_invslg_drp=1.0d0/8.0d0*rhog*cd*ia_invslg_drp*vfg   !next-too small
!
         vfgl_invslg=vfgl_invslg_ann+vfgl_invslg_drp
!                
      ENDIF
!
!.....restore the original value 
!
      IF(flag.eq.1)THEN
         alphag=alphag_store
         alphal=1.0d0-alphag
         invann=invann_store
         invslg=invslg_store
      ENDIF
!             
      RETURN
      ENDSUBROUTINE invann_invslg_for_invchn
