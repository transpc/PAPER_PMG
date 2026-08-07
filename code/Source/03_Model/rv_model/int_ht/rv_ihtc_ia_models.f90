!
      SUBROUTINE rv_ihtc_ia_preCHF(i,reg_idx)
! 
      USE VOL_DATA 
      USE Zparam          , ONLY: ndim,pi   
      USE Zconst2         , ONLY: hydraulicd,ggc
      USE Zvector         , ONLY: ul_o,ug_o,vrel_o
      USE Zvector         , ONLY: vl_o,vg_o
      USE Zrv_flowmap     , ONLY: alp_tb
      USE Zrv_ihtc_models , ONLY: F9
      USE Zrv_int_friction, ONLY: xa,dcon,web,drop_max
      USE Zwall_HTC       , ONLY: reflod,rey_reflod
      USE Zrv_model       , ONLY: ia_option 
!
      IMPLICIT NONE
!
      INTEGER i,reg_idx
!
      REAL(8) We_crit,dia_min,beta,xag,alphabub,dia_max,vfg,vfg2   
      REAL(8) cF9,alphatb   
      REAL(8) alpha_ad,alpha_ef,xal,gamma,gamma_s,rey_f,v_crit,f11,alpha_ff,alpha_ad_s,vfg_s,alpha_fd,vfg_ss,vfg_hat_2,d_drp
      REAL(8) alphadrp
      REAL(8) vg,vl,del_rho,jg
      REAL(8) lapno,dbbl,vsupf,dbbm,eta,lo,nlo,jj,nref,rhol,rhog,dh,mu,dpdx,b1,eps,reb,dbbs,dbb_sb,agf0,sig
      REAL(8) f99,alpha_gs,alpha_tb,lambda,cf10,cf11,d_star,agf
      REAL(8) dr
      REAL(8) dbubl,dbubm 
!
! reflod models
      REAL(8) alpanf,centr,vlf,vlg,fluxm,rhofg,vcritl,rvcrit,lamcr,scrchh,vfgbb,alphad,alphed,voidb,alp,voidx, &
               alpdrp,velgdf
!                     
!     REAL(8),PARAMETER :: pcrit=22.06d6      !critical pressure of H2O 
      REAL(8),PARAMETER :: alpha_am=0.9999d0          
!
      LOGICAL, SAVE :: lsj_coding,ljr_mars
!      DATA lsj_coding,ljr_mars /.false.,.true./   !pik-bug-rv-ihtc    
!      DATA lsj_coding,ljr_mars /.true.,.true./      
!
!      IF(liquid.eq.1.or.gas.eq.1) THEN
!         RETURN
!      ENDIF
!
      IF(ia_option.eq.1) THEN
         lsj_coding=.true.
         ljr_mars=.true.
      ELSEIF(ia_option.eq.2) THEN
         lsj_coding=.false.
         ljr_mars=.true.          
      ENDIF      
!
      ggc=9.81d0
      vg=ug_o(i)
      vl=ul_o(i)
      rhog=cell%rhog(i)
      rhol=cell%rhol(i)      
      del_rho=cell%rhol(i)-cell%rhog(i)
      del_rho=DMAX1(1.d-5,del_rho) !LSJ
      dh=hydraulicd(i)
      mu=cell%lviscosl(i)
      sig=cell%sigma(i)
      drop_max=dh
!      cell%length(i)=dh     
!
      IF(reflod)THEN
         dcon(1)=5.0d-3
!         dcon(2)=1.5d-3
!         dcon(2)=5.0d-4  ! MARS source, FLECHT SEASET min drp size: 2.5 mm / original: 1.5 mm / 31805 
         dcon(2)=2.5d-4  ! MARS source interphaseDrag Line 1962
         dcon(3)=2.0d-4
         web(1)=5.0d0
         web(2)=1.5d0
!         web(3)=1.5d0
!
!         use TRACE 5.0 We number of reflood We (based on dmax) =12.0
!         sauter mean diam = 1/3 * dmax
!         We = 3/1*12.0 = 4.0
!
         web(3)=4.0d0
      ELSE
         dcon(1)=5.0d-3
         dcon(2)=2.5d-3
         dcon(3)=2.0d-4
         web(1)=5.0d0
         web(2)=1.5d0
         web(3)=6.0d0          
      ENDIF    
      xa(:)=DMAX1(web(:)*sig,1.0d-10)  
!
!.....Bubbly
!     
      IF(reg_idx.eq.2) THEN
!       
         IF(lsj_coding) THEN  !LSJ coding
            xag=DMAX1(0.0,DMIN1(cell%alphag(i),cell%alpha_bs(i)))
!            
            We_crit=5.d0
            dia_min=0.005d0
            beta=1.d0
            alphabub=DMAX1(xag,1.0D-5)
            dia_max=DMIN1(dia_min*alphabub**(1.0/3.0),hydraulicd(i))
            vfg=(ug_o(i)-ul_o(i))*beta*beta                                                    
            vfg2=DMAX1(vfg*vfg,We_crit*cell%sigma(i)/cell%rhol(i)/dia_max)
            vfg=DSQRT(vfg2)
!
            cell%dbb(i)=We_crit*cell%sigma(i)/cell%rhol(i)/vfg2 
            cell%ia_bubbly(i)=3.6d0*alphabub/cell%dbb(i)               
!       
         ELSE  !LJR coding
!                   
!           dbb_mars: dbbl, dbbm         
!            xag=cell%alphag(i)      
            xag=DMAX1(0.0,DMIN1(cell%alphag(i),cell%alpha_bs(i)))
            xal=1.d0-xag
!            
            alphabub=DMAX1(xag,1.d-5)
            lapno=DSQRT(sig/ggc/del_rho)
            vsupf=DMAX1(0.01d0,xal*vl)
!            
            dbbl=DMIN1(0.9d0*dh,2.0d0*lapno)            
            dbbm=2.96d0*sig**0.6*dh**0.48*mu**0.08d0/(rhol**0.68*vsupf**1.12)
!
!           dbb_space: dbbs
            eta=1.d0
            lo=DSQRT(sig/ggc/del_rho)
            nlo=lo/dh
            jj=xag*vg+xal*vl
            jg=jj/DMAX1(1.d-5,xag*rhog*vg+xal*rhol*vl)
            jg=DABS(jg)
!
            nref=DMAX1(1.d-5,xal*rhol*vl*dh/mu)
            nref=DMIN1(nref,1.d6)  !LSJ
            dpdx=10.d0
            b1=0.25d0*cell%length(i)*pi*dh*dh*(1.d0-DEXP(-.0005839*nref))
!
            eps=ggc*xag*vg/DEXP(0.0005839*nref)+jg*dpdx*b1
            reb=eps**0.33333*lo**1.333*rhol/mu
!            
            reb=DMAX1(1.0d0,reb) !pik-apr1400-rec 
            dbbs=1.99*lo/eta*nlo**(-0.335)*reb**(-0.239)
!             
            IF(ljr_mars) THEN
               dbb_sb=DMAX1(1.d-4,DMIN1(dbbm,dbbl)) 
!               agf0=3.6d0*alphabub/dbb_sb
               agf0=6.0d0*alphabub/dbb_sb  !modified by LSJ, since CUPID has used 6*alphag/D1 for the IA
            ELSE
               dbb_sb=dbbs
               agf0=3.02d0*xag*eta*nlo**0.335*reb*0.239/lo
            ENDIF   
!               
            cell%dbb(i)=dbb_sb
            cell%ia_bubbly(i)=agf0
         ENDIF   
      ENDIF   
!
!.....Slug
!      
      IF(reg_idx.eq.3) THEN
!       
         IF(lsj_coding) THEN   !LSJ coding   
            xag=DMAX1(cell%alpha_bs(i),DMIN1(cell%alphag(i),cell%alpha_de(i))) 
!            
            We_crit=5.d0
            dia_min=0.005d0
            cF9=F9(xag,cell%alpha_bs(i),cell%alpha_sa(i)) 
            beta=cF9
            alphabub=cell%alpha_bs(i)*cF9
            alphatb=(xag-alphabub)/(1.d0-alphabub)
!               
!            dia_max=DMIN1(dia_min*alphabub**(1.0/3.0),hydraulicd(i))
!            vfg=(ug_o(i)-ul_o(i))*beta*beta                                                    
!            vfg2=DMAX1(vfg*vfg,We_crit*cell%sigma(i)/cell%rhol(i)/dia_max)
!            vfg=DSQRT(vfg2)
!
            lapno=DSQRT(cell%sigma(i)/(9.81d0*DMAX1(1.0d0,cell%rhol(i)*cell%rhog(i))))
            dbubl=DMIN1(0.9d0*hydraulicd(i),2.0d0*lapno)
!
            vsupf=cell%alphal(i)*vl
            vsupf=DMAX1(0.01d0,vsupf)
!
            dbubm=2.96d0*cell%sigma(i)**0.6d0*hydraulicd(i)**0.48d0*cell%lviscosl(i)**0.08d0/(cell%rhol(i)**0.68d0*vsupf**1.12d0)
            dbubl=DMAX1(1.0d-4,DMIN1(dbubm,dbubl))
!
!            cell%dsb(i)=We_crit*cell%sigma(i)/cell%rhol(i)/vfg2 
            cell%dsb(i)=dbubl  ! yjm
!            
            cell%dbb(i)=cell%dsb(i)
!            cell%ia_slug_sb(i)=3.6d0*alphabub*(1.d0-alphatb)*cF9/cell%dsb(i)  
            cell%ia_slug_sb(i)=6.0d0*alphabub*(1.d0-alphatb)/cell%dsb(i)  ! yjm ref. MARS coding / interphase line 1069
            cell%ia_slug_tb(i)=4.5d0*alphatb*2.d0/hydraulicd(i)  
            alp_tb(i)=alphatb
!        
         ELSE  !LJR coding
            xag=DMAX1(cell%alpha_bs(i),DMIN1(cell%alphag(i),cell%alpha_de(i))) 
            xal=1.d0-xag
!            
            lapno=DSQRT(sig/ggc/del_rho)
            vsupf=DMAX1(0.01d0,xal*vl)
!                        
            dbbl=DMIN1(0.9d0*dh,2.0d0*lapno)
            dbbm=2.96d0*sig**0.6*dh**0.48*mu**0.08/(rhol**0.68d0*vsupf**1.12)
!            
            cell%dsb(i)=DMAX1(1.d-4,DMIN1(dbbm,dbbl))
            cell%dbb(i)=cell%dsb(i)
!           
!           alpha_gs
            f99=(xag-cell%alpha_bs(i))/(cell%alpha_sa(i)-cell%alpha_bs(i))
            cf9=DEXP(-8.d0*f99)
            alpha_gs=cell%alpha_bs(i)*cf9
            alphabub=cell%alpha_bs(i)*cF9
!
!           taylor/small bubble portion
!            alpha_tb =(xag-alpha_gs)/DMAX1(1.d-8,1.d0-alpha_gs)                 ! MARS manual (alpha_gs)
            alpha_tb=(xag-cell%alpha_bs(i))/DMAX1(1.d-8,1.d0-cell%alpha_bs(i))  ! MARS code   (alpha_bs)
!            
            cell%ia_slug_tb(i)=4.5d0/dh*alpha_tb*2.0d0
            cell%ia_slug_sb(i)=3.6d0*alphabub/cell%dbb(i)*(1.d0-alpha_tb)
!            cell%ia_slug_sb(i)=6.0d0*alphabub/cell%dbb(i)*(1.d0-alpha_tb)            
            alp_tb(i)=alpha_tb
            
         ENDIF   
!
      ENDIF  
!
!.....Churn
!
      IF(reg_idx.eq.4) THEN  
!       
         IF(lsj_coding) THEN  !LSJ coding      
!           slug         
            xag=cell%alpha_de(i)
!         
            We_crit = 5.d0
            dia_min = 0.005d0
            cF9 = F9(xag,cell%alpha_bs(i),cell%alpha_sa(i)) 
            beta = cF9
            alphabub=cell%alpha_bs(i)*cF9
            alphatb=(xag-alphabub)/(1.d0-alphabub)
!               
            dia_max = DMIN1( dia_min * alphabub**(1.0/3.0), hydraulicd(i) )
            vfg = ( ug_o(i) - ul_o(i) ) * beta * beta                                                    
            vfg2 = DMAX1( vfg*vfg, We_crit * cell%sigma(i) / cell%rhol(i) / dia_max )
            vfg = DSQRT(vfg2)
!
            cell%dsb(i)=We_crit * cell%sigma(i) / cell%rhol(i) / vfg2 
            cell%dbb(i)=cell%dsb(i)         
            cell%ia_slug_sb(i)=3.6d0*alphabub*(1.d0-alphatb)*cF9/cell%dsb(i)  
!            cell%ia_slug_sb(i)=6.0d0*alphabub*(1.d0-alphatb)*cF9/cell%dsb(i)              
            cell%ia_slug_tb(i)=4.5d0*alphatb*2.d0/hydraulicd(i)  
            alp_tb(i)=alphatb         
!
!           anular: annular film(a_an)
            xag=cell%alpha_sa(i)
!          
            xal=1.0d0-xag
            alpha_ad=1.0d-4 
            alpha_ef=DMAX1(alpha_ad+alpha_ad,DMIN1(2.0d-3*cell%rhog(i)/cell%rhol(i),2.0d-4))
            gamma=(xal-alpha_ad)/(alpha_ef-alpha_ad)
            IF(xag>cell%alpha_sa(i).AND.xal<alpha_ef) THEN
               gamma_s=gamma
            ELSE
               gamma_s=1.0d0
            ENDIF
            rey_f=xal*cell%rhol(i)*DABS(vl_o(i,ndim))*hydraulicd(i)/cell%lviscosl(i)    
            dr=DMAX1(1.d-5,cell%rhol(i)-cell%rhog(i)) 
            v_crit=3.2d0*(DMAX1(cell%sigma(i),1.0d-7)*ggc*dr/(cell%rhog(i)*cell%rhog(i)))**0.25 
            f11=gamma_s*DMAX1(0.0d0,1.0d0-rey_f**0.25*1.d-4)*DEXP(-7.5d-5*(xag*vg_o(i,ndim)/v_crit)**6) 
            alpha_ff=DMAX1(0.0d0,xal*f11)  
            cell%ia_annular_ann(i)=9.0d0/hydraulicd(i)*DSQRT(1.0d0-alpha_ff)*(30.0d0*alpha_ff)**0.125
!
!           anular: drops(a_drp, d2)
            vfg=vg_o(i,NDIM)-vl_o(i,NDIM)
            IF(xag>cell%alpha_sa(i).AND.xal<alpha_ef) THEN
               alpha_ad_s=alpha_ad*gamma+1.d-5*(1.d0-gamma)
               vfg_s=vfg*(1.d0-f11*gamma)
            ELSE
               alpha_ad_s=alpha_ad
               vfg_s=vfg*(1.d0-f11)
            ENDIF
            alpha_fd=DMAX1((xal-alpha_ff)/(1.d0-alpha_ff),alpha_ad_s)
            IF(xal<1.d-6) THEN
               vfg_ss=vfg_s*xal*1.d6
            ELSE
               vfg_ss=vfg_s
            ENDIF
!            vfg_hat_2=DMAX1(1.5d0*cell%sigma(i),1.d-10)/(cell%rhog(i)*DMIN1(0.0025d0*alpha_fd**(1.d0/3.d0),hydraulicd(i)))
            vfg_hat_2=xa(2)/(cell%rhog(i)*DMIN1(dcon(2)*alpha_fd**(1.d0/3.d0),hydraulicd(i)))
            vfg_hat_2=DMAX1(vfg_ss*vfg_ss,vfg_hat_2)
!            d_drp=1.5d0*cell%sigma(i)/(cell%rhog(i)*vfg_hat_2)
            d_drp=web(2)*cell%sigma(i)/(cell%rhog(i)*vfg_hat_2)
            d_drp=DMIN1(dh,drop_max,DMAX1(84.0d-6,d_drp))
            cell%ia_annular_drp(i)=3.6d0*alpha_fd*(1.d0-alpha_ff)/d_drp
!            cell%ia_annular_drp(i)=6.0d0*alpha_fd*(1.d0-alpha_ff)/d_drp            
            cell%ddrp(i)=d_drp          
!       
         ELSE   !LJR coding
!           slug             
            xag=cell%alpha_de(i)
            xal=1.d0-xag           
!            
            lapno=DSQRT(sig/ggc/del_rho)
            vsupf=DMAX1(0.01d0,xal*vl)
!                        
            dbbl=DMIN1(0.9d0*dh,2.0d0*lapno)
            dbbm=2.96d0*sig**0.6*dh**0.48*mu**0.08/(rhol**0.68d0*vsupf**1.12)
!            
            cell%dsb(i)=DMAX1(1.d-4,DMIN1(dbbm,dbbl))
            cell%dbb(i)=cell%dsb(i)
!           
!           alpha_gs
            f99=(xag-cell%alpha_bs(i))/(cell%alpha_sa(i)-cell%alpha_bs(i))
            cf9=DEXP(-8.d0*f99)
            alpha_gs=cell%alpha_bs(i)*cf9
            alphabub=cell%alpha_bs(i)*cF9
!
!           taylor/small bubble portion
!            alpha_tb =(xag-alpha_gs)/DMAX1(1.d-8,1.d0-alpha_gs)                 ! MARS manual (alpha_gs)
            alpha_tb=(xag-cell%alpha_bs(i))/DMAX1(1.d-8,1.d0-cell%alpha_bs(i))  ! MARS code   (alpha_bs)
!            
            cell%ia_slug_tb(i)=4.5d0/dh*alpha_tb*2.0d0
            cell%ia_slug_sb(i)=3.6d0*alphabub/cell%dbb(i)*(1.d0-alpha_tb)
!            cell%ia_slug_sb(i)=6.0d0*alphabub/cell%dbb(i)*(1.d0-alpha_tb)           
            alp_tb(i)=alpha_tb         
!
!           anular: annular film(a_an)
            xag=cell%alpha_sa(i)
            xal=1.d0-xag
!
            v_crit=3.2d0*(sig*ggc*del_rho)**0.25/DSQRT(rhog)
            rey_f=xal*rhol*vl*dh/mu
            lambda=xag*vg/DSQRT(v_crit)            
            cf10=DMIN1(1.0d0+DSQRT(lambda)+0.05d0*DABS(lambda), 6.d0)
            alpha_ef=DMAX1(2.d-4,DMIN1(2.d-3*rhog/rhol,2.d-4))
            gamma=(xal-1.d-4)/DMAX1(1.d-8,alpha_ef-1.d-4)
            IF(xag.gt.cell%alpha_sa(i) .and. xal.lt.alpha_ef)then
               gamma_s=gamma
            ELSE
               gamma_s=1.d0
            ENDIF
            cf11=gamma_s*DMAX1(0.d0,(1.d0-1.d-4*rey_f**0.25))*DEXP(-7.5d-5*lambda**6)
            alpha_ff=DMAX1(0.d0,xal*cf11)
!            
            cell%ia_annular_ann(i)=4.d0*2.5d0*(30.d0*alpha_ff)**0.125/dh*DSQRT(1.d0-alpha_ff)
!
!           anular: drop
            IF(xag.gt.cell%alpha_sa(i) .and. xal.lt.alpha_ef)then
               alpha_ad_s=1.d-4*gamma+1.d-5*(1.d0-gamma)
            ELSE
               alpha_ad_s=1.d-4
            ENDIF
            alpha_fd=dmax1((xal-alpha_ff)/(1.d0-alpha_ff),alpha_ad_s)
            vfg=vg-vl
            IF(xag.gt.cell%alpha_sa(i) .and. xal.lt.alpha_ef)then
               vfg_s=vfg*(1.d0-cf11*lambda)
            ELSE
               vfg_s=vfg*(1.d0-cf11)
            ENDIF
            IF(xal.lt.1.d-6)then
               vfg_ss=vfg_s*xal*1.d6
            ELSE
               vfg_ss=vfg_s
            ENDIF
            vfg_hat_2=1.5d0*sig/(rhog*dmin1(0.0025d0*alpha_fd**(1.d0/3.d0),dh))
            vfg_hat_2=dmax1(vfg_ss*vfg_ss,vfg_hat_2)
!            
            cell%ddrp(i)=dmax1(1.5d0*sig,1.d-10)/(rhog*vfg_hat_2)
            cell%ddrp(i)=DMIN1(dh,DMAX1(84.0d-6,cell%ddrp(i)))
            cell%ia_annular_drp(i)=3.6d0*alpha_fd*(1.d0-alpha_ff)/cell%ddrp(i) 
!            cell%ia_annular_drp(i)=6.0d0*alpha_fd*(1.d0-alpha_ff)/cell%ddrp(i)  
!
!reflod model by LSJ 180123
            IF(reflod) THEN 
               alpanf=1.d0-xag
               centr=7.5d-5
               vlf=vl_o(i,ndim)
               vlg=vg_o(i,ndim)
               fluxm=1.0d-4*dsqrt(dsqrt(dabs(vlf*alpanf*rhol*dh/mu)))
               rhofg=dmax1(1.0d0,rhol-rhog) 
               vcritl=5.66d0*dsqrt(dsqrt(dmax1(sig,1.0d-7)*rhofg)/rhog) 
               rvcrit=xag*vlg/vcritl 
               lamcr=rvcrit
               scrchh=centr*lamcr**6
               IF(scrchh.le.200.d0) THEN
                  vfgbb=dexp(-scrchh)*dmax1(0.0d0,(1.0d0-fluxm))
               ELSE
                  vfgbb=0.d0
               ENDIF
               alphad=1.d-4
               alphed=1.d-3
               scrchh=(alpanf-alphad)/(alphed-alphad) 
               alphad=alphad*scrchh+1.0d-5*(1.0d0-scrchh)
               vfgbb=vfgbb*scrchh
               voidb=DMAX1(0.d0,alpanf*vfgbb)               
!            
               vfg=scrchh
               alp=dmax1((alpanf-voidb)/(1.d0-voidb),alphad)
               vfg_hat_2=max(vfg*vfg,1.5d0*sig/(rhog*dmin1(0.0015d0*alp**0.333333333d0,dh)))
               vfg_hat_2=dmax1(vfg_ss*vfg_ss,vfg_hat_2)
!
               voidx=1.d0-alp               
               rey_reflod=1.5d0*sig*voidx**3/(cell%lviscosg(i)*dsqrt(vfg_hat_2*voidx))
!
               cell%ddrp(i)=1.5d0*sig/(rhog*vfg_hat_2)
               cell%ddrp(i)=dmax1(cell%ddrp(i),84.d-6) !Limit drop diameter to > 84 microns (TRAC)
               cell%ddrp(i)=DMIN1(cell%ddrp(i),dh)
!
               cell%ia_annular_drp(i)=3.6d0*alp/cell%ddrp(i) 
            ENDIF   
!                                       
         ENDIF
      ENDIF      
!
!.....Annular
!
      IF(reg_idx.eq.5) THEN
!       
         IF(lsj_coding) THEN   !LSJ coding           
!           annular film: a_an 
            xag=DMIN1(alpha_am,DMAX1(cell%alphag(i),cell%alpha_sa(i)))
            xal=1.0d0-xag
!            
            alpha_ad=1.0d-4 
            alpha_ef=DMAX1(alpha_ad+alpha_ad,DMIN1(2.0d-3*cell%rhog(i)/cell%rhol(i),2.0d-4))
            gamma=(xal-alpha_ad)/(alpha_ef-alpha_ad)
            IF(xag>cell%alpha_sa(i).AND.xal<alpha_ef) THEN
               gamma_s=gamma
            ELSE
               gamma_s=1.0d0
            ENDIF
            rey_f=xal*cell%rhol(i)*DABS(vl_o(i,ndim))*hydraulicd(i)/cell%lviscosl(i)     
            dr=DMAX1(1.d-5,cell%rhol(i)-cell%rhog(i))
            v_crit=3.2d0*(DMAX1(cell%sigma(i),1.0d-7)*ggc*dr/(cell%rhog(i)*cell%rhog(i)))**0.25 
            f11=gamma_s*DMAX1(0.0d0,1.0d0-rey_f**0.25*1.d-4)*DEXP(-7.5d-5*(xag*vg_o(i,ndim)/v_crit)**6) 
            alpha_ff=DMAX1(0.0d0,xal*f11)  
            cell%ia_annular_ann(i)=9.0d0/hydraulicd(i)*DSQRT(1.0d0-alpha_ff)*(30.0d0*alpha_ff)**0.125
!
!           drops: a_drp, d2
            vfg=vg_o(i,NDIM)-vl_o(i,NDIM)
            IF(xag>cell%alpha_sa(i).AND.xal<alpha_ef) THEN
               alpha_ad_s=alpha_ad*gamma+1.d-5*(1.d0-gamma)
               vfg_s=vfg*(1.d0-f11*gamma)
            ELSE
               alpha_ad_s=alpha_ad
               vfg_s=vfg*(1.d0-f11)
            ENDIF
            alpha_fd=DMAX1((xal-alpha_ff)/(1.d0-alpha_ff),alpha_ad_s)
            IF(xal<1.d-6) THEN
               vfg_ss=vfg_s*xal*1.d6
            ELSE
               vfg_ss=vfg_s
            ENDIF
!            vfg_hat_2=DMAX1(1.5d0*cell%sigma(i),1.d-10)/(cell%rhog(i)*DMIN1(0.0025d0*alpha_fd**(1.d0/3.d0),hydraulicd(i)))
            vfg_hat_2=xa(2)/(cell%rhog(i)*DMIN1(dcon(2)*alpha_fd**(1.d0/3.d0),hydraulicd(i)))
            vfg_hat_2=DMAX1(vfg_ss*vfg_ss,vfg_hat_2)
!            d_drp=1.5d0*cell%sigma(i)/(cell%rhog(i)*vfg_hat_2)
            d_drp=web(2)*cell%sigma(i)/(cell%rhog(i)*vfg_hat_2)
            d_drp=DMIN1(dh,drop_max,DMAX1(84.0d-6,d_drp))
!            cell%ia_annular_drp(i)=3.6d0*alpha_fd*(1.d0-alpha_ff)/d_drp
!            cell%ia_annular_drp(i)=6.0d0*alpha_fd*(1.d0-alpha_ff)/d_drp            
            cell%ddrp(i)=d_drp 
            cell%ia_annular_drp(i)=3.6d0*alpha_fd*(1.d0-alpha_ff)/cell%ddrp(i) 
!       
         ELSE  !LJR coding
!            xag=cell%alphag(i)
            xag=DMIN1(alpha_am,DMAX1(cell%alphag(i),cell%alpha_sa(i)))
            xal=1.d0-xag
!
!           annular portion
            v_crit=3.2d0*(sig*ggc*del_rho)**0.25/DSQRT(rhog)
            rey_f=xal*rhol*vl*dh/mu
            lambda=xag*vg/DSQRT(v_crit)            
            cf10=DMIN1(1.0d0+DSQRT(lambda)+0.05d0*DABS(lambda), 6.d0)
            alpha_ef=DMAX1(2.d-4,DMIN1(2.d-3*rhog/rhol,2.d-4))
            gamma=DMAX1(0.d0,DMIN1(1.d0,(xal-1.d-4)/DMAX1(1.d-8,alpha_ef-1.d-4)))
            IF(xag.gt.cell%alpha_sa(i) .and. xal.lt.alpha_ef)then
               gamma_s=gamma
            ELSE
               gamma_s=1.d0
            ENDIF
            cf11=gamma_s*DMAX1(0.d0,(1.d0-1.d-4*rey_f**0.25))*DEXP(-7.5d-5*lambda**6)
            alpha_ff=DMAX1(0.d0,xal*cf11)
!            
            cell%ia_annular_ann(i)=4.d0*2.5d0*(30.d0*alpha_ff)**0.125/dh*DSQRT(1.d0-alpha_ff)
!
!           mist (drop) portion
            IF(xag.gt.cell%alpha_sa(i) .and. xal.lt.alpha_ef)then
               alpha_ad_s=1.d-4*gamma+1.d-5*(1.d0-gamma)
            ELSE
               alpha_ad_s=1.d-4
            ENDIF
            alpha_fd=dmax1((xal-alpha_ff)/(1.d0-alpha_ff),alpha_ad_s)
            vfg=vg-vl
            IF(xag.gt.cell%alpha_sa(i) .and. xal.lt.alpha_ef)then
               vfg_s=vfg*(1.d0-cf11*lambda)
            ELSE
               vfg_s=vfg*(1.d0-cf11)
            ENDIF
            IF(xal.lt.1.d-6)then
               vfg_ss=vfg_s*xal*1.d6
            ELSE
               vfg_ss=vfg_s
            ENDIF
            vfg_hat_2=1.5d0*sig/(rhog*dmin1(0.0025d0*alpha_fd**(1.d0/3.d0),dh))
            vfg_hat_2=dmax1(vfg_ss*vfg_ss,vfg_hat_2)
!            
            cell%ddrp(i)=dmax1(1.5d0*sig,1.d-10)/(rhog*vfg_hat_2)
            cell%ddrp(i)=DMIN1(dh,DMAX1(84.0d-6,cell%ddrp(i)))
            cell%ia_annular_drp(i)=3.6d0*alpha_fd*(1.d0-alpha_ff)/cell%ddrp(i)
!            cell%ia_annular_drp(i)=6.0d0*alpha_fd*(1.d0-alpha_ff)/cell%ddrp(i)  
!
!reflod model by LSJ 180123
            IF(reflod) THEN
               alpanf=1.d0-xag
               centr=7.5d-5
               vlf=vl_o(i,ndim)
               vlg=vg_o(i,ndim)
               fluxm=1.0d-4*dsqrt(dsqrt(dabs(vlf*alpanf*rhol*dh/mu)))
               rhofg=dmax1(1.0d0,rhol-rhog) 
               vcritl=5.66d0*dsqrt(dsqrt(dmax1(sig,1.0d-7)*rhofg)/rhog) 
               rvcrit=xag*vlg/vcritl 
               lamcr=rvcrit
               scrchh=centr*lamcr**6
               IF(scrchh.le.200.d0) THEN
                  vfgbb=dexp(-scrchh)*dmax1(0.0d0,(1.0d0-fluxm))
               ELSE
                  vfgbb=0.d0
               ENDIF
               alphad=1.d-4
               alphed=1.d-3
               scrchh=(alpanf-alphad)/(alphed-alphad) 
               alphad=alphad*scrchh+1.0d-5*(1.0d0-scrchh)
               vfgbb=vfgbb*scrchh
               voidb=DMAX1(0.d0,alpanf*vfgbb)               
!            
               vfg=scrchh
               alp=dmax1((alpanf-voidb)/(1.d0-voidb),alphad)
               vfg_hat_2=max(vfg*vfg,1.5d0*sig/(rhog*dmin1(0.0015d0*alp**0.333333333d0,dh)))
               vfg_hat_2=dmax1(vfg_ss*vfg_ss,vfg_hat_2)
!
               voidx=1.d0-alp               
               rey_reflod=1.5d0*sig*voidx**3/(cell%lviscosg(i)*dsqrt(vfg_hat_2*voidx))
!
               cell%ddrp(i)=1.5d0*sig/(rhog*vfg_hat_2)
               cell%ddrp(i)=dmax1(cell%ddrp(i),84.d-6) !Limit drop diameter to > 84 microns (TRAC)
               cell%ddrp(i)=DMIN1(cell%ddrp(i),dh)
!
               cell%ia_annular_drp(i)=3.6d0*alp/cell%ddrp(i) 
            ENDIF   
                  
         ENDIF   
      ENDIF   
!
!.....Mist
!
      IF(reg_idx.eq.6.or.reg_idx.eq.11.or.reg_idx.eq.12) THEN   
!       
         IF(lsj_coding) THEN   !LSJ coding
!            We_crit = 6.d0
!            dia_min = 0.0002d0
            We_crit = web(2)
            dia_min = dcon(2)

!
            xag=DMAX1(cell%alpha_sa(i),DMIN1(cell%alphag(i),alpha_am))
            xal=1.d0-xag      
!         
            alphadrp=DMAX1(xal,1.d-4)         
!          
            vfg  = ( ug_o(i) - ul_o(i) )
!            vfg2 = DMAX1( vfg*vfg, 1.d-6 )
            vfg2=DMAX1(vfg*vfg,xa(2)/(cell%rhog(i)*dmin1(dcon(2)*alphadrp**0.333333333d0,dh)))
!
!           small bubble
            cell%ddrp(i) = We_crit*cell%sigma(i)/cell%rhog(i)/vfg2
            cell%ddrp(i) = DMIN1(dh,DMAX1(dia_min,cell%ddrp(i)))
            cell%ddrp(i) = DMIN1(drop_max,DMAX1(84.0d-6,cell%ddrp(i))) ! yjm                             
! 
!           droplets
            cell%ia_mist(i)=3.6d0*alphadrp/cell%ddrp(i)
            cell%ia_mpo(i)=3.6d0*alphadrp/cell%ddrp(i)
            cell%ia_mpr(i)=3.6d0*alphadrp/cell%ddrp(i)
!            cell%ia_mist(i)=6.0d0*alphadrp/cell%ddrp(i)
!            cell%ia_mpo(i)=6.0d0*alphadrp/cell%ddrp(i)
!            cell%ia_mpr(i)=6.0d0*alphadrp/cell%ddrp(i)    
!       
         ELSE   !LJR coding
!            xag=cell%alphag(i)
            xag=DMAX1(cell%alpha_sa(i),DMIN1(cell%alphag(i),alpha_am))
            xal=1.d0-xag
!
            d_star=0.0025d0
            IF(xag.eq.1.d0 .and. cell%quala(i).ne.0.d0) alphadrp=DMAX1(xal,1.d-3)
            IF(xag.ne.1.d0 .or.  cell%quala(i).eq.0.d0) alphadrp=DMAX1(xal,1.d-4)
            We_crit=1.5d0
            vfg=vrel_o(i)*vrel_o(i)
            vfg2=DMAX1(vfg,We_crit*sig/(rhog*DMIN1(d_star*alphadrp**(1.d0/3.d0),dh)))
!
            cell%ddrp(i)=DMAX1(We_crit*cell%sigma(i),1.d-10)/(rhog*vfg2)
            cell%ddrp(i)=DMIN1(dh,DMAX1(84.0d-6,cell%ddrp(i)))
!
            agf=3.6d0*alphadrp/cell%ddrp(i)
!            agf=6.0d0*alphadrp/cell%ddrp(i)            
            cell%ia_mpr(i)=agf
!
!reflod model by LSJ 180123
            IF(reflod) THEN
               IF(cell%quala(i).ne.0.0d0.and.cell%alphag(i).eq.1.0d0)then
                  alpdrp=max(cell%alphal(i),1.0d-3) 
               ELSE 
                  alpdrp=max(cell%alphal(i),1.0d-4)
               ENDIF  
               velgdf=dabs(vg_o(i,ndim)-vl_o(i,ndim))
               scrchh=velgdf
!            
               vfg=scrchh
               alp=alpdrp
!               
               vfg_hat_2=max(vfg*vfg,1.5d0*sig/(rhog*dmin1(0.0015d0*alp**0.333333333d0,dh)))
               vfg_hat_2=dmax1(vfg_ss*vfg_ss,vfg_hat_2)
!
               voidx=1.d0-alp               
               rey_reflod=1.5d0*sig*voidx**3/(cell%lviscosg(i)*dsqrt(vfg_hat_2*voidx))
!
               cell%ddrp(i)=1.5d0*sig/(rhog*vfg_hat_2)
               cell%ddrp(i)=dmax1(cell%ddrp(i),84.d-6) !Limit drop diameter to > 84 microns (TRAC)
               cell%ddrp(i)=DMIN1(cell%ddrp(i),dh)
!
               cell%ia_mpr(i)=3.6d0*alp/cell%ddrp(i) 
            ENDIF            
         ENDIF         
!
      ENDIF             
!      
      RETURN  
      END SUBROUTINE rv_ihtc_ia_preCHF
!
!---------------------------------------------------------------------
!
      SUBROUTINE rv_ihtc_ia_posDRY(i,reg_idx)
! 
      USE VOL_DATA 
      USE Zparam          , ONLY: ndim,pi   
      USE Zconst2         , ONLY: hydraulicd,grav,ggc
      USE Zvector         , ONLY: ul_o,ug_o,vrel_o
      USE Zvector         , ONLY: vl_o,vg_o
      USE STM_TBL_cupid   , ONLY: pcrit
      USE Zrv_ihtc_models , ONLY: F17
      USE Zwall_HTC       , ONLY: reflod,rey_reflod     
      USE Zrv_int_friction, ONLY: xa,dcon,web,drop_max      
!
      IMPLICIT NONE
!
      INTEGER i,reg_idx
!
      REAL(8) We_crit,dia_min,beta,xag,alphabub,dia_max,vfg,vfg2   
      REAL(8) xal,vfg_ss,vfg_hat_2
      REAL(8) alphaian,cF15,cF16,cF17,cF18,alphab
      REAL(8) cF21,alphadrp,dd,pstar,dmin
      REAL(8) vg,vl,del_rho
      REAL(8) vfg0
      REAL(8) rhol,rhog,dh,mu,sig
      REAL(8) d_star,agf
!
! reflod models
      REAL(8) vfgbb,alpian,scrchh,voidb,alpbub,alphbc,velgdf,alp,voidx,alphac,void,alpdrp
!      
!.....pcrit is constant read in stread.f90 why change?
!     REAL(8),PARAMETER :: pcrit=22.06d6 !critical pressure of H2O 
      REAL(8),PARAMETER :: alpha_am=0.9999d0                  
!
      LOGICAL, SAVE :: lsj_coding=.true.
!      DATA lsj_coding,ljr_mars /.false.,.true./      
!
!      IF(liquid.eq.1.or.gas.eq.1) THEN
!         RETURN
!      ENDIF
!
      ggc=DSQRT(DOT_PRODUCT(grav,grav))
      vg=ug_o(i)
      vl=ul_o(i)
      rhog=cell%rhog(i)
      rhol=cell%rhol(i)      
      del_rho=cell%rhol(i)-cell%rhog(i)
      del_rho=DMAX1(1.d-5,del_rho)
      dh=hydraulicd(i)
      drop_max=dh
      mu=cell%lviscosl(i)
      sig=cell%sigma(i)
!      cell%length(i)=dh
!
      IF(reflod)THEN
         dcon(1)=5.0d-3
!         dcon(2)=1.5d-3
!         dcon(2)=5.0d-4  ! MARS source, FLECHT SEASET min drp size: 2.5 mm / original: 1.5 mm / 31805 
         dcon(2)=2.5d-4  ! MARS source interphaseDrag Line 1962
         dcon(3)=2.0d-4
         web(1)=5.0d0
         web(2)=1.5d0
!         web(3)=1.5d0
!
!         use TRACE 5.0 We number of reflood We (based on dmax) =12.0
!         sauter mean diam = 1/3 * dmax
!         We = 3/1*12.0 = 4.0
!
         web(3)=4.0d0
      ELSE
         dcon(1)=5.0d-3
         dcon(2)=2.5d-3
         dcon(3)=2.0d-4
         web(1)=5.0d0
         web(2)=1.5d0
         web(3)=6.0d0          
      ENDIF    
      xa(:)=DMAX1(web(:)*sig,1.0d-10)
!
!.....Inverted annular
!
      IF(reg_idx.eq.8) THEN
!       
         IF(lsj_coding) THEN   !LSJ coding                
            xag=DMAX1(0.0,DMIN1(cell%alphag(i),cell%alpha_bs(i)))
!            
            We_crit = 5.d0
            dia_min = 0.005d0
!         
            alphaian=xag         
            cF17=F17(i,alphaian)
            alphab=cF17*alphaian
            alphabub=DMAX1((alphaian-alphab)/(1.d0-alphab),1.0D-7)
            cF16=DMIN1(1.d0-cF17,1.d0)
            cF15=DSQRT(1.d0-alphab)         
            beta=cF16
!         
            dia_max=DMIN1(dia_min*alphabub**(1.0/3.0),hydraulicd(i))
            vfg=(ug_o(i)-ul_o(i))*beta*beta                                                    
            vfg2=DMAX1(vfg*vfg,We_crit*cell%sigma(i)/cell%rhol(i)/dia_max)
            vfg=DSQRT(vfg2)
!
!           small bubble
            cell%dbb(i)=We_crit*cell%sigma(i)/cell%rhol(i)/vfg2 
            cell%ia_invann_sb(i)=3.6d0*alphabub*(1.d0-alphab)/cell%dbb(i)*beta
!
!-          annular film
            cell%ia_invann_ann(i)=9.0d0*cF15/hydraulicd(i)
!       
         ELSE   !LJR coding
!            xag=cell%alphag(i)
            xag=DMAX1(0.0,DMIN1(cell%alphag(i),cell%alpha_bs(i)))
!
!           bubble portion
            IF(xag.lt.cell%alpha_bs(i))then
               alphaian=xag
            ELSEIF(xag.lt.cell%alpha_cd(i))then
               alphaian=cell%alpha_bs(i)
            ENDIF
            cf18=dmin1(xag/0.05d0,0.99999d0)
            cf17=dexp(-8.d0*(cell%alpha_bs(i)-alphaian)/cell%alpha_bs(i))*cf18
            cf16=1.d0-cf17
            alphab=alphaian*cf17
            alphabub=dmax1((alphaian-alphab)/(1.d0-alphab),1.d-7)
!
            IF(xag.ge.1.d-5)then
               vfg0=vg-vl
            ELSE
               vfg0=(vg-vl)*xag*1.d-5
            ENDIF
            We_crit=5.d0
            vfg=We_crit*sig/(rhol*dmin1(0.005d0*alphabub**(1.d0/3.d0),dh))
            vfg2=dmax1(vfg0*vfg0,vfg)
!
            cell%dbb(i)=We_crit*sig/dmax1(1.d-8,(rhol*vfg2))
            cell%ia_invann_sb(i)=3.6d0*alphabub/cell%dbb(i)*(1.d0-alphab)*cf16
!
!           annular portion
            cf15=dsqrt(1.d0-alphab)
            cell%ia_invann_ann(i)=4.d0*cf15*2.5d0/dh
!
!reflod model by LSJ 180123
            IF(reflod) THEN
               alphab=cell%alpha_bs(i)
               IF(reflod) alphab=dmax1(alphab,0.25d0)
               alpian=alphaian
               alphbc=0.05d0
               vfgbb=dexp(-(alphab-alpian)*8.0d0/alphab)*dmin1(xag/alphbc,0.9999999d0) 
               voidb=alpian*vfgbb  
               scrchh=dsqrt(voidb) 
               alpbub=dmax1((alpian-voidb)/(1.0d0-voidb),1.0d-7) 
               vfgbb=1.0d0-vfgbb
               velgdf=dabs(vg_o(i,ndim)-vl_o(i,ndim))
               scrchh=velgdf*vfgbb*vfgbb 
!            
               vfg=scrchh
               alp=alpbub
               vfg_hat_2=dmax1(vfg*vfg,5.d0*sig/(rhol*dmin1(5.d-03*alp**0.333333333d0,dh)))
               vfg_hat_2=dmax1(vfg_ss*vfg_ss,vfg_hat_2)
!
               voidx=1.d0-alp   
               rey_reflod=5.d0*sig*voidx/(cell%lviscosl(i)*dsqrt(vfg_hat_2))
!
               cell%dbb(i)=5.d0*sig/dmax1(1.d-8,(rhol*vfg_hat_2))
               cell%ia_invann_sb(i)=3.6d0*alp/cell%dbb(i)
            ENDIF   
!                        
         ENDIF
      ENDIF   
!
!.....Inverted churn
!
      IF(reg_idx.eq.9) THEN 
!       
         IF(lsj_coding) THEN  !LSJ coding             
!           inverted annular
            We_crit = 5.d0
            dia_min = 0.005d0
            xag=cell%alpha_bs(i)
            alphaian=xag
!
            cF17=F17(i,alphaian)
            alphab=cF17*alphaian
            alphabub=DMAX1( (alphaian-alphab)/(1.d0-alphab),1.0D-7 )                  
            cF16=DMIN1(1.d0-cF17,1.d0)
            cF15=DSQRT(1.d0-alphab)         
            beta = cF16
            dia_max = DMIN1( dia_min * alphabub**(1.0/3.0), hydraulicd(i) )
            vfg = ( ug_o(i) - ul_o(i) ) * beta * beta                                                    
            vfg2 = DMAX1( vfg*vfg, We_crit * cell%sigma(i) / cell%rhol(i) / dia_max )
            vfg = DSQRT(vfg2)
!
!           small bubble
            cell%dbb(i)=We_crit * cell%sigma(i) / cell%rhol(i) / vfg2 
            cell%ia_invann_sb(i)=3.6d0*alphabub*(1.d0-alphab)/cell%dbb(i)*beta
!
!           annular film
            cell%ia_invann_ann(i)=9.0d0*cF15/hydraulicd(i)
!
!           inverted slug
!            We_crit = 6.d0
!            dia_min = 0.0025d0
            We_crit = web(3)
            dia_min = dcon(3)
!
            xag=cell%alpha_cd(i)
            xal=1.d0-xag
            cF21=DEXP( (cell%alpha_sa(i)-xag)/(cell%alpha_sa(i)-cell%alpha_bs(i)) )
            alphadrp=(1.d0-cell%alpha_sa(i))*cF21
            alphab=(xal-alphadrp)/(1.d0-alphadrp)
            beta = cF21
            vfg = ( ug_o(i) - ul_o(i) ) * beta * beta
            vfg2= vfg*vfg
!         
!           droplet diameter 
            dd=We_crit*cell%sigma(i)/cell%rhog(i)/vfg2 
            pstar=cell%p(i)/pcrit
!            IF(pstar.lt.0.025d0) THEN
!               dmin=0.0025d0
!            ELSEIF(pstar.gt.0.25d0) THEN
!               dmin=0.0002d0
!            ELSE
!               dmin=0.0025d0+(0.0002d0-0.0025d0)/(0.25d0-0.025d0)*(pstar-0.025d0)
!            ENDIF
!            cell%ddrp(i)=DMIN1(DMAX1(dd,dmin),hydraulicd(i),0.0025d0)
            IF(pstar.lt.0.25d0)THEN
               IF(pstar.lt.0.025d0)THEN
                  dmin=dcon(2)
               ELSE
                  dmin=dcon(2)+(dcon(3)-dcon(2))*4.44444d0*(pstar-0.025d0)
               ENDIF   
            ELSE
               dmin=dcon(3)
            ENDIF
            cell%ddrp(i)=DMIN1(DMAX1(dd,dmin),hydraulicd(i),drop_max) ! yjm
! 
!           droplets
            cell%ia_invslg_drp(i)=3.6d0*alphadrp*(1.d0-alphab)/cell%ddrp(i)
!
!           annular film
            cell%ia_invslg_ann(i)=11.25d0*alphab/hydraulicd(i)     
!       
         ELSE  !LJR coding
!           inverted annular         
            xag=cell%alpha_bs(i)
!
!           bubble portion
            IF(xag.lt.cell%alpha_bs(i))then
               alphaian=xag
            ELSEIF(xag.lt.cell%alpha_cd(i))then
               alphaian=cell%alpha_bs(i)
            ENDIF
            cf18=dmin1(xag/0.05d0,0.99999d0)
            cf17=dexp(-8.d0*(cell%alpha_bs(i)-alphaian)/cell%alpha_bs(i))*cf18
            cf16=1.d0-cf17
            alphab=alphaian*cf17
            alphabub=dmax1((alphaian-alphab)/(1.d0-alphab),1.d-7)
!
            IF(xag.ge.1.d-5)then
               vfg0=vg-vl
            ELSE
               vfg0=(vg-vl)*xag*1.d-5
            ENDIF
            We_crit=5.d0
            vfg=We_crit*sig/(rhol*dmin1(0.005d0*alphabub**(1.d0/3.d0),dh))
            vfg2=dmax1(vfg0*vfg0,vfg)
!
            cell%dbb(i)=We_crit*sig/dmax1(1.d-8,(rhol*vfg2))
            cell%ia_invann_sb(i)=3.6d0*alphabub/cell%dbb(i)*(1.d0-alphab)*cf16
!
!           annular portion
            cf15=dsqrt(1.d0-alphab)
            cell%ia_invann_ann(i)=4.d0*cf15*2.5d0/dh 
!
!reflod model by LSJ 180123
            IF(reflod) THEN
               alphab=cell%alpha_bs(i)
               IF(reflod) alphab=dmax1(alphab,0.25d0)
               alpian=alphaian
               alphbc=0.05d0
               vfgbb=dexp(-(alphab-alpian)*8.0d0/alphab)*dmin1(xag/alphbc,0.9999999d0) 
               voidb=alpian*vfgbb  
               scrchh=dsqrt(voidb) 
               alpbub=dmax1((alpian-voidb)/(1.0d0-voidb),1.0d-7) 
               vfgbb=1.0d0-vfgbb
               velgdf=dabs(vg_o(i,ndim)-vl_o(i,ndim))
               scrchh=velgdf*vfgbb*vfgbb 
!            
               vfg=scrchh
               alp=alpbub
               vfg_hat_2=dmax1(vfg*vfg,5.d0*sig/(rhol*dmin1(5.d-03*alp**0.333333333d0,dh)))
               vfg_hat_2=dmax1(vfg_ss*vfg_ss,vfg_hat_2)
!
               voidx=1.d0-alp   
               rey_reflod=5.d0*sig*voidx/(cell%lviscosl(i)*dsqrt(vfg_hat_2))
!
               cell%dbb(i)=5.d0*sig/dmax1(1.d-8,(rhol*vfg_hat_2))
               cell%ia_invann_sb(i)=3.6d0*alp/cell%dbb(i)
            ENDIF               
!
!           inverted slug
            xag=cell%alpha_cd(i)
            xal=1.d0-xag       
!
!           annular portion
            cf21=DEXP((xag-cell%alpha_sa(i))/(cell%alpha_sa(i)-cell%alpha_bs(i)))
            alphadrp=(1.d0-cell%alpha_sa(i))*cf21
            alphab=(xal-alphadrp)/(1.d0-alphadrp)
            cell%ia_invslg_ann(i)=4.5d0*alphab*2.5d0/dh
!
!            drop portion
!            vfg=vrel_o(i)*f21*f21                  !MARS manual
            vfg=DMAX1(vrel_o(i)*cf21*cf21,0.001d0) !MARS code
            dd=DMAX1(6.d0*sig,1.d-10)/(rhog*vfg*vfg)
!            
            pstar=cell%p(i)/pcrit
            IF(pstar.lt.0.025d0)then
               dmin=0.0025d0
            ELSEIF(pstar.gt.0.25d0)then
               dmin=0.0002d0
            ELSE
               dmin=0.0025d0-(0.0025d0-0.0002d0)/(0.25d0-0.025d0)*(pstar-0.025d0)
            ENDIF
            dd=DMAX1(dd,dmin)                           
            cell%ddrp(i)=DMIN1(dd,dh)                   ! < MARS code (bubbleDropSurface)
!            
            cell%ia_invslg_drp(i)=3.6d0*alphadrp*(1.d0-alphab)/cell%ddrp(i)  
!
!reflod model by LSJ 180123
            IF(reflod) THEN
               alphab=cell%alpha_bs(i)
               IF(reflod) alphab=dmax1(alphab,0.25d0)
               alphac=cell%alpha_sa(i)
               IF(reflod) alphac=dmax1(alphac,0.8d0)
               void=exp(-(alphac-xag)/(alphac-alphab)) 
               alpdrp=(1.0d0-alphac)*void 
               voidb=(xal-alpdrp)/(1.0d0-alpdrp) 
               velgdf=dabs(vg_o(i,ndim)-vl_o(i,ndim))
               scrchh=velgdf*void*void 
!            
               vfg=dmax1(dabs(scrchh),0.001d0) 
               alp=alpdrp
               vfg_hat_2=dmax1(vfg*vfg,1.5d0*sig/(rhol*dmin1(5.d-03*alp**0.333333333d0,dh)))
               vfg_hat_2=dmax1(vfg_ss*vfg_ss,vfg_hat_2)
!
               dd=1.5d0*sig/(rhog*vfg*vfg)
               pstar=cell%p(i)/pcrit
               IF(pstar.lt.0.025d0)then
                  dmin=0.0015d0
               ELSEIF(pstar.gt.0.25d0)then
                  dmin=0.0002d0
               ELSE
                  dmin=0.0015d0-(0.0015d0-0.0002d0)*4.44444d0*(pstar-0.025d0)
               ENDIF
               dd=DMAX1(dd,dmin)
               dd=dmin1(dh,0.0015d0,dd)                   
!
               rey_reflod=rhog*vfg*dd*((1.d0-alp)**2.5d0)/cell%lviscosg(i)
!
               cell%ddrp(i)=dd                              ! < MARS code (bubbleDropSurface)
               cell%ia_invslg_drp(i)=3.6d0*alp/cell%ddrp(i)
            ENDIF  
!                               
         ENDIF   
      ENDIF      
!
!.....Inverted slug
!
      IF(reg_idx.eq.10) THEN
!       
         IF(lsj_coding) THEN   !LSJ coding               
!            We_crit = 6.d0
!            dia_min = 0.0025d0
            We_crit = web(3)
            dia_min = dcon(3)
!            
            xag=DMAX1(cell%alpha_cd(i),DMIN1(cell%alphag(i),cell%alpha_sa(i)))
            xal=1.d0-xag
            cF21=DEXP( (cell%alpha_sa(i)-xag)/(cell%alpha_sa(i)-cell%alpha_bs(i)) )
            alphadrp=(1.d0-cell%alpha_sa(i))*cF21
            alphab=(xal-alphadrp)/(1.d0-alphadrp)
            beta = cF21
            vfg = ( ug_o(i) - ul_o(i) ) * beta * beta
            vfg2= vfg*vfg
!         
!           droplet diameter 
!            dd=We_crit*cell%sigma(i)/cell%rhog(i)/vfg2 
            dd=We_crit*cell%sigma(i)/cell%rhog(i)/DMAX1(1.0d-5,vfg2) 
            pstar=cell%p(i)/pcrit
!            IF(pstar.lt.0.025d0) THEN
!               dmin=0.0025d0
!            ELSEIF(pstar.gt.0.25d0) THEN
!               dmin=0.0002d0
!            ELSE
!               dmin=0.0025d0+(0.0002d0-0.0025d0)/(0.25d0-0.025d0)*(pstar-0.025d0)
!            ENDIF
!            cell%ddrp(i)=DMIN1(DMAX1(dd,dmin),hydraulicd(i),0.0025d0)
            IF(pstar.lt.0.25d0)THEN
               IF(pstar.lt.0.025d0)THEN
                  dmin=dcon(2)
               ELSE
                  dmin=dcon(2)+(dcon(3)-dcon(2))*4.44444d0*(pstar-0.025d0)
               ENDIF   
            ELSE
               dmin=dcon(3)
            ENDIF
            cell%ddrp(i)=DMIN1(DMAX1(dd,dmin),hydraulicd(i),drop_max)
! 
!           droplets
            cell%ia_invslg_drp(i)=3.6d0*alphadrp*(1.d0-alphab)/cell%ddrp(i)
!
!           annular film
            cell%ia_invslg_ann(i)=11.25d0*alphab/hydraulicd(i)
!       
         ELSE  !LJR coding
!            xag=cell%alphag(i)         
            xag=DMAX1(cell%alpha_cd(i),DMIN1(cell%alphag(i),cell%alpha_sa(i)))
            xal=1.d0-xag
!
!           annular portion
            cf21=DEXP((xag-cell%alpha_sa(i))/(cell%alpha_sa(i)-cell%alpha_bs(i)))
            alphadrp=(1.d0-cell%alpha_sa(i))*cf21
            alphab=(xal-alphadrp)/(1.d0-alphadrp)
            cell%ia_invslg_ann(i)=4.5d0*alphab*2.5d0/dh
!
!            drop portion
            We_crit = 6.d0
!            vfg=vrel_o(i)*f21*f21                  !MARS manual
            vfg=DMAX1(vrel_o(i)*cf21*cf21,0.001d0) !MARS code
            dd=DMAX1(We_crit*sig,1.d-10)/(rhog*vfg*vfg)
!            
            pstar=cell%p(i)/pcrit
            IF(pstar.lt.0.025d0)then
               dmin=0.0025d0
            ELSEIF(pstar.gt.0.25d0)then
               dmin=0.0002d0
            ELSE
               dmin=0.0025d0-(0.0025d0-0.0002d0)/(0.25d0-0.025d0)*(pstar-0.025d0)
            ENDIF
            dd=DMAX1(dd,dmin)         
            cell%ddrp(i)=DMIN1(dd,dh)                   ! < MARS code (bubbleDropSurface)
!            
            cell%ia_invslg_drp(i)=3.6d0*alphadrp*(1.d0-alphab)/cell%ddrp(i)
!
!reflod model by LSJ 180123
            IF(reflod) THEN
               alphab=cell%alpha_bs(i)
               IF(reflod) alphab=dmax1(alphab,0.25d0)
               alphac=cell%alpha_sa(i)
               IF(reflod) alphac=dmax1(alphac,0.8d0)
               void=exp(-(alphac-xag)/(alphac-alphab)) 
               alpdrp=(1.0d0-alphac)*void 
               voidb=(xal-alpdrp)/(1.0d0-alpdrp) 
               velgdf=dabs(vg_o(i,ndim)-vl_o(i,ndim))
               scrchh=velgdf*void*void 
!            
               vfg=dmax1(dabs(scrchh),0.001d0) 
               alp=alpdrp
               vfg_hat_2=dmax1(vfg*vfg,1.5d0*sig/(rhol*dmin1(5.d-03*alp**0.333333333d0,dh)))
               vfg_hat_2=dmax1(vfg_ss*vfg_ss,vfg_hat_2)
!
               dd=1.5d0*sig/(rhog*vfg*vfg)
               pstar=cell%p(i)/pcrit
               IF(pstar.lt.0.025d0)then
                  dmin=0.0015d0
               ELSEIF(pstar.gt.0.25d0)then
                  dmin=0.0002d0
               ELSE
                  dmin=0.0015d0-(0.0015d0-0.0002d0)*4.44444d0*(pstar-0.025d0)
               ENDIF
               dd=DMAX1(dd,dmin)
               dd=dmin1(dh,0.0015d0,dd)                   
!
               rey_reflod=rhog*vfg*dd*((1.d0-alp)**2.5d0)/cell%lviscosg(i)
!
               cell%ddrp(i)=dd                              ! < MARS code (bubbleDropSurface)
               cell%ia_invslg_drp(i)=3.6d0*alp/cell%ddrp(i)
            ENDIF               
         ENDIF
      ENDIF    

!
!.....MPO
!
      IF(reg_idx.eq.12) THEN   
!       
         IF(lsj_coding) THEN   !LSJ coding
!            We_crit = 6.d0
!            dia_min = 0.0002d0
!
            xag=DMAX1(cell%alpha_sa(i),DMIN1(cell%alphag(i),alpha_am))
            xal=1.d0-xag      
!         
            alphadrp=DMAX1(xal,1.d-4)         
!          
            vfg  = ( ug_o(i) - ul_o(i) )
            vfg2 = DMAX1( vfg*vfg, 1.d-6 )
!
!           small bubble
!            cell%ddrp(i) = We_crit*cell%sigma(i)/cell%rhog(i)/vfg2                      
!
!           droplets diameter
            cell%ddrp(i)=web(3)*cell%sigma(i)/cell%rhog(i)/vfg2    
            cell%ddrp(i)=DMIN1(dh,drop_max,DMAX1(84.0d-6,cell%ddrp(i)))  
! 
!           droplets
            cell%ia_mpo(i)=3.6d0*alphadrp/cell%ddrp(i)
!       
         ELSE   !LJR coding
!            xag=cell%alphag(i)
            xag=DMAX1(cell%alpha_sa(i),DMIN1(cell%alphag(i),alpha_am))
            xal=1.d0-xag
!
            d_star=0.0025d0
            IF(xag.eq.1.d0 .and. cell%quala(i).ne.0.d0) alphadrp=DMAX1(xal,1.d-3)
            IF(xag.ne.1.d0 .or.  cell%quala(i).eq.0.d0) alphadrp=DMAX1(xal,1.d-4)
            We_crit=1.5d0
            vfg=vrel_o(i)*vrel_o(i)
            vfg2=DMAX1(vfg,We_crit*sig/(rhog*DMIN1(d_star*alphadrp**(1.d0/3.d0),dh)))
!
            cell%ddrp(i)=DMAX1(We_crit*cell%sigma(i),1.d-10)/(rhog*vfg2)
            cell%ddrp(i)=DMIN1(dh,DMAX1(84.0d-6,cell%ddrp(i)))
!
            agf=3.6d0*alphadrp/cell%ddrp(i)
!            agf=6.0d0*alphadrp/cell%ddrp(i)            
            cell%ia_mpo(i)=agf
!
!reflod model by LSJ 180123
            IF(reflod) THEN
               alpdrp=max(cell%alphal(i),1.0d-4)
               velgdf=dabs(vg_o(i,ndim)-vl_o(i,ndim))
               scrchh=velgdf
!            
               vfg=scrchh
               alp=alpdrp
!               
!               dd=1.5d0*sig/(rhog*vfg*vfg)
               dd=1.5d0*sig/DMAX1(1.0d-4,rhog*vfg*vfg)
!               
               pstar=cell%p(i)/pcrit
               IF(pstar.lt.0.025d0)then
                  dmin=0.0015d0
               ELSEIF(pstar.gt.0.25d0)then
                  dmin=0.0002d0
               ELSE
                  dmin=0.0015d0-(0.0015d0-0.0002d0)*4.44444d0*(pstar-0.025d0)
               ENDIF
               dd=DMAX1(dd,dmin)
               dd=dmin1(dh,0.0015d0,dd)                   
!
               rey_reflod=rhog*vfg*dd*((1.d0-alp)**2.5d0)/cell%lviscosg(i)
!
               cell%ddrp(i)=dd                              ! < MARS code (bubbleDropSurface)
               cell%ia_mpo(i)=3.6d0*alp/cell%ddrp(i)
            ENDIF
!            
         ENDIF             
      ENDIF   
!
!
      RETURN  
      END SUBROUTINE rv_ihtc_ia_posDRY
!
