!
      SUBROUTINE condensation_ncg(i)
!
!     This routine calculates wall heat transfer coefficient for condensation near wall
!
      USE VOL_DATA      , ONLY: cell
      USE STM_TBL_cupid , ONLY: st_tbl,nfluid,ndxstd,nt
      USE Zncg          , ONLY: advn_cell,ra_cell,wmole_gas
      USE Zwall_HTC     , ONLY: tw,dia_rod,incnd,h_bundle,HTC_d,pvblk
      USE Zwall_HTC     , ONLY: f_direc,tw
      USE Zio_unit      , ONLY: unit_log
      USE Zvector       , ONLY: vg_n
!
      IMPLICIT NONE        
!                                                                                                                            
      REAL(8) advh2o,dc,dcnst,hdl,hdt,molh2o,pa,pb,pres_intf,prop(36)
      REAL(8) pw,reg,rggb,rmolg,sc,tol,x,hnat,gr,ray                  
      REAL(8) rhoai,rhosi,rhomi 
      REAL(8) gravity
      INTEGER i
      INTEGER ier,iprop(2,36),iter,nsig 
      LOGICAL erx,lprop(2,36) 
      EQUIVALENCE(iprop(1,1),lprop(1,1),prop(1)) 
      EXTERNAL  pres_intf
!                                                                       
      DATA erx/.false./ 
      DATA gravity/9.81d0/ 
!                                                                       
!     molh20 = mol.weight of water                                         
!     advh2o = atomic diffusion volume of water, ref.given in rnoncn       
      DATA molh2o,advh2o/18.0d0,12.7d0/ 
!                                                                       
!.....Calculate properties at tsat.      
!     
      pvblk=cell%pps(i)
      rmolg=8314.3d0     
      !wmole_gas(i),advn_cell(i),ra_cell(i): calculated in ncg_cell.f90
      rggb=(cell%p(i)-pvblk)/(ra_cell(i)*cell%tst(i))
!                                                                       
!.....Mass diffusivity is calculated using eqn.11-4.4 (Fuller) in the properties of gases and liquid by reid,praudnitz,sherwood,third ed. mc-graw-hill book co.,1977.                                
!                                                                       
      dcnst=(sqrt(1.0d0/molh2o+1.0d0/wmole_gas(i)))/((advh2o**0.33333d0+advn_cell(i)**0.33333d0)**2)
      dc=0.0101325d0*dcnst*cell%tg(i)**1.75d0/cell%p(i)
!                                                                       
!.....Calculate Reynolds number
!      
      reg=cell%rhog(i)*abs(vg_n(i,f_direc))*dia_rod/cell%eviscosg(i) 
!                                                                       
!.....For high flow - use Gilliland correlation.                           
!     sc = Schmidt number                                               
      sc=cell%lviscosg(i)/(cell%rhog(i)*dc)
      hdt=0.023d0*(dc/dia_rod)*(reg**0.83d0)*(sc**0.44d0) 
!                                                                       
!.....For laminar: Rohsenow-Choi correlation.(Heat and mass transfer analogy)     
!      
      hdl=4.0d0*dc/dia_rod 
!                                                                       
!.....Calculate saturation pressure at tw
!      
      prop(1)=tw 
      prop(9)=1.d0 
      IF(nfluid.eq.1)THEN 
         CALL sth2x1_cupid(st_tbl(ndxstd),st_tbl(ndxstd+nt),prop,erx) 
      ELSEIF(nfluid.eq.2)THEN 
         CALL std2x1_cupid(st_tbl(ndxstd),prop,erx) 
      ELSEIF(nfluid.eq.15)THEN
         CALL nth2x1_cupid(st_tbl(ndxstd),prop,erx) 
      ELSE 
         WRITE(unit_log,3010) 
 3010 FORMAT  ('0******** System fluid is not h2o, d2o, or h2on')
         RETURN 
      ENDIF 
      pw=prop(2) 
      rhosi=1/prop(3) 
      rhoai=(cell%p(i)-pw)/(ra_cell(i)*tw) 
      rhomi=rhosi+rhoai 
!
!.... Use the density difference as the driving potential. See Eckart & Drake 2nd edition page 474.
!      
      gr=abs(rhomi-cell%rhog(i))/rhomi
      gr=gravity*gr*cell%rhog(i)**2*h_bundle**3/cell%lviscosg(i)**2 
      ray=gr*sc 
!
!.....Churchhill-Chu from Incropera & DeWitt page 501 vertical bodies; eq.(926)
!      
      hnat=0.387d0*ray**0.1666667d0 
      hnat=hnat/(1.0d0+(0.492d0/sc)**0.5625d0)**0.296296d0 
      hnat=(0.825d0+hnat)**2 
      hnat=dc*hnat/h_bundle
      HTC_d=max(hdt,hdl,hnat)  !Calculate max of natural, laminar and turbulent values. 
!                                                                       
      IF(erx) GOTO 1000 
!                                                                       
!.....Set convergence criteria and bound values for Brent_bi   
!
      incnd=0 
      tol=0.0005d0 
      nsig=8 
      iter=20 
      pa=pw 
      pb=cell%pps(i)
      CALL Brent_bi(pres_intf,tol,nsig,pa,pb,iter,ier,i) 
      IF(ier.ne.0) GOTO 1010 
!                                                                       
      incnd=1 
      x=pres_intf(pb,i) 
      RETURN 
!                                                                       
 1000 CONTINUE 
      WRITE(unit_log,2020)tw,i      
 2020 FORMAT  ('0******** Subroutine sth2x1, std2x1, or stpu2t returned an error flag for temperature = ',1pg13.5,' in volume = ',i10,'.')
      RETURN 
!                                                                       
 1010 CONTINUE 
      incnd=2       
      
      RETURN
      END SUBROUTINE condensation_ncg      
!
!
!   
      FUNCTION pres_intf(pvi,i) 
!                                                                       
!.....Calculate the difference between heat flux (qfluxo) from vapor-gas mixture to liquid film         
!     and heat fluxv (phiv) from liquid film to wall for a given pressure pvi. (jjj 10/31/1997)                            
!                                                                       
!
      USE VOL_DATA      , ONLY: cell
      USE STM_TBL_cupid , ONLY: st_tbl,pcrit,nfluid,ndxstd,nt
      USE Zwall_HTC     , ONLY: incnd,HTC_tg,HTC_cond,HTC_d,qflux_t,dia_rod,f_direc,qual_eq,tw,mflux_liqa,mflux_tota
      USE Zio_unit      , ONLY: unit_log    
      USE Zvector       , ONLY: vg_n
!
      IMPLICIT none 
!                                                                       
      REAL(8) pres_intf,pvi,gravity,Re_film,film_thick
      REAL(8) cpfi,ftr,hdb,hf,hfgi,pratio,pr, &
      qual_max,reyi,rhofi,rhogi,rhovb,s(36),thcofi,tsati,viscfi,viscgi,z
      REAL(8) tsati1(1),pvi1(1),rhofi1(1),rhogi1(1),thcofi1(1),viscfi1(1),viscgi1(1)
      EQUIVALENCE(tsati,tsati1)
      EQUIVALENCE(rhofi,rhofi1)
      EQUIVALENCE(rhogi,rhogi1)
      EQUIVALENCE(thcofi,thcofi1)
      EQUIVALENCE(viscfi,viscfi1)
      EQUIVALENCE(viscgi,viscgi1)
      REAL(8) htcnon,hturb,hshah,qfluxv
      INTEGER iones(2),i
      LOGICAL ls(2,36) 
      EQUIVALENCE(s(1),ls(1,1)) 
!                                                                       
      LOGICAL erx 
      DATA erx/.false./,iones/1,1/ 
      DATA gravity/9.81d0/    
!
!.....Calculate properties at pvi(liquiq film - vapor interface)
!                                                                       
      s(2)=pvi 
      s(9)=cell%quals(i)
      IF(nfluid.eq.1)then 
         CALL sth2x2_cupid(st_tbl(ndxstd),st_tbl(ndxstd+nt),s,erx) 
      ELSEIF(nfluid.eq.2)then 
         CALL std2x2_cupid(st_tbl(ndxstd),s,erx)
      ELSEIF(nfluid.eq.15)then 
         CALL nth2x2_cupid(st_tbl(ndxstd),s,erx) 
      ENDIF 
      IF(erx) GOTO 1000 
!                                                                       
      rhofi=1.0d0/s(11) 
      rhogi=1.0d0/s(12) 
      tsati=s(1) 
      hfgi=s(16)-s(15) 

      pvi1(1)=pvi
      CALL thcond_cupid(i,nfluid,iones,iones,1,1,tsati1,pvi1,rhofi1,'liquid',thcofi1,erx)
      CALL viscos_cupid(i,nfluid,iones,iones,1,1,tsati1,rhofi1,'liquid',viscfi1,erx)
      CALL viscos_cupid(i,nfluid,iones,iones,1,1,tsati1,rhogi1,'vapor',viscgi1,erx)
      cpfi=s(21) 
!                                                                       
!.....condensation cal. within Brent_bi in sub. noncnd                        
!     
      HTC_tg=0.0d0 
      qflux_t=0.0d0
      htcnon=0.0d0       
      hturb=0.0d0       
      qfluxv=0.0d0 
      IF(tw.ge.tsati) GOTO 350
!
!.....Condnsation inside a vertical pipe or plate.                                   
!                                                                       
      Re_film=DABS(mflux_liqa)*dia_rod/viscfi
      film_thick=0.9086d0*(Re_film*(viscfi/rhofi)**2/gravity)**0.333333d0
      film_thick=DMAX1(film_thick,1.0d-5)
      HTC_tg=thcofi/film_thick   !Nusselt HTC.                                                                   
      HTC_tg=DMAX1(HTC_tg,4.36d0*thcofi/dia_rod) 
!
!.....Horizontal stratified condensation --- !!!cyj: temporary Inactivation because CUPID cannot automatically determine the angle of each cell
!   
      !HTC_tg=0.296d0*((rhofi*DMAX1(rhofi-rhogi),0.0d0)*gravity*hfgi*thcofi**3)/(dia_rod*viscfi*DMAX1((tsati-tw),1.0d0))**0.25d0      
!                                                                       
      htcnon=HTC_tg 
!   
  100 IF(abs(vg_n(i,f_direc).le.1.d-3)) GOTO 200 
!                                                                       
!.....turbulent condensation heat transfer correlation
!                                   
      qual_max=DMAX1(1.0d-9,DMIN1(1.0d0,qual_eq))
      reyi=mflux_tota*dia_rod/viscfi 
      pr=viscfi*cpfi/thcofi 
      hdb=0.023d0*thcofi*reyi**0.8d0*pr**0.4d0/dia_rod 
      hf=hdb*(1.0d0-qual_max)**0.8d0 
      z=(pvi/pcrit)**0.4d0*(1.0d0/qual_max-1.0d0)**0.8d0 
      ftr=1.0d0 
      IF(z.ne.0.0d0)ftr=1.0d0+3.8d0/z**0.95d0 
      hshah=hf*ftr 
      hturb=hshah 
!                                                                                                     
      HTC_tg=DMAX1(HTC_tg,hturb) 
      htcnon=HTC_tg 
!                                                                       
  200 CONTINUE 
      qflux_t=htcnon*(tw-tsati) 
!                                                                       
      IF(incnd.eq.1) GOTO 300 
      GOTO 350 
!                                                                       
!.....final call to pres_intf from noncnd    
!      
  300 CONTINUE 
      HTC_cond=htcnon 
      qflux_t=HTC_cond*(tw-tsati) 
!                                                                       
  350 CONTINUE
!                                                                       
!.....vapor calculation       
!    
      pratio=(cell%p(i)-pvi)/(cell%p(i)-cell%pps(i)) 
      pratio=DMAX1(pratio,1.0d-9) 
      rhovb=(1.0d0-cell%quala(i))*cell%rhog(i) 
      qfluxv=-HTC_d*(cell%hgsat(i)-cell%hlsat(i))*rhovb*log(pratio) 
!                                                                       
  400 pres_intf=qflux_t-qfluxv 
!                                                                                           
      IF(tw.gt.tsati)incnd=2 
      RETURN 
!                                                                       
 1000 CONTINUE 
      WRITE(unit_log,2000)pvi 
 2000 FORMAT  ('0******** Subroutine sth2x2 returned an error flag for pvi =',1pe14.7,'.  Called from pres_intf.')
      pres_intf=0.0d0 
      RETURN 
      END FUNCTION pres_intf     
!
!
!   
      SUBROUTINE Brent_bi(func,tol,nsig,x1,x2,itmax,ier,i) 
!                                                                   
!.....Using Brent's method, find root of a function func known to lie between x1,and x2.
!     The root, return as x2, will be refined until its accuracy is tol.   
!      -itmax .. maximum allowed number of iteration.                        
!      -nsig  .. significant digial of machine floading point prescision.    
!     Ref: 'Numerical Recipes - The Art of Scientific Computing' by W. H. Press, et. el., Cambrige University Press. 1986                
!     The subroutine is modified from the function zbrent in the book. (jjj 10/31/1997)       
!                                                                       
      IMPLICIT none 
!      
      REAL(8) a,b,c,d,e,eps,fa,fb,fc,func,p,q,r,s,tol,tol1,xm,x1,x2 
      INTEGER ier,iter,itmax,nsig,i
      EXTERNAL  func
!                                                                       
      ier=0 
      eps=3.0d0 * 10.0d0**(-nsig) 
      a=x1 
      b=x2 
      fa=func(a,i) 
      fb=func(b,i) 
      IF(fb*fa.gt.0.d0)then 
         ier=1 
         RETURN 
      ENDIF 
      fc=fb 
      DO 11 iter=1,itmax 
         IF(fb*fc.gt.0.d0)then 
            c=a 
            fc=fa 
            d=b-a 
            e=d 
         ENDIF 
         IF(abs(fc).lt.abs(fb))then 
            a=b 
            b=c 
            c=a 
            fa=fb 
            fb=fc 
            fc=fa 
         ENDIF 
         tol1=2.d0*eps*abs(b)+0.5d0*tol 
         xm=.5d0*(c-b) 
         IF(abs(xm).le.tol1.or.fb.eq.0.d0)then 
            x2=b 
            RETURN 
         ENDIF 
         IF(abs(e).ge.tol1.and.abs(fa).gt.abs(fb))then 
            s=fb/fa 
            IF(a.eq.c)then 
               p=2.d0*xm*s 
               q=1.d0-s 
            ELSE 
               q=fa/fc 
               r=fb/fc 
               p=s*(2.d0*xm*q*(q-r)-(b-a)*(r-1.d0)) 
               q=(q-1.d0)*(r-1.d0)*(s-1.d0) 
            ENDIF 
            IF(p.gt.0.d0)q=-q 
            p=abs(p) 
            IF(2.d0*p.lt.min(3.d0*xm*q-abs(tol1*q),abs(e*q)))then 
               e=d 
               d=p/q 
            ELSE 
               d=xm 
               e=d 
            ENDIF 
         ELSE 
            d=xm 
            e=d 
         ENDIF 
         a=b 
         fa=fb 
         IF(abs(d).gt.tol1)then 
            b=b+d 
         ELSE 
            b=b+sign(tol1,xm) 
         ENDIF 
         fb=func(b,i) 
   11 END DO 
      ier=2     
      x2=b 
      RETURN 
      END SUBROUTINE Brent_bi        
