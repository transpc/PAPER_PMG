!
      SUBROUTINE thcond_cupid(icell,ifluid,ihld,ixhld,nh,iscskp,temp,pres,rho,state,thcon,err)
!
!     bsw-transp : copy from MARS. 2006.4.5                                                                       
!     thcond_cupid  - compute thermal conductivity for given fluid           
!                                                                       
!     Author:   J. E. Tolli, EG&G Idaho, Inc.                          
!     Date:     7/89                                                   
!     Modified: 2/90  (helium,nitrogen,sodium,NaK)                     
!     Modified: 10/90 (nitrogen vapor table look-up)                   
!     Modified: 3/91  (lithium-lead)                                   
!     Language: FORTRAN77                                              
!                                                                       
!     Modified by Won-Jae Lee for inclusion of D2O, CO2 and He Properties                                                                
!                                                                       
!     Parameters:                                                      
!        IFluid = fluid number (input)                                  
!        ihld   = primary list vector holding subscripts for required   
!                 arrays (input)                                        
!        ixhld  = secondary list vector (USEs values stored in ihld for 
!                 subscripts) containing pointers to needed values in   
!                 alternate DATA base (input)                           
!        ishld  = scratch array holding pointers for table interpolation
!                 (USEs same indexing as ihld) (output)                 
!        nh     = vector length (input)                                 
!        iscskp = skip factor for ihld and ixhld (input)                
!        temp   = array containing temperatures (input)                 
!        pres   = array containing pressures (input)                    
!        rho    = array containing fluid densities (input)              
!        state  = fluid state (input)                                   
!                 'liquid' = liquid state                               
!                 'vapor'  = vapor state                                
!        thcon  = array containing thermal conductivities (output)      
!        err    = error flag (output)                                   
!                                                                       
! --------------------------------------------------------------  
!   Modified by Won-Jae Lee on Oct '02 for CO2 calculations
!
      USE VOL_DATA       , ONLY: cell
      USE STM_TBL_cupid  , ONLY: st_tbl,       &
                                 ndxstd,nfluid
!      
      IMPLICIT NONE 
!
      REAL(8) p_stm(36),vis
      REAL(8) vis1(1)
      EQUIVALENCE(vis,vis1)
      INTEGER idum,iones(2),ip,icell
      DATA iones/1,1/
!      
! --------------------------------------------------------------  
!                                                                       
!  Arguments:                                                           
      REAL(8) pres(*),rho(*),temp(*),thcon(*),rhoc  
      INTEGER ifluid,ihld(2,*),iscskp,ixhld(2,*),nh 
      LOGICAL err,erx 
      CHARACTER*(*)state 
!                                                                       
!  Local variables:                                                     
      INTEGER i,ih,ix,m 
      REAL(8) ag0,ag1,ag2,ag3,bg0,bg1,bg2,cc,x1,x2
      INTEGER nptct 
      PARAMETER (nptct=13) 
!     REAL(8) c(4),thcont(2,nptct) 
      REAL(8) c(4)
!  Loca variable defined by LWJ
      REAL(8) convt,term,tk
      INTEGER j
!                                                                       
!  Data statements:                                                     
!     DATA thcont/0.000d0,0.566989d0,273.150d0,0.566989d0,310.928d0,    &
!     0.629295d0,366.483d0,0.679140d0,394.261d0,0.685371d0,422.039d0,   &
!     0.685371d0,449.817d0,0.679140d0,477.594d0,0.660448d0,533.150d0,   &
!     0.604373d0,588.706d0,0.517144d0,616.483d0,0.448606d0,647.039d0,   &
!     0.274148d0,1.0d+75,0.274148d0/                                    
!
!   Data for D2O Thermal Conductivity
!
      DATA c/0.56340135d0,0.14504443d-2,-0.79650470d-5,0.71584948d-8/ 
!                                                                       
      DATA cc/2.1482d+5/,ag0/1.76d-2/,ag1/5.87d-5/,ag2/1.04d-7/,ag3/-   &
      4.51d-11/,bg0/1.0351d-4/,bg1/.4198d-6/,bg2/-2.771d-11/            
!       
!--------------------------------------------------------------------------
!  CO2 Thermal Conductivity - Modified by Won-Jae Lee on Oct. '02
!                           
      REAL(8) bc(8),dc(4),ci(5),t,p
      REAL(8) t1(1),p1(1)
      EQUIVALENCE(t,t1)
      EQUIVALENCE(p,p1)
      REAL(8) cp,cv,ro,beta,capa,tong,cint0,cint,r2,tstar
      REAL(8) ro1(1)
      EQUIVALENCE(ro,ro1)
!     REAL(8) ttcrit,rorocrit,ppcrit,ez,emew,egamma,r,eta,eta0,gam
      REAL(8) rorocrit,ppcrit,emew,egamma,r,eta,eta0,gam
      REAL(8) tref,caparef,qd,qdinv,pi,boltzmann,chi,chiref,delchi
      REAL(8) omega,omega0,thcon0,delthcon,critthcon
      DATA bc /0.4226159d+0,0.6280115d+0,-0.5387661d+0,0.6735941d+0,    &
               0.0d0,0.0d0,-0.4362677d+0,0.2255388d+0/
      DATA ci /2.387869d-2,4.350794d0,-10.33404d0,7.98159d0,-1.940558d0/
      DATA dc /2.447164d-02,8.705605d-05,-6.547950d-08,6.594919d-11/
!     definition of constants (specified for correlation)\
!     critical parameters
!      DATA ttcrit,rorocrit,ppcrit /304.107d0,467.69d0,7372100.d0/
!     DATA ttcrit,rorocrit,ppcrit /304.1282d0,467.6d0,7377300.d0/
      DATA rorocrit,ppcrit /467.6d0,7377300.d0/
      ! critical exponents
!     DATA ez,emew,egamma /0.06d0,0.630d0,1.2415d0/
      DATA emew,egamma /0.630d0,1.2415d0/
      ! critical amplitudes
      DATA r,eta0,gam,boltzmann /1.01d0,1.5d-10,0.052d0,1.380662d-23/
      ! cutoff wavenumber
      DATA qdinv /4.0d-10/
      ! Reference temperature for Delchi
      DATA tref /450.d0/
!---------------------------------------------------------------------------
!  Helium Thermal Conductivity - Upgraded by Won-Jae Lee on Jun. '03
!
      REAL(8) cd(4),cex(11),sum
      INTEGER jj
      DATA cd/3.739232544d0,-26.20316969d0,59.82252246d0,-49.26397634d0/
      DATA cex/ 0.186297053d-3,-0.1427549651d-3,0.3290833592d-4,   &
              -0.7275964435d-6,-0.5213335363d-7,0.4492659933d-7,   &
              -0.5924416513d-8,0.7087321137d-5,-0.6013335678d-5,   &  
               0.8067145814d-6,0.3995125013d-6/
!
!---------------------------------------------------------------------------
!  Coefficients of Collision Integral for N2, O2, Ar and Air
      REAL*8 bi(5),delt,ohm,tstar1,tau,visc0,delvis,tvisc   ! coefficients of collision integral (N2, O2, Ar, Air)
      DATA bi/0.431d0,-0.4623d0,0.08406d0,0.005341d0,-0.00331d0/
      INTEGER ii
!  
!  Nitrogen Viscosity
!     REAL*8 moln2
!     DATA moln2/28.013482377d0/
      REAL*8 tcn2,rocn2,pcn2,mn2,ekn2,sign2,ksin2,gamn2,qdn2,trefn2
      DATA tcn2,rocn2,pcn2,mn2,ekn2,sign2,ksin2,gamn2,qdn2,trefn2/126.192d0,11.1839d0,3.3958d0,28.01348d0,98.94d0,0.3656d0,0.17d0,0.055d0,0.4d0,252.384d0/
!
      REAL*8 vnin2(5),vtin2(5),vgamin2(5)
      INTEGER vdin2(5),vlin2(5)
      DATA vnin2,vtin2,vdin2,vgamin2,vlin2/10.72d0,0.03989d0,0.001208d0,-7.402d0,4.62d0,0.1d0,0.25d0,3.2d0,0.9d0,0.3d0,&
      2,10,12,2,1,0.d0,1.d0,1.d0,1.d0,1.d0,0,1,1,2,3/
!
      REAL*8 nin2(9),tin2(9),gamin2(9)
      INTEGER din2(9),lin2(9)
      DATA nin2,tin2,din2,gamin2,lin2/1.511d0,2.117d0,-3.332d0,8.862d0,31.11d0,-73.13d0,20.03d0,-0.7096d0,0.2672d0,&
      0.0d0,-1.0d0,-0.7d0,0.0d0,0.03d0,0.2d0,0.8d0,0.6d0,1.9d0,0,0,0,1,2,3,4,8,10,0.d0,0.d0,0.d0,0.d0,0.d0,1.d0,1.d0,1.d0,1.d0,0,0,0,0,0,1,2,2,2/
!  Oxygen Viscosity
      REAL*8 tco2,roco2,pco2,mo2,eko2,sigo2,ksio2,gamo2,qdo2,trefo2
      DATA tco2,roco2,pco2,mo2,eko2,sigo2,ksio2,gamo2,qdo2,trefo2/154.581d0,13.63d0,5.043d0,31.9988d0,118.5d0,0.3428d0,0.24d0,0.055d0,0.51d0,309.162d0/
!     REAL*8 molo2
!     DATA molo2/31.9988d0/
!
      REAL*8 vnio2(5),vtio2(5),vgamio2(5)
      INTEGER vdio2(5),vlio2(5)
      DATA vnio2,vtio2,vdio2,vgamio2,vlio2/17.67d0,0.4042d0,0.0001077d0,0.3510d0,-13.67d0,0.05d0,0.0d0,2.1d0,0.0d0,0.5d0,1,5,12,8,1,0.d0,0.d0,0.d0,1.d0,1.d0,0,0,0,1,2/
!
      REAL*8 nio2(9),tio2(9),gamio2(9)
      INTEGER dio2(9),lio2(9)
      DATA nio2,tio2,dio2,gamio2,lio2/1.036d0,6.283d0,-4.262d0,15.31d0,8.898d0,-0.7336d0,6.728d0,-4.374d0,-0.4747d0,0.0d0,-0.9d0,-0.6d0,0.d0,0.d0,0.3d0,4.3d0,0.5d0,1.8d0,&
      0,0,0,1,3,4,5,7,10,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,1.d0,1.d0,1.d0,0,0,0,0,0,0,2,2,2/
!  Argon Viscosity
!     REAL*8 molar
!     DATA molar/39.948d0/
      REAL*8 tcar,rocar,pcar,mar,ekar,sigar,ksiar,gamar,qdar,trefar
      DATA tcar,rocar,pcar,mar,ekar,sigar,ksiar,gamar,qdar,trefar/150.687d0,13.40743d0,4.863d0,39.948d0,143.2d0,0.335d0,0.13d0,0.055d0,0.32d0,301.374d0/
!
      REAL*8 vniar(6),vtiar(6),vgamiar(6)
      INTEGER vdiar(6),vliar(6)
      DATA vniar,vtiar,vdiar,vgamiar,vliar/12.19d0,13.99d0,0.005027d0,-18.93d0,-6.698d0,-3.827d0,0.42d0,0.0d0,0.95d0,0.5d0,0.9d0,0.8d0,&
      1,2,10,5,1,2,0.d0,0.d0,0.d0,1.d0,1.d0,1.d0,0,0,0,2,4,4/ 
!
      REAL*8 niar(10),tiar(10),gamiar(10)
      INTEGER diar(10),liar(10)
      DATA niar,tiar,diar,gamiar,liar/0.8158d0,-0.4320d0,0.d0,13.73d0,10.07d0,0.7375d0,-33.96d0,20.47d0,-2.274d0,-3.973d0,&
      0.0d0,-0.77d0,-1.0d0,0.d0,0.d0,0.d0,0.8d0,1.2d0,0.8d0,0.5d0,3*0,1,2,4,5,6,9,1,6*0.d0,4*1.d0,6*0,2,2,2,4/ 
!  Air Viscosity
!     REAL*8 molair
!     DATA molair/28.013482377d0/
      REAL*8 tcair,rocair,pcair,mair,ekair,sigair,ksiair,gamair,qdair,trefair
      DATA tcair,rocair,pcair,mair,ekair,sigair,ksiair,gamair,qdair,trefair/132.6312d0,10.4477d0,3.78502d0,28.9586d0,103.3d0,0.360d0,0.11d0,0.055d0,0.31d0,265.262d0/
!
      REAL*8 vniair(5),vtiair(5),vgamiair(5)
      INTEGER vdiair(5),vliair(5)
      DATA vniair,vtiair,vdiair,vgamiair,vliair/10.72d0,1.122d0,0.002019d0,-8.876d0,-0.02916d0,0.2d0,0.05d0,2.4d0,0.6d0,3.6d0,1,4,9,1,8,0.d0,0.d0,0.d0,1.d0,1.d0,0,0,0,1,1/
!
      REAL*8 niair(9),tiair(9),gamiair(9)
      INTEGER diair(9),liair(9)
      DATA niair,tiair,diair,gamiair,liair/1.308d0,1.405d0,-1.036d0,8.743d0,14.76d0,-16.62d0,3.793d0,-6.412d0,-0.3778d0,0.d0,-1.1d0,&
      -0.3d0,0.1d0,0.d0,0.5d0,2.7d0,0.3d0,1.3d0,3*0,1,2,3,7,7,11,5*0.d0,4*1.d0,5*0,2,2,2,2/
!
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!                                                                       
!  Execution:                                                           
!--initialize error flag                                                
      err=.false. 
!                                                                       
!--check for valid state specifier                                      
      IF(state.ne.'liquid'.and.state.ne.'vapor')then 
         err=.true. 
         RETURN 
      ENDIF 
!                                                                       
!--initialize pointer index                                             
      ih=1 
!                                                                       
!--branch to correlation for given fluid    
!    Light Water: 10
!    Heavy Water: 20
!    CO2        : 30
!    He         : 40
!    H2         : 50
!    O2         : 60
!    N2         : 70
!    Na         : 80
!    Ar         : 90
!    Air        : 100
!    LBE        : 110
!                            
      IF    (ifluid.eq. 1 .or. ifluid.eq.15) THEN
!
!--light water                                                          
!
         i=ihld(2,ih) 
         ix=ixhld(2,i) 
         IF (state.eq.'liquid') then
             tk=temp(i)
             rhoc=rho(ix)
         ELSE
            tk=temp(ix)
            rhoc=rho(ix)
         ENDIF
         CALL cond_lw_single(tk,rhoc,convt)
         thcon(i)=dmax1(convt,0.01d0) 
         ih=ih+iscskp 
      ELSEIF(ifluid.eq. 2) THEN
!
!--heavy water
!
         IF(state.eq.'liquid')then 
!                                                                       
!--saturated or subcooled liquid;  correlation from Flowtran program    
!--(Savannah River);  FORTRAN coding by R. J. Wagner, C. S. Miller,     
!--and J. E. Tolli, EG&G Idaho, Inc.                                    
!  Factor of 0.001 moved from correlation into coefficient (c).         
!
!--saturated or subcooled liquid;  correlation from Flowtran program    
!--(Savannah River);  FORTRAN coding by R. J. Wagner, C. S. Miller,     
!--and J. E. Tolli, EG&G Idaho, Inc.                                    
!  Factor of 0.001 moved from correlation into coefficient (c).         
            DO m=1,nh 
               i=ihld(2,ih) 
               term=max(temp(i)-273.15d0,1.0d-6) 
               thcon(i)=(((c(4)*term+c(3))*term+c(2))*term+c(1)) 
               ih=ih+iscskp 
            END DO 
         ELSE 
!                                                                       
!--saturated or superheated vapor;  FORTRAN coding by R. J. Wagner,     
!--C. S. Miller, and J. E. Tolli, EG&G Idaho, Inc.                      
            DO m=1,nh 
               i=ihld(2,ih) 
               ix=ixhld(2,i) 
               term=max(temp(ix)-273.15d0,1.0d-6) 
               x1=((ag3*term+ag2)*term+ag1)*term+ag0 
               x2=(bg2*term+bg1)*term+bg0 
               thcon(i)=x1+rho(ix)*(x2+cc*rho(ix)*term**(-4.2d0)) 
               ih=ih+iscskp 
            END DO
         ENDIF 
      ELSEIF(ifluid.eq. 3) THEN
!---------------------------------------------------------------------------
!
!   CO2 Thermal Conductivity
!       written by Won-Jae Lee on Nov. '02
!
!
!  reference
!    Vesovic, V. et. al.
!    "The Transport Properties of Carbon Dioxide"
!    J. of Physical and Chemistry Reference Data
!    Vol. 19: 763-808, 1990
!
!    Simplified Form
!        
!         thcon = thcon0 + delthcon + critthcon
!                thcon0:    conductivity in zero-density limit eq(29)
!                delthcon:  excess conductivity outside critical region eq(63)
!		 critthcon: conductivity at critical region eq(58)-(61)
!
!         specific heat (isobaric and isochoric)
!                Cp : from steam table
!                Cv : Cv=Cp-(beta)**2*T*v/(kapa) (Gas Dynamics, Ernst Becker (1968))
!
!  units: t(K), ro(kg/m3), thcon*(W/m/K), p(Pa)                                               
!
         DO m=1,nh
            i=ihld(2,ih)
            ix=ixhld(2,i)
            t=temp(i)
            p=pres(i)
            ro=rho(ix)
            tstar=t/251.196d0 
!     Thermal Conductivity in Zero-density Limit
      ! calculate r**2
            cint0=0.0d0
            DO j=1,5
               cint0=cint0+ci(j)*(t/100.d0)**(2-j)
            END DO
            cint=1.d0+dexp(-183.5d0/t)*cint0
            r2=0.4d0*cint
      ! calculate tong        	
            tong=0.0d0
            DO j=1,8
               tong=tong+bc(j)/tstar**(j-1)
            END DO
            thcon0=0.475598d0*dsqrt(t)*(1.d0+r2)/tong
!     Excess Viscosoty
            delthcon=0.0d0
            DO j=1,4
               delthcon=delthcon+dc(j)*ro**j
            END DO
!     Critical Region Conductivity
            IF(ro.lt.1.d-6) then
               critthcon=0.0d0
               GOTO 35
            ENDIF
      ! steam table search
            p_stm(1)=t
            P_stm(3)=1.d0/ro
!
!----------------------------------------------------------------------- 
!    Find beta, capa, Cp for the calculation of Cv vs. T, ro
!
!!    For superheat or subcooled fluid over 0.1 K (arbitrary)
!        IF(abs((t-v_da(i)%satt)).ge.0.1d0) then
            call strtro_cupid(st_tbl(ndxstd),p_stm,ip,erx)
            beta=p_stm(6)
            capa=p_stm(7)
            cp=p_stm(8)
      ! calculate cv
            cv=cp-beta**2*t/(capa*ro)
cv=max(cv,0.0d0)  ! prevent negative values *ARBITRARY*
      ! calculate delchi
        ! steam table search at Tref and given specific volume
            p_stm(1)=tref
           p_stm(3)=1.d0/ro  
!          CALL std2x4(st_tbl(ndxstd),p_stm,ip,erx)  ! replace std2o
           CALL strtro_cupid(st_tbl(ndxstd),p_stm,ip,erx)
           caparef=p_stm(7)
!  calculate delchi
!       VESOVIC Form
!        chi=ppcrit/(ttcrit*rorocrit**2)*t*capa    ! Vesovic Form
!        chiref=ppcrit/(ttcrit*rorocrit**2)*tref*caparef  ! Vesovic Form
!
!      NIST Model
!         modified form : t/ttcrit term is removed
           chi=ppcrit/(rorocrit**2)*ro*ro*capa
           chiref=ppcrit/(rorocrit**2)*ro*ro*caparef
           delchi=chi-chiref*tref/t
           IF(delchi.le.0.0d0) then
              critthcon=0.0d0
              GOTO 35
           ENDIF
      ! calculate eta
           eta=eta0*(delchi/gam)**(emew/egamma)
      ! calculate omega
!        pi=3.14159265359d0
           qd=1.d0/qdinv
           omega=2.d0/pi*((cp-cv)/cp*atan(eta*qd)+cv/cp*eta*qd)
      ! calculate omega0
           omega0=2.d0/pi*(1.d0-                                    &
           dexp(-1.d0/((qd*eta)**(-1.d0)+((qd*eta*rorocrit/ro)**2)/3.d0)))
      ! calculate viscosity
           CALL viscos_cupid(icell,nfluid,iones,iones,1,idum,t1,ro1, &
           state,vis1,err)
      ! calcuate conductivity enhancement at critical region
           critthcon=ro*cp*r*boltzmann*t/(6.d0*pi*vis*eta)*(omega-omega0)
!
! total thermal conductivity
!
   35      CONTINUE
           thcon(i)=(thcon0+delthcon)*1.d-3+critthcon   ! unit conversion
!  
           ih=ih+iscskp
         END DO
      ELSEIF(ifluid.eq. 4) THEN
!
!---------------------------------------------------------------------------
!
!
!  He-4 Th. Conductivity by Won-Jae Lee
!  Ref: Chemical Properties Handbook 
!      k = A + B*T + C*T**2 (for GAS)
!  http://www.knovel.com/knovel2/Toc.jsp?SpaceID=10093&BookID=49
!    Properties for only gas are modeled, since
!      Tcrit = 5.2 K
!      Application Ranges: 100 K < T(K) < 2000 K
!    units: t (K), thcon (W/m/K)                                             
!
!      DO 41 m=1,nh
!        i=ihld(2,ih)
!        ix=ixhld(2,i)
!        t=temp(i)
!		thcon0=0.05516d0+3.254d-4*t-2.2723d-8*t**2
!		thcon(i)=thcon0
!      ih=ih+iscskp
!   41 CONTINUE
!
!--------------------------------------------------------------------------
!
!  He-4 Th. Conductivity by Won-Jae Lee
!
!  Ref: B.A. Hands and V.D. Arp
!       "A correlation of thermal conductivity data for Helium"
!       Cryogenics, 21 (12): 697-703, 1981
!
!      Simplified Form - Thermal Conductivity at critical region is not modeled, 
!                        since operating range is far from critical point
!        
!         thcon = thcon0 + delthcon + critthcon
!                 thcon0:  dilute gas data - eq(3)
!                 delthcon:  excess conductivity outside critical region eq(13)
!		  critthcon: conductivity at critical region eq(7) - set to 0.0
! 
!    Application Ranges: 12 K < Temp < 1500 K
!                               Pres < 100 MPa
!    units: t (K), thcon (W/m/K)                                             
!
         DO m=1,nh
            i=ihld(2,ih)
            ix=ixhld(2,i)
            t=temp(i)
            ro=rho(ix)
!  dilute gas data
            sum=0.0d0
            do jj=1,4
               sum=sum+cd(jj)/t**jj
            enddo
            thcon0=2.7870034d-3*t**0.7034007057d0*dexp(sum)
!  excess term
            sum=0.0d0
            do jj=1,11
               if(jj.le.4) then
                  sum=sum+cex(jj)*t**((jj-1)/3.d0)*ro
               elseif(jj.le.7) then
                  sum=sum+cex(jj)*t**((jj-5)/3.d0)*ro**3
               elseif(jj.le.10) then
                  sum=sum+cex(jj)*t**((jj-8)/3.d0)*ro**2*dlog(ro/68.d0)
               else
                  sum=sum+cex(jj)/t*ro**2*dlog(ro/68.d0)
               endif
            enddo
            delthcon=sum
!  critical term (3.5 ~ 12 K) - currently set to ZERO (far out of range)
            critthcon=0.0d0
!  total thermal conductivity in W/m/K
            thcon(i)=thcon0+delthcon+critthcon
            thcon(i)=max(thcon(i),1.d-2)
            ih=ih+iscskp
         ENDDO
      ELSEIF(ifluid.eq. 5) THEN
!                                                                       
!  H2 Thermal Conductivity         
!  NOT IMPLEMENTED YET - 
!
         err=.true. 
      ELSEIF(ifluid.eq. 6) THEN
!
!---------------------------------------------------------------------------
!
!
!  O2 Thermal Conductivity         
!
!  Ref: E.W. Lemmon and R.T. Jacobsen
!       "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argoen and Air"
!       Internaltional Journal of Thermophysics, Vol. 25, No. 1, Jan, 2004
!
!         thcon = thcon0 + delthcon + critthcon
!                 thcon0:  dilute gas thermal conductivity
!                 delthcon:  residual fluid thermal conductivity
!		  critthcon: thermal conductivity in critical enhancement
! 
!    Application Ranges: 12 K < Temp < 1500 K
!                               Pres < 100 MPa
!    units: t (K), thcon (mW/m/K)                                             
!
!
         DO m=1,nh
            i=ihld(2,ih)
            ix=ixhld(2,i)
            t=temp(i)
            p=pres(i)
            ro=rho(ix)
!
!  Calculate Viscosity in the equation
!
            tau=tco2/t  ! Tc/T
            delt=ro/roco2/mo2   ! ro/roc/conversion to mol/dm3
            tstar=t/eko2  ! T/(e/K)
            tstar1=dlog(tstar)
            ohm=0.0d0
            DO ii=1,5
               ohm=ohm+bi(ii)*tstar1**(ii-1)
            ENDDO
            ohm=dexp(ohm)
            visc0=0.0266958d0*dsqrt(mo2*t)/ohm/sigo2**2   ! dilute gas viscosity (micoPa.s)
!     residual fluid viscosity
            delvis=0.0d0
            DO ii=1,5
               delvis=delvis+vnio2(ii)*tau**vtio2(ii)*delt**vdio2(ii)*     &
               dexp(-vgamio2(ii)*delt**vlio2(ii))
            ENDDO
!		
            tvisc=(visc0+delvis)   ! total viscosity in (microPa.s)
!
! Calculate Conductivity
!
! dilite gas
            thcon0=nio2(1)*visc0+nio2(2)*tau**tio2(2)+nio2(3)*tau**tio2(3)
! residual 
            delthcon=0.d0
            DO ii=4,9
               delthcon=delthcon+nio2(ii)*tau**tio2(ii)*delt**dio2(ii)*     &
               dexp(-gamio2(ii)*delt**lio2(ii))
            ENDDO
!     Critical Region Conductivity
            IF(ro.lt.1.d-6) then
               critthcon=0.0d0
               GOTO 65
            ENDIF
      ! steam table search
            p_stm(1)=t
            p_stm(3)=1.d0/ro
!
!----------------------------------------------------------------------- 
!    Find beta, capa, Cp for the calculation of Cv vs. T, ro
!
!!    For superheat or subcooled fluid over 0.1 K (arbitrary)
!        IF(abs((t-v_da(i)%satt)).ge.0.1d0) then
            call strtro_cupid(st_tbl(ndxstd),p_stm,ip,erx)
            beta=p_stm(6)
            capa=p_stm(7)
            cp=p_stm(8)
!!    For nearly saturated fluid within 0.1 K (arbitrary)
!!    Properties calculated at saturation (overcome peaking uncertainty)
!        ELSE
!          call strtx(st_tbl(ndxstd),p_stm,erx)
!		  IF(erx) then
!	        call strtro(st_tbl(ndxstd),p_stm,ip,erx)
!            beta=p_stm(6)
!            capa=p_stm(7)
!            cp=p_stm(8)
!		  ELSE
!!         call strtro(st_tbl(ndxstd),p_stm,ip,erx)
!            IF(state.eq.'liquid') then
!              beta=p_stm(17)
!              capa=p_stm(19)
!              cp=p_stm(21)
!            ELSE
!              beta=p_stm(18)
!              capa=p_stm(20)
!              cp=p_stm(22)
!            ENDIF
!		  ENDIF
!        ENDIF
!     calculate cv
            cv=cp-beta**2*t/(capa*ro)
            cv=max(cv,0.0d0)  ! prevent negative values *ARBITRARY*
! calculate delchi
! steam table search at Tref and given specific volume
            p_stm(1)=trefo2
            p_stm(3)=1.d0/ro  
!           CALL std2x4(st_tbl(ndxstd),p_stm,ip,erx)  ! replace std2o
            CALL strtro_cupid(st_tbl(ndxstd),p_stm,ip,erx)
            caparef=p_stm(7)
!  calculate delchi
!           chi=pco2/(roco2**2)*capa    
!           chiref=pco2/(roco2**2)*caparef  
!
            chi=pco2/(roco2**2)*capa
            chiref=pco2/(roco2**2)*caparef
            delchi=chi-chiref*trefo2/t
            IF(delchi.le.0.0d0) then
               critthcon=0.0d0
               GOTO 65
            ENDIF
! calculate eta
            eta=ksio2*(delchi/gamo2)**(0.63d0/1.2415d0)
! calculate omega
!           pi=3.14159265359d0
            qd=qdo2
            omega=2.d0/pi*((cp-cv)/cp*atan(eta/qd)+cv/cp*eta/qd)
! calculate omega0
            omega0=2.d0/pi*(1.d0-                                    &
            dexp(-1.d0/((eta/qd)**(-1.d0)+((eta/qd*roco2/ro)**2)/3.d0)))
      ! calculate viscosity
!        CALL viscos(nfluid,iones,iones,iones,idum,t,p,ro,dum,      &
!        state,vis,err)
      ! calcuate conductivity enhancement at critical region
            critthcon=ro*cp*r*boltzmann*t/(6.d0*pi*tvisc*eta)*(omega-omega0)
! check unit of eta???
! total thermal conductivity
!
   65       CONTINUE
            thcon(i)=(thcon0+delthcon)*1.d-3+critthcon   ! unit conversion
!           thcon(i)=thcon0+delthcon+critthcon   ! unit conversion   !!!  need to check wjl
!  
            ih=ih+iscskp
         ENDDO
!
      ELSEIF(ifluid.eq. 7) THEN
!
!---------------------------------------------------------------------------
!
!
!  N2 Thermal Conductivity         
!
!  Ref: E.W. Lemmon and R.T. Jacobsen
!       "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argoen and Air"
!       Internaltional Journal of Thermophysics, Vol. 25, No. 1, Jan, 2004
!
!         thcon = thcon0 + delthcon + critthcon
!                 thcon0:  dilute gas thermal conductivity
!                 delthcon:  residual fluid thermal conductivity
!		  critthcon: thermal conductivity in critical enhancement
! 
!    Application Ranges: 12 K < Temp < 1500 K
!                               Pres < 100 MPa
!    units: t (K), thcon (mW/m/K)                                             
!
!
         DO m=1,nh
            i=ihld(2,ih)
            ix=ixhld(2,i)
            t=temp(i)
            p=pres(i)
            ro=rho(ix)
!
!  Calculate Viscosity in the equation
!
            tau=tcn2/t  ! Tc/T
            delt=ro/rocn2/mn2   ! ro/roc/conversion to mol/dm3
            tstar=t/ekn2  ! T/(e/K)
            tstar1=dlog(tstar)
            ohm=0.0d0
            DO ii=1,5
               ohm=ohm+bi(ii)*tstar1**(ii-1)
            ENDDO
            ohm=dexp(ohm)
            visc0=0.0266958d0*dsqrt(mn2*t)/ohm/sign2**2   ! dilute gas viscosity (micoPa.s)
!     residual fluid viscosity
            delvis=0.0d0
            DO ii=1,5
               delvis=delvis+vnin2(ii)*tau**vtin2(ii)*delt**vdin2(ii)*     &
               dexp(-vgamin2(ii)*delt**vlin2(ii))
            ENDDO
!		
            tvisc=(visc0+delvis)   ! total viscosity in (microPa.s)
!
! Calculate Conductivity
!
! dilite gas
            thcon0=nin2(1)*visc0+nin2(2)*tau**tin2(2)+nin2(3)*tau**tin2(3)
! residual 
            delthcon=0.d0
            DO ii=4,9
               delthcon=delthcon+nin2(ii)*tau**tin2(ii)*delt**din2(ii)*     &
               dexp(-gamin2(ii)*delt**lin2(ii))
            ENDDO
!     Critical Region Conductivity
            IF(ro.lt.1.d-6) then
               critthcon=0.0d0
               GOTO 75
            ENDIF
      ! steam table search
            p_stm(1)=t
            p_stm(3)=1.d0/ro
!
!----------------------------------------------------------------------- 
!    Find beta, capa, Cp for the calculation of Cv vs. T, ro
!
!!    For superheat or subcooled fluid over 0.1 K (arbitrary)
!        IF(abs((t-v_da(i)%satt)).ge.0.1d0) then
            call strtro_cupid(st_tbl(ndxstd),p_stm,ip,erx)
            beta=p_stm(6)
            capa=p_stm(7)
            cp=p_stm(8)
!!    For nearly saturated fluid within 0.1 K (arbitrary)
!!    Properties calculated at saturation (overcome peaking uncertainty)
!        ELSE
!          call strtx(st_tbl(ndxstd),p_stm,erx)
!		  IF(erx) then
!	        call strtro(st_tbl(ndxstd),p_stm,ip,erx)
!            beta=p_stm(6)
!            capa=p_stm(7)
!            cp=p_stm(8)
!		  ELSE
!!         call strtro(st_tbl(ndxstd),p_stm,ip,erx)
!            IF(state.eq.'liquid') then
!              beta=p_stm(17)
!              capa=p_stm(19)
!              cp=p_stm(21)
!            ELSE
!              beta=p_stm(18)
!              capa=p_stm(20)
!              cp=p_stm(22)
!            ENDIF
!		  ENDIF
!        ENDIF
!
!     calculate cv
            cv=cp-beta**2*t/(capa*ro)
            cv=max(cv,0.0d0)  ! prevent negative values *ARBITRARY*
! calculate delchi
! steam table search at Tref and given specific volume
            p_stm(1)=trefn2
            p_stm(3)=1.d0/ro  
!        CALL std2x4(st_tbl(ndxstd),p_stm,ip,erx)  ! replace std2o
            CALL strtro_cupid(st_tbl(ndxstd),p_stm,ip,erx)
            caparef=p_stm(7)
!  calculate delchi
!           chi=pcn2/(rocn2**2)*capa    
!           chiref=pcn2/(rocn2**2)*caparef  
!
            chi=pcn2/(rocn2**2)*capa
            chiref=pcn2/(rocn2**2)*caparef
            delchi=chi-chiref*trefn2/t
            IF(delchi.le.0.0d0) then
               critthcon=0.0d0
               GOTO 75
            ENDIF
! calculate eta
            eta=ksin2*(delchi/gamn2)**(0.63d0/1.2415d0)
! calculate omega
!           pi=3.14159265359d0
            qd=qdn2
            omega=2.d0/pi*((cp-cv)/cp*atan(eta/qd)+cv/cp*eta/qd)
! calculate omega0
            omega0=2.d0/pi*(1.d0-                                    &
            dexp(-1.d0/((eta/qd)**(-1.d0)+((eta/qd*rocn2/ro)**2)/3.d0)))
      ! calculate viscosity
!        CALL viscos(nfluid,iones,iones,iones,idum,t,p,ro,dum,      &
!        state,vis,err)
      ! calcuate conductivity enhancement at critical region
             critthcon=ro*cp*r*boltzmann*t/(6.d0*pi*tvisc*eta)*(omega-omega0)
! check unit of eta???
! total thermal conductivity
!
   75       CONTINUE
            thcon(i)=(thcon0+delthcon)*1.d-3+critthcon   ! unit conversion
!           thcon(i)=thcon0+delthcon+critthcon   ! unit conversion   !!!  need to check wjl
!  
            ih=ih+iscskp
         ENDDO
!
      ELSEIF(ifluid.eq. 8) THEN
!
!---------------------------------------------------------------------------
!  D. LMR-K.S. Ha for liquid metal properties - Na
!  Na Thermal Conductivity         
!  "Thermophysical properties of Sodium" by G.H.Golden and J.V.Tokar
!  Conversion factor : 1.73073467d0 
      IF(state.eq.'liquid')then 
         DO m=1,nh 
            i=ihld(2,ih) 
            term=max((temp(i)-273.15d0)*1.8d0+32.d0,100.d0) 
            thcon(i)=(54.306d0-(1.878d-2-2.0914d-6*term)*term)*1.73073467d0
            ih=ih+iscskp 
         ENDDO 
      ELSE 
         DO m=1,nh 
            i=ihld(2,ih) 
            ix=ixhld(2,i) 
            term=max((temp(ix)-273.15d0)*1.8d0+32.d0,100.d0) 
            thcon(i)=(0.1639d-2+(0.3977d-4-0.9697d-8*term)*term)*1.73073467d0
            ih=ih+iscskp 
         ENDDO 
      ENDIF 
!
      ELSEIF(ifluid.eq. 9) THEN
!---------------------------------------------------------------------------
!
!  Ar Thermal Conductivity         
!
!  Ref: E.W. Lemmon and R.T. Jacobsen
!       "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argoen and Air"
!       Internaltional Journal of Thermophysics, Vol. 25, No. 1, Jan, 2004
!
!         thcon = thcon0 + delthcon + critthcon
!                 thcon0:  dilute gas thermal conductivity
!                 delthcon:  residual fluid thermal conductivity
!		  critthcon: thermal conductivity in critical enhancement
! 
!    Application Ranges: 12 K < Temp < 1500 K
!                               Pres < 100 MPa
!    units: t (K), thcon (mW/m/K)                                             
!
!
         DO m=1,nh
            i=ihld(2,ih)
            ix=ixhld(2,i)
            t=temp(i)
            p=pres(i)
            ro=rho(ix)
!
!  Calculate Viscosity in the equation
!
            tau=tcar/t  ! Tc/T
            delt=ro/rocar/mar   ! ro/roc/conversion to mol/dm3
            tstar=t/ekar  ! T/(e/K)
            tstar1=dlog(tstar)
            ohm=0.0d0
            DO ii=1,5
               ohm=ohm+bi(ii)*tstar1**(ii-1)
            ENDDO
            ohm=dexp(ohm)
            visc0=0.0266958d0*dsqrt(mar*t)/ohm/sigar**2   ! dilute gas viscosity (micoPa.s)
!     residual fluid viscosity
            delvis=0.0d0
            DO ii=1,5
               delvis=delvis+vniar(ii)*tau**vtiar(ii)*delt**vdiar(ii)* &
               dexp(-vgamiar(ii)*delt**vliar(ii))
            ENDDO
!		
            tvisc=(visc0+delvis)   ! total viscosity in (microPa.s)
!
! Calculate Conductivity
!
! dilite gas
            thcon0=niar(1)*visc0+niar(2)*tau**tiar(2)+niar(3)*tau**tiar(3)
! residual 
            delthcon=0.d0
            DO ii=4,9
               delthcon=delthcon+niar(ii)*tau**tiar(ii)*delt**diar(ii)*     &
               dexp(-gamiar(ii)*delt**liar(ii))
            ENDDO
!     Critical Region Conductivity
            IF(ro.lt.1.d-6) then
               critthcon=0.0d0
               GOTO 95
            ENDIF
      ! steam table search
            p_stm(1)=t
            p_stm(3)=1.d0/ro
!
!----------------------------------------------------------------------- 
!    Find beta, capa, Cp for the calculation of Cv vs. T, ro
!
!!    For superheat or subcooled fluid over 0.1 K (arbitrary)
!        IF(abs((t-v_da(i)%satt)).ge.0.1d0) then
            call strtro_cupid(st_tbl(ndxstd),p_stm,ip,erx)
            beta=p_stm(6)
            capa=p_stm(7)
            cp=p_stm(8)
!!    For nearly saturated fluid within 0.1 K (arbitrary)
!!    Properties calculated at saturation (overcome peaking uncertainty)
!        ELSE
!          call strtx(st_tbl(ndxstd),p_stm,erx)
!		  IF(erx) then
!	        call strtro(st_tbl(ndxstd),p_stm,ip,erx)
!            beta=p_stm(6)
!            capa=p_stm(7)
!            cp=p_stm(8)
!		  ELSE
!!         call strtro(st_tbl(ndxstd),p_stm,ip,erx)
!            IF(state.eq.'liquid') then
!              beta=p_stm(17)
!              capa=p_stm(19)
!              cp=p_stm(21)
!            ELSE
!              beta=p_stm(18)
!              capa=p_stm(20)
!              cp=p_stm(22)
!            ENDIF
!		  ENDIF
!        ENDIF
!     calculate cv
            cv=cp-beta**2*t/(capa*ro)
            cv=max(cv,0.0d0)  ! prevent negative values *ARBITRARY*
! calculate delchi
! steam table search at Tref and given specific volume
            p_stm(1)=trefar
            p_stm(3)=1.d0/ro  
!           CALL std2x4(st_tbl(ndxstd),p_stm,ip,erx)  ! replace std2o
            CALL strtro_cupid(st_tbl(ndxstd),p_stm,ip,erx)
            caparef=p_stm(7)
!  calculate delchi
!           chi=pcar/(rocar**2)*capa    
!           chiref=pcar/(rocar**2)*caparef  
!
            chi=pcar/(rocar**2)*capa
            chiref=pcar/(rocar**2)*caparef
            delchi=chi-chiref*trefar/t
            IF(delchi.le.0.0d0) then
               critthcon=0.0d0
               GOTO 95
            ENDIF
! calculate eta
           eta=ksiar*(delchi/gamar)**(0.63d0/1.2415d0)
! calculate omega
!           pi=3.14159265359d0
            qd=qdar
            omega=2.d0/pi*((cp-cv)/cp*atan(eta/qd)+cv/cp*eta/qd)
! calculate omega0
            omega0=2.d0/pi*(1.d0-                                    &
            dexp(-1.d0/((eta/qd)**(-1.d0)+((eta/qd*rocar/ro)**2)/3.d0)))
      ! calculate viscosity
!       CALL viscos(nfluid,iones,iones,iones,idum,t,p,ro,dum,      &
!       state,vis,err)
      ! calcuate conductivity enhancement at critical region
            critthcon=ro*cp*r*boltzmann*t/(6.d0*pi*tvisc*eta)*(omega-omega0)
! check unit of eta???
! total thermal conductivity
!
   95       CONTINUE
            thcon(i)=(thcon0+delthcon)*1.d-3+critthcon   ! unit conversion
!           thcon(i)=thcon0+delthcon+critthcon   ! unit conversion   !!!  need to check wjl
!  
            ih=ih+iscskp
         ENDDO
!
      ELSEIF(ifluid.eq.10) THEN
!
!---------------------------------------------------------------------------
!
!
!  Air Thermal Conductivity         
!
!  Ref: E.W. Lemmon and R.T. Jacobsen
!       "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argoen and Air"
!       Internaltional Journal of Thermophysics, Vol. 25, No. 1, Jan, 2004
!
!         thcon = thcon0 + delthcon + critthcon
!                 thcon0:  dilute gas thermal conductivity
!                 delthcon:  residual fluid thermal conductivity
!		  critthcon: thermal conductivity in critical enhancement
! 
!    Application Ranges: 12 K < Temp < 1500 K
!                               Pres < 100 MPa
!    units: t (K), thcon (mW/m/K)                                             
!
!
         DO m=1,nh
            i=ihld(2,ih)
            ix=ixhld(2,i)
            t=temp(i)
            p=pres(i)
            ro=rho(ix)
!
!  Calculate Viscosity in the equation
!
            tau=tcair/t  ! Tc/T
            delt=ro/rocair/mair   ! ro/roc/conversion to mol/dm3
            tstar=t/ekair  ! T/(e/K)
            tstar1=dlog(tstar)
            ohm=0.0d0
            DO ii=1,5
               ohm=ohm+bi(ii)*tstar1**(ii-1)
            ENDDO
            ohm=dexp(ohm)
            visc0=0.0266958d0*dsqrt(mair*t)/ohm/sigair**2   ! dilute gas viscosity (micoPa.s)
!     residual fluid viscosity
            delvis=0.0d0
            DO ii=1,5
               delvis=delvis+vniair(ii)*tau**vtiair(ii)*delt**vdiair(ii)*     &
               dexp(-vgamiair(ii)*delt**vliair(ii))
            ENDDO
!		
            tvisc=(visc0+delvis)   ! total viscosity in (microPa.s)
!
! Calculate Conductivity
!
! dilite gas
            thcon0=niair(1)*visc0+niair(2)*tau**tiair(2)+niair(3)*tau**tiair(3)
! residual 
            delthcon=0.d0
            DO ii=4,9
               delthcon=delthcon+niair(ii)*tau**tiair(ii)*delt**diair(ii)*     &
               dexp(-gamiair(ii)*delt**liair(ii))
            ENDDO
!     Critical Region Conductivity
            IF(ro.lt.1.d-6) then
               critthcon=0.0d0
               GOTO 105
            ENDIF
      ! steam table search
            p_stm(1)=t
            p_stm(3)=1.d0/ro
!
!----------------------------------------------------------------------- 
!    Find beta, capa, Cp for the calculation of Cv vs. T, ro
!
!!    For superheat or subcooled fluid over 0.1 K (arbitrary)
!        IF(abs((t-v_da(i)%satt)).ge.0.1d0) then
            call strtro_cupid(st_tbl(ndxstd),p_stm,ip,erx)
            beta=p_stm(6)
            capa=p_stm(7)
            cp=p_stm(8)
!!    For nearly saturated fluid within 0.1 K (arbitrary)
!!    Properties calculated at saturation (overcome peaking uncertainty)
!        ELSE
!          call strtx(st_tbl(ndxstd),p_stm,erx)
!		  IF(erx) then
!	        call strtro(st_tbl(ndxstd),p_stm,ip,erx)
!            beta=p_stm(6)
!            capa=p_stm(7)
!            cp=p_stm(8)
!		  ELSE
!!         call strtro(st_tbl(ndxstd),p_stm,ip,erx)
!            IF(state.eq.'liquid') then
!              beta=p_stm(17)
!              capa=p_stm(19)
!              cp=p_stm(21)
!            ELSE
!              beta=p_stm(18)
!              capa=p_stm(20)
!              cp=p_stm(22)
!            ENDIF
!		  ENDIF
!        ENDIF
!     calculate cv
            cv=cp-beta**2*t/(capa*ro)
            cv=max(cv,0.0d0)  ! prevent negative values *ARBITRARY*
! calculate delchi
! steam table search at Tref and given specific volume
            p_stm(1)=trefair
            p_stm(3)=1.d0/ro  
!           CALL std2x4(st_tbl(ndxstd),p_stm,ip,erx)  ! replace std2o
            CALL strtro_cupid(st_tbl(ndxstd),p_stm,ip,erx)
            caparef=p_stm(7)
!  calculate delchi
!           chi=pcair/(rocair**2)*capa    
!           chiref=pcair/(rocair**2)*caparef  
!
            chi=pcair/(rocair**2)*capa
            chiref=pcair/(rocair**2)*caparef
            delchi=chi-chiref*trefair/t
            IF(delchi.le.0.0d0) then
               critthcon=0.0d0
               GOTO 105
            ENDIF
! calculate eta
            eta=ksiair*(delchi/gamair)**(0.63d0/1.2415d0)
! calculate omega
!           pi=3.14159265359d0
            qd=qdair
            omega=2.d0/pi*((cp-cv)/cp*atan(eta/qd)+cv/cp*eta/qd)
! calculate omega0
            omega0=2.d0/pi*(1.d0-                                    &
            dexp(-1.d0/((eta/qd)**(-1.d0)+((eta/qd*rocair/ro)**2)/3.d0)))
      ! calculate viscosity
!        CALL viscos(nfluid,iones,iones,iones,idum,t,p,ro,dum,      &
!        state,vis,err)
      ! calcuate conductivity enhancement at critical region
            critthcon=ro*cp*r*boltzmann*t/(6.d0*pi*tvisc*eta)*(omega-omega0)
! check unit of eta???
! total thermal conductivity
!
  105       CONTINUE
            thcon(i)=(thcon0+delthcon)*1.d-3+critthcon   ! unit conversion
!           thcon(i)=thcon0+delthcon+critthcon   ! unit conversion   !!!  need to check wjl
!  
            ih=ih+iscskp
         END DO
!
      ELSEIF(ifluid.eq.11) THEN
!---------------------------------------------------------------------------
!  D. LMR-K.S. Ha for liquid metal properties - LBE
!  LBE Thermal Conductivity         
!  HLMC Handbook "CHAPTER II, Thermophysical and electric properties"
!  ENEA : by G. Benamati ...
!  SCK-CEN : by V. Sobolev, H. Ait Abderrahim
      IF(state.eq.'liquid')then 
         DO m=1,nh 
            i=ihld(2,ih) 
            term=max(temp(i),397.7d0) 
            thcon(i)=3.61d0+1.517d-2*term-1.741d-6*term**2
            ih=ih+iscskp 
         END DO 
!
!  LBE viscosity for vapor phase is not known, thus following 
!  correlation from the ATHENA program by K. E. Carlson has not meaning
      ELSE 
         DO m=1,nh 
            i=ihld(2,ih) 
            ix=ixhld(2,i) 
            term=max(temp(ix),397.7d0) 
            thcon(i)=7.56596736566d-3+3.72286713291d-5*term-8.97435897576d-10*term**2
            ih=ih+iscskp 
         END DO 
      ENDIF 
!
      ELSEIF(ifluid.eq.16) THEN
!
!  r410a
!      
         IF(state.eq.'liquid')then 
            DO m=1,nh 
               i=ihld(2,ih) 
               term=max(temp(i),120.0d0) 
               thcon(i)=0.213d0-6.33525d-4*term+4.73928d-7*term**2              !!!CYJ 20.01.07 [W/m-k], R^2=0.99939
               ih=ih+iscskp 
            END DO 
         ELSE 
            DO m=1,nh 
               i=ihld(2,ih) 
               term=max(temp(i),cell%ts(icell),120.0d0) 
               thcon(i)=-0.00362d0+5.69827e-7*DEXP((term-261.48373d0)/11.79588d0)+0.01186d0*DEXP((term-261.48373d0)/214.58719d0)     !!!CYJ 20.01.07 [W/m-k], R2=0.9995
               ih=ih+iscskp 
            END DO     
         ENDIF 
!
      ELSEIF(ifluid.eq.17) THEN
!
!.....r134a
!      
         IF(state.eq.'liquid')then 
            DO m=1,nh 
               i=ihld(2,ih) 
               term=max(temp(i),170.0d0) 
               thcon(i)=0.29386d0-0.00121d0*term+2.392d-6*term**2-2.42026d-9*term**3              !!!CYJ 20.01.07 [W/m-k], R^2=0.99989
               ih=ih+iscskp 
            END DO 
         ELSE 
            DO m=1,nh 
               i=ihld(2,ih) 
               term=max(temp(i),cell%ts(icell),170.0d0) 
               thcon(i)=-0.0179d0+8.56712e-7*DEXP((term-268.9396d0)/10.17695d0)+0.02896d0*DEXP((term-268.9396d0)/313.02497d0)     !!!CYJ 20.01.07 [W/m-k], R2=0.9998
               ih=ih+iscskp 
            END DO        
         ENDIF 
!
      ELSE
         err=.true. 
      ENDIF
!
      END SUBROUTINE thcond_cupid               
      

!---------------------------------------------------------------------------
!
      Subroutine cond_lw_single (tk,rhoc,convt)
!
!   Subroutine cond (tk,rhoc,convt) written by WJL on Dec. '97
!   This Subroutine calculates the LIGHT WATER THERMAL CONDUCTIVITY
!   based on ASME '93 Steam Table, Appendix 7
!   Input quantities are t(K) and rhoc(kg/m3)
!   Output quantities is convt(W/K.m)
!
      IMPLICIT NONE 
!
      REAL(8) tk,rhoc,convt
      REAL(8) A0,A1,A2,A3,B0,B1,B2,BB1,BB2,C1,C2,C3,C4,C5,C6,D1,D2,D3,D4,Rhostr,Tstar
      REAL(8) tratio,dratio,dtstar,dtstar06,q,r,s,ylamb,xlamb,dlamb,barlam,baslam
      REAL(8) t0,t1
!
!   Constants for Light Water Property Generation
!
      DATA A0, A1   / +1.02811D-2, +2.99621D-2 /
      DATA A2, A3   / +1.56146D-2, -4.22464D-3 /
      DATA B0, B1   / -3.97070D-1, +4.00302D-1 /
      DATA B2       / +1.06000D0 /
      DATA BB1, BB2 / -1.71587D-1, +2.39219D0  /
      DATA C1, C2   / +6.42857D-1, -4.11717D0  /
      DATA C3, C4   / -6.17937D0, +3.08976D-3  /
      DATA C5, C6   / +8.22994D-2, +1.00932D+1 /
      DATA D1, D2   / +7.01309D-2, +1.18520D-2 /
      DATA D3, D4   / +1.69937D-3, -1.02000D0  /
!      DATA NAME     / 'CONDV' /
      DATA RHOSTR   / +317.7D0 /
!      DATA TACOR    /  32.D0 /
!      DATA TBCOR    / 1.8D0 /
      DATA TSTAR    / +647.3D0 /
!
!
!      IF (TE .LT. TMIN .OR. TE .GT. TMAX) CALL STER(LEV5, 12, PE, TE)
!
      TRATIO = tk/TSTAR
      DRATIO = rhoc/RHOSTR
      DTSTAR = ABS(TRATIO - 1.D0) + C4
      DTSTAR06=DTSTAR**(-0.6D0)
      Q = 2.D0 + C5*DTSTAR06
      R = Q + 1.D0
      S = 1.D0/DTSTAR
      IF (tk .LT. TSTAR) S=C6*DTSTAR06
      YLAMB = C2*sqrt(TRATIO*TRATIO*TRATIO)+C3*(RHOSTR/rhoc)**5
!
!       THE TEST ON YLAMB IS TO PREVENT ERROR RETURNS FROM THE
!     THE SYSTEM EXPONENTIAL SUBROUTINE.  THE VALUE OF -690 REPRESENTS
!     LN(10E-300).  THE LIMIT MAY BE DIFFERENT ON THE USERS COMPUTER.
!
      IF (YLAMB .LT. -690.D0) YLAMB = -690.D0
      XLAMB = D4*EXP(YLAMB)
      t0=DRATIO**Q
      t1=DRATIO**1.8D0
      DLAMB = (D1*(TSTAR/tk)**10 + D2)*t1*EXP(C1*(1.D0 - DRATIO**2.8D0)) &
             +                    D3*S*t0*EXP((Q/R)*(1.D0 - DRATIO**R))  &
             + XLAMB
      BARLAM = B0 + B1*DRATIO + B2*EXP(BB1*(DRATIO+BB2)*(DRATIO+BB2))
      BASLAM = SQRT(TRATIO)*(A0 + TRATIO*(A1 + TRATIO*(A2 + TRATIO*A3)))
      convt = BASLAM + BARLAM + DLAMB
!
      end subroutine cond_lw_single
      
   
