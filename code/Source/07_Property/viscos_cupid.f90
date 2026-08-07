!
!     SUBROUTINE viscos_cupid(icell,ifluid,ihld,ixhld,nh,iscskp,temp,pres,rho,satt,state,visc,err)
      SUBROUTINE viscos_cupid(icell,ifluid,ihld,ixhld,nh,iscskp,temp,rho,state,visc,err)
!jjj  10/31/1997
!                                                                       
!   Light Water
!
!   Correlation from IAPS '75 Release
!   Reference:  ASME '93 Steam Table, Appendix 6, Appendix C
!   Written by Won-Jae Lee on Dec. '97
!      viscos  - compute dynamic viscosity for given fluid              
!                                                                       
!      Author:   J. E. Tolli, EG&G Idaho, Inc.                          
!      Date:     7/89                                                   
!                                                                       
!      Modified by Won-Jae Lee for inclusion of D2O, CO2 and He Properties                                                                
!      Ar, N2, O2, Dry Air Added by WJ Lee                                                                 
!                                                                       
!      Calling sequence:                                                
!                                                                       
!                call  viscos (ip1,ip2,ip3,ip4,ip5,rp6,rp7,rp8,rp9,cp10,
!                              rp11,lp12)                               
!                                                                       
!      Parameters:                                                      
!                                                                       
!                ip1  = ifluid = fluid number (input)                   
!                                                                       
!                ip2  = ihld   = primary list vector holding subscripts 
!                                for required arrays (input)            
!                                                                       
!                ip3  = ixhld  = secondary list vector (uses values     
!                                stored in ip2 for subscripts)          
!                                containing pointers to needed values in
!                                alternate data base (input)            
!                                                                       
!                ip4  = nh     = vector length (input)                  
!                                                                       
!                ip5  = iscskp = skip factor for ip2 and ip3 (input)    
!                                                                       
!                rp6  = temp   = array containing temperatures (input)  
!                                                                       
!                rp7  = pres   = array containing pressures (input)     
!                                                                       
!                rp8  = rho    = array containing fluid densities       
!                                (input)                                
!                                                                       
!                rp9  = satt   = array containing saturation            
!                                temperatures corresponding to rp7      
!                                (input)                                
!                                                                       
!                cp10 = state  = fluid state (input)                    
!                                'liquid' = liquid state                
!                                'vapor'  = vapor state                 
!                                                                       
!                rp11 = visc   = array containing viscosities (output)  
!                                                                       
!                lp12 = err    = error flag (output)                    
!                  
!if def,impnon,1 
      USE VOL_DATA   , ONLY: cell
!
      IMPLICIT NONE
      SAVE 
!                                                                       
!  Arguments.                                                           
!.....Input
      INTEGER :: icell,ifluid,ihld(2,*),ixhld(2,*),iscskp,nh
      CHARACTER*(*)state 
      REAL(8) :: temp(*),rho(*)
!.....Output
      LOGICAL :: err 
      REAL(8) :: visc(*) 
!                                                                       
!.....Local variables.                                                     
      REAL(8) term,tm6,tm7 
      INTEGER i,ih,ix,m 
!     REAL(8) crtr,a0g,a1g,a2g,b1g,c1g,d1g,e1g,g1g,g2g,g3g,g4g,t1,t2, &
      REAL(8) a0g,a1g,a2g,b1g,c1g,d1g,e1g,g1g,g2g,g3g,g4g,t1,t2, &
              f1g,f2g,f3g,f4g,th,tt         
!     REAL(8) a(9),b(6) 
      REAL(8) b(6) 
!                                                                       
!  Common blocks:                                                       
!      INCLUDE 'ufiles.h' 
! Local Variable defined by LWJ
      REAL(8) rhoc,tk
!     INCLUDE 'contrl.h'
!
!   Data  statements:                                                   
!   For Heavy Water Properties
!
      DATA tm6/1.0d-6/,tm7/1.0d-7/ 
!     DATA a/-7.691234564d0,-26.08023696d0,-168.1706546d0,6.423285504d1,&
!            -1.189646225d2,4.167117320d0,2.097506760d1,1.d9,6.d0/
!                                                                       
      DATA b/0.42325045d+01,0.28363668d+01,-0.32995982d+01,             &
             0.54697999d+01,-0.39023042d+01,0.11637261d+01/
!                                                                       
      DATA a0g/3.53d-8/,a1g/6.765d-11/,a2g/1.021d-14/,b1g/0.407d-7/,c1g/&
               10.4d-6/,d1g/1.858d-7/,e1g/5.9d-10/,g1g/176.d0/,g2g/-1.6d0/,g3g/  &
               .0048d0/,g4g/-.474074074d-5/t1/573.15d0/,t2/648.15d0/
      DATA f1g/-.2885d-5/,f2g/.2427d-7/,f3g/-.6789333d-10/,f4g/         &
               0.6317037037d-13/
!                                                                       
!  note:  1.544878727e-3 = 1.0/647.3                                    
!     DATA crtr/1.544878727d-3/ 
!       
!--------------------------------------------------------------------------
!  CO2 Viscosoty - by Won-Jae Lee on Oct. '02
!                           
      REAL*8 ac(5),dc(5),t,ro,visc0,delvis,tstar,tau,tong,tongg
      INTEGER j
      DATA ac /0.235156d0,-0.491266d0,5.211155d-2,5.347906d-2,          &
              -1.537102d-2/
      DATA dc /0.4071119d-2,0.7198037d-4,0.2411697d-16,0.2971072d-22,    &
              -0.1627888d-22/
!---------------------------------------------------------------------------
!  He Viscosity
      REAL*8 x,dd,bh,ch,dh,eta0,eta0a,eta0b,etae
!
!---------------------------------------------------------------------------
!  Coefficients of Collision Integral for N2, O2, Ar and Air
      REAL*8 bi(5),delt,ohm,tstar1   ! coefficients of collision integral (N2, O2, Ar, Air)
      DATA bi/0.431d0,-0.4623d0,0.08406d0,0.005341d0,-0.00331d0/
      INTEGER ii
!  
!  Nitrogen Viscosity
!     REAL*8 moln2
!     DATA moln2/28.013482377d0/
!     REAL*8 tcn2,rocn2,pcn2,mn2,ekn2,sign2,ksin2,gamn2,qdn2,trefn2
      REAL*8 tcn2,rocn2,mn2,ekn2,sign2
!     DATA tcn2,rocn2,pcn2,mn2,ekn2,sign2,ksin2,gamn2,qdn2,trefn2/126.192d0,11.1839d0,3.3958d0,28.01348d0,98.94d0,0.3656d0,0.17d0,0.055d0,0.4d0,252.384d0/
      DATA tcn2,rocn2,mn2,ekn2,sign2/126.192d0,11.1839d0,28.01348d0,98.94d0,0.3656d0/
      REAL*8 nin2(5),tin2(5),gamin2(5)
      INTEGER din2(5),lin2(5)
      DATA nin2,tin2,din2,gamin2,lin2/10.72d0,0.03989d0,0.001208d0,-7.402d0,4.62d0,0.1d0,0.25d0,3.2d0,0.9d0,0.3d0,2,10,12,2,1,0.d0,1.d0,1.d0,1.d0,1.d0,0,1,1,2,3/
!  Oxygen Viscosity
!     REAL*8 tco2,roco2,pco2,mo2,eko2,sigo2,ksio2,gamo2,qdo2,trefo2
      REAL*8 tco2,roco2,mo2,eko2,sigo2
!     DATA tco2,roco2,pco2,mo2,eko2,sigo2,ksio2,gamo2,qdo2,trefo2/154.581d0,13.63d0,5.043d0,31.9988d0,118.5d0,0.3428d0,0.24d0,0.055d0,0.51d0,309.162d0/
      DATA tco2,roco2,mo2,eko2,sigo2/154.581d0,13.63d0,31.9988d0,118.5d0,0.3428d0/
!     REAL*8 molo2
!     DATA molo2/31.9988d0/
      REAL*8 nio2(5),tio2(5),gamio2(5)
      INTEGER dio2(5),lio2(5)
      DATA nio2,tio2,dio2,gamio2,lio2/17.67d0,0.4042d0,0.0001077d0,0.3510d0,-13.67d0,0.05d0,0.0d0,2.1d0,0.0d0,0.5d0,1,5,12,8,1,0.d0,0.d0,0.d0,1.d0,1.d0,0,0,0,1,2/
!  Argon Viscosity
!     REAL*8 molar
!     DATA molar/39.948d0/
!     REAL*8 tcar,rocar,pcar,mar,ekar,sigar,ksiar,gamar,qdar,trefar
      REAL*8 tcar,rocar,mar,ekar,sigar
!     DATA tcar,rocar,pcar,mar,ekar,sigar,ksiar,gamar,qdar,trefar/150.687d0,13.40743d0,4.863d0,39.948d0,143.2d0,0.335d0,0.13d0,0.055d0,0.32d0,301.374d0/
      DATA tcar,rocar,mar,ekar,sigar/150.687d0,13.40743d0,39.948d0,143.2d0,0.335d0/
      REAL*8 niar(6),tiar(6),gamiar(6)
      INTEGER diar(6),liar(6)
      DATA niar,tiar,diar,gamiar,liar/12.19d0,13.99d0,0.005027d0,-18.93d0,-6.698d0,-3.827d0,0.42d0,0.0d0,0.95d0,0.5d0,0.9d0,0.8d0,1,2,10,5,1,2,0.d0,0.d0,0.d0,1.d0,1.d0,1.d0,0,0,0,2,4,4/ 
!  Air Viscosity
!     REAL*8 molair
!     DATA molair/28.013482377d0/
!     REAL*8 tcair,rocair,pcair,mair,ekair,sigair,ksiair,gamair,qdair,trefair
      REAL*8 tcair,rocair,mair,ekair,sigair
!     DATA tcair,rocair,pcair,mair,ekair,sigair,ksiair,gamair,qdair,trefair/132.6312d0,10.4477d0,3.78502d0,28.9586d0,103.3d0,0.360d0,0.11d0,0.055d0,0.31d0,265.262d0/
      DATA tcair,rocair,mair,ekair,sigair/132.6312d0,10.4477d0,28.9586d0,103.3d0,0.360d0/
      REAL*8 niair(5),tiair(5),gamiair(5)
      INTEGER diair(5),liair(5)
      DATA niair,tiair,diair,gamiair,liair/10.72d0,1.122d0,0.002019d0,-8.876d0,-0.02916d0,0.2d0,0.05d0,2.4d0,0.6d0,3.6d0, &
      1,4,9,1,8,0.d0,0.d0,0.d0,1.d0,1.d0,0,0,0,1,1/
!
!---------------------------------------------------------------------------
!                                                                       
!  Execution:                                                           
!                                                                       
!--initialize error flag                                                
!                                                                       
      err=.false. 
!                                                                       
!--check for valid state specifier                                      
!                                                                       
      IF(state.ne.'liquid'.and.state.ne.'vapor')then 
         err=.true. 
         RETURN
      ENDIF 
!                                                                       
!--initialize pointer index                                             
!                                                                       
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
      IF(    ifluid.eq. 1 .or. ifluid.eq.15) THEN
!
!--light water                                                          
!
!   Correlation from IAPS '75 Release
!   Reference:  ASME '93 Steam Table, Appendix 6, Appendix C
!   Modified by Won-Jae Lee on Dec. '97
!
!   Comment
!   IAPS '75 Correlation applies both liquid and vapor
!   in the ranges of 0 C < T < 800 C
!                    0 kg/m3 < rho < 1050 kg/m3
!                    0 MPa < P < 100 MPa
!
         DO  m=1,nh 
             i=ihld(2,ih) 
             ix=ixhld(2,i) 
             if(state.eq.'liquid') then
                tk=temp(i)
                tk=max(273.16d0,min(tk,cell%ts(icell)))
                rhoc=rho(ix)
             else
                tk=temp(i)
                tk=max(273.16d0,max(tk,cell%ts(icell)))
                rhoc=rho(ix)
             end if
             call viscos_lw_single(tk,rhoc,term)
             visc(i)=term
             ih=ih+iscskp 
         END DO 
      ELSEIF(ifluid.eq. 2) THEN
!
!--heavy water                                                          
!                                                                       
         IF(state.eq.'liquid')then 
!                                                                       
!--liquid dynamic viscosity;  correlation from Flowtran program         
!--(Savannah River);  FORTRAN coding by R. J. Wagner, C. S. Miller,     
!--and J. E. Tolli, EG&G Idaho, Inc.                                    
!                                                                       
            DO m=1,nh 
               i=ihld(2,ih) 
               th=643.89d0/temp(i)-1.d0 
               term=b(1)+(b(2)+(b(3)+(b(4)+(b(5)+b(6)*th)*th)*th)*th)*th 
               visc(i)=exp(term)*tm6 
               ih=ih+iscskp 
            END DO 
         ELSE 
!                                                                       
!--saturated or superheated vapor;  correlation from J. M. Sicilian     
!--and R. P. Harper, "Heavy Water Properties for the Transient Reactor  
!--Analysis Code (TRAC)", FSI-85-14-Q6-1, Appendix A, Flow Science Inc.,
!--December 1985;  FORTRAN coding by R. J. Wagner, C. S. Miller, and    
!--J. E. Tolli, EG&G Idaho, Inc.                                        
!                                                                       
!--note:  c1g = 8.04e-6 + 2.36e-6                                       
!                                                                       
            DO m=1,nh 
               i=ihld(2,ih) 
               ix=ixhld(2,i) 
               tt=max(temp(i)-273.15d0,tm6) 
               IF(temp(i).le.t1)then 
                  term=(b1g*tt+c1g)-rho(ix)*(d1g-e1g*(tt)) 
               ELSEIF(temp(i).lt.t2)then 
                  term=b1g*tt+c1g+(f1g+f2g*tt+f3g*tt**2+f4g*tt**3)*rho(ix)+   &
                  rho(ix)*(g1g+g2g*tt+g3g*tt**2+g4g*tt**3)*(a0g+a1g*rho(ix)+  &
                  a2g*rho(ix)**2)
               ELSE 
                  term=b1g*tt+c1g+rho(ix)*(a0g+a1g*rho(ix)+a2g*rho(ix)**2) 
               ENDIF 
               visc(i)=max(term,tm7) 
               ih=ih+iscskp 
            END DO 
         ENDIF 
      ELSEIF(ifluid.eq. 3) THEN
!---------------------------------------------------------------------------
!
!   CO2 Viscosity - by Won-Jae Lee on Oct. '02
!
!  reference
!    A. Fenghour et. al.
!    "The Viscosity of Carbon Dioxide"
!    J. of Physical and Chemistry Reference Data
!    Vol. 27: 31-44, 1998
!        
!         visc = visc0 + delvisc + critvisc
!                visc0:   viscosity in zero-density limit
!                delvis:  excess viscosity outside critical region
!		 critvis: viscosity at critical region ( ~ 0.0)
!
!  units: t(K), ro(kg/m3), visc*(microPa.s)                                               
!
         DO m=1,nh
            i=ihld(2,ih)
            ix=ixhld(2,i)
            t=temp(i)
            tstar=t/251.196d0
            ro=rho(ix)
!     Viscosity in Zero-density Limit
            tau=log(tstar)
            tongg=0.d0
            DO j=1,5
               tongg=tongg+ac(j)*tau**(j-1)
            END DO
            tong=exp(tongg)
            visc0=1.00697d0*sqrt(t)/tong
!     Excess Viscosoty
            delvis= dc(1)*ro+dc(2)*ro**2+dc(3)*ro**6/tstar**3 &
                   +dc(4)*ro**8+dc(5)*ro**8/tstar
!     Total Viscosity
            visc(i)=visc0+delvis
            visc(i)=visc(i)/1.d6! conversion to Pa.s
            ih=ih+iscskp
         END DO
      ELSEIF(ifluid.eq. 4) THEN
!---------------------------------------------------------------------------
!
!  He-4 Viscosity by Won-Jae Lee
!  Ref: Chemical Properties Handbook 
!      visc = A + B*T + C*T**2 (for GAS)
!  http://www.knovel.com/knovel2/Toc.jsp?SpaceID=10093&BookID=49
!    Properties for only gas are modeled, since
!      Tcrit = 5.2 K
!      Application Ranges: 150 K < T(K) < 2000 K
!    units: t (K), visc (Pa.s)                                             
!
!      DO 41 m=1,nh
!        i=ihld(2,ih)
!        ix=ixhld(2,i)
!        t=temp(i)
!        visc0=71.094d0+4.43d-1*t-5.18d-5*t**2  ! in microP
!		visc(i)=visc0/1.d7  ! conversion to Pa.s
!        ih=ih+iscskp
!   41 END DO
!
!  REPLACED WITH NEW NIST MODEL
!  Ref: Arp, V.D., McCarty, R.D. and Friend, D.G.
!  "Thermophysical Properties of Helium-4 from 0.8 to 1500 K with Pressures of 2000 MPa",
!  NIST Technical Note 1334, Bolder, CO 1998
!
         DO m=1,nh
            i=ihld(2,ih)
            ix=ixhld(2,i)
            tt=temp(i)
!
            x=5.7037825d0
            if(tt.le.300.d0) x=log(tt)
            ro=rho(ix)
            dd=ro/1000.d0
            bh=-47.5295259d0/x+87.6799309d0-42.0741589d0*x   &
               +8.33128289d0*x**2-0.589252385d0*x**3
            ch= 547.309267d0/x-904.870586d0+431.404928d0*x   &
               -81.4504854d0*x**2+5.37008433d0*x**3
            dh=-1684.39324d0/x+3331.08630d0-1632.19172d0*x   &
               +308.804413d0*x**2-20.2936367d0*x**3
!
            eta0a=exp(-0.135311743d0/x+1.00347841d0+1.20654649d0*x   &
                      -0.149564551d0*x**2+0.0125208416d0*x**3)
            etae=exp(bh*dd+ch*dd**2+dh*dd**3)
            if(tt.gt.100.d0) then
               eta0b=196.d0*tt**0.71938d0                             &
                     *exp(12.451d0/tt-295.67d0/tt**2-4.1249d0)
               if(tt.lt.110.d0) then
                  eta0=eta0a+(eta0b-eta0a)*(tt-100.d0)/10.d0
               else
                  eta0=eta0b
               endif
               visc0=eta0a*etae+eta0-eta0a
            else
               visc0=eta0a*etae
            endif
!
            visc(i)=visc0/1.d7
            ih=ih+iscskp
         END DO
      ELSEIF(ifluid.eq. 5) THEN 
!
!---------------------------------------------------------------------------
!                                                                       
!  H2 Dynamic Viscosity            
!  NOT IMPLEMENTED YET - 
!
         err=.true. 
      ELSEIF(ifluid.eq. 6) THEN 
!
!---------------------------------------------------------------------------
!
!  O2 Dynamic Viscosity         
!
!  Ref: E.W. Lemmon and R.T. Jacobsen
!       "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argoen and Air"
!       Internaltional Journal of Thermophysics, Vol. 25, No. 1, Jan, 2004
!
!         visc = visc0 + delvisc
!                visc0:   viscosity of dilute gas
!                delvis:  residual fluid viscosity 
!
!  units: t(K), ro(kg/m3), visc*(microPa.s)                                               
!
         DO m=1,nh
            i=ihld(2,ih)
            ix=ixhld(2,i)
            tt=temp(i)
            ro=rho(ix)
!
            tau=tco2/tt  ! Tc/T
            delt=ro/roco2/mo2   ! ro/roc/conversion to mol/dm3
            tstar=tt/eko2  ! T/(e/K)
            tstar1=log(tstar)
            ohm=0.d0
            DO ii=1,5
               ohm=ohm+bi(ii)*tstar1**(ii-1)
            ENDDO
            ohm=exp(ohm)
            visc0=0.0266958d0*sqrt(mo2*tt)/ohm/sigo2**2   ! dilute gas viscosity (micoPa.s)
!     residual fluid viscosity
            delvis=0.d0
            DO ii=1,5
               delvis=delvis+nio2(ii)*tau**tio2(ii)*delt**dio2(ii)* &
                      exp(-gamio2(ii)*delt**lio2(ii))
            ENDDO
!		
            visc(i)=(visc0+delvis)/1.d6
            ih=ih+iscskp
         END DO
      ELSEIF(ifluid.eq. 7) THEN 
!
!---------------------------------------------------------------------------
!
!  N2 Dynamic Viscosity        
!
!  Ref: E.W. Lemmon and R.T. Jacobsen
!       "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argoen and Air"
!       Internaltional Journal of Thermophysics, Vol. 25, No. 1, Jan, 2004
!
!         visc = visc0 + delvisc
!                visc0:   viscosity of dilute gas
!                delvis:  residual fluid viscosity 
!
!  units: t(K), ro(kg/m3), visc*(microPa.s)                                               
!
         DO m=1,nh
            i=ihld(2,ih)
            ix=ixhld(2,i)
            tt=temp(i)
            ro=rho(ix)
!
            tau=tcn2/tt  ! Tc/T
!        delt=ro/rocn2/moln2   ! ro/roc/conversion to mol/dm3
            delt=ro/rocn2/mn2   ! ro/roc/conversion to mol/dm3
            tstar=tt/ekn2  ! T/(e/K)
            tstar1=log(tstar)
            ohm=0.d0
            DO ii=1,5
               ohm=ohm+bi(ii)*tstar1**(ii-1)
            ENDDO
            ohm=exp(ohm)
            visc0=0.0266958d0*sqrt(mn2*tt)/ohm/sign2**2   ! dilute gas viscosity (micoPa.s)
!     residual fluid viscosity
            delvis=0.d0
            DO ii=1,5
               delvis=delvis+nin2(ii)*tau**tin2(ii)*delt**din2(ii)*     &
                      exp(-gamin2(ii)*delt**lin2(ii))
            ENDDO
!		
            visc(i)=(visc0+delvis)/1.d6
            ih=ih+iscskp
         END DO
      ELSEIF(ifluid.eq. 8) THEN 
!
!---------------------------------------------------------------------------
!>>
!  D. LMR-K.S. Ha for liquid metal properties - Na
!--sodium 
!  "Thermophysical properties of Sodium" by G.H.Golden and J.V.Tokar
!  Conversion factor : 1.73073467d0 
!
         IF(state.eq.'liquid')then 
!
!  Shpilrain et al. equation in data book of Golden and Tokar
            DO m=1,nh 
               i=ihld(2,ih) 
               th=temp(i)
               term=-1.6814d0-0.4296d0*log10(th)+234.65/th
               visc(i)=10.d0**term*1.d-3 
               ih=ih+iscskp 
            END DO 
         ELSE 
!  term:Fahrenheit, Conversion factor:4.13379d-4
            DO m=1,nh 
               i=ihld(2,ih) 
               term=max((temp(i)-273.15d0)*1.8d0+32.d0,100.d0) 
               visc(i)=(0.019d0+(0.1375d-4+0.1709d-9*term)*term)*4.13379d-4
               ih=ih+iscskp 
            END DO 
         ENDIF 
!
      ELSEIF(ifluid.eq. 9) THEN 
!
!---------------------------------------------------------------------------
!
!  Argon Dynamic Viscosity         
!
!  Ref: E.W. Lemmon and R.T. Jacobsen
!       "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argoen and Air"
!       Internaltional Journal of Thermophysics, Vol. 25, No. 1, Jan, 2004
!
!         visc = visc0 + delvisc
!                visc0:   viscosity of dilute gas
!                delvis:  residual fluid viscosity 
!
!  units: t(K), ro(kg/m3), visc*(microPa.s)                                               
!
         DO m=1,nh
            i=ihld(2,ih)
            ix=ixhld(2,i)
            tt=temp(i)
            ro=rho(ix)
!
            tau=tcar/tt  ! Tc/T
            delt=ro/rocar/mar   ! ro/roc/conversion to mol/dm3
            tstar=tt/ekar  ! T/(e/K)
            tstar1=log(tstar)
            ohm=0.d0
            DO ii=1,5
               ohm=ohm+bi(ii)*tstar1**(ii-1)
            ENDDO
            ohm=exp(ohm)
            visc0=0.0266958d0*sqrt(mar*tt)/ohm/sigar**2   ! dilute gas viscosity (micoPa.s)
!     residual fluid viscosity
            delvis=0.d0
            DO ii=1,5
               delvis=delvis+niar(ii)*tau**tiar(ii)*delt**diar(ii)* &
                      exp(-gamiar(ii)*delt**liar(ii))
            ENDDO
!		
            visc(i)=(visc0+delvis)/1.d6
            ih=ih+iscskp
         END DO
      ELSEIF(ifluid.eq.10) THEN 
!
!---------------------------------------------------------------------------
!
!  Air Dynamic Viscosity        
!
!  Ref: E.W. Lemmon and R.T. Jacobsen
!       "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argoen and Air"
!       Internaltional Journal of Thermophysics, Vol. 25, No. 1, Jan, 2004
!
!         visc = visc0 + delvisc
!                visc0:   viscosity of dilute gas
!                delvis:  residual fluid viscosity 
!
!  units: t(K), ro(kg/m3), visc*(microPa.s)                                               
!
         DO m=1,nh
            i=ihld(2,ih)
            ix=ixhld(2,i)
            tt=temp(i)
            ro=rho(ix)
!
            tau=tcair/tt  ! Tc/T
            delt=ro/rocair/mair   ! ro/roc/conversion to mol/dm3
            tstar=tt/ekair  ! T/(e/K)
            tstar1=log(tstar)
            ohm=0.d0
            DO ii=1,5
               ohm=ohm+bi(ii)*tstar1**(ii-1)
            ENDDO
            ohm=exp(ohm)
            visc0=0.0266958d0*sqrt(mair*tt)/ohm/sigair**2   ! dilute gas viscosity (micoPa.s)
!     residual fluid viscosity
            delvis=0.d0
            DO ii=1,5
               delvis=delvis+niair(ii)*tau**tiair(ii)*delt**diair(ii)* &
                      exp(-gamiair(ii)*delt**liair(ii))
            ENDDO
!		
            visc(i)=(visc0+delvis)/1.d6
            ih=ih+iscskp
         END DO
      ELSEIF(ifluid.eq.11) THEN 
!
!---------------------------------------------------------------------------
!>>
!  D. LMR-K.S. Ha for liquid metal properties - LBE
!--Lead Bismuth Eutetic
!  2nd meeting NEA Nuclear Science Committee ....
!  "Thermophysical properties of LBE Alloy for use in reactor safety analysis" 
!  by Koji Morita et al.
!
         IF(state.eq.'liquid')then 
            DO m=1,nh 
               i=ihld(2,ih) 
               th=temp(i)
               visc(i)=0.490d-3*exp(760.1d0/th)
               ih=ih+iscskp 
            END DO 
         ELSE 
!
!  LBE viscosity for vapor phase is not known, thus following 
!  correlation from the ATHENA program by K. E. Carlson has not meaning
            DO m=1,nh 
               i=ihld(2,ih) 
               th=temp(i)
               term=log(th)
               visc(i)=2.16867957522d-4-6.72194914063d-5*th+5.39196834628d-6*th**2
               ih=ih+iscskp 
            END DO 
         ENDIF 
!
      ELSEIF(ifluid.eq.16) THEN 
!
!.....r12
!    
         IF(state.eq.'liquid')then 
            DO m=1,nh 
               i=ihld(2,ih) 
               th=max(temp(i),120.0d0) 
               visc(i)=1.26494d-4+0.00605d0*EXP(-(th-124.51017d0)/6.13061d0)+0.00331d0*EXP(-(th-124.51017d0)/39.95463d0)     !!!CYJ 20.01.07 [Pa-s], R2=0.99898
               ih=ih+iscskp 
            END DO 
         ELSE 
            DO m=1,nh 
               i=ihld(2,ih) 
               th=max(temp(i),120.0d0) 
               visc(i)=-6.26376d-5+5.50892d-1*EXP((th-50.10553d0)/15.66878d0)+6.48982d-5*EXP((th-50.10553d0)/1833.67125d0)     !!!CYJ 20.01.07 [Pa-s], R2=0.99965
               ih=ih+iscskp 
            END DO 
         ENDIF 
      ELSEIF(ifluid.eq.17) THEN 
!
!.....r134a
!     
         IF(state.eq.'liquid')then 
            DO m=1,nh 
               i=ihld(2,ih) 
               th=max(temp(i),170.0d0) 
               visc(i)=-8.19201d-6+0.00116d0*EXP(-(th-169.30137d0)/15.87371d0)+0.00104d0*EXP(-(th-169.30137d0)/78.06749d0)     !!!CYJ 20.01.07 [Pa-s], R2=0.99987
               ih=ih+iscskp 
            END DO 
         ELSE 
            DO m=1,nh 
               i=ihld(2,ih) 
               th=max(temp(i),170.0d0) 
               visc(i)=-3.72684d-5+7.26365d-1*EXP((th-134.2855d0)/12.87191d0)+4.29147d-5*EXP((th-134.2855d0)/1234.7172d0)     !!!CYJ 20.01.07 [Pa-s], R2=0.99954
               ih=ih+iscskp 
            END DO 
         ENDIF 
!
      ELSE
         err=.true. 
      ENDIF
!
      END SUBROUTINE viscos_cupid   
!
!---------------------------------------------------------------------------
!
      Subroutine viscos_lw_single (tk,rhoc,term)
!
!   Subroutine cond (tk,rhoc,convt) written by WJL on Dec. '97
!   This Subroutine calculates the LIGHT WATER THERMAL CONDUCTIVITY
!   based on ASME '93 Steam Table, Appendix 7
!   Input quantities are t(K) and rhoc(kg/m3)
!   Output quantities is term(Pa.s)

      IMPLICIT NONE
!.....Input
      REAL(8) :: tk,rhoc
!.....Output
      REAL(8) :: term
!.....Local variables
      INTEGER :: k,ii,jj
      REAL(8) :: ts1,dt1,r1,dr,sum,temp_,dt1i
!
!   Data Statements:
!   For Light Water Properties
!
      REAL(8) a(4),bb(6,5) 
      DATA a/0.0181583d0,0.0177624d0,0.0105287d0,-0.0036744d0/ 
      DATA bb/0.501938d0,0.162888d0,-0.130356d0,0.907919d0,-0.551119d0,      &
              0.146543d0,0.235622d0,0.789393d0,0.673665d0,1.207552d0,        &
              0.0670665d0,-0.0843370d0,-0.274637d0,-0.743539d0,-0.959456d0,- &
              0.687343d0,-0.497089d0,0.195286d0,0.145831d0,0.263129d0,       &
              0.347247d0,0.213486d0,0.100754d0,-0.032932d0,-0.0270448d0,-    &
              0.0253093d0,-0.0267758d0,-0.0822904d0,0.0602253d0,-0.0202595d0/   
!
!   Test ranges of validity
!
      !if ((tk.ge.1073.15d0).or.(rhoc.ge.1050.d0).or.(ptest.ge.1.d8))  &
      !   print *, tk, 'Ranges Invalid - Subroutine Viscos'
!
        ts1=647.27d0/tk
        dt1=ts1-1.d0
        r1=rhoc/317.763d0
        dr=r1-1.d0
!
        sum=0.d0 
        temp_=1.d0
        DO  k=1,4 
          sum=sum+a(k)*temp_
          temp_=temp_*ts1
        ENDDO 
        term=1.d-06*sqrt(1.d0/ts1)/sum
!
        if(tk.eq.647.27d0)  return
!
        sum=0.d0
        dt1i=1.d0
        DO ii=1,6 
           temp_=1.d0
           DO jj=1,5 
              sum=sum+bb(ii,jj)*dt1i*temp_
              temp_=temp_*dr
           ENDDO
           dt1i=dt1i*dt1 
        ENDDO
! 
        term=term*exp(r1*sum) 
!
        END SUBROUTINE viscos_lw_single
!!
!      SUBROUTINE viscos_cupid(ifluid,ihld,ixhld,nh,iscskp,temp,pres,rho,satt, &
!                        state,visc,err)
!!jjj  10/31/1997
!!                                                                       
!!   Light Water
!!
!!   Correlation from IAPS '75 Release
!!   Reference:  ASME '93 Steam Table, Appendix 6, Appendix C
!!   Written by Won-Jae Lee on Dec. '97
!!      viscos  - compute dynamic viscosity for given fluid              
!!                                                                       
!!      Author:   J. E. Tolli, EG&G Idaho, Inc.                          
!!      Date:     7/89                                                   
!!                                                                       
!!      Modified by Won-Jae Lee for inclusion of D2O, CO2 and He Properties                                                                
!!      Ar, N2, O2, Dry Air Added by WJ Lee                                                                 
!!                                                                       
!!      Calling sequence:                                                
!!                                                                       
!!                call  viscos (ip1,ip2,ip3,ip4,ip5,rp6,rp7,rp8,rp9,cp10,
!!                              rp11,lp12)                               
!!                                                                       
!!      Parameters:                                                      
!!                                                                       
!!                ip1  = ifluid = fluid number (input)                   
!!                                                                       
!!                ip2  = ihld   = primary list vector holding subscripts 
!!                                for required arrays (input)            
!!                                                                       
!!                ip3  = ixhld  = secondary list vector (uses values     
!!                                stored in ip2 for subscripts)          
!!                                containing pointers to needed values in
!!                                alternate data base (input)            
!!                                                                       
!!                ip4  = nh     = vector length (input)                  
!!                                                                       
!!                ip5  = iscskp = skip factor for ip2 and ip3 (input)    
!!                                                                       
!!                rp6  = temp   = array containing temperatures (input)  
!!                                                                       
!!                rp7  = pres   = array containing pressures (input)     
!!                                                                       
!!                rp8  = rho    = array containing fluid densities       
!!                                (input)                                
!!                                                                       
!!                rp9  = satt   = array containing saturation            
!!                                temperatures corresponding to rp7      
!!                                (input)                                
!!                                                                       
!!                cp10 = state  = fluid state (input)                    
!!                                'liquid' = liquid state                
!!                                'vapor'  = vapor state                 
!!                                                                       
!!                rp11 = visc   = array containing viscosities (output)  
!!                                                                       
!!                lp12 = err    = error flag (output)                    
!!                  
!!if def,impnon,1 
!      IMPLICIT none 
!      SAVE 
!!                                                                       
!!  Arguments.                                                           
!      REAL(8) pres(*),rho(*),satt(*),temp(*),visc(*) 
!      INTEGER ifluid,ihld(2,*),iscskp,ixhld(2,*),nh 
!      LOGICAL err 
!      CHARACTER*(*)state 
!!                                                                       
!!  Local variables.                                                     
!      REAL(8) term,tm6,tm7 
!      INTEGER i,ih,ix,m 
!      REAL(8) crtr,a0g,a1g,a2g,b1g,c1g,d1g,e1g,g1g,g2g,g3g,g4g,t1,t2,   &
!      f1g,f2g,f3g,f4g,psat,templ,fr,fr1,tflim,tempx,tglim,th,tt         
!      REAL(8) a(9),b(6) 
!!                                                                       
!!  Common blocks:                                                       
!!      INCLUDE 'ufiles.h' 
!! Local Variable defined by LWJ
!      REAL(8) rhoc,tk
!!      INCLUDE 'contrl.h'
!!
!!   Data  statements:                                                   
!!   For Heavy Water Properties
!!
!      DATA tm6/1.0d-6/,tm7/1.0d-7/ 
!      DATA a/-7.691234564d0,-26.08023696d0,-168.1706546d0,6.423285504d1,&
!      -1.189646225d2,4.167117320d0,2.097506760d1,1.d9,6.0d0/            
!!                                                                       
!      DATA b/0.42325045d+01,0.28363668d+01,-0.32995982d+01,             &
!      0.54697999d+01,-0.39023042d+01,0.11637261d+01/                    
!!                                                                       
!      DATA a0g/3.53d-8/,a1g/6.765d-11/,a2g/1.021d-14/,b1g/0.407d-7/,c1g/&
!      10.4d-6/,d1g/1.858d-7/,e1g/5.9d-10/,g1g/176.d0/,g2g/-1.6d0/,g3g/  &
!      .0048d0/,g4g/-.474074074d-5/t1/573.15d0/,t2/648.15d0/             
!      DATA f1g/-.2885d-5/,f2g/.2427d-7/,f3g/-.6789333d-10/,f4g/         &
!      0.6317037037d-13/                                                 
!!                                                                       
!!  note:  1.544878727e-3 = 1.0/647.3                                    
!      DATA crtr/1.544878727d-3/ 
!!       
!!--------------------------------------------------------------------------
!!  CO2 Viscosoty - by Won-Jae Lee on Oct. '02
!!                           
!      REAL*8 ac(5),dc(5),t,ro,visc0,delvis,tstar,tau,tong,wm,tongg
!  INTEGER j
!      DATA ac /0.235156d0,-0.491266d0,5.211155d-2,5.347906d-2,          &
!              -1.537102d-2/
!      DATA dc /0.4071119d-2,0.7198037d-4,0.2411697d-16,0.2971072d-22,    &
!              -0.1627888d-22/
!!---------------------------------------------------------------------------
!!  He Viscosity
!      REAL*8 x,dd,bh,ch,dh,eta0,eta0a,eta0b,etae
!!
!!---------------------------------------------------------------------------
!!  Coefficients of Collision Integral for N2, O2, Ar and Air
!      REAL*8 bi(5),delt,ohm,tstar1   ! coefficients of collision integral (N2, O2, Ar, Air)
!  DATA bi/0.431d0,-0.4623d0,0.08406d0,0.005341d0,-0.00331d0/
!  INTEGER ii
!!  
!!  Nitrogen Viscosity
!  REAL*8 moln2
!  DATA moln2/28.013482377d0/
!      REAL*8 tcn2,rocn2,pcn2,mn2,ekn2,sign2,ksin2,gamn2,qdn2,trefn2
!  DATA tcn2,rocn2,pcn2,mn2,ekn2,sign2,ksin2,gamn2,qdn2,trefn2/       &
!  126.192d0,11.1839d0,3.3958d0,28.01348d0,98.94d0,0.3656d0,0.17d0,   &
!  0.055d0,0.4d0,252.384d0/
!  REAL*8 nin2(5),tin2(5),gamin2(5)
!  INTEGER din2(5),lin2(5)
!  DATA nin2,tin2,din2,gamin2,lin2/                    &
!  10.72d0,0.03989d0,0.001208d0,-7.402d0,4.62d0,                      &
!  0.1d0,0.25d0,3.2d0,0.9d0,0.3d0,                                    &
!  2,10,12,2,1,                                                       &
!  0.d0,1.d0,1.d0,1.d0,1.d0,                                          &
!  0,1,1,2,3/
!!  Oxygen Viscosity
!      REAL*8 tco2,roco2,pco2,mo2,eko2,sigo2,ksio2,gamo2,qdo2,trefo2
!  DATA tco2,roco2,pco2,mo2,eko2,sigo2,ksio2,gamo2,qdo2,trefo2/       &
!  154.581d0,13.63d0,5.043d0,31.9988d0,118.5d0,0.3428d0,0.24d0,       &
!  0.055d0,0.51d0,309.162d0/
!  REAL*8 molo2
!  DATA molo2/31.9988d0/
!  REAL*8 nio2(5),tio2(5),gamio2(5)
!  INTEGER dio2(5),lio2(5)
!  DATA nio2,tio2,dio2,gamio2,lio2/                    &
!  17.67d0,0.4042d0,0.0001077d0,0.3510d0,-13.67d0,                    &
!  0.05d0,0.0d0,2.1d0,0.0d0,0.5d0,                                    &
!  1,5,12,8,1,                                                        &
!  0.d0,0.d0,0.d0,1.d0,1.d0,                                          &
!  0,0,0,1,2/
!!  Argon Viscosity
!  REAL*8 molar
!  DATA molar/39.948d0/
!      REAL*8 tcar,rocar,pcar,mar,ekar,sigar,ksiar,gamar,qdar,trefar
!  DATA tcar,rocar,pcar,mar,ekar,sigar,ksiar,gamar,qdar,trefar/       &
!  150.687d0,13.40743d0,4.863d0,39.948d0,143.2d0,0.335d0,0.13d0,      &
!  0.055d0,0.32d0,301.374d0/
!  REAL*8 niar(6),tiar(6),gamiar(6)
!  INTEGER diar(6),liar(6)
!  DATA niar,tiar,diar,gamiar,liar/                    &
!  12.19d0,13.99d0,0.005027d0,-18.93d0,-6.698d0,-3.827d0,             &
!  0.42d0,0.0d0,0.95d0,0.5d0,0.9d0,0.8d0,                             &
!  1,2,10,5,1,2,                                                      &
!  0.d0,0.d0,0.d0,1.d0,1.d0,1.d0,                                     &
!  0,0,0,2,4,4/ 
!!  Air Viscosity
!  REAL*8 molair
!  DATA molair/28.013482377d0/
!      REAL*8 tcair,rocair,pcair,mair,ekair,sigair,ksiair,gamair,qdair    &
!         ,trefair
!  DATA tcair,rocair,pcair,mair,ekair,sigair,ksiair,gamair,qdair      &
!       ,trefair/       &
!  132.6312d0,10.4477d0,3.78502d0,28.9586d0,103.3d0,0.360d0,0.11d0,   &
!  0.055d0,0.31d0,265.262d0/
!  REAL*8 niair(5),tiair(5),gamiair(5)
!  INTEGER diair(5),liair(5)
!  DATA niair,tiair,diair,gamiair,liair/               &
!  10.72d0,1.122d0,0.002019d0,-8.876d0,-0.02916d0,                    &
!  0.2d0,0.05d0,2.4d0,0.6d0,3.6d0,                                    &
!  1,4,9,1,8,                                                         &
!  0.d0,0.d0,0.d0,1.d0,1.d0,                                          &
!  0,0,0,1,1/
!!
!!---------------------------------------------------------------------------
!!                                                                       
!!  Execution:                                                           
!!                                                                       
!!--initialize error flag                                                
!!                                                                       
!      err=.false. 
!!                                                                       
!!--check for valid state specifier                                      
!!                                                                       
!      IF(state.ne.'liquid'.and.state.ne.'vapor')then 
!         err=.true. 
!         GOTO 999 
!      ENDIF 
!!                                                                       
!!--initialize pointer index                                             
!!                                                                       
!      ih=1 
!!                                                                       
!!--branch to correlation for given fluid    
!!    Light Water: 10
!!    Heavy Water: 20
!!    CO2        : 30
!!    He         : 40
!!    H2         : 50
!!    O2         : 60
!!    N2         : 70
!!    Na         : 80
!!    Ar         : 90
!!    Air        : 100
!!    LBE        : 110
!!                            
!      GOTO(10,20,30,40,50,60,70,80,90,100,110),ifluid 
!      IF(ifluid.eq.15) GOTO 10 
!      err=.true. 
!      GOTO 999 
!!
!!--light water                                                          
!!
!!   Correlation from IAPS '75 Release
!!   Reference:  ASME '93 Steam Table, Appendix 6, Appendix C
!!   Modified by Won-Jae Lee on Dec. '97
!!
!   10 CONTINUE 
!!
!!   Comment
!!   IAPS '75 Correlation applies both liquid and vapor
!!   in the ranges of 0 C < T < 800 C
!!                    0 kg/m3 < rho < 1050 kg/m3
!!                    0 MPa < P < 100 MPa
!!
!      IF(.TRUE.) then   !Dynamic property model
!!
!!   1-D part array calculation
!!
!         IF(state.eq.'liquid')then 
!            GOTO 11 
!         ELSE 
!            GOTO 13 
!         ENDIF 
!!
!!--saturated or subcooled liquid;  correlation from                     
!!--The American Society Of Mechanical Engineers,                        
!!--Thermodynamic And Transport Properties Of Steam,                     
!!--United Engineering Center, 345 East 47-th Street,                    
!!--New York, N.Y., 10017, (1967);                                       
!!--Springer-Verlag New York Inc., 1969                                  
!!--Properties of water and steam in SI-units                            
!!--Prepared by Ernst Schmidt                                            
!!--FORTRAN coding by R. J. Wagner and J. E. Tolli, EG&G Idaho, Inc.     
!!                                                                       
!   11    CONTINUE 
!         DO 12 m=1,nh 
!            i=ihld(2,ih) 
!            ix=ixhld(2,i) 
!            psat=pres(i) 
!            templ=temp(i) 
!            templ=max(273.16d0,min(templ,satt(i))) 
!            IF(templ.lt.647.3d0)then 
!               fr=templ*crtr 
!               fr1=1.0d0-fr 
!               psat=2.212d+7*dexp((((((a(5)*fr1+a(4))*fr1+a(3))*fr1+a(2))*  &
!               fr1+a(1))*fr1)/(((a(7)*fr1+a(6))*fr1+1.0d0)*fr)-fr1/(a(8)*  &
!               fr1*fr1+a(9)))                                              
!            ENDIF 
!            tflim=max(573.15d0,min(647.3d0,templ)) 
!            visc(i)=(647.3d0-tflim)*dexp(570.58059d0/(templ-140.d0))*(      &
!            3.2555630478758d-7+3.4088115981118d-18*max((pres(i)-psat),0.d0)&
!            *(templ-305.d0))+(tflim-573.15d0)*(1.0842886041807d-7+         &
!            5.488873904248d-10*(templ-273.15d0)+rho(ix)*(                  &
!            4.7606203631266d-10+rho(ix)*(9.1233984665203d-13+rho(ix) *     &
!            1.3769386378961d-16)))                                         
!            ih=ih+iscskp 
!   12    END DO 
!         GOTO 999 
!!                                                                       
!!--saturated or superheated vapor;  correlation from                    
!!--The American Society Of Mechanical Engineers,                        
!!--Thermodynamic And Transport Properties Of Steam,                     
!!--United Engineering Center, 345 East 47-th Street,                    
!!--New York, N.Y., 10017, (1967);                                       
!!--Springer-Verlag New York Inc., 1969                                  
!!--Properties of water and steam in SI-units                            
!!--Prepared by Ernst Schmidt                                            
!!--FORTRAN coding by R. J. Wagner and J. E. Tolli, EGG Idaho, Inc.      
!!                                                                       
!   13    DO 14 m=1,nh 
!            i=ihld(2,ih) 
!            ix=ixhld(2,i) 
!!  Check for sub-critical state.                                        
!            IF(pres(i).ge.2.212d+7)then 
!               psat=pres(i) 
!               templ=temp(i) 
!               templ=max(273.16d0,min(templ,satt(i))) 
!               IF(templ.lt.647.3d0)then 
!                  fr=templ*crtr 
!                  fr1=1.0d0-fr 
!                  psat=2.212d+7*dexp((((((a(5)*fr1+a(4))*fr1+a(3))*fr1+a(2))&
!                  *fr1+a(1))*fr1)/(((a(7)*fr1+a(6))*fr1+1.0d0)*fr)-fr1/(a( &
!                  8)*fr1*fr1+a(9)))                                        
!               ENDIF 
!               tflim=max(573.15d0,min(647.3d0,templ)) 
!               visc(i)=(647.3d0-tflim)*dexp(570.58059d0/(templ-140.d0))*(   &
!               3.2555630478758d-7+3.4088115981118d-18*max((pres(i)-psat),  &
!               0.d0)*(templ-305.d0))+(tflim-573.15d0)*(1.0842886041807d-7+ &
!               5.488873904248d-10*(templ-273.15d0)+rho(ix)*(               &
!               4.7606203631266d-10+rho(ix)*(9.1233984665203d-13+rho(ix) *  &
!               1.3769386378961d-16)))                                      
!            ELSE 
!               tempx=max(temp(i),satt(i))-273.15d0 
!!  Vapor viscosity at pres = 1 bar.                                     
!               term=4.07d-8*tempx+8.04d-6 
!               IF(pres(i).ne.1.0d+5.or.tempx.lt.100.0d0)then 
!!  Saturated and superheated vapor viscosity.                           
!                  tglim=max(340.0d0,min(365.0d0,tempx)) 
!                  term=term+rho(ix)*((tglim-340.0d0)*(1.412d-9+rho(ix)*(   &
!                  2.706d-12+rho(ix) * 4.084d-16))-(365.0d0-tglim)*(        &
!                  7.432d-9-2.36d-11*tempx))                                
!               ENDIF 
!               visc(i)=term 
!            ENDIF 
!            ih=ih+iscskp 
!   14    END DO 
!         GOTO 999 
!!
!      ELSE
!         DO  m=1,nh 
!           i=ihld(2,ih) 
!           ix=ixhld(2,i) 
!           if (state.eq.'liquid') then
!              tk=temp(i)
!              tk=max(273.16d0,min(tk,satt(i)))
!              rhoc=rho(ix)
!           else
!              tk=temp(i)
!              tk=max(273.16d0,max(tk,satt(i)))
!              rhoc=rho(ix)
!           end if
!           call viscos_lw(tk,rhoc,term)
!           visc(i)=term
!           ih=ih+iscskp 
!         END DO 
!         GOTO 999 
!
!! end of LWJ 1D Update
!      ENDIF
!!
!!--heavy water                                                          
!!                                                                       
!   20 CONTINUE
!      IF(state.eq.'liquid')then 
!         GOTO 21 
!      ELSE 
!         GOTO 23 
!      ENDIF 
!!                                                                       
!!--liquid dynamic viscosity;  correlation from Flowtran program         
!!--(Savannah River);  FORTRAN coding by R. J. Wagner, C. S. Miller,     
!!--and J. E. Tolli, EG&G Idaho, Inc.                                    
!!                                                                       
!   21 DO 22 m=1,nh 
!         i=ihld(2,ih) 
!         th=643.89d0/temp(i)-1.0d0 
!         term=b(1)+(b(2)+(b(3)+(b(4)+(b(5)+b(6)*th)*th)*th)*th)*th 
!         visc(i)=dexp(term)*tm6 
!         ih=ih+iscskp 
!   22 END DO 
!      GOTO 999 
!!                                                                       
!!--saturated or superheated vapor;  correlation from J. M. Sicilian     
!!--and R. P. Harper, "Heavy Water Properties for the Transient Reactor  
!!--Analysis Code (TRAC)", FSI-85-14-Q6-1, Appendix A, Flow Science Inc.,
!!--December 1985;  FORTRAN coding by R. J. Wagner, C. S. Miller, and    
!!--J. E. Tolli, EG&G Idaho, Inc.                                        
!!                                                                       
!!--note:  c1g = 8.04e-6 + 2.36e-6                                       
!!                                                                       
!   23 DO 24 m=1,nh 
!         i=ihld(2,ih) 
!         ix=ixhld(2,i) 
!         tt=max(temp(i)-273.15d0,tm6) 
!         IF(temp(i).le.t1)then 
!            term=(b1g*tt+c1g)-rho(ix)*(d1g-e1g*(tt)) 
!         ELSEIF(temp(i).lt.t2)then 
!            term=b1g*tt+c1g+(f1g+f2g*tt+f3g*tt**2+f4g*tt**3)*rho(ix)+   &
!            rho(ix)*(g1g+g2g*tt+g3g*tt**2+g4g*tt**3)*(a0g+a1g*rho(ix)+  &
!            a2g*rho(ix)**2)                                             
!         ELSE 
!            term=b1g*tt+c1g+rho(ix)*(a0g+a1g*rho(ix)+a2g*rho(ix)**2) 
!         ENDIF 
!         visc(i)=max(term,tm7) 
!         ih=ih+iscskp 
!   24 END DO 
!      GOTO 999 
!!---------------------------------------------------------------------------
!!
!!   CO2 Viscosity - by Won-Jae Lee on Oct. '02
!!
!   30 CONTINUE
!!
!!  reference
!!    A. Fenghour et. al.
!!    "The Viscosity of Carbon Dioxide"
!!    J. of Physical and Chemistry Reference Data
!!    Vol. 27: 31-44, 1998
!!        
!!         visc = visc0 + delvisc + critvisc
!!                visc0:   viscosity in zero-density limit
!!                delvis:  excess viscosity outside critical region
!!		 critvis: viscosity at critical region ( ~ 0.0)
!!
!!  units: t(K), ro(kg/m3), visc*(microPa.s)                                               
!!
!      DO 31 m=1,nh
!        i=ihld(2,ih)
!        ix=ixhld(2,i)
!        t=temp(i)
!        tstar=t/251.196d0
!        ro=rho(ix)
!!     Viscosity in Zero-density Limit
!        tau=dlog(tstar)
!        tongg=0.0
!        DO 32 j=1,5
!          tongg=tongg+ac(j)*tau**(j-1)
!   32   END DO
!        tong=dexp(tongg)
!        visc0=1.00697d0*dsqrt(t)/tong
!!     Excess Viscosoty
!        delvis=dc(1)*ro+dc(2)*ro**2+dc(3)*ro**6/tstar**3              &
!     &         +dc(4)*ro**8+dc(5)*ro**8/tstar
!!     Total Viscosity
!        visc(i)=visc0+delvis
!        visc(i)=visc(i)/1.d6! conversion to Pa.s
!        ih=ih+iscskp
!   31 END DO
!      RETURN
!!---------------------------------------------------------------------------
!!
!   40 CONTINUE
!!
!!  He-4 Viscosity by Won-Jae Lee
!!  Ref: Chemical Properties Handbook 
!!      visc = A + B*T + C*T**2 (for GAS)
!!  http://www.knovel.com/knovel2/Toc.jsp?SpaceID=10093&BookID=49
!!    Properties for only gas are modeled, since
!!      Tcrit = 5.2 K
!!      Application Ranges: 150 K < T(K) < 2000 K
!!    units: t (K), visc (Pa.s)                                             
!!
!!      DO 41 m=1,nh
!!        i=ihld(2,ih)
!!        ix=ixhld(2,i)
!!        t=temp(i)
!!        visc0=71.094d0+4.43d-1*t-5.18d-5*t**2  ! in microP
!!		visc(i)=visc0/1.d7  ! conversion to Pa.s
!!        ih=ih+iscskp
!!   41 END DO
!!
!!  REPLACED WITH NEW NIST MODEL
!!  Ref: Arp, V.D., McCarty, R.D. and Friend, D.G.
!!  "Thermophysical Properties of Helium-4 from 0.8 to 1500 K with Pressures of 2000 MPa",
!!  NIST Technical Note 1334, Bolder, CO 1998
!!
!      DO 41 m=1,nh
!        i=ihld(2,ih)
!        ix=ixhld(2,i)
!        tt=temp(i)
!!
!        x=5.7037825d0
!        if (tt.le.300.0d0) x=dlog(tt)
!        ro=rho(ix)
!        dd=ro/1000.0d0
!        bh=-47.5295259d0/x+87.6799309d0-42.0741589d0*x   &
!           +8.33128289d0*x**2-0.589252385d0*x**3
!        ch= 547.309267d0/x-904.870586d0+431.404928d0*x   &
!            -81.4504854d0*x**2+5.37008433d0*x**3
!        dh=-1684.39324d0/x+3331.08630d0-1632.19172d0*x   &
!            +308.804413d0*x**2-20.2936367d0*x**3
!!
!        eta0a=dexp(-0.135311743d0/x+1.00347841d0+1.20654649d0*x   &
!              -0.149564551d0*x**2+0.0125208416d0*x**3)
!        etae=dexp(bh*dd+ch*dd**2+dh*dd**3)
!        if (tt.gt.100.0d0) then
!          eta0b=196.d0*tt**0.71938d0                             &
!                *dexp(12.451d0/tt-295.67d0/tt**2-4.1249d0)
!          if (tt.lt.110.0d0) then
!            eta0=eta0a+(eta0b-eta0a)*(tt-100.0d0)/10.0d0
!          else
!            eta0=eta0b
!          endif
!          visc0=eta0a*etae+eta0-eta0a
!        else
!          visc0=eta0a*etae
!        endif
!!
!        visc(i)=visc0/1.d7
!        ih=ih+iscskp
!   41 END DO
!      RETURN
!!
!!---------------------------------------------------------------------------
!!
!   50 CONTINUE
!!                                                                       
!!  H2 Dynamic Viscosity            
!!  NOT IMPLEMENTED YET - 
!!
!!
!         err=.true. 
!      RETURN
!!
!!---------------------------------------------------------------------------
!!
!   60 CONTINUE
!!
!!  O2 Dynamic Viscosity         
!!
!!  Ref: E.W. Lemmon and R.T. Jacobsen
!!       "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argoen and Air"
!!       Internaltional Journal of Thermophysics, Vol. 25, No. 1, Jan, 2004
!!
!!         visc = visc0 + delvisc
!!                visc0:   viscosity of dilute gas
!!                delvis:  residual fluid viscosity 
!!
!!  units: t(K), ro(kg/m3), visc*(microPa.s)                                               
!!
!!
!!
!      DO 61 m=1,nh
!        i=ihld(2,ih)
!        ix=ixhld(2,i)
!        tt=temp(i)
!        ro=rho(ix)
!!
!        tau=tco2/tt  ! Tc/T
!        delt=ro/roco2/mo2   ! ro/roc/conversion to mol/dm3
!        tstar=tt/eko2  ! T/(e/K)
!        tstar1=dlog(tstar)
!ohm=0.0d0
!DO ii=1,5
!  ohm=ohm+bi(ii)*tstar1**(ii-1)
!ENDDO
!ohm=dexp(ohm)
!        visc0=0.0266958d0*dsqrt(mo2*tt)/ohm/sigo2**2   ! dilute gas viscosity (micoPa.s)
!!     residual fluid viscosity
!        delvis=0.0d0
!DO ii=1,5
!  delvis=delvis+nio2(ii)*tau**tio2(ii)*delt**dio2(ii)*     &
!  dexp(-gamio2(ii)*delt**lio2(ii))
!ENDDO
!!		
!        visc(i)=(visc0+delvis)/1.d6
!        ih=ih+iscskp
!   61 END DO
!      RETURN
!!
!!---------------------------------------------------------------------------
!!
!   70 CONTINUE
!!
!!  N2 Dynamic Viscosity        
!!
!!  Ref: E.W. Lemmon and R.T. Jacobsen
!!       "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argoen and Air"
!!       Internaltional Journal of Thermophysics, Vol. 25, No. 1, Jan, 2004
!!
!!         visc = visc0 + delvisc
!!                visc0:   viscosity of dilute gas
!!                delvis:  residual fluid viscosity 
!!
!!  units: t(K), ro(kg/m3), visc*(microPa.s)                                               
!!
!!
!!
!      DO 71 m=1,nh
!        i=ihld(2,ih)
!        ix=ixhld(2,i)
!        tt=temp(i)
!        ro=rho(ix)
!!
!        tau=tcn2/tt  ! Tc/T
!!        delt=ro/rocn2/moln2   ! ro/roc/conversion to mol/dm3
!        delt=ro/rocn2/mn2   ! ro/roc/conversion to mol/dm3
!        tstar=tt/ekn2  ! T/(e/K)
!        tstar1=dlog(tstar)
!ohm=0.0d0
!DO ii=1,5
!  ohm=ohm+bi(ii)*tstar1**(ii-1)
!ENDDO
!ohm=dexp(ohm)
!        visc0=0.0266958d0*dsqrt(mn2*tt)/ohm/sign2**2   ! dilute gas viscosity (micoPa.s)
!!     residual fluid viscosity
!        delvis=0.0d0
!DO ii=1,5
!  delvis=delvis+nin2(ii)*tau**tin2(ii)*delt**din2(ii)*     &
!  dexp(-gamin2(ii)*delt**lin2(ii))
!ENDDO
!!		
!        visc(i)=(visc0+delvis)/1.d6
!        ih=ih+iscskp
!   71 END DO
!      RETURN
!!
!!---------------------------------------------------------------------------
!!>>
!!  D. LMR-K.S. Ha for liquid metal properties - Na
!!--sodium 
!!  "Thermophysical properties of Sodium" by G.H.Golden and J.V.Tokar
!!  Conversion factor : 1.73073467d0 
!!
!   80 IF(state.eq.'liquid')then 
!         GOTO 81 
!      ELSE 
!         GOTO 83 
!      ENDIF 
!!
!!
!!  Shpilrain et al. equation in data book of Golden and Tokar
!   81 DO m=1,nh 
!         i=ihld(2,ih) 
!         th=temp(i)
!         term=-1.6814d0-0.4296d0*log10(th)+234.65/th
!         visc(i)=10.d0**term*1.d-3 
!         ih=ih+iscskp 
!      END DO 
!      GOTO 999 
!!  term:Fahrenheit, Conversion factor:4.13379d-4
!   83 DO m=1,nh 
!         i=ihld(2,ih) 
!         term=max((temp(i)-273.15d0)*1.8d0+32.d0,100.d0) 
!         visc(i)=(0.019d0+(0.1375d-4+0.1709d-9*term)*term)*4.13379d-4
!         ih=ih+iscskp 
!      END DO 
!!
!      RETURN
!!
!!---------------------------------------------------------------------------
!!
!   90 CONTINUE
!!
!!
!!  Argon Dynamic Viscosity         
!!
!!  Ref: E.W. Lemmon and R.T. Jacobsen
!!       "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argoen and Air"
!!       Internaltional Journal of Thermophysics, Vol. 25, No. 1, Jan, 2004
!!
!!         visc = visc0 + delvisc
!!                visc0:   viscosity of dilute gas
!!                delvis:  residual fluid viscosity 
!!
!!  units: t(K), ro(kg/m3), visc*(microPa.s)                                               
!!
!!
!!
!      DO 91 m=1,nh
!        i=ihld(2,ih)
!        ix=ixhld(2,i)
!        tt=temp(i)
!        ro=rho(ix)
!!
!        tau=tcar/tt  ! Tc/T
!        delt=ro/rocar/mar   ! ro/roc/conversion to mol/dm3
!        tstar=tt/ekar  ! T/(e/K)
!        tstar1=dlog(tstar)
!ohm=0.0d0
!DO ii=1,5
!  ohm=ohm+bi(ii)*tstar1**(ii-1)
!ENDDO
!ohm=dexp(ohm)
!        visc0=0.0266958d0*dsqrt(mar*tt)/ohm/sigar**2   ! dilute gas viscosity (micoPa.s)
!!     residual fluid viscosity
!        delvis=0.0d0
!DO ii=1,5
!  delvis=delvis+niar(ii)*tau**tiar(ii)*delt**diar(ii)*     &
!  dexp(-gamiar(ii)*delt**liar(ii))
!ENDDO
!!		
!        visc(i)=(visc0+delvis)/1.d6
!        ih=ih+iscskp
!   91 END DO
!      RETURN
!!
!!---------------------------------------------------------------------------
!!
!  100 CONTINUE
!!
!!  Air Dynamic Viscosity        
!!
!!  Ref: E.W. Lemmon and R.T. Jacobsen
!!       "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argoen and Air"
!!       Internaltional Journal of Thermophysics, Vol. 25, No. 1, Jan, 2004
!!
!!         visc = visc0 + delvisc
!!                visc0:   viscosity of dilute gas
!!                delvis:  residual fluid viscosity 
!!
!!  units: t(K), ro(kg/m3), visc*(microPa.s)                                               
!!
!!
!!
!      DO 101 m=1,nh
!        i=ihld(2,ih)
!        ix=ixhld(2,i)
!        tt=temp(i)
!        ro=rho(ix)
!!
!        tau=tcair/tt  ! Tc/T
!        delt=ro/rocair/mair   ! ro/roc/conversion to mol/dm3
!        tstar=tt/ekair  ! T/(e/K)
!        tstar1=dlog(tstar)
!ohm=0.0d0
!DO ii=1,5
!  ohm=ohm+bi(ii)*tstar1**(ii-1)
!ENDDO
!ohm=dexp(ohm)
!        visc0=0.0266958d0*dsqrt(mair*tt)/ohm/sigair**2   ! dilute gas viscosity (micoPa.s)
!!     residual fluid viscosity
!        delvis=0.0d0
!DO ii=1,5
!  delvis=delvis+niair(ii)*tau**tiair(ii)*delt**diair(ii)*     &
!  dexp(-gamiair(ii)*delt**liair(ii))
!ENDDO
!!		
!        visc(i)=(visc0+delvis)/1.d6
!        ih=ih+iscskp
!  101 END DO
!      RETURN
!!
!!---------------------------------------------------------------------------
!!>>
!!  D. LMR-K.S. Ha for liquid metal properties - LBE
!!--Lead Bismuth Eutetic
!!  2nd meeting NEA Nuclear Science Committee ....
!!  "Thermophysical properties of LBE Alloy for use in reactor safety analysis" 
!!  by Koji Morita et al.
!!
!  110 IF(state.eq.'liquid')then 
!         GOTO 111 
!      ELSE 
!         GOTO 113 
!      ENDIF 
!!
!  111 DO m=1,nh 
!         i=ihld(2,ih) 
!         th=temp(i)
!         visc(i)=0.490d-3*exp(760.1d0/th)
!         ih=ih+iscskp 
!      END DO 
!      GOTO 999 
!!
!!  LBE viscosity for vapor phase is not known, thus following 
!!  correlation from the ATHENA program by K. E. Carlson has not meaning
!  113 DO m=1,nh 
!         i=ihld(2,ih) 
!         th=temp(i)
!         term=log(th)
!         visc(i)=2.16867957522d-4-6.72194914063d-5*th+5.39196834628d-6*th**2
!         ih=ih+iscskp 
!      END DO 
!      RETURN
!!<<
!!
!!
!  999 CONTINUE 
!      RETURN 
!      END SUBROUTINE viscos_cupid        
