!
MODULE Zrv_ihtc_models
!
      USE Vol_Data
      USE Zrv_htc
!
      IMPLICIT NONE
!
CONTAINS     
! 
!--------------------------------------------------------------
!  Function: Plesset_Zwick
!--------------------------------------------------------------
   FUNCTION Plesset_Zwick(dia,DelTsl,rhog,rhol,hfg,cpl,condl) RESULT(htc)
!
      IMPLICIT NONE
!      
      REAL(8) dia,DelTsl,rhog,rhol,hfg,cpl,condl
      REAL(8) math_pi
      REAL(8) htc
!
      math_pi=4.0*DATAN(1.0)
!      
      htc=-condl/dia*12.0/math_pi*DelTsl*rhol*cpl/(rhog*hfg)
!
   END FUNCTION Plesset_Zwick
!      
!--------------------------------------------------------------
!  Function: Lee_Ryley
!--------------------------------------------------------------
   FUNCTION Lee_Ryley(dia,kcond,RE,a,b) RESULT(htc)
!      
      IMPLICIT NONE
!      
      REAL(8) dia,kcond,RE,a,b
      REAL(8) htc
!
      htc=kcond/dia*(a+b*RE**0.5)
!
   END FUNCTION Lee_Ryley      
!      
!--------------------------------------------------------------
!  Function: Unal
!--------------------------------------------------------------
!  FUNCTION Unal(dia,Dhyd,rhog,rhol,hfg,press,vl,alphab) RESULT(htc)
   FUNCTION Unal(dia,rhog,rhol,hfg,press,vl,alphab) RESULT(htc)
!
      IMPLICIT NONE
!      
!     REAL(8) dia,Dhyd,rhog,rhol,hfg,press,vl,alphab
      REAL(8) dia,rhog,rhol,hfg,press,vl,alphab
      REAL(8) C,phi,cF5
      REAL(8) htc
!
      IF(press.le.1.1272D6) THEN
         C=65.0-5.69D-5*(press-1.0D5)
      ELSE
         C=2.5D9*press**(-1.418)
      ENDIF
!
      phi=DMAX1(1.d0,DABS(vl/0.61d0)**0.47)
      IF(alphab<0.25d0) THEN
         cF5=1.8d0*C*phi*DEXP(-45.d0*alphab)+0.075d0
      ELSE
         cF5=0.075d0
      ENDIF
!
      htc=cF5*hfg*rhog*rhol*dia/3.6d0/DMAX1(rhol-rhog,1.0D-7)
!
   END FUNCTION Unal   
!      
!--------------------------------------------------------------
!  Function: Sieder_Tate
!--------------------------------------------------------------
!  FUNCTION Sieder_Tate(condl,rhol,mul,cpl,dia,Dhyd,vg,vl) RESULT(htc)
   FUNCTION Sieder_Tate(condl,rhol,mul,cpl,Dhyd,vg,vl) RESULT(htc)
!   
      IMPLICIT NONE
!      
!     REAL(8) condl,rhol,mul,cpl,dia,Dhyd,vg,vl
      REAL(8) condl,rhol,mul,cpl,Dhyd,vg,vl
      REAL(8) Ref,Prf,vfg
      REAL(8) htc
!
      vfg=DMIN1(DABS(vg-vl),0.8)
      Ref=rhol*Dhyd*vfg/mul
      Prf=cpl*mul/condl
!
      htc=1.18942*(Ref*Prf)**0.5*condl/Dhyd
!
   END FUNCTION Sieder_Tate   
!      
!--------------------------------------------------------------
!  Function: Theofanous
!--------------------------------------------------------------   
   FUNCTION Theofanous(rhol,cpl,vl) RESULT(htc)
!
      IMPLICIT NONE
!      
      REAL(8) rhol,cpl,vl
      REAL(8) htc
!
      htc=1.0D-3*rhol*cpl*DABS(vl)
!
   END FUNCTION Theofanous   
!      
!--------------------------------------------------------------
!  Function: Dittus_Boelter
!--------------------------------------------------------------   
   FUNCTION Dittus_Boelter(Dhyd,cond,rho,mu,vfg,alpha) RESULT(htc)
!
      IMPLICIT NONE
!      
      REAL(8) Dhyd,cond,rho,mu,vfg,alpha
      REAL(8) Re
      REAL(8) htc
!
!reflod model by LSJ 180123
!      IF(.not.reflod) THEN
         Re=alpha*rho*DABS(vfg)*Dhyd/mu
!      ELSE
!         Re=rey_reflod
!      ENDIF
!
      htc=0.023*Re**0.8*cond/Dhyd
!
   END FUNCTION Dittus_Boelter   
!      
!--------------------------------------------------------------
!  Function: Umbrella_cap for SCL condition
!--------------------------------------------------------------
   FUNCTION umbrella_cap(Hil,alphag,pres) RESULT(h)
!
      IMPLICIT NONE
!      
      REAL(8) Hil,alphag,pres
      REAL(8) HilMax,xp
      REAL(8) h
      REAL(8), PARAMETER :: plow=8.618D6
      REAL(8), PARAMETER :: dp=10.342D6-plow
!
      HilMax=17539.0*DMAX1(4.724,472.4*alphag*(1-alphag)) &
                    *DMAX1(0.0,DMIN1((alphag-1.0D-10)/(0.1-1.0D-10),1.0))
      xp=DMAX1(0.0,DMIN1((pres-plow)/dp,1.0))
!
      h=DMIN1(Hil,Hil*xp+HilMax*(1.0-xp))
!
   END FUNCTION umbrella_cap   
!      
!--------------------------------------------------------------
!  Function: Power-interpolation between the superheated and subcooled HTCs
!--------------------------------------------------------------
   FUNCTION pow_interp(DelTs,Hsup,Hsub) RESULT(Hi)
!
      IMPLICIT NONE
!      
      REAL(8) DelTs,Hsup,Hsub
      REAL(8) x1,x
      REAL(8) Hi
!
      x1=DMAX1(0.d0,DMIN1(0.5d0*(DelTs+1.d0),1.d0))
      x=x1**2*(3.d0-2.d0*x1)
!
      Hi=Hsub**x*Hsup**(1.d0-x)
!      
   END FUNCTION pow_interp   
!
!------------------------------------------------------------------------
!  FUNCTION: F1
!------------------------------------------------------------------------
   FUNCTION F1(alphab)
!   
      IMPLICIT NONE
!      
      REAL(8) alphab
      REAL(8) F1
      REAL(8), PARAMETER :: alpha_trig=0.01  !< alpha_trig = 0.001 in RELAP5/MOD3.3
!
      F1=DMIN1(alpha_trig,alphab)/alphab
!      
   END FUNCTION F1
!
!------------------------------------------------------------------------
!  FUNCTION: F2
!------------------------------------------------------------------------
   FUNCTION F2(alphab)
!
      IMPLICIT NONE
!      
      REAL(8) alphab
      REAL(8) F2
      REAL(8), PARAMETER :: alpha_trig=0.25
!
      F2=DMIN1(alpha_trig,alphab)/alphab
!      
   END FUNCTION F2
!
!------------------------------------------------------------------------
!  FUNCTION: F3
!------------------------------------------------------------------------
   FUNCTION F3(alphag,Xn,DelTs)
!   
      IMPLICIT NONE
!      
      REAL(8) alphag,Xn,DelTs
      REAL(8) F3,F4
!
!.....Effect of NC at low void fraction
!
      F4=DMIN1(1.0D-5,alphag*(1.d0-Xn))*1.0D5
!
      IF(DelTs<=-1.0)THEN
         F3=1.0
      ELSEIF(0.0<=DelTs)THEN
         F3=DMAX1(0.0,F4)
      ELSE
         F3=DMAX1(0.0,F4*(1+DelTs)-DelTs) ! interpolation between F4 and 1 according to  (-DelTs).
      ENDIF
!      
   END FUNCTION F3
!
!------------------------------------------------------------------------
!  FUNCTION: F6
!------------------------------------------------------------------------
   FUNCTION F6(DelTs)
!   
      IMPLICIT NONE
!      
      REAL(8) DelTs
      REAL(8) F6,eta
!
!.....Effective up to 2 K degree of superheat.
!
      eta=DABS(DMAX1(-2.d0,DelTs))
!
      F6=1.d0+eta*(100.d0+25.d0*eta)
!
   END FUNCTION F6
!
!------------------------------------------------------------------------
!  FUNCTION: F7
!------------------------------------------------------------------------
   FUNCTION F7(alphag)
!   
      IMPLICIT NONE
!      
      REAL(8) alphag
      REAL(8) F7
!
      F7=DMAX1(alphag,1.0D-5)/DMAX1(alphag,1.0D-9)
!   
   END FUNCTION F7   
!
!------------------------------------------------------------------------
!  FUNCTION: F9
!------------------------------------------------------------------------   
   FUNCTION F9(alphag,alpha_bs,alpha_sa)
!
      IMPLICIT NONE
!      
      REAL(8) alphag,alpha_bs,alpha_sa
      REAL(8) F9
!
      F9=DEXP(-8.0*(alphag-alpha_bs)/     &
         DMAX1(1.0D-6,(alpha_sa-alpha_bs)))
!         
   END FUNCTION F9   
!
!------------------------------------------------------------------------
!  FUNCTION: F10
!------------------------------------------------------------------------      
!  FUNCTION F10(alphag,rhog,rhol,sigma,vg,vl,ggc,lambda)
   FUNCTION F10(alphag,rhog,rhol,sigma,vg,ggc,lambda)
!
      IMPLICIT NONE
!      
!     REAL(8) alphag,rhog,rhol,sigma,vg,vl,ggc,lambda
      REAL(8) alphag,rhog,rhol,sigma,vg,ggc,lambda
      REAL(8) sigma_st,vcrit,vcrit_tmp
      REAL(8) F10
      REAL(8) dr
!
      sigma_st=DMAX1(sigma,1.0D-7)
      dr=DMAX1(1.d-5,rhol-rhog)
      vcrit_tmp=DSQRT(sigma_st*ggc*dr)
      vcrit=3.2d0*DSQRT(vcrit_tmp/rhog)
      lambda=DABS(alphag*vg/vcrit)
!
      F10=DMIN1(1.d0+DSQRT(lambda)+0.05d0*lambda,6.d0)
!      
   END FUNCTION F10   
!
!------------------------------------------------------------------------
!  FUNCTION: F12
!------------------------------------------------------------------------      
   FUNCTION F12(delTs)
!
      IMPLICIT NONE
!      
      REAL(8) delTs
      REAL(8) F12
      REAL(8) ksi
!
      ksi=DMAX1(0.d0,-DelTs)
!      
!      F12=1.d0+ksi*(250.d0+50.d0*ksi)    
      F12=1.0d0+0.25d0*ksi*ksi  ! yjm  ref. MARS source code
!      
   END FUNCTION F12
!
!------------------------------------------------------------------------
!  FUNCTION: F13
!------------------------------------------------------------------------      
   FUNCTION F13(delTs,cpf,hfg)
!
      IMPLICIT NONE
!      
      REAL(8) delTs,cpf,hfg
      REAL(8) F13
!
      F13=2.d0+7.0*DMIN1(1.d0+cpf*DMAX1(DelTs,0.d0)/hfg,8.d0)
!      
   END FUNCTION F13  
!
!------------------------------------------------------------------------
!  FUNCTION: F17
!------------------------------------------------------------------------      
   FUNCTION F17(i,xag)
!
      IMPLICIT NONE
!      
      INTEGER i
!      
      REAL(8) xag
      REAL(8) F17,cF18,alphaian
 !
      alphaian=xag
      cF18=DMIN1(xag/0.05d0,0.999999d0)
!      
      F17=DEXP(-8.d0*(cell%alpha_bs(i)-alphaian)/cell%alpha_bs(i))*cF18
!      
   END FUNCTION F17    
!
!------------------------------------------------------------------------
!  FUNCTION: F19
!------------------------------------------------------------------------      
   FUNCTION F19(delTs)
!
      IMPLICIT NONE
!      
      REAL(8) delTs
      REAL(8) F19
!
      F19=2.5d0-DelTs*(0.2d0-0.1d0*DelTs)
!      
   END FUNCTION F19
!
!------------------------------------------------------------------------
!  FUNCTION: F21
!------------------------------------------------------------------------      
   FUNCTION F21(i,alphag)
!
      IMPLICIT NONE
!      
      INTEGER i
!      
      REAL(8) alphag
      REAL(8) F21
!
      F21=DEXP((cell%alpha_sa(i)-alphag)/(cell%alpha_sa(i)-cell%alpha_bs(i)))
!      
   END FUNCTION F21
!
!------------------------------------------------------------------------
!  FUNCTION: F22
!------------------------------------------------------------------------      
   FUNCTION F22(alphag)
!
      IMPLICIT NONE
!      
      REAL(8) alphag
      REAL(8) F22,ag4
!
      ag4=alphag/4.d0
!      
      F22=DMAX1(0.02d0,DMIN1(ag4*(1.d0-ag4),0.2d0))
!      
   END FUNCTION F22
!
!------------------------------------------------------------------------
!  FUNCTION: F23
!------------------------------------------------------------------------      
   FUNCTION F23(alphal)
!      
      IMPLICIT NONE
!      
      REAL(8) alphal
      REAL(8) F23,alphadrp
!
      alphadrp=DMAX1(alphal,1.d-4)
!      
      F23=alphadrp/DMAX1(alphal,1.d-12)
!      
   END FUNCTION F23
!
!------------------------------------------------------------------------
!  FUNCTION: F24
!------------------------------------------------------------------------      
   FUNCTION F24(delts,alphal)
!   
      IMPLICIT NONE
!      
      REAL(8) delts,alphal
      REAL(8) F24
      REAL(8) cF25,cF26
!
      cF25=10*10*10*10*10*DMIN1(alphal,1.d-5) 
      cF26=1.d0-5.d0*DMIN1(0.2d0,DMAX1(0.d0,delts))
!
      F24=DMAX1(0.d0,cF26*(cF25-1.d0)+1.d0)
!      
   END FUNCTION F24   
!
!------------------------------------------------------------------------
!  FUNCTION: F30
!------------------------------------------------------------------------      
   FUNCTION F30(i,alphal,vg,vl,delts)
!
      USE Zconst2      , ONLY: hydraulicd,ggc 
!            
      IMPLICIT NONE
!      
      INTEGER i
!      
      REAL(8) delts,alphal,vg,vl
      REAL(8) F30
      REAL(8) cF32,cF33,cF34
      REAL(8) vm,vtb,gm,rm,dr
!
      gm=cell%alphal(i)*cell%rhol(i)*DABS(vl)+cell%alphag(i)*cell%rhog(i)*DABS(vg)
      rm=cell%alphal(i)*cell%rhol(i)+cell%alphag(i)*cell%rhog(i)
      vm=gm/rm
!      dr=DABS(cell%rhol(i)-cell%rhog(i))
      dr=DMAX1(1.d-5,cell%rhol(i)-cell%rhog(i))
      vtb=0.35d0*DSQRT(ggc*hydraulicd(i)*dr/cell%rhol(i))
!
      cF32=1.d0-DMIN1(1.d0,100.d0*alphal)
      cF33=DMAX1(0.d0,2.d0*DMIN1(1.d0,vm/vtb)-1.d0)
      cF34=DMAX1(0.d0,DMIN1(1.d0,-0.5d0*delts))
!
      F30=DMAX1(cF32,cF33,cF34)
!      
   END FUNCTION F30
!
!------------------------------------------------------------------------
!  FUNCTION: F35
!------------------------------------------------------------------------      
!  FUNCTION F35(i,alphal,vg,vl,delts)
   FUNCTION F35(i,vg,vl,delts)
!      
      USE Zconst2         , ONLY: hydraulicd,ggc 
!            
      IMPLICIT NONE
!      
      INTEGER i
!     REAL(8) delts,alphal,vg,vl
      REAL(8) delts,vg,vl
      REAL(8) F35
      REAL(8) cF33,cF36
      REAL(8) vm,vtb,gm,rm,dr
!
      gm=cell%alphal(i)*cell%rhol(i)*DABS(vl)+cell%alphag(i)*cell%rhog(i)*DABS(vg)
      rm=cell%alphal(i)*cell%rhol(i)+cell%alphag(i)*cell%rhog(i)
      vm=gm/rm
!      dr=DABS(cell%rhol(i)-cell%rhog(i))
      dr=DMAX1(1.d-5,cell%rhol(i)-cell%rhog(i))
      vtb=0.35d0*DSQRT(ggc*hydraulicd(i)*dr/cell%rhol(i))
!
      cF33=DMAX1(0.d0,2.d0*DMIN1(1.d0,vm/vtb)-1.d0)
      cF36=DMIN1(1.d0,0.5d0*delts )
!
      F35=DMAX1(cF33,cF36)
!      
   END FUNCTION F35   
!
!------------------------------------------------------------------------
! FUNCTION: F16
!------------------------------------------------------------------------      
!   FUNCTION F16( alphag,alpha_bs )
!      IMPLICIT NONE
!      REAL(8), INTENT(IN) :: alphag
!      REAL(8) :: xag,alphaian,alpha_bs
!      REAL(8) :: F16,F17,F18
!!
!      xag=DMAX1(0.0,DMIN1(alphag,alpha_bs))
!      alphaian=xag
!!      
!      F18=DMIN1( xag/0.05d0, 0.999999d0 )
!      F17=DEXP( -8.d0*(alpha_bs-alphaian)/alpha_bs )*F18 
!      F16=DMIN1(1.d0-F17,1.d0)
!!      
!   END FUNCTION F16    
!
END MODULE Zrv_ihtc_models
