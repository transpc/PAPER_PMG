!
      SUBROUTINE thcond_cupid_rv (tk,rhoc,convt)
!
!     SUBROUTINE cond (tk,rhoc,convt) written by WJL on Dec. '97
!     This SUBROUTINE calculates the LIGHT WATER THERMAL CONDUCTIVITY
!     based on ASME '93 Steam Table, AppENDix 7
!     Input quantities are t(K) and rhoc(kg/m3)
!     Output quantities is convt(W/K.m)
!
      IMPLICIT NONE
!  
!     Constants for Light Water Property Generation
!
      REAL(8) A0,A1,A2,A3,B0,B1,B2,BB1,BB2,C1,C2,C3,C4,C5,C6,D1,D2,D3,D4,Rhostr,Tstar
      REAL(8) tratio,dratio,dtstar,dtstar06,q,r,s,ylamb,xlamb,dlamb,barlam,baslam
!      REAL(8) tk(*),rhoc(*),convt(*)
      REAL(8) tk,rhoc,convt
      REAL(8) t0,t1

      DATA A0, A1   / +1.02811D-2, +2.99621D-2 /
      DATA A2, A3   / +1.56146D-2, -4.22464D-3 /
      DATA B0, B1   / -3.97070D-1, +4.00302D-1 /
      DATA B2       / +1.06000D+0 /
      DATA BB1, BB2 / -1.71587D-1, +2.39219D+0 /
      DATA C1, C2   / +6.42857D-1, -4.11717D+0 /
      DATA C3, C4   / -6.17937D+0, +3.08976D-3 /
      DATA C5, C6   / +8.22994D-2, +1.00932D+1 /
      DATA D1, D2   / +7.01309D-2, +1.18520D-2 /
      DATA D3, D4   / +1.69937D-3, -1.02000D+0 /
      DATA RHOSTR   / +317.7D+0 /
      DATA TSTAR    / +647.3D+0 /
!
!     IF (TE .LT. TMIN .OR. TE .GT. TMAX) CALL STER(LEV5, 12, PE, TE)
!
!
!      TRATIO = tk(i)/TSTAR
      TRATIO = tk/TSTAR
!      DRATIO = rhoc(i)/RHOSTR
      DRATIO = rhoc/RHOSTR
      DTSTAR = DABS(TRATIO - 1.0D+0) + C4
      DTSTAR06=DTSTAR**(-0.6D+0)
      Q = 2.0D+0 + C5*DTSTAR06
      R = Q + 1.0D+0
      S = 1.0D+0/DTSTAR
!      IF (tk(i) .LT. TSTAR) S=C6*DTSTAR06
      IF (tk .LT. TSTAR) S=C6*DTSTAR06
!      YLAMB = C2*dsqrt(TRATIO*TRATIO*TRATIO)+C3*dsqrt(RHOSTR/rhoc(i))
      YLAMB = C2*dsqrt(TRATIO*TRATIO*TRATIO)+C3*dsqrt(RHOSTR/rhoc)
!
!.....THE TEST ON YLAMB IS TO PREVENT ERROR RETURNS FROM THE
!     THE SYSTEM EXPONENTIAL SUBROUTINE.  THE VALUE OF -690 REPRESENTS
!     LN(10E-300).  THE LIMIT MAY BE DIFFERENT ON THE USERS COMPUTER.
!
      IF (YLAMB .LT. -690.0D+0) YLAMB = -690.0D+0
      XLAMB = D4*DEXP(YLAMB)
      t0=DRATIO**Q
      t1=DRATIO**1.8D+0
!      DLAMB = (D1*(TSTAR/tk(i))**10 + D2)*t1*DEXP(C1*       &
!       (1.0D+0 - DRATIO*t1))+D3*S*t0*DEXP((Q/R)*(1.0D+0  &
!        - DRATIO*t0)) + XLAMB
      DLAMB = (D1*(TSTAR/tk)**10 + D2)*t1*DEXP(C1*       &
       (1.0D+0 - DRATIO*t1))+D3*S*t0*DEXP((Q/R)*(1.0D+0  &
        - DRATIO*t0)) + XLAMB
      BARLAM = B0 + B1*DRATIO + B2*DEXP(BB1*(DRATIO+BB2)*(DRATIO+BB2))
      BASLAM = DSQRT(TRATIO)*(A0 + TRATIO*(A1 + TRATIO*(A2+ TRATIO*A3)))
!      convt(i) = BASLAM + BARLAM + DLAMB
      convt = BASLAM + BARLAM + DLAMB
!
!      convt(i)= dmax1(convt(i),0.01d0)
      convt= dmax1(convt,0.01d0)      
!
!
      RETURN
      END SUBROUTINE thcond_cupid_rv
!

!
      SUBROUTINE viscos_rv (temp,rhoc,satt,term)
!
!     SUBROUTINE cond (tk,rhoc,convt) written by WJL on Dec. '97
!     This SUBROUTINE calculates the LIGHT WATER THERMAL CONDUCTIVITY
!     based on ASME '93 Steam Table, AppENDix 7
!     Input quantities are t(K) and rhoc(kg/m3)
!     Output quantities is term(Pa.s)
!
      IMPLICIT NONE
!
!.....Data Statements:
!     For Light Water Properties
!
      REAL(8) tk
      REAL(8) ts1,ts2,ts3
      REAL(8) r1
      REAL(8) dr1,dr2,dr3,dr4
      REAL(8) dt1,dt2,dt3,dt4,dt5
      REAL(8) sum,sum1,sum2,sum3,sum4,sum5,sum6
      REAL(8) temp,rhoc,satt,term 
      REAL(8) a(4),bb(6,5) 
      DATA a  /0.0181583d0,0.0177624d0,0.0105287d0,-0.0036744d0/       
      DATA bb /0.501938d0,0.162888d0,-0.130356d0,0.907919d0,-0.551119d0,        &
               0.146543d0,0.235622d0,0.789393d0,0.673665d0,1.207552d0,           &
               0.0670665d0,-0.0843370d0,-0.274637d0,-0.743539d0,-0.959456d0,     &
              -0.687343d0,-0.497089d0,0.195286d0,0.145831d0,0.263129d0,          &
               0.347247d0,0.213486d0,0.100754d0,-0.032932d0,-0.0270448d0,        &
              -0.0253093d0,-0.0267758d0,-0.0822904d0,0.0602253d0,-0.0202595d0/   
!
!.....Test ranges of validity
!
      tk=dmax1(273.16d0,dmax1(temp,satt))
      ts1=647.27d0/tk
      ts2=ts1*ts1
!      ts2=418958.4529d0/tk/tk
      ts3=ts2*ts1
!      ts3=271179237.8d0/tk/tk/tk
      sum= a(1)+a(2)*ts1+a(3)*ts2 +a(4)*ts3
      term=1.d-06*dsqrt(1.d0/ts1)/sum
!
      IF(tk.eq.647.27d0)  return
!
      r1=rhoc/317.763d0
      dr1=r1-1.0d0
      dr2=dr1*dr1
      dr3=dr2*dr1
      dr4=dr3*dr1
      dt1=ts1-1.0d0
!      dt1=(647.27d0-tk)/tk
      dt2=dt1*dt1
      dt3=dt2*dt1
      dt4=dt3*dt1
      dt5=dt4*dt1
      sum1=bb(1,1)+bb(1,2)*dr1+bb(1,3)*dr2+bb(1,4)*dr3+bb(1,5)*dr4
      sum2=bb(2,1)+bb(2,2)*dr1+bb(2,3)*dr2+bb(2,4)*dr3+bb(2,5)*dr4
      sum3=bb(3,1)+bb(3,2)*dr1+bb(3,3)*dr2+bb(3,4)*dr3+bb(3,5)*dr4
      sum4=bb(4,1)+bb(4,2)*dr1+bb(4,3)*dr2+bb(4,4)*dr3+bb(4,5)*dr4
      sum5=bb(5,1)+bb(5,2)*dr1+bb(5,3)*dr2+bb(5,4)*dr3+bb(5,5)*dr4
      sum6=bb(6,1)+bb(6,2)*dr1+bb(6,3)*dr2+bb(6,4)*dr3+bb(6,5)*dr4
      sum=sum1+sum2*dt1+sum3*dt2+sum4*dt3+sum5*dt4+sum6*dt5
!
      term=term*dexp(r1*sum)
!
      RETURN
      END SUBROUTINE viscos_rv      
