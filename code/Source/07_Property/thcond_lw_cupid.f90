!!
!      SUBROUTINE cond_lw_cupid (tk,rhoc,convt,n)
!!
!!     SUBROUTINE cond (tk,rhoc,convt) written by WJL on Dec. '97
!!     This SUBROUTINE calculates the LIGHT WATER THERMAL CONDUCTIVITY
!!     based on ASME '93 Steam Table, AppENDix 7
!!     Input quantities are t(K) and rhoc(kg/m3)
!!     Output quantities is convt(W/K.m)
!!
!      IMPLICIT NONE
!!  
!!     Constants for Light Water Property Generation
!!
!      INTEGER(4) i,n
!      REAL(8) A0,A1,A2,A3,B0,B1,B2,BB1,BB2,C1,C2,C3,C4,C5,C6,D1,D2,D3,D4,Rhostr,Tstar
!      REAL(8) tratio,dratio,dtstar,dtstar06,q,r,s,ylamb,xlamb,dlamb,barlam,baslam
!!      REAL(8) tk(*),rhoc(*),convt(*)
!      REAL(8) tk(n),rhoc(n),convt(n)
!      REAL(8) t0,t1
!!
!      DATA A0, A1   / +1.02811D-2, +2.99621D-2 /
!      DATA A2, A3   / +1.56146D-2, -4.22464D-3 /
!      DATA B0, B1   / -3.97070D-1, +4.00302D-1 /
!      DATA B2       / +1.06000D+0 /
!      DATA BB1, BB2 / -1.71587D-1, +2.39219D+0 /
!      DATA C1, C2   / +6.42857D-1, -4.11717D+0 /
!      DATA C3, C4   / -6.17937D+0, +3.08976D-3 /
!      DATA C5, C6   / +8.22994D-2, +1.00932D+1 /
!      DATA D1, D2   / +7.01309D-2, +1.18520D-2 /
!      DATA D3, D4   / +1.69937D-3, -1.02000D+0 /
!      DATA RHOSTR   / +317.7D+0 /
!      DATA TSTAR    / +647.3D+0 /
!!
!!     IF (TE .LT. TMIN .OR. TE .GT. TMAX) CALL STER(LEV5, 12, PE, TE)
!!
!      DO i=1,n
!!
!      TRATIO = tk(i)/TSTAR
!      DRATIO = rhoc(i)/RHOSTR
!      DTSTAR = DABS(TRATIO - 1.0D+0) + C4
!      DTSTAR06=DTSTAR**(-0.6D+0)
!      Q = 2.0D+0 + C5*DTSTAR06
!      R = Q + 1.0D+0
!      S = 1.0D+0/DTSTAR
!      IF (tk(i) .LT. TSTAR) S=C6*DTSTAR06
!      YLAMB = C2*dsqrt(TRATIO*TRATIO*TRATIO)+C3*(RHOSTR/rhoc(i))**5 
!!
!!.....THE TEST ON YLAMB IS TO PREVENT ERROR RETURNS FROM THE
!!     THE SYSTEM EXPONENTIAL SUBROUTINE.  THE VALUE OF -690 REPRESENTS
!!     LN(10E-300).  THE LIMIT MAY BE DIFFERENT ON THE USERS COMPUTER.
!!
!      IF (YLAMB .LT. -690.0D+0) YLAMB = -690.0D+0
!      XLAMB = D4*DEXP(YLAMB)
!      t0=DRATIO**Q
!      t1=DRATIO**1.8D+0
!      DLAMB = (D1*(TSTAR/tk(i))**10 + D2)*t1*DEXP(C1*       &
!       (1.0D+0 - DRATIO*t1))+D3*S*t0*DEXP((Q/R)*(1.0D+0  &
!        - DRATIO*t0)) + XLAMB
!      BARLAM = B0 + B1*DRATIO + B2*DEXP(BB1*(DRATIO+BB2)*(DRATIO+BB2))
!      BASLAM = DSQRT(TRATIO)*(A0 + TRATIO*(A1 + TRATIO*(A2+ TRATIO*A3)))
!      convt(i) = BASLAM + BARLAM + DLAMB
!!
!      convt(i)= dmax1(convt(i),0.01d0)
!!
!      ENDDO
!!
!      RETURN
!      END SUBROUTINE cond_lw_cupid
!
!
      SUBROUTINE cond_lw_cupid (tk,rhoc,convt,n)
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
      INTEGER(4) i,n
      REAL(8) A0,A1,A2,A3,B0,B1,B2,BB1,BB2,C1,C2,C3,C4,C5,C6,D1,D2,D3,D4,Rhostr,Tstar
      REAL(8) Rhostrr,Tstarr
      REAL(8) tratio,dratio,dtstar,dtstar06,q,r,s,ylamb,xlamb,dlamb,barlam,baslam
!      REAL(8) tk(*),rhoc(*),convt(*)
      REAL(8) tk(n),rhoc(n),convt(n)
      REAL(8) t0,t1
!
      DATA A0, A1   / +1.02811D-2, +2.99621D-2 /
      DATA A2, A3   / +1.56146D-2, -4.22464D-3 /
      DATA B0, B1   / -3.97070D-1, +4.00302D-1 /
      DATA B2       / +1.06000D0 /
      DATA BB1, BB2 / -1.71587D-1, +2.39219D0  /
      DATA C1, C2   / +6.42857D-1, -4.11717D0  /
      DATA C3, C4   / -6.17937D0 , +3.08976D-3 /
      DATA C5, C6   / +8.22994D-2, +1.00932D1  /
      DATA D1, D2   / +7.01309D-2, +1.18520D-2 /
      DATA D3, D4   / +1.69937D-3, -1.02000D0  /
      DATA RHOSTR   / +317.7D0 /
      DATA TSTAR    / +647.3D0 /
!
!     IF (TE .LT. TMIN .OR. TE .GT. TMAX) CALL STER(LEV5, 12, PE, TE)
!
        TSTARR=1.d0/TSTAR
        RHOSTRR=1.d0/RHOSTR
      DO i=1,n
         TRATIO = tk(i)*TSTARR
         DRATIO = rhoc(i)*RHOSTRR
!        TRATIO = tk(i)/TSTAR
!        DRATIO = rhoc(i)/RHOSTR
         DTSTAR = ABS(TRATIO - 1.D0) + C4
         DTSTAR06=DTSTAR**(-0.6D0)
         Q = 2.D0 + C5*DTSTAR06
         R = Q + 1.D0
         IF (tk(i) .LT. TSTAR) THEN
            S=C6*DTSTAR06
         ELSE
            S = 1.D0/DTSTAR
         ENDIF
         YLAMB = C2*sqrt(TRATIO*TRATIO*TRATIO)+C3*(RHOSTR/rhoc(i))**5 
!
!.....THE TEST ON YLAMB IS TO PREVENT ERROR RETURNS FROM THE
!     THE SYSTEM EXPONENTIAL SUBROUTINE.  THE VALUE OF -690 REPRESENTS
!     LN(10E-300).  THE LIMIT MAY BE DIFFERENT ON THE USERS COMPUTER.
!
!        IF (YLAMB .LT. -690.D0) YLAMB = -690.D0
         YLAMB=MAX(YLAMB,-690.D0)
         XLAMB = D4*EXP(YLAMB)
         t0=DRATIO**Q
         t1=DRATIO**1.8D0
         DLAMB = (D1*(TSTAR/tk(i))**10 + D2)*t1*EXP(C1   *(1.D0 - DRATIO*t1))  &
                +                       D3*S*t0*EXP((Q/R)*(1.D0 - DRATIO*t0))  &
                + XLAMB
         BARLAM = B0 + B1*DRATIO + B2*EXP(BB1*(DRATIO+BB2)*(DRATIO+BB2))
         BASLAM = SQRT(TRATIO)*(A0 + TRATIO*(A1 + TRATIO*(A2+ TRATIO*A3)))
         convt(i) = BASLAM + BARLAM + DLAMB
!
         convt(i)= MAX(convt(i),0.01d0)
      ENDDO   
!
      END SUBROUTINE cond_lw_cupid
