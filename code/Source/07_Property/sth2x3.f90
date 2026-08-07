!
      SUBROUTINE sth2x3_cupid(s,it,err,  &
                              a1,a2,a5)
!
!     compute water thermodynamic properties as a function of temperature and pressure
!
      USE STM_TBL_cupid  , ONLY: crt,crp,      &
                                 k,            &
                                 a31,a3,a41,a4, &
                                 nt,np,ns,ns2
!
      IMPLICIT NONE
!
!.....Input
      REAL(8) :: a1(nt)
      REAL(8) :: a2(np)
      REAL(8) :: a5(6,nt,np)
!.....Output
      INTEGER :: it
      LOGICAL :: err
      REAL(8) :: s(26)
!.....Local variables
      INTEGER :: i
      INTEGER :: ip,jp
      LOGICAL :: s1,s2 
      REAL(8) :: unp
      REAL(8) :: hfg1,hfg2
      REAL(8) :: dpdt1,dpdt2,cv,ren
      REAL(8) :: fr,fr1,frn,frc,frc2,frd,frd2,px
      REAL(8) :: f1,f2,c0,c1,c2,c3,d1,d2
      REAL(8) :: ta,tb,pa,pb
      REAL(8) :: a0(6,2),b0(6,2)
      REAL(8) :: ac(6),ad(6)
!
      INTEGER iunp(2) 
      EQUIVALENCE(unp,iunp(1)) 
!                                                                       
!.....check for valid input                                               
!
      IF(s(1).lt.a1(1).or.s(1).gt.5000.d0) GOTO 1001 
      IF(s(2).le.0.d0 .or.s(2).gt.a2(np) ) GOTO 1001 
!
      unp=s(23)
      ip=iunp(1)
      jp=iunp(2)
      IF(ip.le.0.or.ip.gt.nt) ip=1
      IF(jp.le.0.or.jp.ge.np) jp=1
!
      IF(s(1).ge.a1(nt)) THEN
         ip=nt 
         it=4 
         s(9)=1.d0 
         GOTO 50
      ENDIF
!
!.....set indexes in temperature and pressure tables for saturation computations
!
   11 IF(s(1).ge.a1(ip)) GOTO 10
      ip=ip-1
      GOTO 11
   10 IF(s(1).lt.a1(ip+1)) GOTO 12
      ip=ip+1
      GOTO 10
   12 IF(ip.ge.ns) THEN
         it=4
         s(9)=1.d0
         GOTO 50
      ENDIF
!
      s1=.false. 
      s2=.false. 
      fr=s(1)/crt 
      fr1=1.d0-fr 
      s(10)=crp*exp((((((k(5)*fr1+k(4))*fr1+k(3))*fr1+k(2))*fr1+k(1))*  &
      fr1)/(((k(7)*fr1+k(6))*fr1+1.d0)*fr)-fr1/(k(8)*fr1*fr1+k(9)))    
!
   13 IF(s(10).ge.a2(jp)) GOTO 14
      jp=jp-1
      IF(jp.gt.0) GOTO 13
      jp=1
      s2=.true.
      GOTO 15
   14 IF(s(10).lt.a2(jp+1)) GOTO 15
      jp=jp+1
      IF(jp.lt.np) GOTO 14
      s1=.true.
!
   15 CONTINUE
      IF(s2.or.a2(jp).le.a31(ip)) THEN
         ta=a1(ip)
         pa=a31(ip)
         DO i=1,6
            a0(i,1)=a3(i,1,ip)
            a0(i,2)=a3(i,2,ip)
         ENDDO
      ELSE
         ta=a41(jp)
         pa=a2(jp)
         DO i=1,6
            a0(i,1)=a4(i,1,jp)
            a0(i,2)=a4(i,2,jp)
         ENDDO
      ENDIF
      IF(s1.or.a2(jp+1).ge.a31(ip+1)) THEN
         tb=a1(ip+1)
         pb=a31(ip+1)
         DO i=1,6
            b0(i,1)=a3(i,1,ip+1)
            b0(i,2)=a3(i,2,ip+1)
         ENDDO
      ELSE
         pb=a2(jp+1)
         tb=a41(jp+1)
         DO i=1,6
            b0(i,1)=a4(i,1,jp+1)
            b0(i,2)=a4(i,2,jp+1)
         ENDDO
      ENDIF
!
!.....compare input pressure to saturation pressure to determine vapor or liquid
!
      fr1=s(1)-ta 
      fr=fr1/(tb-ta) 
      IF(s(2).le.s(10)) THEN
         it=3 
         s(9)=1.d0 
      ELSE
         it=1 
         s(9)=0.d0 
      ENDIF
   50 CONTINUE
!
!.....search for single phase indexes.                                    
!
   51 IF(s(2).ge.a2(jp)) GOTO 53
      jp=jp-1
      IF(jp.gt.0) GOTO 51
      GOTO 90
   53 IF(s(2).lt.a2(jp+1)) GOTO 54
      jp=jp+1
      IF(jp.lt.np) GOTO 53
      GOTO 1001
   54 CONTINUE
!
      frn=s(1)-a1(ip) 
      IF(it.eq.3) THEN
         GOTO 70
      ELSEIF(it.gt.3) THEN
         GOTO 73
      ELSEIF(it.lt.3) THEN
!
!........liquid phase.                                                       
!
         IF(a2(jp).ge.s(10)) THEN
            px=a2(jp) 
            s2=.true. 
            IF(jp.gt.ns2) GOTO 63
            IF(a41(jp).gt.a1(ip+1)) GOTO 63
               frc=frn/(a41(jp)-a1(ip)) 
               frc2=s(1)-a41(jp)
               DO i=1,6
                  ac(i)=a4(i,1,jp)
               ENDDO
               GOTO 62 
   63          frc=frn/(a1(ip+1)-a1(ip)) 
               frc2=s(1)-a1(ip+1) 
               frd=frc 
               frd2=frc2 
               DO i=1,6
                  ac(i)=a5(i,ip+1,jp)
                  ad(i)=a5(i,ip+1,jp+1)
               ENDDO
               GOTO 66 
         ELSE
            px=s(10) 
            s2=.false. 
         ENDIF
   62    IF(jp.ge.ns2) GOTO 67
         IF(a41(jp+1).gt.a1(ip+1)) GOTO 67
         frd=frn/(a41(jp+1)-a1(ip)) 
         frd2=s(1)-a41(jp+1) 
         DO i=1,6
            ad(i)=a4(i,1,jp+1)
         ENDDO
         GOTO 66
   67    frd=frn/(a1(ip+1)-a1(ip)) 
         frd2=s(1)-a1(ip+1) 
         DO i=1,6
            ad(i)=a5(i,ip+1,jp+1)
         ENDDO
!
   66    IF(s2) THEN
            d1=a5(3,ip,jp)*a5(1,ip,jp)*(frn-frc2) 
            d2=ac(3)*ac(1)*(frn-frc2) 
            c0=a5(1,ip,jp)
            c1=d1 
            c2=3.d0*(ac(1)-a5(1,ip,jp))-d2-2.d0*d1 
            c3=d2+d1-2.d0*(ac(1)-a5(1,ip,jp)) 
            s(11)=c0+frc*(c1+frc*(c2+frc*c3)) 
            s(13)=a5(2,ip,jp)+(ac(2)-a5(2,ip,jp))*frc 
            s(17)=a5(3,ip,jp)+(ac(3)-a5(3,ip,jp))*frc 
            s(19)=a5(4,ip,jp)+(ac(4)-a5(4,ip,jp))*frc 
            s(21)=a5(5,ip,jp)+(ac(5)-a5(5,ip,jp))*frc 
            s(25)=a5(6,ip,jp)+(ac(6)-a5(6,ip,jp))*frc 
         ELSE
            hfg1=a0(2,2)-a0(2,1)+pa*(a0(1,2)-a0(1,1))
            hfg2=b0(2,2)-b0(2,1)+pb*(b0(1,2)-b0(1,1)) 
            dpdt1=hfg1/(ta*(a0(1,2)-a0(1,1))) 
            dpdt2=hfg2/(tb*(b0(1,2)-b0(1,1)))
            f1=a0(1,1)*(a0(3,1)-a0(4,1)*dpdt1) 
            f2=b0(1,1)*(b0(3,1)-b0(4,1)*dpdt2) 
            d1=f1*(tb-ta) 
            d2=f2*(tb-ta) 
            c0=a0(1,1)
            c1=d1 
            c2=3.d0*(b0(1,1)-a0(1,1))-d2-2.d0*d1 
            c3=d2+d1-2.d0*(b0(1,1)-a0(1,1))
            s(11)=c0+fr*(c1+fr*(c2+fr*c3)) 
            s(13)=a0(2,1)+(b0(2,1)-a0(2,1))*fr 
            s(17)=a0(3,1)+(b0(3,1)-a0(3,1))*fr 
            s(19)=a0(4,1)+(b0(4,1)-a0(4,1))*fr 
            s(21)=a0(5,1)+(b0(5,1)-a0(5,1))*fr 
            s(25)=a0(6,1)+(b0(6,1)-a0(6,1))*fr 
         ENDIF
         d1=a5(3,ip,jp+1)*a5(1,ip,jp+1)*(frn-frd2) 
         d2=ad(3)*ad(1)*(frn-frd2) 
         c0=a5(1,ip,jp+1)
         c1=d1 
         c2=3.d0*(ad(1)-a5(1,ip,jp+1))-d2-2.d0*d1 
         c3=d2+d1-2.d0*(ad(1)-a5(1,ip,jp+1)) 
         s(12)=c0+frd*(c1+frd*(c2+frd*c3)) 
         s(14)=a5(2,ip,jp+1)+(ad(2)-a5(2,ip,jp+1))*frd 
         s(18)=a5(3,ip,jp+1)+(ad(3)-a5(3,ip,jp+1))*frd 
         s(20)=a5(4,ip,jp+1)+(ad(4)-a5(4,ip,jp+1))*frd 
         s(22)=a5(5,ip,jp+1)+(ad(5)-a5(5,ip,jp+1))*frd 
         s(26)=a5(6,ip,jp+1)+(ad(6)-a5(6,ip,jp+1))*frd 
!
         IF(s(11).gt.s(12)) THEN
            s(15)=s(19)*s(11) 
            s(16)=s(20)*s(12) 
            fr1=s(16)-s(15) 
            IF(abs(fr1).lt.1.d-10) THEN
               fr=(s(2)-px)/(a2(jp+1)-px) 
               s(3)=s(11)*(1.d0-fr)+s(12)*fr 
               fr=(s(3)-s(11))/(s(12)-s(11)) 
            ELSE
               fr=s(11)+s(12)-(a2(jp+1)-px)*s(15)*s(16)/fr1 
               fr1=sqrt(fr*fr-4.d0*s(11)*s(12)*(s(16)*(1.d0-s(19)*(s(2)-px))-s(&
               15)*(1.d0-s(20)*(s(2)-a2(jp+1))))/fr1)
               s(3)=0.5d0*(fr+fr1) 
               IF(s(3).gt.s(11))s(3)=0.5d0*(fr-fr1) 
               IF(s(3).ge.s(12)) THEN
                  fr=(s(3)-s(11))/(s(12)-s(11)) 
               ELSE
                  fr=(s(2)-px)/(a2(jp+1)-px) 
                  s(3)=s(11)*(1.d0-fr)+s(12)*fr 
                  fr=(s(3)-s(11))/(s(12)-s(11)) 
               ENDIF
            ENDIF
         ELSE
            s(3)=s(11) 
            fr=0.d0 
         ENDIF
         s( 4)=s(13)+(s(14)-s(13))*fr
         s( 5)=s( 4)+ s( 2)*s( 3)
         s( 6)=s(17)+(s(18)-s(17))*fr
         s( 7)=s(19)+(s(20)-s(19))*fr
         s( 8)=s(21)+(s(22)-s(21))*fr
         s(24)=s(25)+(s(26)-s(25))*fr
         iunp(1)=ip 
         iunp(2)=jp 
         s(23)=unp 
         err=.false. 
         RETURN 
      ENDIF
!
!.....vapor phase.                                                        
!
   70 s1=.false. 
      IF(a2(jp+1).le.s(10)) GOTO 71 
      frc=s(10) 
      hfg1=a0(2,2)-a0(2,1)+pa*(a0(1,2)-a0(1,1))
      hfg2=b0(2,2)-b0(2,1)+pb*(b0(1,2)-b0(1,1))
      dpdt1=hfg1/(ta*(a0(1,2)-a0(1,1)))
      dpdt2=hfg2/(tb*(b0(1,2)-b0(1,1)))
      f1=a0(1,2)*(a0(3,2)-a0(4,2)*dpdt1)
      f2=b0(1,2)*(b0(3,2)-b0(4,2)*dpdt2)
      d1=f1*(tb-ta) 
      d2=f2*(tb-ta) 
      c0=a0(1,2)
      c1=d1 
      c2=3.d0*(b0(1,2)-a0(1,2))-d2-2.d0*d1
      c3=d2+d1-2.d0*(b0(1,2)-a0(1,2))
      s(12)=c0+fr*(c1+fr*(c2+fr*c3)) 
      s(14)=a0(2,2)+(b0(2,2)-a0(2,2))*fr
      s(18)=a0(3,2)+fr*tb/s(1)*(b0(3,2)-a0(3,2))
      s(20)=a0(4,2)+(s(10)-pa)/(pb-pa)*pb/s(10)*(b0(4,2)-a0(4,2))
      s(22)=a0(5,2)+(b0(5,2)-a0(5,2))*fr
      s(26)=a0(6,2)+(b0(6,2)-a0(6,2))*fr
      GOTO 72 
   71 IF(a41(jp+1).lt.a1(ip)) GOTO 73 
      frd=(s(1)-a41(jp+1))/(a1(ip+1)-a41(jp+1)) 
      DO i=1,6
         ac(i)=a4(i,2,jp+1)
      ENDDO
      GOTO 74 
!
   73 IF(ip.eq.nt) THEN
!
!........vapor phase, temperature greater than highest table temperature.    
!
         fr=a5(1,ip,jp+1)*(a2(jp+1)-a2(jp)) 
         s( 3)=a5(1,ip,jp)*fr/(fr+(a5(1,ip,jp)-a5(1,ip,jp+1))*(s(2)-a2(jp))) 
         fr=(s(3)-a5(1,ip,jp))/(a5(1,ip,jp+1)-a5(1,ip,jp)) 
         frc=fr*a5(1,ip,jp+1)/s(3) 
         s( 5)=a5(2,ip,jp)+(a5(2,ip,jp+1)-a5(2,ip,jp))*fr+s(2)*s(3) 
         s( 8)=a5(5,ip,jp)+(a5(5,ip,jp+1)-a5(5,ip,jp))*frc
         frd=s(1)/a1(nt) 
         s( 3)=s(3)*frd 
         s( 5)=s(5)+s(8)*frn 
         s( 4)=s(5)-s(2)*s(3) 
         s( 6)=(a5(3,ip,jp)+(a5(3,ip,jp+1)-a5(3,ip,jp))*frc)/frd 
         s( 7)= a5(4,ip,jp)+(a5(4,ip,jp+1)-a5(4,ip,jp))*fr 
         s(24)= a5(6,ip,jp)+(a5(6,ip,jp+1)-a5(6,ip,jp))*fr 
         cv=s(8)-a1(nt)*s(6)*s(6)*s(3)/s(7) 
         s(24)=s(24)+cv*log(frd**(s(8)/cv)) 
         iunp(1)=ip 
         iunp(2)=jp 
         s(23)=unp 
         err=.false. 
         RETURN 
      ENDIF
      frd=frn/(a1(ip+1)-a1(ip)) 
      s1=.true. 
      DO i=1,6
         ac(i)=a5(i,ip,jp+1)
      ENDDO
   74 frc=a2(jp+1) 
      c0=ac(1) 
      d1=ac(1)*ac(3)*(a1(ip+1)-a1(ip)) 
      d2=a5(3,ip+1,jp+1)*a5(1,ip+1,jp+1)*(a1(ip+1)-a1(ip)) 
      c1=d1 
      c2=3.d0*(a5(1,ip+1,jp+1)-ac(1))-d2-2.d0*d1 
      c3=d2+d1-2.d0*(a5(1,ip+1,jp+1)-ac(1)) 
      s(12)=c0+frd*(c1+frd*(c2+frd*c3)) 
      s(14)=ac(2)+(a5(2,ip+1,jp+1)-ac(2))*frd 
      s(18)=ac(3)+frd*a1(ip+1)/s(1)*(a5(3,ip+1,jp+1)-ac(3)) 
      s(20)=ac(4)+(a5(4,ip+1,jp+1)-ac(4))*frd 
      s(22)=ac(5)+(a5(5,ip+1,jp+1)-ac(5))*frd 
      s(26)=ac(6)+(a5(6,ip+1,jp+1)-ac(6))*frd 
   72 IF(s1) THEN
         c0=a5(1,ip,jp)
         d1=a5(1,ip,jp)*a5(3,ip,jp)*(a1(ip+1)-a1(ip)) 
         d2=a5(3,ip+1,jp)*a5(1,ip+1,jp)*(a1(ip+1)-a1(ip)) 
         c1=d1 
         c2=3.d0*(a5(1,ip+1,jp)-a5(1,ip,jp))-d2-2.d0*d1 
         c3=d2+d1-2.d0*(a5(1,ip+1,jp)-a5(1,ip,jp)) 
         s(11)=c0+frd*(c1+frd*(c2+frd*c3)) 
         s(13)=a5(2,ip,jp)+(a5(2,ip+1,jp)-a5(2,ip,jp))*frd 
         s(17)=a5(3,ip,jp)+frd*a1(ip+1)/s(1)*(a5(3,ip+1,jp)-a5(3,ip,jp))
         s(19)=a5(4,ip,jp)+(a5(4,ip+1,jp)-a5(4,ip,jp))*frd 
         s(21)=a5(5,ip,jp)+(a5(5,ip+1,jp)-a5(5,ip,jp))*frd 
         s(25)=a5(6,ip,jp)+(a5(6,ip+1,jp)-a5(6,ip,jp))*frd 
      ELSE
         IF(a41(jp).lt.a1(ip)) THEN
            frd=frn/(a1(ip+1)-a1(ip)) 
            c0=a5(1,ip,jp)
            d1=a5(1,ip,jp)*a5(3,ip,jp)*(a1(ip+1)-a1(ip)) 
            d2=a5(3,ip+1,jp)*a5(1,ip+1,jp)*(a1(ip+1)-a1(ip)) 
            c1=d1 
            c2=3.d0*(a5(1,ip+1,jp)-a5(1,ip,jp))-d2-2.d0*d1 
            c3=d2+d1-2.d0*(a5(1,ip+1,jp)-a5(1,ip,jp)) 
            s(11)=c0+frd*(c1+frd*(c2+frd*c3)) 
            s(13)=a5(2,ip,jp)+(a5(2,ip+1,jp)-a5(2,ip,jp))*frd 
            s(17)=a5(3,ip,jp)+frd*a1(ip+1)/s(1)*(a5(3,ip+1,jp)-a5(3,ip,jp))
            s(19)=a5(4,ip,jp)+(a5(4,ip+1,jp)-a5(4,ip,jp))*frd 
            s(21)=a5(5,ip,jp)+(a5(5,ip+1,jp)-a5(5,ip,jp))*frd 
            s(25)=a5(6,ip,jp)+(a5(6,ip+1,jp)-a5(6,ip,jp))*frd 
         ELSE
            frd=(s(1)-a41(jp))/(a1(ip+1)-a41(jp)) 
            c0=a4(1,2,jp)
            d1=a4(1,2,jp)*a4(3,2,jp)*(a1(ip+1)-a1(ip)) 
            d2=a5(3,ip+1,jp)*a5(1,ip+1,jp)*(a1(ip+1)-a1(ip)) 
            c1=d1 
            c2=3.d0*(a5(1,ip+1,jp)-a4(1,2,jp))-d2-2.d0*d1 
            c3=d2+d1-2.d0*(a5(1,ip+1,jp)-a4(1,2,jp))
            s(11)=c0+frd*(c1+frd*(c2+frd*c3)) 
            s(13)=a4(2,2,jp)+(a5(2,ip+1,jp)-a4(2,2,jp))*frd 
            s(17)=a4(3,2,jp)+frd*a1(ip+1)/s(1)*(a5(3,ip+1,jp)-a4(3,2,jp))
            s(19)=a4(4,2,jp)+(a5(4,ip+1,jp)-a4(4,2,jp))*frd 
            s(21)=a4(5,2,jp)+(a5(5,ip+1,jp)-a4(5,2,jp))*frd 
            s(25)=a4(6,2,jp)+(a5(6,ip+1,jp)-a4(6,2,jp))*frd 
         ENDIF
      ENDIF
      fr=s(12)*(frc-a2(jp)) 
      s(3)=s(11)*fr/(fr+(s(11)-s(12))*(s(2)-a2(jp))) 
      fr=(s(3)-s(11))/(s(12)-s(11)) 
      frn=fr*s(12)/s(3) 
      s(4)=s(13)+(s(14)-s(13))*fr 
      s(5)=s(4)+s(2)*s(3) 
      s(6)=s(17)+frn*(s(18)-s(17)) 
      s(7)=s(19)+(s(20)-s(19))*fr 
      s(8)=s(21)+frn*(s(22)-s(21)) 
      s(24)=s(25)+(s(26)-s(25))*fr 
!
      iunp(1)=ip 
      iunp(2)=jp 
      s(23)=unp 
      err=.false. 
      RETURN 
!
!.....vapor phase, pressure less than lowest table pressure               
!
   90 IF(it.eq.1) GOTO 1001
      IF(s(1).lt.a41(1)) THEN
         s( 3)=(fr*pb*b0(1,2)/tb+(1.d0-fr)*pa*a0(1,2)/ta)*s(1)/s(2)
         s( 4)=a0(2,2)+(b0(2,2)-a0(2,2))*fr
         s( 6)=1.d0/s(1)
         s( 8)=a0(5,2)+(b0(5,2)-a0(5,2))*fr
         s(24)=a0(6,2)+(b0(6,2)-a0(6,2))*fr
         ren=s(2)*s(3)/s(1)
         s(24)=s(24)-ren*log(s(2)/s(10))
         s( 5)=s(4)+s(2)*s(3)
         s( 7)=1.d0/s(2)
         jp=1
         iunp(1)=ip
         iunp(2)=jp
         s(23)=unp
         err=.false.
         RETURN
      ENDIF
!
      IF(ip.eq.nt) THEN
         frd=s(1)/a1(nt)
         frc=a2(1)*a5(1,ip,1)
         s( 3)=frc*frd/s(2)
         s( 8)=a5(5,ip,1)
         s( 5)=a5(2,ip,1)+frc+s(8)*(s(1)-a1(nt))
         s( 4)=s(5)-s(2)*s(3)
         s( 6)=a5(3,ip,1)/frd
         s( 7)=a5(4,ip,1)*a2(1)/s(2)
         ren=s(2)*s(3)/s(1)
         s(24)=a5(6,ip,1)+s(8)*log(s(1)/a1(nt))-ren*log(s(2)/a2(1))
         iunp(1)=ip
         iunp(2)=jp
         s(23)=unp
         err=.false.
         RETURN
      ENDIF
!
      IF(a1(ip).lt.a41(1)) THEN
         fr=(s(1)-a41(1))/(a1(ip+1)-a41(1))
         s( 3)=(fr*a5(1,ip+1,1)/a1(ip+1)+(1.d0-fr)*a4(1,2,1)/a41(1))*a2(1)*s(1)/s(2)
         s( 4)=a4(2,2,1)+(a5(2,ip+1,1)-a4(2,2,1))*fr
         s( 6)=a4(3,2,1)+(a5(3,ip+1,1)-a4(3,2,1))*fr*a1(ip+1)/s(1)
         s( 8)=a4(5,2,1)+(a5(5,ip+1,1)-a4(5,2,1))*fr
         s(24)=a4(6,2,1)+(a5(6,ip+1,1)-a4(6,2,1))*fr
      ELSE
         fr=(s(1)-a1(ip))/(a1(ip+1)-a1(ip))
         s( 3)=(fr*a5(1,ip+1,1)/a1(ip+1)+(1.d0-fr)*a5(1,ip,1)/a1(ip))*a2(1)*s(1)/s(2)
         s( 4)=a5(2,ip,1)+(a5(2,ip+1,1)-a5(2,ip,1))*fr
         s( 6)=a5(3,ip,1)+(a5(3,ip+1,1)-a5(3,ip,1))*fr*a1(ip+1)/s(1)
         s( 8)=a5(5,ip,1)+(a5(5,ip+1,1)-a5(5,ip,1))*fr
         s(24)=a5(6,ip,1)+(a5(6,ip+1,1)-a5(6,ip,1))*fr
      ENDIF
      ren=s(2)*s(3)/s(1)
      s(24)=s(24)-ren*log(s(2)/a2(1))
      s( 5)=s(4)+s(2)*s(3)
      s( 7)=1.d0/s(2)
      jp=1
      iunp(1)=ip
      iunp(2)=jp
      s(23)=unp
      err=.false.
      RETURN
 1001 err=.true. 
!
      END SUBROUTINE sth2x3_cupid                         
