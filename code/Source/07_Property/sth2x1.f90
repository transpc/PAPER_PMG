!
      SUBROUTINE sth2x1_cupid(a1,a2,s,err)
!
!     compute water thermodynamic properties as a function of temperature and quality
!
      USE STM_TBL_cupid  , ONLY: pxxx,pxxy,pxx1,pxx2, &
                                 crt,crp,             &
                                 b,c,cc,k,            &
                                 a31,a3,a41,a4,       &
                                 nt,np,ns,ns2
!
      IMPLICIT NONE
!
!.....Input
      REAL(8) :: a1(nt)
      REAL(8) :: a2(np)
!.....Output
      LOGICAL :: err
      REAL(8) :: s(26) 
!.....Local variables
      INTEGER :: i
      INTEGER :: ip,jp
      INTEGER :: iunp(2) 
      LOGICAL :: s1,s2,s3 
      REAL(8) :: plow=611.2444d0
      REAL(8) :: unp
      REAL(8) :: pa,pb,ta,tb,hfg1,hfg2,dpdt1,dpdt2
      REAL(8) :: c0,c1,c2,c3,d1,d2,f1,f2,fr,fr1
      REAL(8) :: am1,am2
      REAL(8) :: bm1,bm2
      REAL(8) :: a0(6,2),b0(6,2)
      REAL(8) :: bma(6,2)
!
!     INCLUDE 'machds.h' 
!                                                                       
      EQUIVALENCE(unp,iunp(1)) 
!                                                                       
!.....check for valid input                                               
!
      s1=.false. 
!
!.....Bug ip cannot reach ns, a31(ns+1)out of bounds
!
!  15 IF(s(1).lt.a1(1).or.s(1).gt.a1(ns)) GOTO 101
   15 IF(s(1).lt.a1(1).or.s(1).ge.a1(ns)) GOTO 101 
      IF(s(9).lt.0.0d0.or.s(9).gt.1.0d0) GOTO 101 
      IF(s1) GOTO 16 
      fr=s(1)/crt 
      fr1=1.0d0-fr 
      s(10)=crp*exp((((((k(5)*fr1+k(4))*fr1+k(3))*fr1+k(2))*fr1+k(1))* &
      fr1)/(((k(7)*fr1+k(6))*fr1+1.0d0)*fr)-fr1/(k(8)*fr1*fr1+k(9)))    
      s(2)=s(10) 
      ENTRY sth2xb_cupid(a1,a2,s,err) 
!
   16 unp=s(23) 
      ip=iunp(1) 
      jp=iunp(2) 
      s2=.false. 
      s3=.false. 
      IF(ip.le.0.or.ip.ge.ns)  ip=1 
      IF(jp.le.0.or.jp.ge.ns2) jp=1 
!
!.....set indexes in temperature and pressure tables for saturation computations                                                        
!
   11 IF(s(1).ge.a1(ip)) GOTO 10 
      ip=ip-1 
      GOTO 11 
   10 IF(s(1).lt.a1(ip+1)) GOTO 12 
      ip=ip+1 
      GOTO 10 
!
   12 CONTINUE
  111 IF(s(10).ge.a2(jp)) GOTO 110 
      jp=jp-1 
      IF(jp.gt.0) GOTO 111 
      jp=1 
      s3=.true. 
      GOTO 112 
  110 IF(s(10).lt.a2(jp+1)) GOTO 112 
      jp=jp+1 
      IF(jp.lt.np) GOTO 110 
      s2=.true. 
!
  112 CONTINUE
      IF(s3.or.a2(jp).le.a31(ip)) THEN
         ta=a1(ip) 
         pa=a31(ip) 
         DO i=1,6
            a0(i,1)=a3(i,1,ip)
            a0(i,2)=a3(i,2,ip)
         ENDDO
      ELSE
         pa=a2(jp) 
         ta=a41(jp) 
         DO i=1,6
            a0(i,1)=a4(i,1,jp)
            a0(i,2)=a4(i,2,jp)
         ENDDO
      ENDIF
!
!.....s2 must be checked first to avoid jp=np and a2(np+1) out of bounds
!
      IF(.not.s2) THEN
         IF(a2(jp+1).ge.a31(ip+1)) THEN
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
      ELSE
         tb=a1(ip+1) 
         pb=a31(ip+1) 
         DO i=1,6
            b0(i,1)=a3(i,1,ip+1)
            b0(i,2)=a3(i,2,ip+1)
         ENDDO
      ENDIF
      fr1=s(1)-ta 
      fr=fr1/(tb-ta) 
!
!.....two phase fluid.                                                    
!
      am1=a0(1,2)-a0(1,1)
      am2=a0(2,2)-a0(2,1)
      bm1=b0(1,2)-b0(1,1)
      bm2=b0(2,2)-b0(2,1)
!
      DO i=1,6
         bma(i,1)=b0(i,1)-a0(i,1)
         bma(i,2)=b0(i,2)-a0(i,2)
      ENDDO
!
      hfg1=am2+pa*am1
      hfg2=bm2+pb*bm1
      dpdt1=hfg1/(ta*am1)
      dpdt2=hfg2/(tb*bm1)
!
      f1=a0(1,1)*(a0(3,1)-a0(4,1)*dpdt1)
      f2=b0(1,1)*(b0(3,1)-b0(4,1)*dpdt2)
      d1=f1*(tb-ta) 
      d2=f2*(tb-ta) 
      c0=a0(1,1)
      c1=d1 
      c2= 3.d0*bma(1,1)-d2-2.d0*d1
      c3=d2+     d1-2.d0*bma(1,1)
      s(11)=c0+fr*(c1+fr*(c2+fr*c3)) 
      f1=a0(1,2)*(a0(3,2)-a0(4,2)*dpdt1)
      f2=b0(1,2)*(b0(3,2)-b0(4,2)*dpdt2)
      d1=f1*(tb-ta) 
      d2=f2*(tb-ta) 
      c0=a0(1,2)
      c1=d1 
      c2= 3.d0*bma(1,2)-d2-2.d0*d1
      c3=d2+     d1-2.d0*bma(1,2)
      s(12)=c0+fr*(c1+fr*(c2+fr*c3)) 
!
!.....two phase fluid.                                                    
!
      s(13)=a0(2,1)+bma(2,1)*fr
      s(17)=a0(3,1)+bma(3,1)*fr
      s(19)=a0(4,1)+bma(4,1)*fr
      s(21)=a0(5,1)+bma(5,1)*fr
      s(25)=a0(6,1)+bma(6,1)*fr
      s(14)=a0(2,2)+bma(2,2)*fr
      s(18)=a0(3,2)+bma(3,2)*(fr*tb/s(1))
      s(20)=a0(4,2)+bma(4,2)*((s(10)-pa)/(pb-pa)*pb/s(10))
      s(22)=a0(5,2)+bma(5,2)*fr
      s(26)=a0(6,2)+bma(6,2)*fr
!
      s(15)=s(13)+s(10)*s(11)
      s(16)=s(14)+s(10)*s(12)
      fr=1.0d0-s(9)
      s( 3)=fr*s(11)+s(9)*s(12)
      s( 4)=fr*s(13)+s(9)*s(14)
      s( 5)=fr*s(15)+s(9)*s(16)
      s(24)=fr*s(25)+s(9)*s(26)
!
      iunp(1)=ip 
      iunp(2)=jp 
      s(23)=unp 
      err=.false. 
      RETURN 
      ENTRY sth2x2_cupid(a1,a2,s,err) 
!
!.....compute water thermodynamic properties as a function of pressure and quality                                                         
!
      s1=.true. 
!
!.....temporary patch to be able to DO ice condenser debug runs             
!
      s(2)=max(s(2),plow) 
      s(10)=s(2) 
      IF(s(2).lt.plow.or.s(2).gt.crp) GOTO 101 
      IF(s(2).lt.pxxx) THEN
         fr=log(s(2)) 
         s(1)=(cc(1)*fr+cc(2))*fr+cc(3) 
      ELSEIF(s(2).gt.pxxy) THEN
         fr=log(pxx2*s(2)) 
         s(1)=((((fr*b(6)+b(5))*fr+b(4))*fr+b(3))*fr+b(2))*fr+b(1) 
         s(1)=min(s(1),crt) 
      ELSE
         fr=log(pxx1*s(2)) 
         s(1)=(((((((fr*c(9)+c(8))*fr+c(7))*fr+c(6))*fr+c(5))*fr+c(4))*fr+ &
         c(3))*fr+c(2))*fr+c(1)
      ENDIF
      fr=s(1)/crt 
      fr1=1.0d0-fr 
      d1=((((k(5)*fr1+k(4))*fr1+k(3))*fr1+k(2))*fr1+k(1))*fr1 
      d2=(((5.0d0*k(5)*fr1+4.0d0*k(4))*fr1+3.0d0*k(3))*fr1+2.0d0*k(2))* &
      fr1+k(1)                                                          
      c2=k(7)*fr1 
      c1=(c2+k(6))*fr1+1.0d0 
      c2=2.0d0*c2+k(6) 
      f2=k(8)*fr1 
      f1=1.0d0/(f2*fr1+k(9)) 
      f2=2.0d0*f2 
      hfg1=1.0d0/(fr*c1) 
      hfg2=fr1*f1 
      pa=crp*exp(d1*hfg1-hfg2) 
      s(1)=max(s(1)+(s(2)-pa)*crt/(pa*((d1*hfg1*(fr*c2-c1)-d2)*hfg1+(   &
      1.0d0-hfg2*f2)*f1)),273.16d0)                                     
      GOTO 15 
  101 err=.true. 
      RETURN
!     
      END SUBROUTINE sth2x1_cupid                         
