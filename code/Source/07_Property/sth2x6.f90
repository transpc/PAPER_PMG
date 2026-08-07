!
      SUBROUTINE sth2x6_cupid(s3,s,it,err, &
                              a1,a2,a5)
!
!     Compute water thermodynamic properties as a function of pressure and internal energy.
      USE STM_TBL_cupid  , ONLY: pxxx,pxxy,pxx1,pxx2, &
                                 crt,crp,             &
                                 b,c,cc,k,g,          &
                                 a31,a3,a41,a4,       &
                                 nt,np,ns
!
      IMPLICIT NONE 
!                                                                       
!.....Input
      LOGICAL :: s3
      REAL(8) :: a1(nt)
      REAL(8) :: a2(np)
      REAL(8) :: a5(6,nt,np)
!.....Output
      INTEGER :: it
      LOGICAL :: err
      REAL(8) :: s(26)
!.....Local variables
      INTEGER :: i
      INTEGER :: ic,ip,jp
      INTEGER :: temp,pres,vbar,ubar,hbar,beta,kapa,cp,qual,tsat,vsubf,         &
                 vsubg,usubf,usubg,hsubf,hsubg,betaf,betag,kapaf,kapag,cpf,cpg, &
                 is23,entpy,entpyf,entpyg
      PARAMETER (temp  = 1,  &
                 pres  = 2,  &
                 vbar  = 3,  &
                 ubar  = 4,  &
                 hbar  = 5,  &
                 beta  = 6,  &
                 kapa  = 7,  &
                 cp    = 8,  &
                 qual  = 9,  &
                 tsat  =10,  &
                 vsubf =11,  &
                 vsubg =12,  &
                 usubf =13,  &
                 usubg =14,  &
                 hsubf =15,  &
                 hsubg =16,  &
                 betaf =17,  &
                 betag =18,  &
                 kapaf =19,  &
                 kapag =20,  &
                 cpf   =21,  &
                 cpg   =22,  &
                 is23  =23,  &
                 entpy =24,  &
                 entpyf=25,  &
                 entpyg=26)
      INTEGER iunp(2) 
      EQUIVALENCE(unp,iunp(1)) 
!
      LOGICAL :: oncel,oncev,s1,s2
!
      REAL(8) :: am1,am2
      REAL(8) :: bm1,bm2
      REAL(8) :: a0(6,2),b0(6,2)
      REAL(8) :: bma(6,2)
      REAL(8) :: ac(6),ad(6)
      REAL(8) :: unp 
      REAL(8) :: c0,c1,c2,c3,cv,d1,d2,dpdt1,dpdt2,                           &
                 f1,f2,fr,fr1,fr2,frc,frc2,frc3,frd,frd2,frd3,frn,hfg1,hfg2, &
                 ta0,pa0,pa,pb,pr,px,r2,ren,s11i,s12i,ta,tb,tr,ut,vt   
!                                                                       
      IF(s3) THEN
         vt=1.d0 !apr1400_lbloca_debug
         ut=1.d0 !apr1400_lbloca_debug
      ENDIF
!
!.....Check for valid input.
!
      IF(s(pres).le.0.d0.or.s(pres).gt.a2(np)) GOTO 1001 
      oncel=.false. 
      oncev=.false. 
!
      unp=s(is23) 
      ip=iunp(1) 
      jp=iunp(2) 
      IF(ip.le.0.or.ip.ge.nt) ip=1 
      IF(jp.le.0.or.jp.ge.np) jp=1 
      s1=.false. 
!
!.....Set indexes in temperature and pressure tables for saturation computations.
!
   11 IF(s(pres).lt.a2(jp))THEN 
         jp=jp-1 
         IF(jp.gt.0) GOTO 11 
         jp=1 
         s1=.true. 
         IF(s(pres).lt.a31(1)) THEN
           it=4
           GOTO 50
         ENDIF
      ELSE 
   10    IF(s(pres).ge.a2(jp+1))THEN 
            jp=jp+1 
            GOTO 10 
         ENDIF 
      ENDIF 
!
      IF(s(pres).ge.crp) THEN
         it=4
         GOTO 50
      ENDIF
!
      IF(.not.s3)THEN 
         IF(s(pres).lt.pxxx)THEN 
            fr=log(s(pres)) 
            s(tsat)=(cc(1)*fr+cc(2))*fr+cc(3) 
         ELSEIF(s(pres).gt.pxxy)THEN 
            fr=log(pxx2*s(pres)) 
            s(tsat)=((((fr*b(6)+b(5))*fr+b(4))*fr+b(3))*fr+b(2))*fr+b(1) 
         ELSE 
            fr=log(pxx1*s(pres)) 
            s(tsat)=(((((((fr*c(9)+c(8))*fr+c(7))*fr+c(6))*fr+c(5))*fr+ &
            c(4))*fr+c(3))*fr+c(2))*fr+c(1)                             
         ENDIF 
         fr=s(tsat)/crt 
         fr1=1.d0-fr 
         d1=((((k(5)*fr1+k(4))*fr1+k(3))*fr1+k(2))*fr1+k(1))*fr1 
         d2=(((g(5)*fr1+g(4))*fr1+g(3))*fr1+g(2))*fr1+g(1)
         c2=k(7)*fr1 
         c1=(c2+k(6))*fr1+1.d0 
         c2=2.d0*c2+k(6) 
         f2=k(8)*fr1 
         f1=1.d0/(f2*fr1+k(9)) 
         f2=2.d0*f2 
         hfg1=1.d0/(fr*c1) 
         hfg2=fr1*f1 
         pa0=crp*exp(d1*hfg1-hfg2) 
         s(tsat)=max(s(tsat)+(s(pres)-pa0)*crt/(pa0*((d1*hfg1*(fr*c2-c1)- &
         d2)*hfg1+(1.d0-hfg2*f2)*f1)),273.16d0)                        
      ENDIF 
!
      ic=ip 
   16 IF(s(tsat).lt.a1(ic))THEN 
         ic=ic-1 
         IF(ic.gt.0) GOTO 16 
         ic=1 
      ELSE 
   13    IF(s(tsat).ge.a1(ic+1))THEN 
            ic=ic+1 
            IF(ic.lt.ns) GOTO 13 
            ic=ic-1 
         ENDIF 
      ENDIF 
      DO i=1,6
         ac(i)=a1(ic+i-1)
      ENDDO
!
      IF(s1.or.a2(jp).le.a31(ic)) THEN
         ta=a1(ic) 
         pa=a31(ic)
         DO i=1,6
            a0(i,1)=a3(i,1,ic)
            a0(i,2)=a3(i,2,ic)
         ENDDO
      ELSE 
         pa=a2(jp) 
         ta=a41(jp)
         DO i=1,6
            a0(i,1)=a4(i,1,jp)
            a0(i,2)=a4(i,2,jp)
         ENDDO
      ENDIF 
      IF(a2(jp+1).ge.a31(ic+1)) THEN 
         pb=a31(ic+1)
         tb=a1(ic+1) 
         DO i=1,6
            b0(i,1)=a3(i,1,ic+1)
            b0(i,2)=a3(i,2,ic+1)
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
!.....Compute vsubf and vsubg to determine liquid, two phase, or vapor state.
!
      fr1=s(tsat)-ta 
      fr=fr1/(tb-ta) 
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
      c1=f1*(tb-ta) 
      d2=f2*(tb-ta) 
      c0=a0(1,1)
      c2=3.d0*bma(1,1)-d2-2.d0*c1
      c3=d2+c1-2.d0*bma(1,1)
!
      s(vsubf)=c0+fr*(c1+fr*(c2+fr*c3))
      s(usubf)=a0(2,1)+bma(2,1)*fr
      IF(s(ubar).le.s(usubf)) THEN
         it=1
         IF(s1) GOTO 1001
         GOTO 50
      ENDIF
      f1=a0(1,2)*(a0(3,2)-a0(4,2)*dpdt1)
      f2=b0(1,2)*(b0(3,2)-b0(4,2)*dpdt2)
      c1=f1*(tb-ta)
      d2=f2*(tb-ta)
      c0=a0(1,2)
      c2=3.d0*bma(1,2)-d2-2.d0*c1
      c3=d2+c1-2.d0*bma(1,2)
!
      s(vsubg)=c0+fr*(c1+fr*(c2+fr*c3))
      s(usubg)=a0(2,2)+bma(2,2)*fr
      IF(s(ubar).ge.s(usubg)) THEN
         it=3
         GOTO 50
      ENDIF
!
!.....Two phase fluid.                                                     
!
      it=2
      s(hsubf)=s(usubf)+s(pres)*s(vsubf)
      s(hsubg)=s(usubg)+s(pres)*s(vsubg)
!
      s(betaf) =a0(3,1)+bma(3,1)*fr
      s(kapaf) =a0(4,1)+bma(4,1)*fr
      s(cpf)   =a0(5,1)+bma(5,1)*fr
      s(entpyf)=a0(6,1)+bma(6,1)*fr
      s(betag) =a0(3,2)+bma(3,2)*(fr*tb/s(tsat))
      s(kapag) =a0(4,2)+bma(4,2)*((s(pres)-pa)/(pb-pa)*pb/s(pres))
      s(cpg)   =a0(5,2)+bma(5,2)*fr
      s(entpyg)=a0(6,2)+bma(6,2)*fr
!
      s(qual) =(s(ubar)-s(usubf))/(s(usubg)-s(usubf))
      fr=1.d0-s(qual)
      s(temp) =s(tsat)
      s(vbar) =fr*s(vsubf)+s(qual)*s(vsubg)
      s(hbar) =fr*s(hsubf)+s(qual)*s(hsubg)
      s(entpy)=fr*s(entpyf)+s(qual)*s(entpyg)
      ip=ic
!
      iunp(1)=ip
      iunp(2)=jp
      s(is23)=unp
      err=.false.
      RETURN
!
!.....Single phase fluid, search for single phase indexes.                 
!
   50 CONTINUE
   51 IF(s(ubar).ge.a5(2,ip,jp)) GOTO 52 
      ip=ip-1 
      IF(ip.gt.0) GOTO 51 
      ip=ip+1 
      GOTO 54 
   52 IF(s(ubar).le.a5(2,ip+1,jp)) GOTO 54 
      ip=ip+1 
      IF(ip.lt.nt) GOTO 52 
      IF(s1) THEN
         s(cp)   =a5(5,ip,jp)
         frd=s(cp)-a2(1)*a5(1,ip,jp)/a1(nt) 
         s(temp) =(s(ubar)-a5(2,ip,jp)+frd*a1(nt))/frd 
         frd=s(temp)/a1(nt) 
         frc=a2(1)*a5(1,ip,jp)
         s(vbar) =frc*frd/s(pres) 
         s(beta) =a5(3,ip,jp)/frd 
         s(kapa) =a5(4,ip,jp)*a2(1)/s(pres) 
         ren=s(pres)*s(vbar)/s(temp) 
         IF(s(temp).le.0.d0) GOTO 1001 
         s(entpy)=a5(6,ip,jp)+s(cp)*log(s(temp)/a1(nt))-ren*log(s(pres)/a2(1))
         s(hbar) =s(ubar)+s(pres)*s(vbar) 
         s(qual) =1.d0 
         iunp(1)=ip 
         iunp(2)=jp 
         s(is23)=unp 
         err=.false. 
         RETURN 
      ENDIF
      ip=ip-1 
      GOTO 53 
   54 IF(s1) THEN
!
!........Vapor phase, pressure less than lowest table pressure.
!
         ut=a4(2,2,1)
         IF(s(ubar).lt.ut) THEN
            ip=1
            IF(s(ubar).lt.a3(2,2,ip))THEN
               ren=a3(1,2,ip)*a31(ip)/a1(ip)
               s(cp)   =a3(5,2,ip)
               s(temp) =a1(ip)+(s(ubar)-a3(2,2,ip))/(a3(5,2,ip)-ren)
               fr=ren*s(temp)
               s(hbar) =s(ubar)+fr
               s(beta) =1.d0/s(temp)
               s(kapa) =1.d0/s(pres)
               s(vbar) =fr*s(kapa)
               IF(s(temp).le.0.d0) GOTO 1001
               s(entpy)=a3(6,2,ip)+s(cp)*log(s(temp)/a1(ip))-ren*log(s(pres)/a31(ip))
            ELSE
  202          IF(s(ubar).gt.a3(2,2,ip+1))THEN
                  ip=ip+1
                  ut=a3(2,2,ip+1)
                  GOTO 202
               ENDIF
               fr=(s(ubar)-a3(2,2,ip))/(a3(2,2,ip+1)-a3(2,2,ip))
               s(temp)=a1(ip)+fr*(a1(ip+1)-a1(ip))
               s(vbar)=(fr*a31(ip+1)*a3(1,2,ip+1)/a1(ip+1)+(1.d0-fr)*a31(ip)*a3(1,2,ip)/&
               a1(ip))*s(temp)/s(pres)
               s(beta)=1.d0/s(temp)
               s(cp)=a3(5,2,ip)+(a3(5,2,ip+1)-a3(5,2,ip))*fr
               ren=s(pres)*s(vbar)/s(temp)
               IF(s(temp).le.0.d0) GOTO 1001
               s(entpy)=a3(6,2,ip)+s(cp)*log(s(temp)/a1(ip))-ren*log(s(pres)/a31(ip))
               s(kapa)=1.d0/s(pres)
               s(hbar)=s(ubar)+s(pres)*s(vbar)
            ENDIF
         ELSE
            IF(a1(ip).ge.a41(1))THEN
               ut=a5(2,ip,jp)
               fr=(s(ubar)-ut)/(a5(2,ip+1,jp)-ut)
               fr1=a1(ip+1)-a1(ip)
               fr2=fr1*fr
               s(temp) =a1(ip)+fr2
               fr1=fr2/fr1
               s(vbar) =(fr1*a5(1,ip+1,jp)/a1(ip+1)+(1.d0-fr1)*a5(1,ip,jp)/a1(ip))*a2(1)*s(temp)/s(pres)
               s(beta) =a5(3,ip,jp)+(a5(3,ip+1,jp)-a5(3,ip,jp))*fr*a1(ip+1)/s(temp)
               s(cp)   =a5(5,ip,jp)+(a5(5,ip+1,jp)-a5(5,ip,jp))*fr
               s(entpy)=a5(6,ip,jp)+(a5(6,ip+1,jp)-a5(6,ip,jp))*fr
            ELSE 
               fr=(s(ubar)-ut)/(a5(2,ip+1,jp)-ut)
               fr1=a1(ip+1)-a41(1)
               fr2=fr1*fr
               s(temp) =a41(1)+fr2
               fr1=fr2/fr1
               s(vbar) =(fr1*a5(1,ip+1,jp)/a1(ip+1)+(1.d0-fr1)*a4(1,2,1)/a41(1))*a2(1)*s(temp)/s(pres)
               s(beta) =a4(3,2,1)+(a5(3,ip+1,jp)-a4(3,2,1))*fr*a1(ip+1)/s(temp)
               s(cp)   =a4(5,2,1)+(a5(5,ip+1,jp)-a4(5,2,1))*fr
               s(entpy)=a4(6,2,1)+(a5(6,ip+1,jp)-a4(6,2,1))*fr
            ENDIF
            ren=s(pres)*s(vbar)/s(temp)
            s(entpy)=s(entpy)-ren*log(s(pres)/a2(1))
            s(kapa) =1.d0/s(pres)
            s(hbar) =s(ubar)+s(pres)*s(vbar)
         ENDIF
         s(qual)=1.d0
         iunp(1)=ip
         iunp(2)=jp
         s(is23)=unp
         err=.false.
         RETURN
      ENDIF
   53 CONTINUE
      frn=s(pres)-a2(jp) 
      frc2=s(pres)-a2(jp+1) 
      s3=.false. 
      IF(it.eq.3) THEN
         GOTO 70
      ELSEIF(it.gt.3) THEN
         IF(ip.ge.ns) GOTO 157
      ELSEIF(it.lt.3) THEN
!
!........Liquid phase.                                                        
!
         ut=s(usubf) 
         vt=s(vsubf) 
   60    IF(a1(ip).lt.s(tsat)) GOTO 57 
         ip=ip-1 
         IF(ip.le.0) GOTO 1001 
         GOTO 60 
      ENDIF
   57 s1=.false. 
      IF(it.eq.4  .or. &
         a1(ip+1).le.s(tsat)) THEN
         ta0=a1(ip+1) 
         s2=.true. 
         IF(ip+1.gt.ns .or. &
            a31(ip+1).le.a2(jp)) THEN
            frc3=frn
            frc=a2(jp+1)-a2(jp)
            px=a2(jp)
            DO i=1,6
               ac(i)=a5(i,ip+1,jp)
            ENDDO
            frd3=frc3
            frd=frc
            s1=.true.
         ELSE
            frc3=s(pres)-a31(ip+1)
            frc=a2(jp+1)-a31(ip+1)
            px=a31(ip+1)
            DO i=1,6
               ac(i)=a3(i,1,ip+1)
            ENDDO
         ENDIF
         pr=(s(pres)-px)/(a2(jp+1)-px) 
         c0=1.d0/ac(1) 
         r2=1.d0/a5(1,ip+1,jp+1)
         c1=c0*ac(4)*(a2(jp+1)-px)
         d2=r2*a5(4,ip+1,jp+1)*(a2(jp+1)-px)
         c2=3.d0*(r2-c0)-d2-2.d0*c1
         c3=d2+c1-2.d0*(r2-c0)
         s11i=c0+pr*(c1+pr*(c2+pr*c3))
         s(vsubf)=1.d0/s11i
         fr1=(s(vsubf)-ac(1))/(a5(1,ip+1,jp+1)-ac(1))
         s(usubf)=ac(2)+(a5(2,ip+1,jp+1)-ac(2))*fr1
      ELSE
         ta0=s(tsat)
         s2=.false.
      ENDIF
!
      IF(.not.s3) THEN
   62    IF(s1) THEN
            DO i=1,6
               ad(i)=a5(i,ip,jp)
            ENDDO
         ELSE
            IF(ip.gt.ns .or. &
               a31(ip).le.a2(jp)) THEN
               frd3=frn
               frd=a2(jp+1)-a2(jp)
               px=a2(jp)
               s1=.true.
               DO i=1,6
                  ad(i)=a5(i,ip,jp)
               ENDDO
            ELSE
               frd3=s(pres)-a31(ip)
               frd=a2(jp+1)-a31(ip)
               px=a31(ip)
               DO i=1,6
                  ad(i)=a3(i,1,ip)
               ENDDO
            ENDIF
         ENDIF
         s(beta) =ad(4)*ad(1)
         s(kapa) =a5(4,ip,jp+1)*a5(1,ip,jp+1)
         fr2=s(kapa)-s(beta) 
         pr=(s(pres)-px)/(a2(jp+1)-px) 
         c0=1.d0/ad(1)
         r2=1.d0/a5(1,ip,jp+1)
         c1=c0*ad(4)*(a2(jp+1)-px)
         d2=r2*a5(4,ip,jp+1)*(a2(jp+1)-px) 
         c2=3.d0*(r2-c0)-d2-2.d0*c1 
         c3=d2+c1-2.d0*(r2-c0) 
         s12i=c0+pr*(c1+pr*(c2+pr*c3)) 
         s(vsubg)=1.d0/s12i 
         fr2=(s(vsubg)-ad(1))/(a5(1,ip,jp+1)-ad(1))
         s(usubg)=ad(2)+(a5(2,ip,jp+1)-ad(2))*fr2
         IF(s(usubg).le.s(ubar)) GOTO 68 
         s2=.true. 
         s(vsubf)=s(vsubg)
         s(usubf)=s(usubg)
         fr1=fr2 
         ip=ip-1 
         ta0=a1(ip+1) 
         IF(ip.le.0) GOTO 1001 
         DO i=1,6
            ac(i)=ad(i)
         ENDDO
         GOTO 62 
      ENDIF
   68 IF(s(usubf).ge.s(ubar)) THEN
         IF(s2) THEN
            s(betaf) =ac(3)+(a5(3,ip+1,jp+1)-ac(3))*fr1
            s(kapaf) =ac(4)+(a5(4,ip+1,jp+1)-ac(4))*fr1
            s(cpf)   =ac(5)+(a5(5,ip+1,jp+1)-ac(5))*fr1
            s(entpyf)=ac(6)+(a5(6,ip+1,jp+1)-ac(6))*fr1
         ELSE
            s(betaf) =a0(3,1)+(b0(3,1)-a0(3,1))*fr
            s(kapaf) =a0(4,1)+(b0(4,1)-a0(4,1))*fr
            s(cpf)   =a0(5,1)+(b0(5,1)-a0(5,1))*fr
            s(entpyf)=a0(6,1)+(b0(6,1)-a0(6,1))*fr
         ENDIF
         s(betag) =ad(3)+(a5(3,ip,jp+1)-ad(3))*fr2
         s(kapag) =ad(4)+(a5(4,ip,jp+1)-ad(4))*fr2
         s(cpg)   =ad(5)+(a5(5,ip,jp+1)-ad(5))*fr2
         s(entpyg)=ad(6)+(a5(6,ip,jp+1)-ad(6))*fr2
         fr=(s(ubar)-s(usubg))/(s(usubf)-s(usubg)) 
         frd=ta0-a1(ip) 
         fr2=frd*fr 
         s(temp) =a1(ip)+fr2 
         fr1=fr2/frd 
         tr=(s(temp)-a1(ip))/frd 
         c0=s(vsubg) 
         c1=s(vsubg)*s(betag)*frd 
         d2=s(vsubf)*s(betaf)*frd 
         c2=3.d0*(s(vsubf)-s(vsubg))-d2-2.d0*c1 
         c3=d2+c1-2.d0*(s(vsubf)-s(vsubg)) 
         s(vbar) =c0+tr*(c1+tr*(c2+tr*c3))
         s(hbar) =s(ubar)+s(pres)*s(vbar)
         s(beta) =s(betag)+(s(betaf)-s(betag))*fr1
         s(kapa) =s(kapag)+(s(kapaf)-s(kapag))*fr1
         s(cp)   =s(cpg)+(s(cpf)-s(cpg))*fr1
         s(entpy)=s(entpyg)+(s(entpyf)-s(entpyg))*fr1
         s(qual) =0.d0
         iunp(1)=ip
         iunp(2)=jp
         s(is23)=unp
         err=.false.
         RETURN
      ELSE
         IF(oncev) GOTO 1002 
         oncel=.true. 
         s(vsubg)=s(vsubf) 
         s(vsubf)=vt      
         s(usubg)=s(usubf) 
         s(usubf)=ut 
         fr2=fr1 
         ip=ip+1 
         IF(ip.lt.ns) THEN
            DO i=1,6
               ad(i)=ac(i)
            ENDDO
            s3=.true.
            GOTO 57
         ELSE
            s3=.false.
            GOTO 157
         ENDIF
      ENDIF
!
   70 ut=s(usubg)
      vt=s(vsubg)
  160 IF(a1(ip+1).le.s(tsat))THEN
!
!........Vapor phase.
!
         ip=ip+1 
         GOTO 160 
      ENDIF 
  157 s1=.false. 
      IF(it.ne.4)THEN 
         IF(a1(ip).lt.s(tsat))THEN 
            ta0=s(tsat) 
            s2=.false. 
            IF(s3) GOTO 168
            GOTO 162
         ENDIF 
      ENDIF 
      ta0=a1(ip) 
      s2=.true. 
      IF(ip.lt.ns)THEN 
         IF(a31(ip).lt.a2(jp+1))THEN 
            frc3=a31(ip)-a2(jp) 
            px=a31(ip)
            DO i=1,6
               ac(i)=a3(i,2,ip)
            ENDDO
            GOTO 164 
         ENDIF 
      ENDIF 
      frc3=a2(jp+1)-a2(jp) 
      px=a2(jp+1) 
      DO i=1,6
         ac(i)=a5(i,ip,jp+1)
      ENDDO
      frd3=frc3 
      s1=.true. 
  164 fr1=ac(1)*frc3 
      pr=(s(pres)-a2(jp))/(px-a2(jp)) 
      c0=1.d0/a5(1,ip,jp)
      r2=1.d0/ac(1) 
      c1=c0*a5(4,ip,jp)*(px-a2(jp)) 
      d2=r2*ac(4)*(px-a2(jp)) 
      c2=3.d0*(r2-c0)-d2-2.d0*c1 
      c3=d2+c1-2.d0*(r2-c0) 
      s12i=c0+pr*(c1+pr*(c2+pr*c3)) 
      s(vsubg)=1.d0/s12i 
      frc2=(s(vsubg)-a5(1,ip,jp))/(ac(1)-a5(1,ip,jp)) 
      frc=frc2*ac(1)/s(vsubg) 
      s(usubg)=a5(2,ip,jp)+(ac(2)-a5(2,ip,jp))*frc2 
      IF(s3) GOTO 168
!
  162 IF(.not.s1) THEN
         IF(ip+1.le.ns)THEN 
            IF(a31(ip+1).lt.a2(jp+1))THEN
               frd3=a31(ip+1)-a2(jp)
               px=a31(ip+1)
               DO i=1,6
                  ad(i)=a3(i,2,ip+1)
               ENDDO
               GOTO 166 
            ENDIF 
         ENDIF 
         frd3=a2(jp+1)-a2(jp) 
         px=a2(jp+1) 
         s1=.true. 
      ENDIF
      DO i=1,6
         ad(i)=a5(i,ip+1,jp+1)
      ENDDO
  166 fr1=ad(1)*frd3 
      pr=(s(pres)-a2(jp))/(px-a2(jp)) 
      c0=1.d0/a5(1,ip+1,jp)
      r2=1.d0/ad(1)
      c1=c0*a5(4,ip+1,jp)*(px-a2(jp)) 
      d2=r2*ad(4)*(px-a2(jp)) 
      c2=3.d0*(r2-c0)-d2-2.d0*c1 
      c3=d2+c1-2.d0*(r2-c0) 
      s11i=c0+pr*(c1+pr*(c2+pr*c3)) 
      s(vsubf)=1.d0/s11i 
      frd2=(s(vsubf)-a5(1,ip+1,jp))/(ad(1)-a5(1,ip+1,jp)) 
      frd=frd2*ad(1)/s(vsubf) 
      s(usubf)=a5(2,ip+1,jp)+(ad(2)-a5(2,ip+1,jp))*frd2 
      IF(s(usubf).lt.s(ubar))THEN 
         s2=.true. 
         ip=ip+1 
         IF(ip.eq.nt) THEN
!
!...........Vapor phase, temperature greater than highest table temperature.     
!
            fr=a5(1,ip,jp+1)*(a2(jp+1)-a2(jp)) 
            s(vbar) =a5(1,ip,jp)*fr/(fr+(a5(1,ip,jp)-a5(1,ip,jp+1))*frn)
            fr=(s(vbar)-a5(1,ip,jp))/(a5(1,ip,jp+1)-a5(1,ip,jp)) 
            frc=fr*a5(1,ip,jp+1)/s(vbar) 
            ut=a5(2,ip,jp)+(a5(2,ip,jp+1)-a5(2,ip,jp))*fr 
            s(beta) =a5(3,ip,jp)+frc*(a5(3,ip,jp+1)-a5(3,ip,jp)) 
            s(cp)   =a5(5,ip,jp)+frc*(a5(5,ip,jp+1)-a5(5,ip,jp)) 
            frd=s(cp)-s(pres)*s(vbar)*s(beta) 
            s(temp) =(s(ubar)-ut+frd*a1(nt))/frd
            frd=s(temp)/a1(nt) 
            s(vbar) =s(vbar)*frd
            s(beta) =s(beta)/frd
            s(hbar) =s(ubar)+s(pres)*s(vbar)
            s(kapa) =a5(4,ip,jp)+(a5(4,ip,jp+1)-a5(4,ip,jp))*fr
            s(entpy)=a5(6,ip,jp)+(a5(6,ip,jp+1)-a5(6,ip,jp))*fr
            cv=s(cp)-a1(nt)*s(beta)*s(beta)*s(vbar)/s(kapa)
            IF(s(temp).le.0.d0) GOTO 1001 
            s(entpy)=s(entpy)+cv*log(frd**(s(cp)/cv))
            IF(s(temp).le.1500.d0) THEN
               s(qual)=1.d0 
               iunp(1)=ip 
               iunp(2)=jp 
               s(is23)=unp 
               err=.false. 
               RETURN 
            ELSE
               GOTO 1001
            ENDIF
         ENDIF
         ta0=a1(ip) 
         DO i=1,6
            ac(i)=ad(i)
         ENDDO
         s(vsubg)=s(vsubf) 
         s(usubg)=s(usubf) 
         frc2=frd2 
         frc=frd 
         GOTO 162 
      ENDIF 
  168 IF(s(usubg).le.s(ubar)) THEN
         IF(.not.s2)THEN
            s(betag) =a0(3,2)+bma(3,2)*(fr*tb/s(tsat))
!           s(kapag) =a0(4,2)+bma(4,2)*((s(pres)-pa)/(pb-pa)*pb/s(pres))
            s(kapag) =a0(4,2)+(s(pres)-pa)/(pb-pa)*pb/s(pres)*bma(4,2)
            s(cpg)   =a0(5,2)+bma(5,2)*fr
            s(entpyg)=a0(6,2)+bma(6,2)*fr
         ELSE
            s(betag) =a5(3,ip,jp)+(ac(3)-a5(3,ip,jp))*frc
            s(kapag) =a5(4,ip,jp)+(ac(4)-a5(4,ip,jp))*frc2
            s(cpg)   =a5(5,ip,jp)+(ac(5)-a5(5,ip,jp))*frc
            s(entpyg)=a5(6,ip,jp)+(ac(6)-a5(6,ip,jp))*frc2
         ENDIF
         s(betaf) =a5(3,ip+1,jp)+(ad(3)-a5(3,ip+1,jp))*frd
         s(kapaf) =a5(4,ip+1,jp)+(ad(4)-a5(4,ip+1,jp))*frd2
         s(cpf)   =a5(5,ip+1,jp)+(ad(5)-a5(5,ip+1,jp))*frd
         s(entpyf)=a5(6,ip+1,jp)+(ad(6)-a5(6,ip+1,jp))*frd2
         fr=(s(ubar)-s(usubg))/(s(usubf)-s(usubg)) 
         frd=a1(ip+1)-ta0
         fr2=frd*fr 
         s(temp)=ta0+fr2 
         fr1=fr2/frd 
         tr=(s(temp)-ta0)/frd 
         c0=s(vsubg) 
         c1=s(vsubg)*s(betag)*frd 
         d2=s(vsubf)*s(betaf)*frd 
         c2=3.d0*(s(vsubf)-s(vsubg))-d2-2.d0*c1 
         c3=d2+c1-2.d0*(s(vsubf)-s(vsubg)) 
         s(vbar)=c0+tr*(c1+tr*(c2+tr*c3)) 
         s(hbar)=s(ubar)+s(pres)*s(vbar) 
         s(beta)=s(betag)+(s(betaf)-s(betag))*fr1 
!.....The above equation gave a bad value on a Puma problem.
!.....This lead to bad derivatives (neg dtgdus) and negative Cpg.
!.....It occured when at aroun tsat=330.
!.....tr was 0.48 for cell 64006. s6 was 0.22
!     while the sat value (s18) was about 0.003.
!.....limit beta gas to be less than beta gas sat
!.....12/1/98 i1 mod3.2.2 change
!
         IF(s(pres).lt.crp)s(beta)=min(s(beta),s(betag)) 
         s(kapa) =s(kapag)+(s(kapaf)-s(kapag))*fr1 
         s(cp)   =s(cpg)+(s(cpf)-s(cpg))*fr1 
         s(entpy)=s(entpyg)+(s(entpyf)-s(entpyg))*fr1 
         s(qual) =1.d0 
         iunp(1)=ip 
         iunp(2)=jp 
         s(is23)=unp 
         err=.false. 
         RETURN 
      ELSE
         IF(oncel) GOTO 1002 
         oncev=.true. 
         s3=.true. 
         s(vsubf)=s(vsubg) 
         s(vsubg)=vt 
         s(usubf)=s(usubg) 
         s(usubg)=ut 
         frd2=frc2 
         frd=frc 
         ip=ip-1 
         IF(ip.le.0) GOTO 1001 
         DO i=1,6
            ad(i)=ac(i)
         ENDDO
         IF(ip.ge.ns) GOTO 157 
         s3=.false. 
         GOTO 57 
      ENDIF
!                                                                       
 1002 WRITE(*,2001)s 
 2001 FORMAT  ('Interpolation failure in sth2x6_cupid.'/1p,5(4e20.10/),2e20.10,i20)                                                     
 1001 err=.true. 
!     
      END SUBROUTINE sth2x6_cupid                         
