      SUBROUTINE strtp_cupid(a,s,it,err) 
!deck strtp                                                             
!                                                                       
!      strtp   - compute water thermodynamic properties as a function of
!                temperature and pressure                               
!                                                                       
!      Calling sequence:                                                
!                                                                       
!                call  strtp (rp1,rp2,ip3,lp4)                          
!                                                                       
!      Parameters:                                                      
!                                                                       
!                rp1 = a   = steam tables (input)                       
!                                                                       
!                rp2 = s   = array into which the computed              
!                            thermodynamic properties are stored        
!                            (input,output)                             
!                                                                       
!                ip3 = it  = flag indicating physical state of steam    
!                            table fluid, i.e., liquid, vapor,          
!                            superheated vapor (output)                 
!                                                                       
!                lp4 = err = error flag (output)                        
!                                                                       
!                                                                       
!      This routine adapted from sth2x3 routine written by R. J. Wagner 
!      for light water steam tables                                     
!                                                                       
      USE STM_TBL_cupid  , ONLY: nt,np,pcrit
!
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL(8) a(*),s(*) 
      INTEGER it 
      LOGICAL err 
      REAL(8) unp,dpsdts,pa,ta,f1,f2,d1,d2,c0,c1,c2,c3,cv,ren,dum,      &
      dpdtcr,pb,tb,fr1,fr,frn,px,frc,frc2,frd,frd2,hfg1,hfg2,dpdt1,     &
      dpdt2                                                             
      INTEGER ip,jp,jpp,kp2,kp,ia,ib,jpq,lpp,lqq,ic,id 
      LOGICAL s1,s2 
      INTEGER iunp(2) 
      EQUIVALENCE(unp,iunp(1)) 
!
!------------------------------------------------------------------------------
!
!      GAS Modification by Won-Jae Lee: MUST in case New GAS is added
!                                                                 
      INCLUDE 'stcom.h' 
      INCLUDE 'gastable.h'
!   D. LMR-K.S. Ha for liquid metal properties- lead-bismuth eutetic(nfluid=11)
      INCLUDE 'lmtable.h'
      INTEGER ntg,npg,nstg,nspg,it3bpg,it4bpg,it5bpg,nprpntg,it3p0g 
! ------------------------------------------------------------------------
      ntg=nt 
      npg=np
      nstg=nst
      nspg=nsp
      it3bpg=it3bp
      it4bpg=it4bp
      it5bpg=it5bp
      nprpntg=nprpnt
      it3p0g=it3p0
      !IF(nfluid.eq.3) then
      !   ntg=ntgc
      !   npg=npgc
      !   nstg=nsc
      !   nspg=ns2c
      !   it3bpg=klpc
      !   it4bpg=klp2c
      !   it5bpg=llpc
      !   nprpntg=nt5c
      !   it3p0g=jplc
      !ELSEIF(nfluid.eq.4) then
      !   ntg=ntgh
      !   npg=npgh
      !   nstg=nsh
      !   nspg=ns2h
      !   it3bpg=klph
      !   it4bpg=klp2h
      !   it5bpg=llph
      !   nprpntg=nt5h
      !   it3p0g=jplh
      !ELSEIF(nfluid.eq.5) then
      !   ntg=ntgh2
      !   npg=npgh2
      !   nstg=nsh2
      !   nspg=ns2h2
      !   it3bpg=klph2
      !   it4bpg=klp2h2
      !   it5bpg=llph2
      !   nprpntg=nt5h2
      !   it3p0g=jplh2
      !ELSEIF(nfluid.eq.6) then
      !   ntg=ntgo
      !   npg=npgo
      !   nstg=nso
      !   nspg=ns2o
      !   it3bpg=klpo
      !   it4bpg=klp2o
      !   it5bpg=llpo
      !   nprpntg=nt5o
      !   it3p0g=jplo
      !ELSEIF(nfluid.eq.7) then
      !   ntg=ntgn
      !   npg=npgn
      !   nstg=nsn
      !   nspg=ns2n
      !   it3bpg=klpn
      !   it4bpg=klp2n
      !   it5bpg=llpn
      !   nprpntg=nt5n
      !   it3p0g=jpln
      !ELSEIF(nfluid.eq.8) then
      !   ntg=ntgna
      !   npg=npgna
      !   nstg=nsna
      !   nspg=ns2na
      !   it3bpg=klpna
      !   it4bpg=klp2na
      !   it5bpg=llpna
      !   nprpntg=nt5na
      !   it3p0g=jplna
      !ELSEIF(nfluid.eq.11) then
      !   ntg=ntlbe
      !   npg=nplbe
      !   nstg=nslbe
      !   nspg=ns2lbe
      !   it3bpg=klplbe
      !   it4bpg=klp2lbe
      !   it5bpg=llplbe
      !   nprpntg=nt5lbe
      !   it3p0g=jpllbe
      !ENDIF
! ------------------------------------------------------------------------
!                                                                       
!                                                                       
!--check for valid input                                                
!                                                                       
      IF(s(1).lt.a(1).or.s(1).gt.5000.d0) GOTO 1001 
      IF(s(2).le.0.0d0.or.s(2).gt.a(it3p0g)) GOTO 1001 
!                                                                       
!if -def,in32,5                                                         
!if def,cray,1                                                          
!      ip = shiftr(s(23),30)                                            
!if -def,cray,1                                                         
!      ip = ishft(s(23),-30)                                            
!      jp = and(s(23),msk)                                              
!if def,in32,3                                                          
      unp=s(23) 
      ip=iunp(1) 
      jp=iunp(2) 
      IF(ip.le.0.or.ip.gt.ntg)ip=1 
      IF(jp.le.0.or.jp.ge.npg)jp=1 
      IF(s(1).ge.a(ntg)) GOTO 46 
!                                                                       
!--set indexes in temperature and pressure tables for saturation        
!--computations                                                         
!                                                                       
   11 IF(s(1).ge.a(ip)) GOTO 10 
      ip=ip-1 
      GOTO 11 
   10 IF(s(1).lt.a(ip+1)) GOTO 12 
      ip=ip+1 
      GOTO 10 
   12 IF(ip.ge.nstg) GOTO 44 
      s1=.false. 
      s2=.false. 
!                                                                       
!--get saturation pressure                                              
!                                                                       
      CALL strsat_cupid(a,1,s(1),s(10),dpsdts,err) 
      IF(err) GOTO 1001 
!                   
      jpp=jp+ntg 
   13 IF(s(10).ge.a(jpp)) GOTO 14 
      jpp=jpp-1 
      IF(jpp.gt.ntg) GOTO 13 
      s2=.true. 
      GOTO 15 
   14 IF(s(10).lt.a(jpp+1)) GOTO 15 
      jpp=jpp+1 
      IF(jpp.lt.it3p0g) GOTO 14 
      s1=.true. 
   15 kp2=it4bpg+(jpp-ntg) * 13 
      kp=it3bpg+ip * 13 
      IF(s2.or.a(jpp).le.a(kp)) GOTO 16 
      pa=a(jpp) 
      ta=a(kp2) 
      ia=kp2 
      GOTO 17 
   16 ta=a(ip) 
      pa=a(kp) 
      ia=kp 
   17 IF(s1.or.a(jpp+1).ge.a(kp+13)) GOTO 18 
      pb=a(jpp+1) 
      tb=a(kp2+13) 
      ib=kp2+13 
      GOTO 19 
   18 tb=a(ip+1) 
      pb=a(kp+13) 
      ib=kp+13 
   19 fr1=s(1)-ta 
!                                                                       
!--compare input pressure to saturation pressure to determine vapor or  
!--liquid                                                               
!                                                                       
      fr=fr1/(tb-ta) 
      IF(s(2).le.s(10)) GOTO 43 
      it=1 
      s(9)=0.0d0 
      GOTO 50 
   46 ip=ntg 
   44 it=4 
      GOTO 45 
   43 it=3 
   45 s(9)=1.0d0 
   50 jpq=jp+ntg 
!                                                                       
!--search for single phase indexes                                      
!                                                                       
   51 IF(s(2).ge.a(jpq)) GOTO 53 
      jpq=jpq-1 
      IF(jpq.gt.ntg) GOTO 51 
      GOTO 90 
   53 IF(s(2).lt.a(jpq+1)) GOTO 54 
      jpq=jpq+1 
      IF(jpq.lt.it3p0g) GOTO 53 
      GOTO 1001 
   54 jp=jpq-ntg 
      lpp=it5bpg+jp*nprpntg+ip * 6 
      lqq=lpp+nprpntg 
      kp2=it4bpg+jp * 13 
      frn=s(1)-a(ip) 
      IF(it-3)60,70,73 
!--liquid phase                                                         
!                                                                       
   60 IF(a(jpq).ge.s(10)) GOTO 61 
      px=s(10) 
      s2=.false. 
      GOTO 62 
   61 px=a(jpq) 
      s2=.true. 
      IF(jp.gt.nspg) GOTO 63 
      IF(a(kp2).gt.a(ip+1)) GOTO 63 
      frc=frn/(a(kp2)-a(ip)) 
      frc2=s(1)-a(kp2) 
      ic=kp2+1 
      GOTO 62 
   63 frc=frn/(a(ip+1)-a(ip)) 
      frc2=s(1)-a(ip+1) 
      ic=lpp+6 
      frd=frc 
      frd2=frc2 
      GOTO 65 
   62 IF(jp.ge.nspg) GOTO 67 
      IF(a(kp2+13).gt.a(ip+1)) GOTO 67 
      frd=frn/(a(kp2+13)-a(ip)) 
      frd2=s(1)-a(kp2+13) 
      id=kp2+14 
      GOTO 66 
   67 frd=frn/(a(ip+1)-a(ip)) 
      frd2=s(1)-a(ip+1) 
   65 id=lqq+6 
   66 IF(s2) GOTO 69 
      hfg1=a(ia+8)-a(ia+2)+pa*(a(ia+7)-a(ia+1)) 
      hfg2=a(ib+8)-a(ib+2)+pb*(a(ib+7)-a(ib+1)) 
      dpdt1=hfg1/(ta*(a(ia+7)-a(ia+1))) 
      IF(tb.ne.tcrit.and.pb.ne.pcrit)then 
         dpdt2=hfg2/(tb*(a(ib+7)-a(ib+1))) 
      ELSE 
         IF(tb.eq.tcrit)then 
            CALL strsat_cupid(a,1,tb,dum,dpdtcr,err) 
         ELSE 
            CALL strsat_cupid(a,2,pb,dum,dpdtcr,err) 
         ENDIF 
         dpdt2=dpdtcr 
      ENDIF 
      f1=a(ia+1)*(a(ia+3)-a(ia+4)*dpdt1) 
      f2=a(ib+1)*(a(ib+3)-a(ib+4)*dpdt2) 
      d1=f1*(tb-ta) 
      d2=f2*(tb-ta) 
      c0=a(ia+1) 
      c1=d1 
      c2=3.d0*(a(ib+1)-a(ia+1))-d2-2.d0*d1 
      c3=d2+d1-2.d0*(a(ib+1)-a(ia+1)) 
      s(11)=c0+fr*(c1+fr*(c2+fr*c3)) 
      s(13)=a(ia+2)+(a(ib+2)-a(ia+2))*fr 
      s(17)=a(ia+3)+(a(ib+3)-a(ia+3))*fr 
      s(19)=a(ia+4)+(a(ib+4)-a(ia+4))*fr 
      s(21)=a(ia+5)+(a(ib+5)-a(ia+5))*fr 
      s(25)=a(ia+6)+(a(ib+6)-a(ia+6))*fr 
      GOTO 56 
   69 d1=a(lpp+2)*a(lpp)*(frn-frc2) 
      d2=a(ic+2)*a(ic)*(frn-frc2) 
      c0=a(lpp) 
      c1=d1 
      c2=3.d0*(a(ic)-a(lpp))-d2-2.d0*d1 
      c3=d2+d1-2.d0*(a(ic)-a(lpp)) 
      s(11)=c0+frc*(c1+frc*(c2+frc*c3)) 
      s(13)=a(lpp+1)+(a(ic+1)-a(lpp+1))*frc 
      s(17)=a(lpp+2)+(a(ic+2)-a(lpp+2))*frc 
      s(19)=a(lpp+3)+(a(ic+3)-a(lpp+3))*frc 
      s(21)=a(lpp+4)+(a(ic+4)-a(lpp+4))*frc 
      s(25)=a(lpp+5)+(a(ic+5)-a(lpp+5))*frc 
   56 d1=a(lqq+2)*a(lqq)*(frn-frd2) 
      d2=a(id+2)*a(id)*(frn-frd2) 
      c0=a(lqq) 
      c1=d1 
      c2=3.d0*(a(id)-a(lqq))-d2-2.d0*d1 
      c3=d2+d1-2.d0*(a(id)-a(lqq)) 
      s(12)=c0+frd*(c1+frd*(c2+frd*c3)) 
      s(14)=a(lqq+1)+(a(id+1)-a(lqq+1))*frd 
      s(18)=a(lqq+2)+(a(id+2)-a(lqq+2))*frd 
      s(20)=a(lqq+3)+(a(id+3)-a(lqq+3))*frd 
      s(22)=a(lqq+4)+(a(id+4)-a(lqq+4))*frd 
      s(26)=a(lqq+5)+(a(id+5)-a(lqq+5))*frd 
      IF(s(11).gt.s(12)) GOTO 83 
      s(3)=s(11) 
      fr=0.0d0 
      GOTO 84 
   83 s(15)=s(19)*s(11) 
      s(16)=s(20)*s(12) 
      fr1=s(16)-s(15) 
      IF(abs(fr1).lt.1.0d-10) GOTO 81 
      fr=s(11)+s(12)-(a(jpq+1)-px)*s(15)*s(16)/fr1 
      fr1=sqrt(fr*fr-4.0d0*s(11)*s(12)*(s(16)*(1.0d0-s(19)*(s(2)-px))-s(&
      15)*(1.0d0-s(20)*(s(2)-a(jpq+1))))/fr1)                           
      s(3)=0.5d0*(fr+fr1) 
      IF(s(3).gt.s(11))s(3)=0.5d0*(fr-fr1) 
      IF(s(3).ge.s(12)) GOTO 82 
   81 fr=(s(2)-px)/(a(jpq+1)-px) 
      s(3)=s(11)*(1.0d0-fr)+s(12)*fr 
   82 fr=(s(3)-s(11))/(s(12)-s(11)) 
   84 s(4)=s(13)+(s(14)-s(13))*fr 
      s(5)=s(4)+s(2)*s(3) 
      s(6)=s(17)+(s(18)-s(17))*fr 
      s(7)=s(19)+(s(20)-s(19))*fr 
      s(8)=s(21)+(s(22)-s(21))*fr 
      s(24)=s(25)+(s(26)-s(25))*fr 
      GOTO 20 
!                                                                       
!--vapor phase                                                          
!                                                                       
   70 s1=.false. 
      IF(a(jpq+1).le.s(10)) GOTO 71 
      frc=s(10) 
      hfg1=a(ia+8)-a(ia+2)+pa*(a(ia+7)-a(ia+1)) 
      hfg2=a(ib+8)-a(ib+2)+pb*(a(ib+7)-a(ib+1)) 
      dpdt1=hfg1/(ta*(a(ia+7)-a(ia+1))) 
      IF(tb.ne.tcrit.and.pb.ne.pcrit)then 
         dpdt2=hfg2/(tb*(a(ib+7)-a(ib+1))) 
      ELSE 
         IF(tb.eq.tcrit)then 
            CALL strsat_cupid(a,1,tb,dum,dpdtcr,err) 
         ELSE 
            CALL strsat_cupid(a,2,pb,dum,dpdtcr,err) 
         ENDIF 
         dpdt2=dpdtcr 
      ENDIF 
      f1=a(ia+7)*(a(ia+9)-a(ia+10)*dpdt1) 
      f2=a(ib+7)*(a(ib+9)-a(ib+10)*dpdt2) 
      d1=f1*(tb-ta) 
      d2=f2*(tb-ta) 
      c0=a(ia+7) 
      c1=d1 
      c2=3.d0*(a(ib+7)-a(ia+7))-d2-2.d0*d1 
      c3=d2+d1-2.d0*(a(ib+7)-a(ia+7)) 
      s(12)=c0+fr*(c1+fr*(c2+fr*c3)) 
      s(14)=a(ia+8)+(a(ib+8)-a(ia+8))*fr 
      s(18)=a(ia+9)+fr*tb/s(1)*(a(ib+9)-a(ia+9)) 
      s(20)=a(ia+10)+(s(10)-pa)/(pb-pa)*pb/s(10)*(a(ib+10)-a(ia+10)) 
      s(22)=a(ia+11)+(a(ib+11)-a(ia+11))*fr 
      s(26)=a(ia+12)+(a(ib+12)-a(ia+12))*fr 
      GOTO 72 
   71 IF(a(kp2+13).lt.a(ip)) GOTO 73 
      frd=(s(1)-a(kp2+13))/(a(ip+1)-a(kp2+13)) 
      ic=kp2+20 
      GOTO 74 
   73 IF(ip.eq.ntg) GOTO 80 
      frd=frn/(a(ip+1)-a(ip)) 
      s1=.true. 
      ic=lqq 
   74 frc=a(jpq+1) 
      c0=a(ic) 
      d1=a(ic)*a(ic+2)*(a(ip+1)-a(ip)) 
      d2=a(lqq+8)*a(lqq+6)*(a(ip+1)-a(ip)) 
      c1=d1 
      c2=3.d0*(a(lqq+6)-a(ic))-d2-2.d0*d1 
      c3=d2+d1-2.d0*(a(lqq+6)-a(ic)) 
      s(12)=c0+frd*(c1+frd*(c2+frd*c3)) 
      s(14)=a(ic+1)+(a(lqq+7)-a(ic+1))*frd 
      s(18)=a(ic+2)+frd*a(ip+1)/s(1)*(a(lqq+8)-a(ic+2)) 
      s(20)=a(ic+3)+(a(lqq+9)-a(ic+3))*frd 
      s(22)=a(ic+4)+(a(lqq+10)-a(ic+4))*frd 
      s(26)=a(ic+5)+(a(lqq+11)-a(ic+5))*frd 
   72 IF(s1) GOTO 75 
      IF(a(kp2).lt.a(ip)) GOTO 77 
      frd=(s(1)-a(kp2))/(a(ip+1)-a(kp2)) 
      ia=kp2+7 
      GOTO 76 
   77 frd=frn/(a(ip+1)-a(ip)) 
   75 ia=lpp 
   76 c0=a(ia) 
      d1=a(ia)*a(ia+2)*(a(ip+1)-a(ip)) 
      d2=a(lpp+8)*a(lpp+6)*(a(ip+1)-a(ip)) 
      c1=d1 
      c2=3.d0*(a(lpp+6)-a(ia))-d2-2.d0*d1 
      c3=d2+d1-2.d0*(a(lpp+6)-a(ia)) 
      s(11)=c0+frd*(c1+frd*(c2+frd*c3)) 
      s(13)=a(ia+1)+(a(lpp+7)-a(ia+1))*frd 
      s(17)=a(ia+2)+frd*a(ip+1)/s(1)*(a(lpp+8)-a(ia+2)) 
      s(19)=a(ia+3)+(a(lpp+9)-a(ia+3))*frd 
      s(21)=a(ia+4)+(a(lpp+10)-a(ia+4))*frd 
      s(25)=a(ia+5)+(a(lpp+11)-a(ia+5))*frd 
      fr=s(12)*(frc-a(jpq)) 
      s(3)=s(11)*fr/(fr+(s(11)-s(12))*(s(2)-a(jpq))) 
      fr=(s(3)-s(11))/(s(12)-s(11)) 
      frn=fr*s(12)/s(3) 
      s(4)=s(13)+(s(14)-s(13))*fr 
      s(5)=s(4)+s(2)*s(3) 
      s(6)=s(17)+frn*(s(18)-s(17)) 
      s(7)=s(19)+(s(20)-s(19))*fr 
      s(8)=s(21)+frn*(s(22)-s(21)) 
      s(24)=s(25)+(s(26)-s(25))*fr 
!if -def,in32,1                                                         
!  20  s(23) = or(shift(ip,30),jp)                                      
!if def,in32,3                                                          
   20 iunp(1)=ip 
      iunp(2)=jp 
      s(23)=unp 
      err=.false. 
      RETURN 
!                                                                       
!--vapor phase, temperature greater than highest table temperature      
!                                                                       
   80 fr=a(lqq)*(a(jpq+1)-a(jpq)) 
      s(3)=a(lpp)*fr/(fr+(a(lpp)-a(lqq))*(s(2)-a(jpq))) 
      fr=(s(3)-a(lpp))/(a(lqq)-a(lpp)) 
      frc=fr*a(lqq)/s(3) 
      s(5)=a(lpp+1)+(a(lqq+1)-a(lpp+1))*fr+s(2)*s(3) 
      s(8)=a(lpp+4)+frc*(a(lqq+4)-a(lpp+4)) 
      frd=s(1)/a(ntg) 
      s(3)=s(3)*frd 
      s(5)=s(5)+s(8)*frn 
      s(4)=s(5)-s(2)*s(3) 
      s(6)=(a(lpp+2)+frc*(a(lqq+2)-a(lpp+2)))/frd 
      s(7)=a(lpp+3)+(a(lqq+3)-a(lpp+3))*fr 
      s(24)=a(lpp+5)+(a(lqq+5)-a(lpp+5))*fr 
      cv=s(8)-a(ntg)*s(6)*s(6)*s(3)/s(7) 
      s(24)=s(24)+cv*log(frd**(s(8)/cv)) 
      GOTO 20 
!                                                                       
!--vapor phase, pressure less than lowest table pressure                
!                                                                       
   90 IF(it.eq.1) GOTO 1001 
      IF(s(1).lt.a(it4bpg+13)) GOTO 92 
      lpp=it5bpg+nprpntg+ip * 6 
      IF(ip.eq.ntg) GOTO 95 
      IF(a(ip).lt.a(it4bpg+13)) GOTO 93 
      ia=ip 
      lqq=lpp 
      GOTO 91 
   93 ia=it4bpg+13 
      lqq=ia+7 
   91 fr=(s(1)-a(ia))/(a(ip+1)-a(ia)) 
      s(3)=(fr*a(lpp+6)/a(ip+1)+(1.0d0-fr)*a(lqq)/a(ia))*a(jpq+1)*s(1)/ &
      s(2)                                                              
      s(4)=a(lqq+1)+(a(lpp+7)-a(lqq+1))*fr 
      s(6)=a(lqq+2)+(a(lpp+8)-a(lqq+2))*fr*a(ip+1)/s(1) 
      s(8)=a(lqq+4)+(a(lpp+10)-a(lqq+4))*fr 
      s(24)=a(lqq+5)+(a(lpp+11)-a(lqq+5))*fr 
      ren=s(2)*s(3)/s(1) 
      s(24)=s(24)-ren*log(s(2)/a(jpq+1)) 
   94 s(5)=s(4)+s(2)*s(3) 
      s(7)=1.0d0/s(2) 
      jp=1 
      GOTO 20 
   92 s(3)=(fr*pb*a(ib+7)/tb+(1.0d0-fr)*pa*a(ia+7)/ta)*s(1)/s(2) 
      s(4)=a(ia+8)+(a(ib+8)-a(ia+8))*fr 
      s(6)=1.0d0/s(1) 
      s(8)=a(ia+11)+(a(ib+11)-a(ia+11))*fr 
      s(24)=a(ia+12)+(a(ib+12)-a(ia+12))*fr 
      ren=s(2)*s(3)/s(1) 
      s(24)=s(24)-ren*log(s(2)/s(10)) 
      GOTO 94 
   95 frd=s(1)/a(ntg) 
      frc=a(ntg+1)*a(lpp) 
      s(3)=frc*frd/s(2) 
      s(8)=a(lpp+4) 
      s(5)=a(lpp+1)+frc+s(8)*(s(1)-a(ntg)) 
      s(4)=s(5)-s(2)*s(3) 
      s(6)=a(lpp+2)/frd 
      s(7)=a(lpp+3)*a(ntg+1)/s(2) 
      ren=s(2)*s(3)/s(1) 
      s(24)=a(lpp+5)+s(8)*log(s(1)/a(ntg))-ren*log(s(2)/a(jpq+1)) 
      GOTO 20 
!                                                                       
!--error                                                                
 1001 err=.true. 
      RETURN 
      END SUBROUTINE strtp_cupid                          
