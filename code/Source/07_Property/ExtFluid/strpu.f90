      SUBROUTINE strpu_cupid(a,s,it,err) 
!deck strpu                                                             
!                                                                       
!      strpu   - compute water thermodynamic properties as a function of
!                pressure and internal energy                           
!                                                                       
!      Calling sequence:                                                
!                                                                       
!                call  strpu (rp1,rp2,ip3,lp4)                          
!                                                                       
!      Parameters:                                                      
!                                                                       
!                rp1 = a   = steam tables (input)                       
!                                                                       
!                rp2 = s   = array into which the computed              
!                            thermodynamic properties are stored        
!                            (input,output)                             
!                                                                       
!                ip3 = it  = flag controlling whether saturation        
!                            temperature is to be calculated, i.e.,     
!                            0 => calculate saturation temperature      
!                            (input);                                   
!                            flag indicating physical state of fluid,   
!                            i.e., liquid, vapor, superheated vapor     
!                            (output)                                   
!                                                                       
!                lp4 = err = error flag (output)                        
!                                                                       
!                                                                       
!      This routine adapted from sth2x6 routine written by R. J. Wagner 
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
      LOGICAL s1,s2,s3 
      INTEGER iunp(2) 
      EQUIVALENCE(unp,iunp(1)) 
      REAL(8) fr1,fr,hfg1,hfg2,dpdt1,dpdt2,dum,dpdtcr,f1,f2,pb,tb,unp,  &
      dpsdts,pa,ta,c1,d2,c0,c2,c3,frn,frc2,ut,vt,frc3,frc,px,frd3,frd,  &
      pr,r2,s11i,fr2,s12i,tr,frd2,cv,ren                                
      INTEGER ib,ip,jp,jpp,ic,kp,kp2,ia,lpp,lqq,id 
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
!      IF(nfluid.eq.3) then
!         ntg=ntgc
!         npg=npgc
!         nstg=nsc
!         nspg=ns2c
!         it3bpg=klpc
!         it4bpg=klp2c
!         it5bpg=llpc
!         nprpntg=nt5c
!         it3p0g=jplc
!      ELSEIF(nfluid.eq.4) then
!         ntg=ntgh
!         npg=npgh
!         nstg=nsh
!         nspg=ns2h
!         it3bpg=klph
!         it4bpg=klp2h
!         it5bpg=llph
!         nprpntg=nt5h
!         it3p0g=jplh
!      ELSEIF(nfluid.eq.5) then
!         ntg=ntgh2
!         npg=npgh2
!         nstg=nsh2
!         nspg=ns2h2
!         it3bpg=klph2
!         it4bpg=klp2h2
!         it5bpg=llph2
!         nprpntg=nt5h2
!         it3p0g=jplh2
!      ELSEIF(nfluid.eq.6) then
!         ntg=ntgo
!         npg=npgo
!         nstg=nso
!         nspg=ns2o
!         it3bpg=klpo
!         it4bpg=klp2o
!         it5bpg=llpo
!         nprpntg=nt5o
!         it3p0g=jplo
!      ELSEIF(nfluid.eq.7) then
!         ntg=ntgn
!         npg=npgn
!         nstg=nsn
!         nspg=ns2n
!         it3bpg=klpn
!         it4bpg=klp2n
!         it5bpg=llpn
!         nprpntg=nt5n
!         it3p0g=jpln
!      ELSEIF(nfluid.eq.8) then
!         ntg=ntgna
!         npg=npgna
!         nstg=nsna
!         nspg=ns2na
!         it3bpg=klpna
!         it4bpg=klp2na
!         it5bpg=llpna
!         nprpntg=nt5na
!         it3p0g=jplna
!      ELSEIF(nfluid.eq.11) then
!         ntg=ntlbe
!         npg=nplbe
!         nstg=nslbe
!         nspg=ns2lbe
!         it3bpg=klplbe
!         it4bpg=klp2lbe
!         it5bpg=llplbe
!         nprpntg=nt5lbe
!         it3p0g=jpllbe
!      ENDIF
! ------------------------------------------------------------------------
!                                                                       
!                                                                       
!--set flag for saturation temperature calculation                      
!                                                                       
      s3=it.ne.0 
!                                                                       
!--check for valid input                                                
!                                                                       
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
      IF(ip.le.0.or.ip.ge.ntg)ip=1 
      IF(jp.le.0.or.jp.ge.npg)jp=1 
      jpp=jp+ntg 
      s1=.false. 
!                                                                       
!--set indexes in temperature and pressure tables for saturation        
!--computations                                                         
!                                                                       
   11 IF(s(2).ge.a(jpp)) GOTO 10 
      jpp=jpp-1 
      IF(jpp.gt.ntg) GOTO 11 
      jpp=jpp+1 
      s1=.true. 
      IF(s(2).lt.a(it3bpg+13)) GOTO 44 
      GOTO 12 
   10 IF(s(2).lt.a(jpp+1)) GOTO 12 
      jpp=jpp+1 
      GOTO 10 
   12 jp=jpp-ntg 
      IF(s(2).ge.pcrit) GOTO 44 
      IF(s3) GOTO 15 
!                                                                       
!--get saturation temperature                                           
!                                                                       
      CALL strsat_cupid(a,2,s(2),s(10),dpsdts,err) 
      IF(err) GOTO 1001 
!                                                                       
   15 ic=ip 
   16 IF(s(10).ge.a(ic)) GOTO 13 
      ic=ic-1 
      IF(ic.gt.0) GOTO 16 
      ic=1 
      GOTO 14 
   13 IF(s(10).lt.a(ic+1)) GOTO 14 
      ic=ic+1 
      IF(ic.lt.nstg) GOTO 13 
      ic=ic-1 
   14 kp=it3bpg+ic * 13 
      kp2=it4bpg+jp * 13 
      IF(s1.or.a(jpp).le.a(kp)) GOTO 19 
      pa=a(jpp) 
      ta=a(kp2) 
      ia=kp2 
      GOTO 20 
   19 ta=a(ic) 
      pa=a(kp) 
      ia=kp 
   20 IF(a(jpp+1).ge.a(kp+13)) GOTO 21 
      pb=a(jpp+1) 
      tb=a(kp2+13) 
      ib=kp2+13 
      GOTO 22 
   21 pb=a(kp+13) 
      tb=a(ic+1) 
      ib=kp+13 
   22 fr1=s(10)-ta 
      fr=fr1/(tb-ta) 
!                                                                       
!--compute vsubf and vsubg to determine liquid, two phase, or vapor     
!--state                                                                
!                                                                       
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
      c1=f1*(tb-ta) 
      d2=f2*(tb-ta) 
      c0=a(ia+1) 
      c2=3.d0*(a(ib+1)-c0)-d2-2.d0*c1 
      c3=d2+c1-2.d0*(a(ib+1)-c0) 
      s(11)=c0+fr*(c1+fr*(c2+fr*c3)) 
   24 s(13)=a(ia+2)+(a(ib+2)-a(ia+2))*fr 
      IF(s(4).le.s(13)) GOTO 41 
      f1=a(ia+7)*(a(ia+9)-a(ia+10)*dpdt1) 
      f2=a(ib+7)*(a(ib+9)-a(ib+10)*dpdt2) 
      c1=f1*(tb-ta) 
      d2=f2*(tb-ta) 
      c0=a(ia+7) 
      c2=3.d0*(a(ib+7)-c0)-d2-2.d0*c1 
      c3=d2+c1-2.d0*(a(ib+7)-c0) 
      s(12)=c0+fr*(c1+fr*(c2+fr*c3)) 
      s(14)=a(ia+8)+(a(ib+8)-a(ia+8))*fr 
      IF(s(4).ge.s(14)) GOTO 43 
!                                                                       
!--two phase fluid                                                      
!                                                                       
      it=2 
      s(15)=s(13)+s(2)*s(11) 
      s(16)=s(14)+s(2)*s(12) 
      s(17)=a(ia+3)+(a(ib+3)-a(ia+3))*fr 
      s(18)=a(ia+9)+fr*tb/s(10)*(a(ib+9)-a(ia+9)) 
      s(19)=a(ia+4)+(a(ib+4)-a(ia+4))*fr 
      s(20)=a(ia+10)+(s(2)-pa)/(pb-pa)*pb/s(2)*(a(ib+10)-a(ia+10)) 
      s(21)=a(ia+5)+(a(ib+5)-a(ia+5))*fr 
      s(22)=a(ia+11)+(a(ib+11)-a(ia+11))*fr 
      s(25)=a(ia+6)+(a(ib+6)-a(ia+6))*fr 
      s(26)=a(ia+12)+(a(ib+12)-a(ia+12))*fr 
      s(9)=(s(4)-s(13))/(s(14)-s(13)) 
      fr=1.0d0-s(9) 
      s(1)=s(10) 
      s(3)=fr*s(11)+s(9)*s(12) 
      s(5)=fr*s(15)+s(9)*s(16) 
      s(24)=fr*s(25)+s(9)*s(26) 
      ip=ic 
!if -def,in32,1                                                         
!  25  s(23) = or(shift(ip,30),jp)                                      
!if def,in32,3                                                          
   25 iunp(1)=ip 
      iunp(2)=jp 
      s(23)=unp 
      err=.false. 
      RETURN 
   44 it=4 
      GOTO 50 
!                                                                       
!--single phase fluid, search for single phase indexes                  
!                                                                       
   41 it=1 
      IF(s1) GOTO 1001 
      GOTO 50 
   43 it=3 
   50 lpp=it5bpg+jp*nprpntg+ip * 6 
   51 IF(s(4).ge.a(lpp+1)) GOTO 52 
      lpp=lpp-6 
      ip=ip-1 
      IF(ip.gt.0) GOTO 51 
      lpp=lpp+6 
      ip=ip+1 
      GOTO 54 
   52 IF(s(4).le.a(lpp+7)) GOTO 54 
      lpp=lpp+6 
      ip=ip+1 
      IF(ip.lt.ntg) GOTO 52 
      IF(s1) GOTO 95 
      lpp=lpp-6 
      ip=ip-1 
      GOTO 53 
   54 IF(s1) GOTO 96 
   53 kp=it3bpg+ip * 13 
      lqq=lpp+nprpntg 
      frn=s(2)-a(jpp) 
      frc2=s(2)-a(jpp+1) 
      s3=.false. 
      IF(it-3)58,70,48 
   48 IF(ip.ge.nstg) GOTO 157 
      GOTO 57 
!                                                                       
!--liquid phase                                                         
!                                                                       
   58 ut=s(13) 
      vt=s(11) 
   60 IF(a(ip).lt.s(10)) GOTO 57 
      ip=ip-1 
      IF(ip.le.0) GOTO 1001 
      kp=kp-13 
      lpp=lpp-6 
      lqq=lqq-6 
      GOTO 60 
   57 s1=.false. 
      IF(it.eq.4) GOTO 61 
      IF(a(ip+1).le.s(10)) GOTO 61 
      ta=s(10) 
      s2=.false. 
      GOTO 55 
   61 ta=a(ip+1) 
      s2=.true. 
      IF(ip+1.gt.nstg) GOTO 63 
      IF(a(kp+13).le.a(jpp)) GOTO 63 
      frc3=s(2)-a(kp+13) 
      frc=a(jpp+1)-a(kp+13) 
      px=a(kp+13) 
      ic=kp+14 
      GOTO 64 
   63 frc3=frn 
      frc=a(jpp+1)-a(jpp) 
      px=a(jpp) 
      ic=lpp+6 
      frd3=frc3 
      frd=frc 
      s1=.true. 
   64 pr=(s(2)-px)/(a(jpp+1)-px) 
      c0=1.d0/a(ic) 
      r2=1.d0/a(lqq+6) 
      c1=c0*a(ic+3)*(a(jpp+1)-px) 
      d2=r2*a(lqq+9)*(a(jpp+1)-px) 
      c2=3.d0*(r2-c0)-d2-2.d0*c1 
      c3=d2+c1-2.d0*(r2-c0) 
      s11i=c0+pr*(c1+pr*(c2+pr*c3)) 
      s(11)=1.d0/s11i 
!  82  fr1 = (s(11)-a(ic))/(a(lqq+6)-a(ic))                             
      fr1=pr 
      s(13)=a(ic+1)+(a(lqq+7)-a(ic+1))*fr1 
   55 IF(s3) GOTO 68 
   62 IF(s1) GOTO 65 
      IF(ip.gt.nstg) GOTO 67 
      IF(a(kp).le.a(jpp)) GOTO 67 
      frd3=s(2)-a(kp) 
      frd=a(jpp+1)-a(kp) 
      px=a(kp) 
      id=kp+1 
      GOTO 66 
   67 frd3=frn 
      frd=a(jpp+1)-a(jpp) 
      px=a(jpp) 
      s1=.true. 
   65 id=lpp 
   66 s(6)=a(id+3)*a(id) 
      s(7)=a(lqq+3)*a(lqq) 
      fr2=s(7)-s(6) 
      pr=(s(2)-px)/(a(jpp+1)-px) 
      c0=1.d0/a(id) 
      r2=1.d0/a(lqq) 
      c1=c0*a(id+3)*(a(jpp+1)-px) 
      d2=r2*a(lqq+3)*(a(jpp+1)-px) 
      c2=3.d0*(r2-c0)-d2-2.d0*c1 
      c3=d2+c1-2.d0*(r2-c0) 
      s12i=c0+pr*(c1+pr*(c2+pr*c3)) 
      s(12)=1.d0/s12i 
! 182  fr2 = (s(12)-a(id))/(a(lqq)-a(id))                               
      fr2=pr 
      s(14)=a(id+1)+(a(lqq+1)-a(id+1))*fr2 
      IF(s(14).le.s(4)) GOTO 68 
      s2=.true. 
      s(11)=s(12) 
      s(13)=s(14) 
      fr1=fr2 
      ip=ip-1 
      ta=a(ip+1) 
      IF(ip.le.0) GOTO 1001 
      kp=kp-13 
      lqq=lqq-6 
      lpp=lpp-6 
      ic=id 
      GOTO 62 
   68 IF(s(13).ge.s(4)) GOTO 59 
      s(12)=s(11) 
      s(11)=vt 
      s(14)=s(13) 
      s(13)=ut 
      fr2=fr1 
      lqq=lqq+6 
      lpp=lpp+6 
      kp=kp+13 
      ip=ip+1 
      IF(ip.lt.nstg) GOTO 158 
      s3=.false. 
      GOTO 157 
  158 id=ic 
      s3=.true. 
      GOTO 57 
   59 IF(s2) GOTO 69 
      s(17)=a(ia+3)+(a(ib+3)-a(ia+3))*fr 
      s(19)=a(ia+4)+(a(ib+4)-a(ia+4))*fr 
      s(21)=a(ia+5)+(a(ib+5)-a(ia+5))*fr 
      s(25)=a(ia+6)+(a(ib+6)-a(ia+6))*fr 
      GOTO 56 
   69 s(17)=a(ic+2)+(a(lqq+8)-a(ic+2))*fr1 
      s(19)=a(ic+3)+(a(lqq+9)-a(ic+3))*fr1 
      s(21)=a(ic+4)+(a(lqq+10)-a(ic+4))*fr1 
      s(25)=a(ic+5)+(a(lqq+11)-a(ic+5))*fr1 
   56 s(18)=a(id+2)+(a(lqq+2)-a(id+2))*fr2 
      s(20)=a(id+3)+(a(lqq+3)-a(id+3))*fr2 
      s(22)=a(id+4)+(a(lqq+4)-a(id+4))*fr2 
      s(26)=a(id+5)+(a(lqq+5)-a(id+5))*fr2 
      fr=(s(4)-s(14))/(s(13)-s(14)) 
      frd=ta-a(ip) 
      fr2=frd*fr 
      s(1)=a(ip)+fr2 
      fr1=fr2/frd 
      tr=(s(1)-a(ip))/frd 
      c0=s(12) 
      c1=s(12)*s(18)*frd 
      d2=s(11)*s(17)*frd 
      c2=3.d0*(s(11)-s(12))-d2-2.d0*c1 
      c3=d2+c1-2.d0*(s(11)-s(12)) 
      s(3)=c0+tr*(c1+tr*(c2+tr*c3)) 
      s(5)=s(4)+s(2)*s(3) 
      s(6)=(c1+tr*(2.d0*c2+3.d0*tr*c3))/(frd*s(3)) 
      s(7)=s(20)+(s(19)-s(20))*fr1 
      s(8)=s(22)+(s(21)-s(22))*fr1 
      s(24)=s(26)+(s(25)-s(26))*fr1 
      s(9)=0.d0 
      GOTO 25 
   70 ut=s(14) 
      vt=s(12) 
  160 IF(a(ip+1).gt.s(10)) GOTO 157 
!                                                                       
!--vapor phase                                                          
!                                                                       
      ip=ip+1 
      kp=kp+13 
      lpp=lpp+6 
      lqq=lqq+6 
      GOTO 160 
  157 s1=.false. 
      IF(it.eq.4) GOTO 161 
      IF(a(ip).ge.s(10)) GOTO 161 
      ta=s(10) 
      s2=.false. 
      GOTO 155 
  161 ta=a(ip) 
      s2=.true. 
      IF(ip.ge.nstg) GOTO 163 
      IF(a(kp).ge.a(jpp+1)) GOTO 163 
      frc3=a(kp)-a(jpp) 
      px=a(kp) 
      ic=kp+7 
      GOTO 164 
  163 frc3=a(jpp+1)-a(jpp) 
      px=a(jpp+1) 
      ic=lqq 
      frd3=frc3 
      s1=.true. 
  164 fr1=a(ic)*frc3 
      pr=(s(2)-a(jpp))/(px-a(jpp)) 
      c0=1.d0/a(lpp) 
      r2=1.d0/a(ic) 
      c1=c0*a(lpp+3)*(px-a(jpp)) 
      d2=r2*a(ic+3)*(px-a(jpp)) 
      c2=3.d0*(r2-c0)-d2-2.d0*c1 
      c3=d2+c1-2.d0*(r2-c0) 
      s12i=c0+pr*(c1+pr*(c2+pr*c3)) 
      s(12)=1.d0/s12i 
      frc2=(s(12)-a(lpp))/(a(ic)-a(lpp)) 
      frc=frc2*a(ic)/s(12) 
      s(14)=a(lpp+1)+(a(ic+1)-a(lpp+1))*frc2 
  155 IF(s3) GOTO 168 
  162 IF(s1) GOTO 165 
      IF(ip+1.gt.nstg) GOTO 167 
      IF(a(kp+13).ge.a(jpp+1)) GOTO 167 
      frd3=a(kp+13)-a(jpp) 
      px=a(kp+13) 
      id=kp+20 
      GOTO 166 
  167 frd3=a(jpp+1)-a(jpp) 
      px=a(jpp+1) 
      s1=.true. 
  165 id=lqq+6 
  166 fr1=a(id)*frd3 
      pr=(s(2)-a(jpp))/(px-a(jpp)) 
      c0=1.d0/a(lpp+6) 
      r2=1.d0/a(id) 
      c1=c0*a(lpp+9)*(px-a(jpp)) 
      d2=r2*a(id+3)*(px-a(jpp)) 
      c2=3.d0*(r2-c0)-d2-2.d0*c1 
      c3=d2+c1-2.d0*(r2-c0) 
      s11i=c0+pr*(c1+pr*(c2+pr*c3)) 
      s(11)=1.d0/s11i 
      frd2=(s(11)-a(lpp+6))/(a(id)-a(lpp+6)) 
      frd=frd2*a(id)/s(11) 
      s(13)=a(lpp+7)+(a(id+1)-a(lpp+7))*frd2 
      IF(s(13).ge.s(4)) GOTO 168 
      s2=.true. 
      ip=ip+1 
      lqq=lqq+6 
      lpp=lpp+6 
      IF(ip.eq.ntg) GOTO 80 
      ta=a(ip) 
      kp=kp+13 
      ic=id 
      s(12)=s(11) 
      s(14)=s(13) 
      frc2=frd2 
      frc=frd 
      GOTO 162 
  168 IF(s(14).le.s(4)) GOTO 159 
      s3=.true. 
      s(11)=s(12) 
      s(12)=vt 
      s(13)=s(14) 
      s(14)=ut 
      frd2=frc2 
      frd=frc 
      ip=ip-1 
      IF(ip.le.0) GOTO 1001 
      kp=kp-13 
      lpp=lpp-6 
      lqq=lqq-6 
      id=ic 
      IF(ip.ge.nstg) GOTO 157 
      s3=.false. 
      GOTO 57 
  159 IF(s2) GOTO 169 
      s(18)=a(ia+9)+fr*tb/s(10)*(a(ib+9)-a(ia+9)) 
      s(20)=a(ia+10)+(s(2)-pa)/(pb-pa)*pb/s(2)*(a(ib+10)-a(ia+10)) 
      s(22)=a(ia+11)+(a(ib+11)-a(ia+11))*fr 
      s(26)=a(ia+12)+(a(ib+12)-a(ia+12))*fr 
      GOTO 156 
  169 s(18)=a(lpp+2)+frc*(a(ic+2)-a(lpp+2)) 
      s(20)=a(lpp+3)+(a(ic+3)-a(lpp+3))*frc2 
      s(22)=a(lpp+4)+frc*(a(ic+4)-a(lpp+4)) 
      s(26)=a(lpp+5)+(a(ic+5)-a(lpp+5))*frc2 
  156 s(17)=a(lpp+8)+frd*(a(id+2)-a(lpp+8)) 
      s(19)=a(lpp+9)+(a(id+3)-a(lpp+9))*frd2 
      s(21)=a(lpp+10)+frd*(a(id+4)-a(lpp+10)) 
      s(25)=a(lpp+11)+(a(id+5)-a(lpp+11))*frd2 
      fr=(s(4)-s(14))/(s(13)-s(14)) 
      frd=a(ip+1)-ta 
      fr2=frd*fr 
      s(1)=ta+fr2 
      fr1=fr2/frd 
      tr=(s(1)-ta)/frd 
      c0=s(12) 
      c1=s(12)*s(18)*frd 
      d2=s(11)*s(17)*frd 
      c2=3.d0*(s(11)-s(12))-d2-2.d0*c1 
      c3=d2+c1-2.d0*(s(11)-s(12)) 
      s(3)=c0+tr*(c1+tr*(c2+tr*c3)) 
      s(5)=s(4)+s(2)*s(3) 
      s(6)=(c1+tr*(2.d0*c2+3.d0*tr*c3))/(frd*s(3)) 
      s(7)=s(20)+(s(19)-s(20))*fr1 
      s(8)=s(22)+(s(21)-s(22))*fr1 
      s(24)=s(26)+(s(25)-s(26))*fr1 
   99 s(9)=1.0d0 
      GOTO 25 
!                                                                       
!--vapor phase, temperature greater than highest table temperature      
!                                                                       
   80 fr=a(lqq)*(a(jpp+1)-a(jpp)) 
      s(3)=a(lpp)*fr/(fr+(a(lpp)-a(lqq))*frn) 
      fr=(s(3)-a(lpp))/(a(lqq)-a(lpp)) 
      frc=fr*a(lqq)/s(3) 
      ut=a(lpp+1)+(a(lqq+1)-a(lpp+1))*fr 
      s(6)=a(lpp+2)+frc*(a(lqq+2)-a(lpp+2)) 
      s(8)=a(lpp+4)+frc*(a(lqq+4)-a(lpp+4)) 
      frd=s(8)-s(2)*s(3)*s(6) 
      s(1)=(s(4)-ut+frd*a(ntg))/frd 
      frd=s(1)/a(ntg) 
      s(3)=s(3)*frd 
      s(6)=s(6)/frd 
      s(5)=s(4)+s(2)*s(3) 
      s(7)=a(lpp+3)+(a(lqq+3)-a(lpp+3))*fr 
      s(24)=a(lpp+5)+(a(lqq+5)-a(lpp+5))*fr 
      cv=s(8)-a(ntg)*s(6)*s(6)*s(3)/s(7) 
      s(24)=s(24)+cv*log(frd**(s(8)/cv)) 
      IF(s(1).le.tmax) GOTO 99 
      GOTO 1001 
!                                                                       
!--vapor phase, pressure less than lowest table pressure                
!                                                                       
   96 ut=a(it4bpg+21) 
      IF(s(4).lt.ut) GOTO 90 
      IF(a(ip).lt.a(it4bpg+13)) GOTO 93 
      ia=ip 
      lqq=lpp 
      ut=a(lqq+1) 
      GOTO 91 
   93 ia=it4bpg+13 
      lqq=ia+7 
   91 fr=(s(4)-ut)/(a(lpp+7)-ut) 
      fr1=a(ip+1)-a(ia) 
      fr2=fr1*fr 
      s(1)=a(ia)+fr2 
      fr1=fr2/fr1 
      s(3)=(fr1*a(lpp+6)/a(ip+1)+(1.d0-fr1)*a(lqq)/a(ia))*a(ntg+1)*s(1)/ &
      s(2)                                                              
      s(6)=a(lqq+2)+(a(lpp+8)-a(lqq+2))*fr*a(ip+1)/s(1) 
      s(8)=a(lqq+4)+(a(lpp+10)-a(lqq+4))*fr 
      s(24)=a(lqq+5)+(a(lpp+11)-a(lqq+5))*fr 
      ren=s(2)*s(3)/s(1) 
      s(24)=s(24)-ren*log(s(2)/a(ntg+1)) 
   94 s(7)=1.d0/s(2) 
   98 s(5)=s(4)+s(2)*s(3) 
      GOTO 99 
   95 s(8)=a(lpp+4) 
      frd=s(8)-a(ntg+1)*a(lpp)/a(ntg) 
      s(1)=(s(4)-a(lpp+1)+frd*a(ntg))/frd 
      frd=s(1)/a(ntg) 
      frc=a(ntg+1)*a(lpp) 
      s(3)=frc*frd/s(2) 
      s(6)=a(lpp+2)/frd 
      s(7)=a(lpp+3)*a(ntg+1)/s(2) 
      ren=s(2)*s(3)/s(1) 
      s(24)=a(lpp+5)+s(8)*log(s(1)/a(ntg))-ren*log(s(2)/a(ntg+1)) 
      GOTO 98 
   90 IF(s(4).lt.a(it3bpg+21)) GOTO 1001 
      ip=1 
      kp=it3bpg+13 
  202 IF(s(5).le.a(kp+21)) GOTO 201 
      ip=ip+1 
      kp=kp+13 
      ut=a(kp+21) 
      GOTO 202 
  201 fr=(s(4)-a(kp+8))/(a(kp+21)-a(kp+8)) 
      s(1)=a(ip)+fr*(a(ip+1)-a(ip)) 
      s(3)=(fr*a(kp+13)*a(kp+20)/a(ip+1)+(1.0d0-fr)*a(kp)*a(kp+7)/a(ip))&
      *s(1)/s(2)                                                        
      s(6)=1.d0/s(1) 
      s(8)=a(kp+11)+(a(kp+24)-a(kp+11))*fr 
      ren=s(2)*s(3)/s(1) 
      s(24)=a(kp+12)+s(8)*log(s(1)/a(ip))-ren*log(s(2)/a(kp)) 
      GOTO 94 
!                                                                       
!--error                                                                
 1001 err=.true. 
      RETURN 
      END SUBROUTINE strpu_cupid                          
