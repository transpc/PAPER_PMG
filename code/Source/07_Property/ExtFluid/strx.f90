      SUBROUTINE strx_cupid(a,s,err) 
!deck strx                                                              
!                                                                       
!      strx    - compute thermodynamic properties as a function of      
!                quality (saturation temperature and pressure previously
!                determined)                                            
!                                                                       
!      Calling sequence:                                                
!                                                                       
!                call  strx (rp1,rp2,lp3)                               
!                                                                       
!      Parameters:                                                      
!                                                                       
!                rp1 = a   = steam tables (input)                       
!                                                                       
!                rp2 = s   = array into which the computed              
!                            thermodynamic properties are stored        
!                            (input,output)                             
!                                                                       
!                lp3 = err = error flag (output)                        
!                                                                       
!                                                                       
!      This routine adapted from sth2xb (entry point in sth2x1) written 
!      by R. J. Wagner for light water steam tables                     
!                                                                       
      USE STM_TBL_cupid  , ONLY: nt,np,pcrit
!
      IMPLICIT none 
!                                                                       
      REAL(8) a(*),s(*) 
      LOGICAL err 
      REAL(8) unp,pa,ta,pb,tb,fr1,fr,dpdt1,dpdt2,dum,dpdtcr,f1,f2,hfg1, &
      hfg2,d1,d2,c0,c1,c2,c3                                            
      INTEGER ip,jp,jpp,kp2,kp,ia,ib 
      LOGICAL s2,s3 
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
!--check for valid input                                                
!                                                                       
      IF(s(9).lt.0.0d0.or.s(9).gt.1.0d0)then 
         err=.true. 
         GOTO 20 
      ENDIF 
!                                                                       
!--compute thermodynamic properties                                     
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
      s2=.false. 
      s3=.false. 
      IF(ip.le.0.or.ip.ge.nstg)ip=1 
      IF(jp.le.0.or.jp.ge.nspg)jp=1 
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
   12 jpp=jp+ntg 
  111 IF(s(10).ge.a(jpp)) GOTO 110 
      jpp=jpp-1 
      IF(jpp.gt.ntg) GOTO 111 
      s3=.true. 
      GOTO 112 
  110 IF(s(10).lt.a(jpp+1)) GOTO 112 
      jpp=jpp+1 
      IF(jpp.lt.it3p0g) GOTO 110 
      s2=.true. 
  112 jp=jpp-ntg 
      kp2=it4bpg+jp * 13 
      kp=it3bpg+ip * 13 
      IF(s3.or.a(jpp).le.a(kp)) GOTO 113 
      pa=a(jpp) 
      ta=a(kp2) 
      ia=kp2 
      GOTO 115 
  113 ta=a(ip) 
      pa=a(kp) 
      ia=kp 
  115 IF(s2.or.a(jpp+1).ge.a(kp+13)) GOTO 116 
      pb=a(jpp+1) 
      tb=a(kp2+13) 
      ib=kp2+13 
      GOTO 117 
  116 tb=a(ip+1) 
      pb=a(kp+13) 
      ib=kp+13 
  117 fr1=s(1)-ta 
      fr=fr1/(tb-ta) 
!                                                                       
!--two phase fluid                                                      
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
      d1=f1*(tb-ta) 
      d2=f2*(tb-ta) 
      c0=a(ia+1) 
      c1=d1 
      c2=3.d0*(a(ib+1)-a(ia+1))-d2-2.d0*d1 
      c3=d2+d1-2.d0*(a(ib+1)-a(ia+1)) 
      s(11)=c0+fr*(c1+fr*(c2+fr*c3)) 
      f1=a(ia+7)*(a(ia+9)-a(ia+10)*dpdt1) 
      f2=a(ib+7)*(a(ib+9)-a(ib+10)*dpdt2) 
      d1=f1*(tb-ta) 
      d2=f2*(tb-ta) 
      c0=a(ia+7) 
      c1=d1 
      c2=3.d0*(a(ib+7)-a(ia+7))-d2-2.d0*d1 
      c3=d2+d1-2.d0*(a(ib+7)-a(ia+7)) 
      s(12)=c0+fr*(c1+fr*(c2+fr*c3)) 
!                                                                       
!--two phase fluid                                                      
!                                                                       
      s(13)=a(ia+2)+(a(ib+2)-a(ia+2))*fr 
      s(14)=a(ia+8)+(a(ib+8)-a(ia+8))*fr 
      s(15)=s(13)+s(10)*s(11) 
      s(16)=s(14)+s(10)*s(12) 
      s(17)=a(ia+3)+(a(ib+3)-a(ia+3))*fr 
      s(18)=a(ia+9)+fr*tb/s(1)*(a(ib+9)-a(ia+9)) 
      s(19)=a(ia+4)+(a(ib+4)-a(ia+4))*fr 
      s(20)=a(ia+10)+(s(10)-pa)/(pb-pa)*pb/s(10)*(a(ib+10)-a(ia+10)) 
      s(21)=a(ia+5)+(a(ib+5)-a(ia+5))*fr 
      s(22)=a(ia+11)+(a(ib+11)-a(ia+11))*fr 
      s(25)=a(ia+6)+(a(ib+6)-a(ia+6))*fr 
      s(26)=a(ia+12)+(a(ib+12)-a(ia+12))*fr 
      fr=1.0d0-s(9) 
      s(3)=fr*s(11)+s(9)*s(12) 
      s(4)=fr*s(13)+s(9)*s(14) 
      s(5)=fr*s(15)+s(9)*s(16) 
      s(24)=fr*s(25)+s(9)*s(26) 
!if -def,in32,1                                                         
!      s(23) = or(shift(ip,30),jp)                                      
!if def,in32,3                                                          
      iunp(1)=ip 
      iunp(2)=jp 
      s(23)=unp 
      err=.false. 
!                                                                       
!--done                                                                 
   20 RETURN 
      END SUBROUTINE strx_cupid                           
