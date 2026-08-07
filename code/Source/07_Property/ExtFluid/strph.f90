      SUBROUTINE strph_cupid(a,s,it,err) 
!
!   compute gas thermodynamic properties as a function of pressure    
!   and total enthalpy                                                  
!                 
!   by Won-jae Lee for gas property calculation in MARS 3D module
!                                                                       
!  ARGUMENT LIST                                                        
!                                                                       
!     a         =   general parameter vector                            
!     s(1)      =   T  temperature  deg K                               
!     s(2)      =   P  pressure     Pa                                  
!     s(3)      =  v  specific volume   cubic meter / kg                
!     s(4)      =  u  specific internal energy  J / kg                  
!     s(5)      =  h  specific enthalpy  J / kg                         
!     s(6)      =  beta  thermal coefficient of expansion  vol / (vol - 
!     s(7)      =  kappa pressure coefficient of expansion  vol / (vol -
!     s(8)      =  Csubp  specific heat  J / (kg - deg K)               
!     s(9)      =  quality X                                            
!     s(10)     =  Psat  or  Tsat     Pa or deg K                       
!     s(11)      = vsubf   fluid specific volume  cubic meters / kg     
!     s(12)      = vsubg   vapor specific volume  cubic meters / kg     
!     s(13)      =  usubf  fluid specific internal energy  J / kg       
!     s(14)      =  usubg  vapor specific internal energy  J / kg       
!     s(15)      =  hsubf  fluid specific enthalpy  J / kg              
!     s(16)      =  hsubg  vapor specific enthalpy  J / kg              
!     s(17)      =  betaf liquid  thermal coefficient of expansion  vol 
!     s(18)      =  betag vapor  thermal coefficient of expansion  vol /
!     s(19)      =  kappaf liquid pressure coefficient of expansion  vol
!     s(20)      =  kappag vapor pressure coefficient of expansion  vol 
!     s(21)      =  Csubpf liquid  specific heat  J / (kg - deg K)      
!     s(22)      =  Csubpg vapor  specific heat  J / (kg - deg K)       
!     s(23)      =  indexes                                             
!     s(24)      =  s   specific entropy  J /  ( kg - deg K)            
!     s(25)      =  ssubf liquid  specific entropy  J /  ( kg - deg K)  
!     s(26)      =  ssubg vapor  specific entropy  J /  ( kg - deg K)   
!     it         =  phase flag                                          
!     err        =  error flag                                          
!                                                                       
!  DIRECT OUTPUTS--(Arguments modified by this routine itself)          
!                                                                       
!    s(1), s(3), s(5), s(4) - s(26), it, err                                  
!
      USE STM_TBL_cupid  , ONLY: nt,np,pcrit
!
      IMPLICIT NONE
!
      REAL(8) a(1),s(26) 
      LOGICAL err 
!     INCLUDE 'machds.h'
      LOGICAL s1,s2,s3 
      INTEGER :: ip,jp,jpp,kp,kp2
!
      INTEGER :: it,ia,ib,ic,id,lpp,lqq
      REAL(8) :: fr,fr1,fr2,frc,frc2,frc3,frd,frd2,frd3,frn
      REAL(8) :: c0,c1,c2,c3,cv
      REAL(8) :: d2,f1,f2
      REAL(8) :: r2,ren
      REAL(8) :: ta,tb,pa,pb,hfg1,hfg2,dpdt1,dpdt2,dpsdts
      REAL(8) :: dum,ht,ut,vt,pr,px,tr,s11i,s12i,dpdtcr
!
      INTEGER iunp(2) 
      REAL(8) :: unp
      EQUIVALENCE(unp,iunp(1)) 
!------------------------------------------------------------------------------
!   Modified to add Saturation Properties on Pressure: Won-Jae Lee on Dec. '03      
      REAL(8) vfsat,vgsat,ufsat,ugsat,hfsat,hgsat,efsat,egsat,    &
	      kfsat,kgsat,cpfsat,cpgsat,entfsat,entgsat
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
      s3=.false. 
!   check for valid input                                               
!   23 IF(s(2).le.0.0d0.or.s(2).gt.a(it3p0g)) GOTO 1001 
      IF(s(2).le.0.0d0.or.s(2).gt.a(it3p0g)) GOTO 1001 
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
!   set indexes in temperature and pressure tables for saturation       
!   computations                                                        
   11 IF(s(2).ge.a(jpp)) GOTO 10 
      jpp=jpp-1 
      IF(jpp.gt.ntg) GOTO 11 
      jpp=jpp+1 
      jp=1 
      s1=.true. 
      IF(s(2).lt.a(it3bpg+13)) GOTO 44 
      GOTO 12 
   10 IF(s(2).lt.a(jpp+1)) GOTO 12 
      jpp=jpp+1 
      GOTO 10 
   12 jp=jpp-ntg 
      IF(s(2).ge.pcrit) GOTO 44 
      IF(s3) GOTO 15 
! ------------------------------------------------------------------------
!     GAS Modification by Won-Jae Lee
      CALL strsat_cupid(a,2,s(2),s(10),dpsdts,err) 
      IF(err) GOTO 1001 
! ------------------------------------------------------------------------
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
!   compute vsubf and vsubg to determine liquid, two phase, or vapor    
!   state                                                               
      hfg1=a(ia+8)-a(ia+2)+pa*(a(ia+7)-a(ia+1)) 
      hfg2=a(ib+8)-a(ib+2)+pb*(a(ib+7)-a(ib+1)) 
      dpdt1=hfg1/(ta*(a(ia+7)-a(ia+1))) 
! ------------------------------------------------------------------------
!      dpdt2=hfg2/(tb*(a(ib+7)-a(ib+1))) 
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
! ------------------------------------------------------------------------
      f1=a(ia+1)*(a(ia+3)-a(ia+4)*dpdt1) 
      f2=a(ib+1)*(a(ib+3)-a(ib+4)*dpdt2) 
      c1=f1*(tb-ta) 
      d2=f2*(tb-ta) 
      c0=a(ia+1) 
      c2=3.d0*(a(ib+1)-c0)-d2-2.d0*c1 
      c3=d2+c1-2.d0*(a(ib+1)-c0) 
!------------------------------------------------------------------------------
!  Mod. by Won-Jae Lee
      s(11)=c0+fr*(c1+fr*(c2+fr*c3)) 
   24 s(13)=a(ia+2)+(a(ib+2)-a(ia+2))*fr 
      s(15)=s(13)+s(2)*s(11) 
!
	  vfsat=s(11)   ! wjlee
	  ufsat=s(13)   ! wjlee
	  hfsat=s(15)   ! wjlee
!
!      IF(s(5).le.s(15)) GOTO 41 
      f1=a(ia+7)*(a(ia+9)-a(ia+10)*dpdt1) 
      f2=a(ib+7)*(a(ib+9)-a(ib+10)*dpdt2) 
      c1=f1*(tb-ta) 
      d2=f2*(tb-ta) 
      c0=a(ia+7) 
      c2=3.d0*(a(ib+7)-c0)-d2-2.d0*c1 
      c3=d2+c1-2.d0*(a(ib+7)-c0) 
      s(12)=c0+fr*(c1+fr*(c2+fr*c3)) 
      s(14)=a(ia+8)+(a(ib+8)-a(ia+8))*fr 
      s(16)=s(14)+s(2)*s(12) 
!
      vgsat=s(12)   ! wjlee
	  ugsat=s(14)   ! wjlee
	  hgsat=s(16)   ! wjlee
!
!      IF(s(5).ge.s(16)) GOTO 43 
!   two phase fluid.                                                    
!      it=2 
      s(17)=a(ia+3)+(a(ib+3)-a(ia+3))*fr 
      s(18)=a(ia+9)+fr*tb/s(10)*(a(ib+9)-a(ia+9)) 
      s(19)=a(ia+4)+(a(ib+4)-a(ia+4))*fr 
      s(20)=a(ia+10)+(s(2)-pa)/(pb-pa)*pb/s(2)*(a(ib+10)-a(ia+10)) 
      s(21)=a(ia+5)+(a(ib+5)-a(ia+5))*fr 
      s(22)=a(ia+11)+(a(ib+11)-a(ia+11))*fr 
      s(25)=a(ia+6)+(a(ib+6)-a(ia+6))*fr 
      s(26)=a(ia+12)+(a(ib+12)-a(ia+12))*fr 
!
      efsat=s(17)
	  egsat=s(18)
	  kfsat=s(19)
	  kgsat=s(20)
	  cpfsat=s(21)
	  cpgsat=s(22)
      entfsat=s(25)
	  entgsat=s(26)
!
      IF(s(5).le.s(15)) GOTO 41 
      IF(s(5).ge.s(16)) GOTO 43 
!
!   two phase fluid.                                                    
!------------------------------------------------------------------------------
      it=2 
      s(9)=(s(5)-s(15))/(s(16)-s(15)) 
      fr=1.0d0-s(9) 
      s(1)=s(10) 
      s(3)=fr*s(11)+s(9)*s(12) 
      s(4)=fr*s(13)+s(9)*s(14) 
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
!   single phase fluid, search for single phase indexes.                
   41 it=1 
      IF(s1) GOTO 1001 
      GOTO 50 
   43 it=3 
   50 lpp=it5bpg+jp*nprpntg+ip * 6 
   51 ht=a(lpp+1)+a(jpp)*a(lpp) 
      IF(s(5).ge.ht) GOTO 52 
      lpp=lpp-6 
      ip=ip-1 
      IF(ip.gt.0) GOTO 51 
      lpp=lpp+6 
      ip=ip+1 
      GOTO 54 
   52 ht=a(lpp+7)+a(jpp)*a(lpp+6) 
      IF(s(5).le.ht) GOTO 54 
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
!   liquid phase.                                                       
   58 ht=s(15) 
      ut=s(13) 
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
      GOTO 62 
   61 ta=a(ip+1) 
      s2=.true. 
      IF(ip+1.gt.nstg) GOTO 63 
      IF(a(kp+13).le.a(jpp)) GOTO 63 
      frc3=s(2)-a(kp+13) 
      px=a(kp+13) 
      frc=a(jpp+1)-a(kp+13) 
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
   82 fr1=(s(11)-a(ic))/(a(lqq+6)-a(ic)) 
      s(13)=a(ic+1)+(a(lqq+7)-a(ic+1))*fr1 
      s(15)=s(13)+s(2)*s(11) 
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
   66 IF(s3) GOTO 68 
      pr=(s(2)-px)/(a(jpp+1)-px) 
      c0=1.d0/a(id) 
      r2=1.d0/a(lqq) 
      c1=c0*a(id+3)*(a(jpp+1)-px) 
      d2=r2*a(lqq+3)*(a(jpp+1)-px) 
      c2=3.d0*(r2-c0)-d2-2.d0*c1 
      c3=d2+c1-2.d0*(r2-c0) 
      s12i=c0+pr*(c1+pr*(c2+pr*c3)) 
      s(12)=1.d0/s12i 
  182 fr2=(s(12)-a(id))/(a(lqq)-a(id)) 
      s(14)=a(id+1)+(a(lqq+1)-a(id+1))*fr2 
      s(16)=s(14)+s(2)*s(12) 
      IF(s(16).le.s(5)) GOTO 68 
      s2=.true. 
      s(11)=s(12) 
      s(13)=s(14) 
      s(15)=s(16) 
      fr1=fr2 
      ip=ip-1 
      ta=a(ip+1) 
      IF(ip.le.0) GOTO 1001 
      kp=kp-13 
      lqq=lqq-6 
      lpp=lpp-6 
      ic=id 
      GOTO 62 
   68 IF(s(15).ge.s(5)) GOTO 59 
      s(12)=s(11) 
      s(11)=vt 
      s(14)=s(13) 
      s(13)=ut 
      s(16)=s(15) 
      s(15)=ht 
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
      fr=(s(5)-s(16))/(s(15)-s(16)) 
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
      s(4)=s(5)-s(2)*s(3) 
      s(6)=(c1+tr*(2.d0*c2+3.d0*tr*c3))/(frd*s(3)) 
      s(7)=s(20)+(s(19)-s(20))*fr1 
      s(8)=s(22)+(s(21)-s(22))*fr1 
      s(24)=s(26)+(s(25)-s(26))*fr1 
      s(9)=0.d0 
!      GOTO 25 
      GOTO 29   ! Won-Jae Lee 
   70 ht=s(16) 
      ut=s(14) 
      vt=s(12) 
!   vapor phase.                                                        
  160 IF(a(ip+1).gt.s(10)) GOTO 157 
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
      GOTO 162 
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
      s(16)=s(14)+s(2)*s(12) 
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
  166 IF(s3) GOTO 168 
      fr1=a(id)*frd3 
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
      s(15)=s(13)+s(2)*s(11) 
      IF(s(15).ge.s(5)) GOTO 168 
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
      s(16)=s(15) 
      frc2=frd2 
      frc=frd 
      GOTO 162 
  168 IF(s(16).le.s(5)) GOTO 159 
      s3=.true. 
      s(11)=s(12) 
      s(12)=vt 
      s(13)=s(14) 
      s(14)=ut 
      s(15)=s(16) 
      s(16)=ht 
      frd2=frc2 
      frd=frc 
      ip=ip-1 
      IF(ip.le.0) GOTO 1001 
      kp=kp-13 
      lpp=lpp-6 
      lqq=lqq-6 
      IF(ip.lt.nstg) GOTO 57 
      GOTO 157 
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
      fr=(s(5)-s(16))/(s(15)-s(16)) 
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
      s(4)=s(5)-s(2)*s(3) 
      s(6)=(c1+tr*(2.d0*c2+3.d0*tr*c3))/(frd*s(3)) 
      s(7)=s(20)+(s(19)-s(20))*fr1 
      s(8)=s(22)+(s(21)-s(22))*fr1 
      s(24)=s(26)+(s(25)-s(26))*fr1 
   99 s(9)=1.0d0 
!      GOTO 25 
      GOTO 29 
!------------------------------------------------------------------------------
!  Mod. by Won-Jae Lee
!  Additional Saturated Properties based on Pressure
   29 CONTINUE
      s(11)=vfsat
	  s(12)=vgsat
	  s(13)=ufsat
	  s(14)=ugsat
	  s(15)=hfsat
	  s(16)=hgsat
	  s(17)=efsat
	  s(18)=egsat
	  s(19)=kfsat
	  s(20)=kgsat
	  s(21)=cpfsat
	  s(22)=cpgsat
	  s(25)=entfsat
	  s(26)=entgsat
      GOTO 25 
!------------------------------------------------------------------------------
!   vapor phase, temperature greater than highest table temperature.    
   80 fr=a(lqq)*(a(jpp+1)-a(jpp)) 
      s(3)=a(lpp)*fr/(fr+(a(lpp)-a(lqq))*frn) 
      fr=(s(3)-a(lpp))/(a(lqq)-a(lpp)) 
      frc=fr*a(lqq)/s(3) 
      ht=a(lpp+1)+(a(lqq+1)-a(lpp+1))*fr+s(2)*s(3) 
      s(8)=a(lpp+4)+frc*(a(lqq+4)-a(lpp+4)) 
      s(1)=(s(5)-ht+s(8)*a(ntg))/s(8) 
      frd=s(1)/a(ntg) 
      s(3)=s(3)*frd 
      s(4)=s(5)-s(2)*s(3) 
      s(6)=(a(lpp+2)+frc*(a(lqq+2)-a(lpp+2)))/frd 
      s(7)=a(lpp+3)+(a(lqq+3)-a(lpp+3))*fr 
      s(24)=a(lpp+5)+(a(lqq+5)-a(lpp+5))*fr 
      cv=s(8)-a(ntg)*s(6)*s(6)*s(3)/s(7) 
      s(24)=s(24)+cv*log(frd**(s(8)/cv)) 
      GOTO 99 
!   vapor phase, pressure less than lowest table pressure               
   96 ht=a(it4bpg+21)+a(ntg+1)*a(it4bpg+20) 
      IF(s(5).lt.ht) GOTO 90 
      IF(a(ip).lt.a(it4bpg+13)) GOTO 93 
      ia=ip 
      lqq=lpp 
      ht=a(lqq+1)+a(ntg+1)*a(lqq) 
      GOTO 91 
   93 ia=it4bpg+13 
      lqq=ia+7 
   91 fr=(s(5)-ht)/(a(lpp+7)+a(ntg+1)*a(lpp+6)-ht) 
      fr1=a(ip+1)-a(ia) 
      fr2=fr1*fr 
      s(1)=a(ia)+fr2 
      fr1=fr2/fr1 
      s(3)=(fr1*a(lpp+6)/a(ip+1)+(1.0d0-fr1)*a(lqq)/a(ia))*a(ntg+1)*s(1)/&
      s(2)                                                              
      s(6)=a(lqq+2)+(a(lpp+8)-a(lqq+2))*fr*a(ip+1)/s(1) 
      s(8)=a(lqq+4)+(a(lpp+10)-a(lqq+4))*fr 
      s(24)=a(lqq+5)+(a(lpp+11)-a(lqq+5))*fr 
      ren=s(2)*s(3)/s(1) 
      s(24)=s(24)-ren*log(s(2)/a(ntg+1)) 
   94 s(7)=1.0d0/s(2) 
   98 s(4)=s(5)-s(2)*s(3) 
      GOTO 99 
   95 s(8)=a(lpp+4) 
      s(1)=(s(5)-a(lpp+1)-a(ntg+1)*a(lpp)+s(8)*a(ntg))/s(8) 
      frd=s(1)/a(ntg) 
      frc=a(ntg+1)*a(lpp) 
      s(3)=frc*frd/s(2) 
      s(6)=a(lpp+2)/frd 
      s(7)=a(lpp+3)*a(ntg+1)/s(2) 
      ren=s(2)*s(3)/s(1) 
      IF(s(1).le.0.0d0) GOTO 1001 
      s(24)=a(lpp+5)+s(8)*log(s(1)/a(ntg))-ren*log(s(2)/a(ntg+1)) 
      GOTO 98 
   90 ht=a(it3bpg+21)+a(it3bpg+13)*a(it3bpg+20) 
      IF(s(5).lt.ht) GOTO 1001 
      ip=1 
      kp=it3bpg+13 
  202 frd=a(kp+21)+a(kp+13)*a(kp+20) 
      IF(s(5).le.frd) GOTO 201 
      ip=ip+1 
      kp=kp+13 
      ht=frd 
      GOTO 202 
  201 fr=(s(5)-ht)/(frd-ht) 
      s(1)=a(ip)+fr*(a(ip+1)-a(ip)) 
      s(3)=(fr*a(kp+13)*a(kp+20)/a(ip+1)+(1.0d0-fr)*a(kp)*a(kp+7)/a(ip))&
      *s(1)/s(2)                                                        
      s(6)=1.0d0/s(1) 
      s(8)=a(kp+11)+(a(kp+24)-a(kp+11))*fr 
      ren=s(2)*s(3)/s(1) 
      IF(s(1).le.0.0d0) GOTO 1001 
      s(24)=a(kp+12)+s(8)*log(s(1)/a(ip))-ren*log(s(2)/a(kp)) 
      GOTO 94 
 1001 err=.true. 
      RETURN 
!      ENTRY std2xe(a,s,it,err) 
!      s3=.true. 
!      GOTO 23 
      END SUBROUTINE strph_cupid                         
