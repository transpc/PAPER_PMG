      SUBROUTINE std2x1_cupid(a,s,err) 
!
!deck std2x1                                                            
!                                                                       
!                           PROLOGUE                                    
!                                                                       
!   TITLE= std2x1 -- steam table interpolation on T and X or P and X    
!                                                                       
!   PURPOSE--                                                           
!                                                                       
!    main subroutine std2x1 (a,s,err)                                   
!   Compute water thermodynamic properties as a function of temperature 
!   and quality.                                                        
!      entry std2x2 (a,s,err)                                           
!   Compute water thermodynamic properties as a function of pressure    
!   and quality.                                                        
!                                                                       
!  ARGUMENT LIST                                                        
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
!     s(11)     =  vsubf   fluid specific volume  cubic meters / kg     
!     s(12)     =  vsubg   vapor specific volume  cubic meters / kg     
!     s(13)     =  usubf  fluid specific internal energy  J / kg        
!     s(14)     =  usubg  vapor specific internal energy  J / kg        
!     s(15)     =  hsubf  fluid specific enthalpy  J / kg               
!     s(16)     =  hsubg  vapor specific enthalpy  J / kg               
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
!     err        =  error flag                                          
!                                                                       
!  DIRECT OUTPUTS--(Arguments modified by this routine itself)          
!                                                                       
!      s(2)-s(8), s(10)-s(26), err                                      
!                                                                       
!  Compute water thermodynamic properties as a function of temperature  
!  and quality                                                          
!                                                                       
      IMPLICIT NONE
!
      INTEGER j,ip,jp,kp,jpp,kp2,ia,ib
      REAL(8) t,ta,tb,tc,pa,pb,prat,fr,fr1,hfg1,hfg2,dpdt1,dpdt2, &
              c0,c1,c2,c3,f1,f2,pp,dfdt
!
!  COMMON BLOCKS                                                        
!                                                                       
      INTEGER nt,np,ns,ns2,klp,klp2,llp,nt5,jpl
      COMMON/std2xc/nt,np,ns,ns2,klp,klp2,llp,nt5,jpl 
!
!                                                                       
!  DECLARATIONS                                                         
      REAL(8) a(1),s(26) 
      LOGICAL err 
      REAL(8) c(6),y,yh 
      REAL(8) pc,ps,b(5) 
      REAL(8) d1,d2 
      REAL(8) tsat,pr 
      LOGICAL s1,s2,s3 
!if def,in32,2                                                          
      INTEGER iunp(2) 
      REAL(8) unp 
      EQUIVALENCE(unp,iunp(1)) 
!                                                                       
!  DATA                                                                 
      REAL(8) plow,crp,tlow
      DATA b/-7.81583d0,17.6012d0,-18.1747d0,-3.92488d0,4.19174d0/ 
      DATA c/0.37228924d+03,0.88331901d+02,0.10970708d+02,              &
      0.97251708d+00,0.51713769d-01,0.12129545d-02/                     
      DATA plow/660.114d0/,crp/2.166d+07/ 
      DATA tlow/276.95d0/ 
!                                                                       
!  EXECUTION                                                            
      s(1)=max(s(1),tlow) 
      s1=.false. 
   15 IF(s(1).lt.a(1).or.s(1).gt.a(ns)) GOTO 101 
      IF(s(9).lt.0.d0.or.s(9).gt.1.d0) GOTO 101 
      IF(s1) GOTO 16 
      tc=643.89d0 
      pc=2.166d+07 
      t=s(1) 
      ta=1.d0-t/tc 
      IF(ta.le.0.d0)ta=0.d0 
      prat=(tc/t)*(b(1)*ta+b(2)*ta**1.9d0+b(3)*ta**2.d0+b(4)*ta**5.5d0+ &
      b(5)*ta**10.d0)                                                   
      ps=pc*exp(prat) 
      s(10)=ps 
      s(2)=s(10) 
      ENTRY std2xb_cupid(a,s,err) 
!if -def,in32,5                                                         
!if def,cray,1                                                          
!  16  ip = shiftr(s(23),30)                                            
!if -def,cray,1                                                         
!  16  ip = ishft(s(23),-30)                                            
!      jp = and(s(23),1023)                                             
!if def,in32,3                                                          
   16 unp=s(23) 
      ip=iunp(1) 
      jp=iunp(2) 
      s2=.false. 
      s3=.false. 
      IF(ip.le.0.or.ip.ge.ns)ip=1 
      IF(jp.le.0.or.jp.ge.ns2)jp=1 
!  Set indexes in temperature and pressure tables for saturation        
!  computations                                                         
   11 IF(s(1).ge.a(ip)) GOTO 10 
      ip=ip-1 
      GOTO 11 
   10 IF(s(1).le.a(ip+1)) GOTO 12 
      ip=ip+1 
      GOTO 10 
   12 jpp=jp+nt 
  111 IF(s(10).ge.a(jpp)) GOTO 110 
      jpp=jpp-1 
      IF(jpp.gt.nt) GOTO 111 
      s3=.true. 
      GOTO 112 
  110 IF(s(10).le.a(jpp+1)) GOTO 112 
      jpp=jpp+1 
      IF(jpp.lt.jpl) GOTO 110 
      s2=.true. 
  112 jp=jpp-nt 
      kp2=klp2+jp * 13 
      kp=klp+ip * 13 
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
!   two phase fluid.                                                    
      hfg1=a(ia+8)-a(ia+2)+pa*(a(ia+7)-a(ia+1)) 
      hfg2=a(ib+8)-a(ib+2)+pb*(a(ib+7)-a(ib+1)) 
      dpdt1=hfg1/(ta*(a(ia+7)-a(ia+1))) 
      dpdt2=hfg2/(tb*(a(ib+7)-a(ib+1))) 
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
!   two phase fluid.                                                    
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
      fr=1.d0-s(9) 
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
   20 RETURN 
      ENTRY std2x2_cupid(a,s,err) 
!   compute water thermodynamic properties as a function of pressure    
!   and quality                                                         
      s1=.true. 
! temporary patch to be able to do ice condenser debug runs             
      s(2)=max(s(2),plow) 
      s(10)=s(2) 
      IF(s(2).lt.plow.or.s(2).gt.crp) GOTO 101 
      tc=643.89d0 
      pc=2.166d+07 
      pp=s(2) 
      tsat=0.0d0 
      pr=(pp/21.671d6) 
      yh=log(pr) 
      y=(1.0d0) 
      DO 5 j=1,6 
         tsat=tsat+c(j)*y 
         y=y*yh 
    5 END DO 
      t=(tsat)+273.15d0 
      ta=1.d0-t/tc 
      IF(ta.le.0.d0)ta=0.d0 
      d1=(tc/t)*(b(1)*ta+b(2)*ta**1.9d0+b(3)*ta**2.d0+b(4)*ta**5.5d0+b( &
      5)*ta**10.d0)                                                     
      d2=b(1)+1.9d0*b(2)*ta**.9d0+2.d0*b(3)*ta**1.d0+5.5d0*b(4)*ta**    &
      4.5d0+10.d0*b(5)*ta**9.d0                                         
      pa=pc*exp(d1) 
      dfdt=(-1/t)*d1-(1/t)*d2 
      s(1)=t-(pa-pp)/(pa*dfdt) 
!  Patch because of troubles computing values at triple point.          
      s(1)=max(s(1),a(1)) 
      GOTO 15 
  101 err=.true. 
      GOTO 20 
      END SUBROUTINE std2x1_cupid                         
