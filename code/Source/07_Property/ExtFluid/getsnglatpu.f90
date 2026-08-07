      SUBROUTINE getsnglatpu_cupid(tables,arg,getprops,u,p,s,err) 
!define win32dvf                                                        
!define erf                                                             
!define fourbyt                                                         
!define hconden                                                         
!define impnon                                                          
!define in32                                                            
!define newnrc                                                          
!define ploc                                                            
!define sphaccm                                                         
!define unix                                                            
!define noselap                                                         
!define noextvol                                                        
!define noextv20                                                        
!define noextsys                                                        
!define noextjun                                                        
!define noextj20                                                        
!define noparcs                                                         
!define nonpa                                                           
!define nomap                                                           
!define logp                                                            
!deck getsnglatpu                                                       
!                                                                       
!  $Id: getsnglatpu.ff,v 1.5 2001/05/18 14:28:29 dbarber Exp dbarber $  
!                                                                       
!  Returns the single phase properties for subcritical liquid or        
!  vapor or for a supercritical fluid at a specified (p,u) point.       
!  This routine replaces sth2x6.F and sth2xf.F                          
!                                                                       
!  Cognizant engineer:  rwt                                             
!                                                                       
      IMPLICIT none 
!                                                                       
      INCLUDE 'stcom.h' 
      INCLUDE 'newstcom.h' 
      INCLUDE 'gibbpnt.h' 
      INCLUDE 'sparms.h' 
!mab+ May 30 2001  Add common block uatleft
      INCLUDE 'uatleft.h'
!mab-
!                                                                       
!      REAL(8) a(11,4),b(8),deltp,deltt,g(11),p,s(*),tables(*),tint,     &
!      tleft,tright,tuatleft,tuatright,u,uatleft,uatright,utatleft,      &
!      utatright,plo,phi,uleftlimitplo,uleftlimitphi,urightlimitplo,     &
!      urightlimitphi                                                    
      REAL(8) a(11,4),b(8),deltp,deltt,g(11),p,s(*),tables(*),tint,     &
      tleft,tright,tuatleft,tuatright,u,uatright,utatleft,              &
      utatright,plo,phi,uleftlimitplo,uleftlimitphi,urightlimitplo,     &
      urightlimitphi                                                    
      INTEGER i,iplo,iphi,otleftlimit,otrightlimit,ptableprop,          &
      ntablelimits,ptablelimits,ofirstpres,olastpres,uindex,            &
      itleftlimitplo,itleftplo,itrightplo,itrightlimitplo,              &
      itleftlimitphi,itleftphi,itrightphi,itrightlimitphi,itleft,       &
      itright,offsetu,ominu,omaxu,ntableprop,ptr,stableprop             
      LOGICAL err,getprops(*),leftshift,rightshift,threepointleft 
!rex+ 20 Mar 2001                                                       
      REAL(8) gtop(11),y 
!rex-                                                                   
      CHARACTER*(*)arg 
!gam                                                                    
      REAL(8) tintm5,tintp5,dtint 
      REAL(8) tintsave,uatleftsave 
!                                                                       
      PARAMETER (offsetu=1) 
      PARAMETER (ominu=3) 
      PARAMETER (omaxu=4) 
      PARAMETER (otleftlimit=5) 
      PARAMETER (otrightlimit=6) 
!gam                                                                    
      PARAMETER (dtint=2.0d0) 
!                                                                       
!      call timstart ('getsnglatpu')                                    
!                                                                       
      err=.false. 
      leftshift=.false. 
      rightshift=.false. 
      threepointleft=.false. 
!                                                                       
!  Find the table and set generic pointers                              
!                                                                       
      CALL setpointers_cupid(arg,ofirstpres,olastpres,ptableprop,stableprop,  &
      ntableprop,ptablelimits,ntablelimits)                             
!                                                                       
!  Find plo, phi, iplo, and iphi                                        
!                                                                       
      CALL getindex_cupid(tables,p,ptable2,ofirstpres,olastpres,iplo,err) 
!                                                                       
      IF(err)stop 1 
!                                                                       
!   Load bounding pressure indices                                      
!                                                                       
      iphi=iplo+1 
      plo=tables(ptable2+iplo) 
      phi=tables(ptable2+iphi) 
      deltp=phi-plo 
!                                                                       
!  Find bounding u indices on plo line from Table 8, 9, or 10           
!                                                                       
      uleftlimitplo=tables(ptablelimits+(iplo-ofirstpres)*ntablelimits+ &
      ominu)                                                            
      urightlimitplo=tables(ptablelimits+(iplo-ofirstpres)*ntablelimits+&
      omaxu)                                                            
      itleftlimitplo=int(tables(ptablelimits+(iplo-ofirstpres)*         &
      ntablelimits+otleftlimit))                                        
      itrightlimitplo=int(tables(ptablelimits+(iplo-ofirstpres)*        &
      ntablelimits+otrightlimit))                                       
!                                                                       
      IF(u.lt.uleftlimitplo)then 
         itleftplo=itleftlimitplo 
         itrightplo=itleftplo+1 
      ELSEIF(u.gt.urightlimitplo)then 
         itleftplo=itrightlimitplo-1 
         itrightplo=itleftplo+1 
      ELSE 
         DO 10 i=itleftlimitplo+1,itrightlimitplo 
            uindex=ptableprop+(iplo-ofirstpres)*stableprop+(i-1)*       &
            ntableprop+offsetu                                          
            IF(tables(uindex).ge.u)then 
               itleftplo=i-1 
               itrightplo=i 
               GOTO 20 
            ENDIF 
   10    END DO 
         err=.true. 
         GOTO 999 
   20    CONTINUE 
      ENDIF 
!                                                                       
!  Find bounding u indices on phi line from Table 8, 9, or 10           
!                                                                       
      uleftlimitphi=tables(ptablelimits+(iphi-ofirstpres)*ntablelimits+ &
      ominu)                                                            
      urightlimitphi=tables(ptablelimits+(iphi-ofirstpres)*ntablelimits+&
      omaxu)                                                            
      itleftlimitphi=int(tables(ptablelimits+(iphi-ofirstpres)*         &
      ntablelimits+otleftlimit))                                        
      itrightlimitphi=int(tables(ptablelimits+(iphi-ofirstpres)*        &
      ntablelimits+otrightlimit))                                       
!                                                                       
      IF(u.lt.uleftlimitphi)then 
         itleftphi=itleftlimitphi 
         itrightphi=itleftphi+1 
      ELSEIF(u.gt.urightlimitphi)then 
         itleftphi=itrightlimitphi-1 
         itrightphi=itleftphi+1 
      ELSE 
         DO 30 i=itleftlimitphi+1,itrightlimitphi 
            uindex=ptableprop+(iphi-ofirstpres)*stableprop+(i-1)*       &
            ntableprop+offsetu                                          
            IF(tables(uindex).ge.u)then 
               itleftphi=i-1 
               itrightphi=i 
               GOTO 40 
            ENDIF 
   30    END DO 
         err=.true. 
         GOTO 999 
   40    CONTINUE 
      ENDIF 
!                                                                       
!  Move temp indices to make a rectangle if a rhombus                   
!                                                                       
      CALL setanchor_cupid(itleftphi,itleftplo,itrightplo,itrightphi,         &
      itleftlimitplo,itrightlimitplo,itleftlimitphi,itrightlimitphi,p,  &
      plo,phi,itleft,itright,threepointleft,err)                        
!                                                                       
      IF(err) GOTO 999 
!                                                                       
  888 tleft=tables(ptable1+itleft) 
      tright=tables(ptable1+itright) 
      deltt=tright-tleft 
!                                                                       
!  Calculate u at the endpoints of the input p line                     
!  using the 4-point formula                                            
!                                                                       
      CALL loadcorners_cupid(tables,a,itleft,itright,iplo,ptableprop,         &
      stableprop,ntableprop,ofirstpres,ptr)                             
!                                                                       
!rex+ 9 Feb 2001; if high press and threepoint or liquid                
      IF(p.gt.21.70d6.and.threepointleft)then 
         CALL herm2d_cupid(tables,ptr,a,arg,gtop,tright,phi,err) 
!                                                                       
         y=(p-plo)/deltp 
!        now get left point                                             
!                                                                       
         CALL herm2d_cupid(tables,ptr,a,arg,g,tleft,plo,err) 
!                                                                       
         g(1)=g(1)+(gtop(1)-g(1))*y 
         g(2)=g(2)+(gtop(2)-g(2))*y 
         g(3)=g(3)+(gtop(3)-g(3))*y 
         g(4)=g(4)+(gtop(4)-g(4))*y 
         g(5)=g(5)+(gtop(5)-g(5))*y 
         g(7)=g(7)+(gtop(7)-g(7))*y 
         tint=tleft+(tright-tleft)*y 
!                                                                       
         GOTO 1020 
      ENDIF 
!rex-                                                                   
!  Recalculate left hand side temperature if we have a triangle         
!rex+ 28 Mar 2001 add another restriction on next if test               
!       if (threepointleft .and. (u .gt. uleftlimitplo)) then           
      IF(threepointleft)then 
         tleft=tleft+(tright-tleft)*(p-plo)/(phi-plo) 
      ENDIF 
!rex-                                                                   
!                                                                       
!  Interpolate for the value of u at left point                         
!                                                                       
      CALL herm2dplus_cupid(tables,ptr,a,arg,g,tleft,p,err) 
!                                                                       
!  Calculate single-phase liquid properties using the elements of g     
!                                                                       
      uatleft=g(1)-tleft*g(4)-p*g(2) 
      utatleft=-tleft*g(7)-p*g(5) 
!      uttatleft = -tleft*g(10) - g(7) - p*g(8)                         
      tuatleft=1.0d0/utatleft 
!      tuuatleft = - uttatleft/utatleft**3                              
!                                                                       
!  Interpolate for the value of u at right point                        
!                                                                       
      CALL herm2dplus_cupid(tables,ptr,a,arg,g,tright,p,err) 
!                                                                       
!  Calculate single-phase liquid properties using the elements of g     
!                                                                       
      uatright=g(1)-tright*g(4)-p*g(2) 
      utatright=-tright*g(7)-p*g(5) 
!      uttatright = -tright*g(10) - g(7) - p*g(8)                       
      tuatright=1.0d0/utatright 
!      tuuatright = -uttatright/utatright**3                            
!                                                                       
!  Use the square we have been working in?                              
!                                                                       
!  if at the leftlimit, do not shift over                               
!                                                                       
!rex+ 23 Mar 2001 do not shift if not threepointleft                    
      IF((u.lt.uatleft).and.(.not.rightshift).and.(.not.threepointleft))&
      then                                                              
!rex-                                                                   
         itleft=itleft-1 
         itright=itleft+1 
         leftshift=.true. 
         IF(itleft.lt.itleftlimitplo)then 
            itleft=itleft+1 
            itright=itleft+1 
            GOTO 777 
         ELSEIF(itleft.lt.itleftlimitphi)then 
            threepointleft=.true. 
            itright=itleftlimitphi 
         ENDIF 
         GOTO 888 
      ELSEIF((u.gt.uatright).and.(.not.leftshift))then 
         itleft=itleft+1 
         itright=itleft+1 
         rightshift=.true. 
         IF((itright.gt.itrightlimitplo).or.(itright.gt.itrightlimitphi)&
         )then                                                          
            itleft=itleft-1 
            itright=itleft+1 
            GOTO 777 
         ENDIF 
         IF(threepointleft)then 
            itright=itleftlimitphi 
            IF(itleft.eq.itleftlimitphi)then 
               itright=itleft+1 
               threepointleft=.false. 
            ENDIF 
         ENDIF 
         GOTO 888 
      ENDIF 
!                                                                       
  777 b(1)=uatleft 
      b(2)=tleft 
      b(3)=tuatleft 
!      b(4) = tuuatleft                                                 
      b(5)=uatright 
      b(6)=tright 
      b(7)=tuatright 
!      b(8) = tuuatright                                                
!                                                                       
      CALL cubic_cupid(b,u,tint,err) 
!                                                                       
      IF((u.gt.uatright).and.(tint.lt.tright))then 
!rex+ 1 Mar 2001 TDV for case pcrit.i above 21 MPa had error            
!       but timestep reductions do not occur from tstate                
!       plus answers were fine if it continued.                         
         IF(p.gt.20.0d6)then 
            tint=tright 
            err=.false. 
            CALL herm2d_cupid(tables,ptr,a,arg,g,tint,p,err) 
            GOTO 1020 
         ELSE 
            err=.true. 
            GOTO 999 
         ENDIF 
!rex-                                                                   
      ENDIF 
!                                                                       
!gam  use linear extripolation instead of cubic fit if cubic fails      
      IF(u.lt.uatleft.and.tint.gt.tleft)then 
         tintsave=tint 
         tint=tleft+(u-uatleft)*(tright-tleft)/(uatright-uatleft) 
      ENDIF 
!                                                                       
!   Calculate single-phase properties using interpolated temperature and
!   the same 4 corner points already loaded into the a array            
!                                                                       
  779 CONTINUE 
!rex+ 8 Mar 2001                                                        
      IF(tint.lt.tleft.and.u.gt.uatleft)then 
         tint=tleft+(u-uatleft)/(uatright-uatleft)*(tright-tleft) 
      ENDIF 
      IF(tint.gt.tright.and.u.lt.uatright)then 
         tint=tleft+(u-uatleft)/(uatright-uatleft)*(tright-tleft) 
      ENDIF 
!rex-                                                                   
      IF(tint.lt.tleft-deltt)then 
!gam  get a better estimate for tint                                    
         uatleftsave=uatleft 
         CALL herm2dleft_cupid(tables,ptr,a,arg,g,tint,p,err) 
         uatleft=g(1)-tint*g(4)-p*g(2) 
         utatleft=-tint*g(7)-p*g(5) 
         tuatleft=1.0d0/utatleft 
!  now do a linear interpolation between this value and the tleft value 
         tint=tint+(u-uatleft)*(tleft-tint)/(uatleftsave-uatleft) 
!  now home in on it by getting u at points tint+dtint and tint-dtint   
!rex+ 27 Mar 2001 dtint to big, use deltt                               
!          tintm5 = tint - dtint                                        
         tintm5=tint-deltt 
!rex-                                                                   
         CALL herm2dleft_cupid(tables,ptr,a,arg,g,tintm5,p,err) 
         uatleft=g(1)-tintm5*g(4)-p*g(2) 
         utatleft=-tintm5*g(7)-p*g(5) 
         tuatleft=1.0d0/utatleft 
!                                                                       
!rex+ 27 Mar 2001 dtint to big, use deltt                               
!          tintp5 = tint + dtint                                        
         tintp5=tint+deltt 
!rex-                                                                   
         CALL herm2dleft_cupid(tables,ptr,a,arg,g,tintp5,p,err) 
         uatright=g(1)-tintp5*g(4)-p*g(2) 
         utatright=-tintp5*g(7)-p*g(5) 
         tuatright=1.0d0/utatright 
!                                                                       
         b(1)=uatleft 
         b(2)=tintm5 
         b(3)=tuatleft 
         b(5)=uatright 
         b(6)=tintp5 
         b(7)=tuatright 
         CALL cubic_cupid(b,u,tint,err) 
!gam  now get the properties at this temperature                        
         CALL herm2dleft_cupid(tables,ptr,a,arg,g,tint,p,err) 
      ELSE 
         CALL herm2d_cupid(tables,ptr,a,arg,g,tint,p,err) 
      ENDIF 
!                                                                       
!                                                                       
!  Calculate single-phase properties using the returned                 
!  elements of g                                                        
!                                                                       
!rex+ 9 Feb 2001                                                        
 1020 CONTINUE 
      IF(g(2).lt.0.0d0)then 
         CALL rholine_cupid(tables,ptr,a,arg,g,tint,p,err) 
      ENDIF 
      IF(g(3).gt.0.0d0)then 
         CALL kappaline_cupid(tables,ptr,a,arg,g,tint,p,err) 
      ENDIF 
      IF(g(5).lt.0.0d0)then 
         CALL betaline_cupid(tables,ptr,a,arg,g,tint,p,err) 
         IF(g(5).lt.0.0d0)then 
            IF(tint.gt.tright)then 
               CALL betaline_cupid(tables,ptr,a,arg,g,tright,p,err) 
            ELSEIF(tint.lt.tleft)then 
               CALL betaline_cupid(tables,ptr,a,arg,g,tleft,p,err) 
            ENDIF 
         ENDIF 
         IF(g(5).lt.0.0d0)then 
            CALL beta2d_cupid(tables,ptr,a,arg,g,tint,p,err) 
         ENDIF 
         IF(g(5).lt.0.0d0)then 
            CALL beta2d_cupid(tables,ptr,a,arg,g,tright,phi,err) 
         ENDIF 
      ENDIF 
      IF(g(7).gt.0.0d0)then 
         CALL cpline_cupid(tables,ptr,a,arg,g,tint,p,err) 
      ENDIF 
!rex-                                                                   
      IF(getprops(temp))then 
         s(temp)=tint 
      ENDIF 
      IF(getprops(vbar))then 
         s(vbar)=g(2) 
      ENDIF 
      IF(getprops(entpy))then 
         s(entpy)=-g(4) 
      ENDIF 
      IF(getprops(kapa))then 
         s(kapa)=-g(3)/g(2) 
      ENDIF 
      IF(getprops(beta))then 
         s(beta)=g(5)/g(2) 
         s(beta)=min(s(beta),0.05d0) 
      ENDIF 
      IF(getprops(cp))then 
         s(cp)=-tint*g(7) 
         s(cp)=min(s(cp),30000.0d0) 
      ENDIF 
      IF(getprops(hbar))then 
         s(hbar)=g(1)-tint*g(4) 
      ENDIF 
      IF(getprops(qual))then 
         IF(arg.eq.'vapor')then 
            s(qual)=1.0d0 
         ELSEIF(arg.eq.'liquid')then 
            s(qual)=0.0d0 
!gam  do not return a value for s(qual) when supercritical              
         ELSEIF(arg.eq.'supercritical')then 
!gam  set sat properties equal to single phase properties               
!gam  and tsat equal to the temperature                                 
            s(tsat)=s(temp) 
            s(vsubf)=s(vbar) 
            s(vsubg)=s(vbar) 
            s(usubf)=s(ubar) 
            s(usubg)=s(ubar) 
            s(hsubf)=s(hbar) 
            s(hsubg)=s(hbar) 
            s(betaf)=s(beta) 
            s(betag)=s(beta) 
            s(kapaf)=s(kapa) 
            s(kapag)=s(kapa) 
            s(cpf)=s(cp) 
            s(cpg)=s(cp) 
            s(entpyf)=s(entpy) 
            s(entpyg)=s(entpy) 
            s(qual)=0.0d0 
         ENDIF 
      ENDIF 
!                                                                       
!  Check interpolated values                                            
!                                                                       
      CALL checkvalue_cupid(arg,p,g,s,getprops,tleft,tright,deltt,err) 
!                                                                       
  999 CONTINUE 
!                                                                       
!      call timstop ('getsnglatpu')                                     
!                                                                       
      RETURN 
      END SUBROUTINE getsnglatpu_cupid                    
