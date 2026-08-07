      SUBROUTINE getsnglatpt_cupid(tables,arg,getprops,t,p,s,err) 
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
!deck getsnglatpt                                                       
!                                                                       
!  $Id: getsnglatpt.ff,v 1.4 2001/04/25 14:55:49 dbarber Exp dbarber $  
!                                                                       
!  Returns the single phase properties for subcritical liquid or        
!  vapor or for a supercritical fluid at a specified (p,t) point.       
!  This routine replaces sth2x3.F                                       
!                                                                       
!  Cognizant engineer:  rwt                                             
!                                                                       
      IMPLICIT none 
!                                                                       
      INCLUDE 'stcom.h' 
      INCLUDE 'newstcom.h' 
      INCLUDE 'gibbpnt.h' 
      INCLUDE 'sparms.h' 
!                                                                       
      REAL(8) a(11,4),deltt,g(11),p,plo,phi,s(*),t,tables(*),tleft,     &
      tleftlimitplo,tleftlimitphi,tright,trightlimitplo,trightlimitphi  
      INTEGER i,iplo,iphi,otleftlimit,otrightlimit,ptableprop,          &
      ntablelimits,ptablelimits,ofirstpres,olastpres,itleftlimitplo,    &
      itleftplo,itrightplo,itrightlimitplo,stableprop,offsetu,omint,    &
      omaxt,ntableprop,itleftlimitphi,itleftphi,itrightphi,             &
      itrightlimitphi,itleft,itright,ptr                                
      LOGICAL err,getprops(*),threepointleft 
      CHARACTER*(*)arg 
!                                                                       
      PARAMETER (offsetu=1) 
      PARAMETER (omint=1) 
      PARAMETER (omaxt=2) 
      PARAMETER (otleftlimit=5) 
      PARAMETER (otrightlimit=6) 
!                                                                       
!      call timstart ('getsnglatpt')                                    
!                                                                       
      err=.false. 
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
      IF(err)stop 2 
!                                                                       
!   Load bounding pressure indices                                      
!                                                                       
      iphi=iplo+1 
      plo=tables(ptable2+iplo) 
      phi=tables(ptable2+iphi) 
!                                                                       
!  Find bounding t indices on plo line from Table 8, 9, or 10           
!                                                                       
      tleftlimitplo=tables(ptablelimits+(iplo-ofirstpres)*ntablelimits+ &
      omint)                                                            
      trightlimitplo=tables(ptablelimits+(iplo-ofirstpres)*ntablelimits+&
      omaxt)                                                            
      itleftlimitplo=int(tables(ptablelimits+(iplo-ofirstpres)*         &
      ntablelimits+otleftlimit))                                        
      itrightlimitplo=int(tables(ptablelimits+(iplo-ofirstpres)*        &
      ntablelimits+otrightlimit))                                       
!                                                                       
      IF(t.lt.tleftlimitplo)then 
         itleftplo=itleftlimitplo 
         itrightplo=itleftplo+1 
      ELSEIF(t.gt.trightlimitplo)then 
         itleftplo=itrightlimitplo-1 
         itrightplo=itleftplo+1 
      ELSE 
         DO 10 i=itleftlimitplo+1,itrightlimitplo 
            IF(tables(ptable1+i).ge.t)then 
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
!  Find bounding t indices on phi line from Table 8, 9, or 10           
!                                                                       
      tleftlimitphi=tables(ptablelimits+(iphi-ofirstpres)*ntablelimits+ &
      omint)                                                            
      trightlimitphi=tables(ptablelimits+(iphi-ofirstpres)*ntablelimits+&
      omaxt)                                                            
      itleftlimitphi=int(tables(ptablelimits+(iphi-ofirstpres)*         &
      ntablelimits+otleftlimit))                                        
      itrightlimitphi=int(tables(ptablelimits+(iphi-ofirstpres)*        &
      ntablelimits+otrightlimit))                                       
!                                                                       
      IF(t.lt.tleftlimitphi)then 
         itleftphi=itleftlimitphi 
         itrightphi=itleftphi+1 
      ELSEIF(t.gt.trightlimitphi)then 
         itleftphi=itrightlimitphi-1 
         itrightphi=itleftphi+1 
      ELSE 
         DO 30 i=itleftlimitphi+1,itrightlimitphi 
            IF(tables(ptable1+i).ge.t)then 
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
      tleft=tables(ptable1+itleft) 
      tright=tables(ptable1+itright) 
      deltt=tright-tleft 
!                                                                       
!  Calculate t at the endpoints of the input p line                     
!  using the 4-point formula                                            
!                                                                       
      CALL loadcorners_cupid(tables,a,itleft,itright,iplo,ptableprop,         &
      stableprop,ntableprop,ofirstpres,ptr)                             
!                                                                       
      CALL herm2d_cupid(tables,ptr,a,arg,g,t,p,err) 
!                                                                       
!rex+ 19 Mar 2001                                                       
 1020 CONTINUE 
      IF(g(2).lt.0.0d0)then 
         CALL rholine_cupid(tables,ptr,a,arg,g,t,p,err) 
      ENDIF 
      IF(g(3).gt.0.0d0)then 
         CALL kappaline_cupid(tables,ptr,a,arg,g,t,p,err) 
      ENDIF 
      IF(g(5).lt.0.0d0)then 
         CALL betaline_cupid(tables,ptr,a,arg,g,t,p,err) 
         IF(g(5).lt.0.0d0)then 
            IF(t.gt.tright)then 
               CALL betaline_cupid(tables,ptr,a,arg,g,tright,p,err) 
            ELSEIF(t.lt.tleft)then 
               CALL betaline_cupid(tables,ptr,a,arg,g,tleft,p,err) 
            ENDIF 
         ENDIF 
         IF(g(5).lt.0.0d0)then 
            CALL beta2d_cupid(tables,ptr,a,arg,g,t,p,err) 
         ENDIF 
         IF(g(5).lt.0.0d0)then 
            CALL beta2d_cupid(tables,ptr,a,arg,g,tright,phi,err) 
         ENDIF 
      ENDIF 
      IF(g(7).gt.0.0d0)then 
         CALL cpline_cupid(tables,ptr,a,arg,g,t,p,err) 
      ENDIF 
!rex-                                                                   
!  Calculate single-phase properties using the returned                 
!  elements of g                                                        
!                                                                       
      IF(getprops(ubar))then 
         s(ubar)=g(1)-t*g(4)-p*g(2) 
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
      ENDIF 
      IF(getprops(cp))then 
         s(cp)=-t*g(7) 
      ENDIF 
      IF(getprops(hbar))then 
         s(hbar)=g(1)-t*g(4) 
      ENDIF 
      IF(getprops(qual))then 
         IF(arg.eq.'liquid')then 
            s(qual)=0.0d0 
         ELSEIF(arg.eq.'vapor')then 
            s(qual)=1.0d0 
         ELSEIF(arg.eq.'supercritical')then 
!gam  set sat properties equal to single phase properties               
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
!gam  do not return a value for s(qual) when supercritical              
!gam          else                                                      
!gam            if (s(temp) .le. tcrit) then                            
!gam              s(qual) = 0.0                                         
!gam            s(qual) = 0.5                                           
!gam            else                                                    
!gam              s(qual) = 1.0                                         
!gam            endif                                                   
         ENDIF 
      ENDIF 
!                                                                       
!  Check interpolated values                                            
!                                                                       
      CALL checkvalue_cupid(arg,p,g,s,getprops,tleft,tright,deltt,err) 
!                                                                       
  999 CONTINUE 
!                                                                       
!      call timstop ('getsnglatpt')                                     
!                                                                       
      RETURN 
      END SUBROUTINE getsnglatpt_cupid                    
