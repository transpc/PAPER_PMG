      SUBROUTINE strpx_cupid(a,s,err) 
!deck strpx                                                             
!                                                                       
!                                                                       
!      strpx   - compute thermodynamic properties as a function of      
!                pressure and quality                                   
!                                                                       
!      Calling sequence:                                                
!                                                                       
!                call  strpx (rp1,rp2,lp3)                              
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
!      This routine adapted from sth2x2 (entry point in sth2x1) written 
!      by R. J. Wagner for light water steam tables                     
!                                                                       
      USE STM_TBL_cupid  , ONLY: nt,np
!
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL(8) a(*),s(*) 
      LOGICAL err 
      REAL(8) dpsdts 
!                                                                       
      EXTERNAL strsat_cupid,strx_cupid 
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
!--temporary patch to be able to do ice condenser debug runs            
!                                                                       
      s(2)=max(s(2),ptrip) 
!                                                                       
!--get saturation temperature                                           
!                                                                       
      CALL strsat_cupid(a,2,s(2),s(1),dpsdts,err) 
      IF(err) GOTO 10 
!                                                                       
!--get thermodynamic properties as a function of quality                
!                                                                       
      s(10)=s(2) 
      CALL strx_cupid(a,s,err) 
!                                                                       
!--done                                                                 
!                                                                       
   10 RETURN 
!                                                                       
!                                                                       
      END SUBROUTINE strpx_cupid                          
