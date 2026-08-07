!
      SUBROUTINE steamtable(imode,np1)
! 
!     Computes various thermodynamic properties.
!
!     imode=0 when this subroutine is called first time
!     imode=1 from the second call
!
      USE Zinterface
      USE VOL_DATA       , ONLY: cell
      USE STM_TBL_cupid  , ONLY: st_tbl,                            &
                                 nt,ns,ns2,                         &
                                 ndxstd,nfluid,                     &
                                 pcrit,                             &
                                 wmole,visa,thcax,thcbx,dcvax,cvaox
      USE Zmpi     , ONLY: jperm
      USE Zcore    , ONLY: np
      USE Zdel_scalar , ONLY:err_stm,suspend_iht_opt,suspend_erg_opt,stmtbl_repeat_for_nc
      USE Zncg     , ONLY: n_ncg_sp,tao,cvao_cell,uao_cell,dcva_cell,ra_cell,qn_cell,wmole_gas,ncg_species
      USE Ztimecon , ONLY: time
      USE Zncg     , ONLY: i_ncg_vis
      USE Zzone    , ONLY: ncell_fluid
!
      USE Zncg     , ONLY: tao,cvao_cell,uao_cell,dcva_cell,ra_cell,qn_cell,n_ncg_sp !ST-pik
!      
      IMPLICIT NONE 
      INCLUDE 'stcom.h' 
!
!     prop - array in which state properties are passed to and
!            from steam table subroutines.
!     prop( 1) = temperature
!         ( 2) = pressure
!         ( 3) = specific volume
!         ( 4) = specific internal energy
!         ( 5) = specific enthalpy
!         ( 6) = single phase beta
!         ( 7) = single phase kappa
!         ( 8) = single phase csubp
!         ( 9) = quality if two-phase
!         (10) = saturation presure
!         (11) = liquid specific volume
!         (12) = vapor  specific volume
!         (13) = liquid specific internal energy
!         (14) = vapor  specific internal energy 
!         (15) = liquid specific enthalpy
!         (16) = vapor  specific enthalpy
!         (17) = liquid beta
!         (18) = vapor  beta
!         (19) = liquid kappa
!         (20) = vapor  kappa
!         (21) = liquid csubp
!         (22) = vapor  csubp
!         (23) = indexs
!         (24) = specific entropy
!         (25) = liquid entropy
!         (26) = vapor entropy
!
!.....Input
      INTEGER :: imode,np1
!.....Local variables
      INTEGER :: i,igas,iliq,iq,itype,k
!      INTEGER iones(2),niter 
      INTEGER :: niter 
      INTEGER :: j
      INTEGER :: iones(2),FluidType,ncg 
!
      LOGICAL :: s3
      LOGICAL :: erx,jstop
      LOGICAL,SAVE :: initial=.true.
      LOGICAL :: repeat
!
      REAL(8) :: ua,duadt,df2dp,dpsdp,dusdp,dpsdug
      REAL(8) :: f11,f12,f21,f22,r1,r2,rdet,delps,delus,toler,dvsdus
      REAL(8) :: dvsdps,df1dug,dusdug,df1dxa,df2dxa,dpsdxa,dusdxa
      REAL(8) :: term1,term2,advdp,advdug,advdxa,rvsubg,sdvdp,sdvdug,sdvdxa
      REAL(8) :: xa,xs,plast,ulast,rs,dp,du,dtgdps,dtgdus,pa
      REAL(8) :: tlimit,tttt 
      REAL(8) :: plimit,pppp 
      REAL(8) :: tg,dv,rdp,rdu,tf,vb,vsat,term
      REAL(8) :: dte
      REAL(8) :: tt,pres,vbar,ubar,hbar,beta,kapa,cp,qual,psat,vsubf,vsubg, &
                 usubf,usubg,hsubf,hsubg,betaf,betag,kapaf,kapag,cpf,cpg,   &
                 entpy,entpyf,entpyg
      REAL(8) :: tsat,psats,vsubfs,vsubgs,usubfs,usubgs,hsubfs,hsubgs, &
                 betafs,betags,kapafs,kapags,cpfs,cpgs,entfs,entgs
      REAL(8) :: tsat1(1),pres1(1)
      EQUIVALENCE(tsat,tsat1)
      EQUIVALENCE(pres,pres1)
      REAL(8) :: prop(36),s(36)
      
      REAL(8) :: treff,trat,mwvap           !for mixture viscosity with multi-NCG
      REAL(8) :: Tintf,timin,delt,pvapi,dpvdt
      REAL(8) :: mwgas,rhovs,visvs,condvs,wgas
      REAL(8) :: wvapi,wvapb,mwmixb,mwmixi,wvref
      REAL(8) :: conmix,vismix,xphi,phiij,dugndt,dugdt      
      REAL(8) :: rhovs1(1),condvs1(1),visvs1(1)
      EQUIVALENCE(rhovs,rhovs1)
      EQUIVALENCE(condvs,condvs1)
      EQUIVALENCE(visvs,visvs1)
      REAL(8),DIMENSION(:),ALLOCATABLE,SAVE :: visg,cong,mw,xmf
      PARAMETER(treff=114.d0)
      PARAMETER(mwvap=18.016d0)
      REAL(8) :: epsilon=1.d-8
      REAL(8) :: p_i,tl_i,tg_i,el_i,eg_i !ST-pik
      REAL(8) :: qn_cell0(n_ncg_sp)      !ST-pik
!
      EQUIVALENCE(prop( 1),tt),     &
                 (prop( 2),pres),   &
                 (prop( 3),vbar),   &
                 (prop( 4),ubar),   &
                 (prop( 5),hbar),   &
                 (prop( 6),beta),   &
                 (prop( 7),kapa),   &
                 (prop( 8),cp),     &
                 (prop( 9),qual),   &
                 (prop(10),psat),   &
                 (prop(11),vsubf),  &
                 (prop(12),vsubg),  &
                 (prop(13),usubf),  &
                 (prop(14),usubg),  &
                 (prop(15),hsubf),  &
                 (prop(16),hsubg),  &
                 (prop(17),betaf),  &
                 (prop(18),betag),  &
                 (prop(19),kapaf),  &
                 (prop(20),kapag),  &
                 (prop(21),cpf),    &
                 (prop(22),cpg),    &
                 (prop(24),entpy),  &
                 (prop(25),entpyf), &
                 (prop(26),entpyg)
!
      EQUIVALENCE(s( 1),tsat),    &
                 (s(10),psats),   &
                 (s(11),vsubfs),  &
                 (s(12),vsubgs),  &
                 (s(13),usubfs),  &
                 (s(14),usubgs),  &
                 (s(15),hsubfs),  &
                 (s(16),hsubgs),  &
                 (s(17),betafs),  &
                 (s(18),betags),  &
                 (s(19),kapafs),  &
                 (s(20),kapags),  &
                 (s(21),cpfs),    &
                 (s(22),cpgs),    &
                 (s(25),entfs),   &
                 (s(26),entgs)
!
!.....Data statements.
!
      DATA iones/1,1/ 
      DATA toler/1.d-6/ 
!
      tlimit(tttt)=min(max(273.16d0,tttt),tcrit) 
      plimit(pppp)=min(max(611.24d0,pppp),pcrit) 
      FluidType=nfluid   
!                                                                       
!.....Set up for loops over components and volumes.                        
!
      err_stm(:)=0
!
      DO i=1,ncell_fluid
!
!........save btflag is input option for thermal condition ex) 3=P,T
!
         itype=1
         pres=cell%p(i)
!         
         IF(stmtbl_repeat_for_nc.eq.1) THEN !somaflow.in
            repeat=.true.
            s(2)=pres         
         ELSE   
            s(2)=cell%p(i)         
         ENDIF
102         CONTINUE         
!
         s(9)=0.d0 
!
         IF(nfluid.eq.1)then 
            CALL sth2x2_cupid(st_tbl(ndxstd),st_tbl(ndxstd+nt),s,erx)
         ELSEIF(nfluid.eq.2)then 
            CALL std2x2_cupid(st_tbl(ndxstd),s,erx)
         ELSEIF(nfluid.eq.15)then 
            CALL nth2x2_cupid(st_tbl(ndxstd),s,erx)
         ELSE 
            CALL strpx_cupid(st_tbl(ndxstd),s,erx) 
         ENDIF          
         cell%ts(i)=tsat 
         cell%tst(i)=tsat
!        dTs/dP=Ts*(vgs-vfs)/(Hgs-Hfs) : Clapeyron-Clausius
         cell%dtsdp(i)=tsat*(vsubgs-vsubfs)/(hsubgs-hsubfs) 
         cell%elsat(i)=usubfs
         cell%egsat(i)=usubgs
         cell%hlsat(i)=hsubfs
         cell%hgsat(i)=hsubgs
!
         IF(cell%alphag(i).le.epsilon) THEN
            igas=0
         ELSE
            igas=1
         ENDIF
!
         IF(cell%alphag(i).ge.1.d0-1.d0*epsilon) THEN
            iliq=0
         ELSE
            iliq=1
         ENDIF
!
!........Liquid properties
!
         IF(iliq.eq.0) THEN
!
!...........Saturated liquid state. 
!
            tf=tsat 
            cpf=cpfs 
            betaf=betafs 
            cell%betal(i)=betaf                  
            kapaf=kapafs 
            vsubf=vsubfs 
         ELSE
!
            ubar=cell%el(i) 
            prop(10)=tsat 
!
!...........Make (P,uf) call, already have sat properties 
!
            IF(nfluid.eq.1)then 
               s3=.true.
               CALL sth2x6_cupid(s3,prop,iq,erx,                    &
                                 st_tbl(ndxstd),                    &
                                 st_tbl(ndxstd+nt),                 &
                                 st_tbl(ndxstd+nt+np1+13*ns+13*ns2))
            ELSEIF(nfluid.eq.2)then 
               CALL std2xf_cupid(st_tbl(ndxstd),prop,iq,erx)
            ELSEIF(nfluid.eq.15)then 
               iq=5
               CALL nth2x6f_cupid(st_tbl(ndxstd),prop,iq,erx,'f') 
            ELSE 
               CALL strpu2_cupid(st_tbl(ndxstd),prop,iq,erx) 
            ENDIF           
!
            IF(iq.eq.2) THEN
!
!..............Superheated liquid state.
!..............Extrapolate specific volume and temperature at constant pressure. 
!
               IF(nfluid.ne.15)then          
                  vb=vsubf*betaf 
                  cpf=cpfs
                  vsat=vsubf 
                  term=(ubar-usubf)/(cpf-cell%p(i)*vb) 
                  tf=tt+term 
                  vsubf=vsubf+vb*term 
                  betaf=vb/vsubf 
                  cell%betal(i)=betaf         
                  kapaf=vsat*kapaf/vsubf 
!         
                  IF(vsubf.gt.0.d0.and.term.le.60.d0) THEN
                     GOTO 130
                  ELSE
                     write(*,*) 'WARNING1: steamtable two phase was detected but temperature sounds too high ',term
                     write(*,*) 'cupid will continue processing but at risk'
                     IF(suspend_erg_opt.eq.1)THEN !ST-pik               
                        p_i=pres
                        tl_i=tsat+40.0d0 !60.0d0 
                        CALL convert_temp2erg_liq(p_i,tl_i,el_i)
                        cell%el(i)=el_i
                     ELSE
                        IF(suspend_iht_opt.gt.0) err_stm(i)=1 
                     ENDIF   
                     GOTO 130
                  ENDIF
               ENDIF         
            ENDIF         
!
!...........Subcooled liquid state.
!
            tf=tt 
            cpf=cp 
            betaf=beta 
            cell%betal(i)=betaf         
            kapaf=kapa 
            vsubf=vbar 
         ENDIF
!
!........Liquid related Derivatives.
!
  130    term=-1.d0/(vsubf*vsubf)          ! -1/v^2
         dv=vsubf*betaf                     ! v*beta
         rdu=1.d0/(cpf-dv*cell%p(i))       ! 1/(Cp-v*beta*P)
!
!........dRho/dU=-v*beta/[v^2*(Cp-v*beta*P)]
!
         cell%drholde(i)=term*dv*rdu 
!
!........dT/dU=1/(Cp-v*beta*P)
!
         cell%dtlde(i)=rdu 
!
         dv=cpf*vsubf*kapaf-tf*(vsubf*betaf)**2    !Cp*v*K-T(v*beta)^2
         dte=cell%p(i)*vsubf*kapaf-tf*vsubf*betaf  !P*v*K-T*v*beta
         rdp=-rdu 
!
!........dRho/dP=[Cp*v*K-T*(v*beta)^2]/[v^2*(Cp-v*beta*P)]
!
         cell%drholdp(i)=term*dv*rdp 
!
!........dT/dP=-(P*v*K-T*v*beta)/(Cp-v*beta*P)
!
         cell%dtldp(i)=dte*rdp 
!
         cell%rhol(i)=1.d0/vsubf 
         cell%tl(i)=tf
         cell%cpl(i)=cpf !hkcho 20081020
!
         IF(iliq.eq.0)THEN
            cell%hl(i)=hsubfs
         ELSE
            cell%hl(i)=hbar
         ENDIF
!
!........Vapor properties
!
         ubar=cell%eg(i) 
         prop(10)=tsat 
!
         xa=cell%quala(i) 
         xs=1.d0-xa 
!
!........Pure steam
!
         IF(xa.le.1.0e-8)THEN
            cell%quala(i)=0.d0
            cell%pps(i)=cell%p(i)
!
!...........Make (P,ug) call, already have sat properties 
!
            IF(nfluid.eq.1)then 
               s3=.true.
               CALL sth2x6_cupid(s3,prop,iq,erx,                    &
                                 st_tbl(ndxstd),                    &
                                 st_tbl(ndxstd+nt),                 &
                                 st_tbl(ndxstd+nt+np1+13*ns+13*ns2))
            ELSEIF(nfluid.eq.2)then 
               CALL std2xf_cupid(st_tbl(ndxstd),prop,iq,erx)
            ELSEIF(nfluid.eq.15)then 
               iq=6
               CALL nth2x6f_cupid(st_tbl(ndxstd),prop,iq,erx,'f') 
            ELSE 
               CALL strpu2_cupid(st_tbl(ndxstd),prop,iq,erx) 
            ENDIF              
            IF(igas.eq.0) THEN
!        
!..............Saturated vapor state. 
!        
               tg=tsat 
               cpg=cpgs 
               betag=betags 
               cell%betag(i)=betag            
               kapag=kapags 
               vsubg=vsubgs 
               cell%eg(i)=usubgs 
            ELSE
!        
               IF(iq.eq.2) THEN
!        
!.................Subcooled steam.
!.................Extrapolate specific volume and temperature at constant pressure.
!        
                  IF(nfluid.ne.15)then             
                     vb=vsubg*betag 
                     vsat=vsubg 
                     term=(ubar-usubg)/(cpg-cell%p(i)*vb) 
                     tg=tt+term 
                     vsubg=vsubg+vb*term 
                     betag=vb/vsubg 
                     cell%betag(i)=betag            
                     kapag=vsat*kapag/vsubg 
!        
                     IF(vsubg.gt.0.d0.and.tg.gt.ttrip.and.term.ge.-50.d0)THEN
                        GOTO 230 
                     ELSE
                        write(*,*) 'WARNING2: steamtable two phase was detected but temperature sounds too high ',term
                        write(*,*) 'cupid will continue processing but at risk'
                        IF(suspend_erg_opt.eq.1)THEN !ST-pik
                           qn_cell0(:)=qn_cell(i,:)
                           p_i=cell%p(i)
                           tg_i=tsat-30.d0 !-50.0d0
                           CALL convert_temp2erg_gas(p_i,tg_i,cell%quala(i),eg_i, &
                                cell%rhog(i),cell%pps_o(i),cell%estm_o(i),              &
                                tao,cvao_cell(i),uao_cell(i),dcva_cell(i),ra_cell(i),qn_cell0)
                           cell%eg(i)=eg_i
                        ELSE
                           IF(suspend_iht_opt.gt.0) err_stm(i)=1                    
                        ENDIF
                        GOTO 230 
                     ENDIF   
                     
                  ENDIF            
               ENDIF            
!        
!..............Superheated vapor state. 
!        
               tg=tt 
               cpg=cp 
               betag=beta
               cell%betag(i)=betag             
               kapag=kapa 
               vsubg=vbar 
            ENDIF
!        
!...........Vapor related Derivatives.
!
  230       term=-1.d0/(vsubg*vsubg) 
            dv=vsubg*betag 
            rdu=1.d0/(cpg-vsubg*betag*cell%p(i)) 
!        
!...........dRho/dU
!
            cell%drhogde(i)=term*dv*rdu 
!       
!...........dT/dU
!
            cell%dtgde(i)=rdu 
!        
            dv=cpg*vsubg*kapag-tg*(vsubg*betag)**2 
            dte=cell%p(i)*vsubg*kapag-tg*vsubg*betag 
            rdp=-rdu 
!        
!...........dRho/dP
!
            cell%drhogdp(i)=term*dv*rdp 
!        
!...........dT/dP
!
            cell%dtgdp(i)=dte*rdp 
!        
            cell%rhog(i)=1.d0/vsubg 
            cell%tg(i)=tg
            cell%cpg(i)=cpg !hkcho 20081020
            cell%estm(i)=usubg
!        
            IF(igas.eq.0)THEN
               cell%hg(i)=hsubgs
            ELSE
               cell%hg(i)=hbar
            ENDIF
!
            cell%dtgdx(i)=0.d0
            cell%drhogdx(i)=0.d0
!
            cell%dtsde(i)=0.d0
            cell%dtsdx(i)=0.d0 
            cell%ha(i)=cell%hg(i) !next-pik
!            
            IF(i_ncg_vis.ge.1)THEN    !i_ncg_vis
               IF(initial)THEN
                  ALLOCATE(visg(n_ncg_sp+1),cong(n_ncg_sp+1),mw(n_ncg_sp+1),xmf(n_ncg_sp+1))
                  visg(:)=0.d0
                  cong(:)=0.d0
                  mw(:)=0.d0
                  xmf(:)=0.d0
                  initial=.FALSE.
                  WRITE(*,*) '     NCG Properties are calculating...'
               ENDIF                
               iones=1
               rhovs=1.d0/vsubgs
               tsat=max(273.16d0,tsat)
               CALL thcond_cupid(i,FluidType,iones,iones,1,1,tsat1,pres1,rhovs1,'vapor',condvs1,erx)            !saturation vapor
               CALL viscos_cupid(i,FluidType,iones,iones,1,1,tsat1,rhovs1,'vapor',visvs1,erx)       !saturation vapor

               cell%lviscosg(i)=visvs
               cell%lcondg(i)=condvs
            ENDIF !i_ncg_vis            
!            
!........Steam with noncondensable gases
!
         ELSE
            s(9)=0.d0 
            itype=7
            IF(itype.lt.7)THEN 
               tsat=cell%ts(i)
            ELSE 
               tsat=cell%tl(i)
            ENDIF 
            tsat=tlimit(tsat) 
            IF(nfluid.eq.1)then 
               CALL sth2x1_cupid(st_tbl(ndxstd),st_tbl(ndxstd+nt),s,erx)
            ELSEIF(nfluid.eq.2)then 
               CALL std2x1_cupid(st_tbl(ndxstd),s,erx)
            ELSEIF(nfluid.eq.15)then 
               CALL nth2x1_cupid(st_tbl(ndxstd),s,erx) 
            ELSE 
               CALL strtx_cupid(st_tbl(ndxstd),s,erx) 
            ENDIF              
            IF(erx)THEN 
               PRINT *, '#### ERROR2: Steam property(NCG1) at time= ', time, 'at cell= ', jperm(i)
               IF(suspend_iht_opt.eq.0)THEN
                  STOP
               ELSE
                  err_stm(i)=2
                  CYCLE
               ENDIF
            ENDIF 
!                                                                
            pres=psats 
            ubar=usubgs 
            jstop=.false. 
            niter=0 
!
            IF(imode.gt.0)THEN
               IF(cell%quala_o(i).lt.0.99999999d0)THEN
                  pres=cell%pps_o(i)
                  ubar=cell%estm_o(i)
               ENDIF
            ENDIF      
!
!...........Beginning of iteration for ps and us.                                
!
  795       plast=pres 
            ulast=ubar 
!
!...........Steam properties (use g subscript to make coding parallel to         
!           non-condensable coding).                                             
!
            IF(nfluid.eq.1)then 
               s3=.false.
               CALL sth2x6_cupid(s3,prop,iq,erx,                    &
                                 st_tbl(ndxstd),                    &
                                 st_tbl(ndxstd+nt),                 &
                                 st_tbl(ndxstd+nt+np1+13*ns+13*ns2))
            ELSEIF(nfluid.eq.2)then 
               CALL std2x6_cupid(st_tbl(ndxstd),prop,iq,erx)
            ELSEIF(nfluid.eq.15)then 
               iq=6
               CALL nth2x6f_cupid(st_tbl(ndxstd),prop,iq,erx,'6') 
            ELSE 
               CALL strpu1_cupid(st_tbl(ndxstd),prop,iq,erx) 
            ENDIF              
            IF(erx)THEN 
               PRINT *, '#### ERROR3: Steam property (NCG2) at time= ', time, 'at cell= ', jperm(i)
               IF(suspend_iht_opt.eq.0)THEN
                  STOP
               ELSE
                  err_stm(i)=3
                  CYCLE
               ENDIF
            ENDIF 
!
            IF(iq.eq.2) THEN
!
!..............Subcooled steam.                                                     
!............. Extrapolate specific volume and temperature at constant pressure.    
!
               IF(nfluid.ne.15)THEN             
                  vb=vsubg*betag 
                  vsat=vsubg 
                  term=(ubar-usubg)/(cpg-pres*vb) 
                  tg=tt+term 
                  vsubg=vsubg+vb*term 
                  betag=vb/vsubg 
                  cell%betag(i)=betag            
                  IF(vsubg.gt.0.d0.and.term.ge.-500.d0)GOTO 830      ! yhy 2012-10-15
!
                  PRINT *, '#### ERROR4: Steam property (NCG3) at time= ', time, 'at cell= ', jperm(i)
!
                  IF(suspend_iht_opt.eq.0)THEN
                     STOP
                  ELSE
                     err_stm(i)=4
                     CYCLE
                  ENDIF
               ENDIF            
            ENDIF            
!
!...........Superheated vapor state.  Or subcooled if nfluid = 15 
!
            tg=tt 
            cpg=cp 
            betag=beta 
            cell%betag(i)=betag            
            kapag=kapa 
            entpyg=entpy 
            vsubg=vbar 
  830       rs=pres*vsubg/tg 
            dte=vsubg*(pres*kapag-tg*betag) 
            dp=pres*vsubg*betag-cpg 
            dtgdps=dte/dp 
            du=-dp 
            dtgdus=1.d0/du 
!
!...........Air properties.                                                      
!
            term=max(tg-tao,0.d0) 
            ua=cvao_cell(i)*tg+0.5d0*dcva_cell(i)*term*term+uao_cell(i)
            duadt=cvao_cell(i)+dcva_cell(i)*term 
            pa=cell%p(i)-pres 
!
!...........2*2 matrix inversion for newton iteration with approx. derivatives.  
!
            f11=xa*duadt*dtgdps 
            f12=xs+xa*duadt*dtgdus 
            f21=-xs*rs-xa*ra_cell(i)
            IF(jstop) GOTO 850 
!
            r1=xs*ubar+xa*ua-cell%eg(i)
            r2=xs*rs*pa-xa*ra_cell(i)*pres 

!...........Gas internal energy when the gas fraction is below 1.0d-8 (yhy)

            IF(cell%alphag(i).le.1.0d-8)THEN !ST-pik
                IF(suspend_erg_opt.eq.1)THEN
                   r1=0.0d0
                   cell%eg(i)=xs*ubar+xa*ua
                ENDIF
            ENDIF
!            
            rdet=-1.d0/(f12*f21) 
            delps=f12*r2*rdet 
            delus=-(-f21*r1+f11*r2)*rdet 
            IF(niter.gt.20) GOTO 840 
!
            pres=max(min(pres+delps,cell%p(i)),pmin) 
            ubar=ubar+delus 
!
            IF(cell%alphag(i).le.1.0d-8)THEN !ST-pik
               IF(suspend_erg_opt.eq.1)THEN
                  IF((abs(pres-plast).lt.toler*cell%p(i))) jstop=.true.   
               ELSE
                  IF((abs(ubar-ulast).lt.toler*cell%eg(i)).and.         &
                     (abs(pres-plast).lt.toler*cell%p(i))) jstop=.true.                  
               ENDIF   
            ELSE
               IF((abs(ubar-ulast).lt.toler*cell%eg(i)).and.         &
                  (abs(pres-plast).lt.toler*cell%p(i))) jstop=.true.                  
            endif
!        
            niter=niter+1 
            GOTO 795 
!
!...........Iteration failed.                                                    
!
  840       CONTINUE 
            PRINT *, '#### ERROR5: NCG was not converged at time= ', time, 'at cell= ', jperm(i)
            IF(suspend_iht_opt.eq.0)THEN
               STOP
            ELSE
               err_stm(i)=5
               CYCLE
            ENDIF
!
!...........Iteration converged, get final mixture properties.                   
!...........Convert approx. derivatives to full derivatives.                     
!...........dvsdus and dvsdps include 1/vsubg factor.                            
!
  850       CONTINUE
!
!...........Repeat liquid property calculations with new pps(i)
            IF(stmtbl_repeat_for_nc.eq.1.and.repeat) THEN 
               repeat=.false.
               s(2)=pres         
               GOTO 102
            ENDIF
            dv=betag             
!
            du=cpg-vsubg*betag*pres 
            dvsdus=dv/du 
            dv=cpg*kapag-tg*vsubg*betag*betag 
            dp=-du 
            dvsdps=dv/dp 
            term=1.d0/pres+dvsdps-dtgdps/tg 
            f21=f21+xs*pa*rs*term 
            term=dvsdus-dtgdus/tg 
            f22=xs*pa*rs*term 
            rdet=1.d0/(f11*f22-f12*f21) 
!
!...........p derivatives.                                                       
!
            df2dp=-xs*rs 
            dpsdp=-f12*df2dp*rdet 
            dusdp=f11*df2dp*rdet 
!
!...........ug derivatives.                                                      
!
            df1dug=1.d0 
            dpsdug=f22*df1dug*rdet 
            dusdug=-f21*df1dug*rdet 
!
!...........xa derivatives.                                                      
!
            df1dxa=ubar-ua 
            df2dxa=ra_cell(i)*pres+rs*pa 
            dpsdxa=(f22*df1dxa-f12*df2dxa)*rdet 
            dusdxa=(-f21*df1dxa+f11*df2dxa)*rdet 
!
!...........Final temperature derivatvies.                                       
!
            cell%dtgdp(i)=dtgdps*dpsdp+dtgdus*dusdp 
            cell%dtgde(i)=dtgdps*dpsdug+dtgdus*dusdug 
            cell%dtgdx(i)=dtgdps*dpsdxa+dtgdus*dusdxa 
!
!...........Final density derivatives calculated using                           
!...........vg = (vs*va)/(vs+va), which gives for the derivative of rhog       
!...........drhog = -(dva/(va*va) + dvs/(vs*vs)).                              
!...........Air formulas for dva/(va*va).                                        
!
            term2=1.d0/(ra_cell(i)*tg) 
            term1=pa*term2/tg 
            advdp=term1*cell%dtgdp(i)+term2*(dpsdp-1.d0) 
            advdug=term1*cell%dtgde(i)+term2*dpsdug 
            advdxa=term1*cell%dtgdx(i)+term2*dpsdxa 
!
!...........Steam formulas for dvs/(vs*vs).                                      
!...........Remember dvsdps and dvsdus include a 1/vsubg factor.                 
!
            rvsubg=1.d0/vsubg 
            sdvdp=(dvsdps*dpsdp+dvsdus*dusdp)*rvsubg 
            sdvdug=(dvsdps*dpsdug+dvsdus*dusdug)*rvsubg 
            sdvdxa=(dvsdps*dpsdxa+dvsdus*dusdxa)*rvsubg 
!
!...........Steam/gas mixture density derivatives.                                         
!
            cell%drhogdp(i)=-(advdp+sdvdp) 
            cell%drhogde(i)=-(advdug+sdvdug) 
            cell%drhogdx(i)=-(advdxa+sdvdxa) 
!
            cell%rhoa(i)=pa/(ra_cell(i)*tg)             ! ***COPAIN
            cell%rhog(i)=1.d0/vsubg+pa/(ra_cell(i)*tg) 
            cell%tg(i)=tg 
            cell%pps(i)=pres 
            cell%estm(i)=usubg 
            cell%hg(i)=ubar+pres*vsubg 
            cell%cpg(i)=cpg !hkcho 20081020
!
!...........Saturation properties.                                               
!
            s(2)=plimit(cell%pps(i)) 
            s(9)=0.d0 
            IF(nfluid.eq.1)then 
               CALL sth2x2_cupid(st_tbl(ndxstd),st_tbl(ndxstd+nt),s,erx)
            ELSEIF(nfluid.eq.2)then 
               CALL std2x2_cupid(st_tbl(ndxstd),s,erx)
            ELSEIF(nfluid.eq.15)then 
               CALL nth2x2_cupid(st_tbl(ndxstd),s,erx) 
            ELSE 
               CALL strpx_cupid(st_tbl(ndxstd),s,erx) 
            ENDIF                         
            IF(erx)then 
               print *, '#### ERROR6: Steam property (NCG3) at time= ', time, 'at cell= ', jperm(i)
               IF(suspend_iht_opt.eq.0)THEN
                  STOP
               ELSE
                  err_stm(i)=6
                  CYCLE
               ENDIF
            ENDIF 
!
            cell%ts(i)=tsat 
            cell%hgsat(i)=hsubgs 
!            
            term=tsat*(vsubgs-vsubfs)/(hsubgs-hsubfs) 
            cell%dtsdp(i)=term*dpsdp 
            cell%dtsde(i)=term*dpsdug 
            cell%dtsdx(i)=term*dpsdxa 
            cell%elsat(i)=cell%hlsat(i)-cell%pps(i)/cell%rhol(i)      
!...........NCG enthalpy 2015.07.28 JHLee (SNU)
            cell%ha(i)=ua+ra_cell(i)*cell%tg(i)
!...........Saturation properties 2022.05.11. IKPARK !ST-pik
            cell%elsat(i)=usubfs
            cell%hlsat(i)=hsubfs
            cell%egsat(i)=usubgs !*(1.0d0-cell%quala(i))+ua*cell%quala(i)        
            cell%hgsat(i)=hsubgs !*(1.0d0-cell%quala(i))+cell%ha(i)*cell%quala(i) 
!...........Viscosity & conductivity with NCG
!
            IF(i_ncg_vis.ge.1)THEN    !i_ncg_vis
               IF(initial)THEN
                  ALLOCATE(visg(n_ncg_sp+1),cong(n_ncg_sp+1),mw(n_ncg_sp+1),xmf(n_ncg_sp+1))
                  visg(:)=0.d0
                  cong(:)=0.d0
                  mw(:)=0.d0
                  xmf(:)=0.d0
                  initial=.FALSE.
                  WRITE(*,*) '     NCG Properties are calculating...'
               ENDIF                
               iones=1
               rhovs=1.d0/vsubgs
               tsat=max(273.16d0,tsat)
               CALL thcond_cupid(i,FluidType,iones,iones,1,1,tsat1,pres1,rhovs1,'vapor',condvs1,erx)            !saturation vapor
               CALL viscos_cupid(i,FluidType,iones,iones,1,1,tsat1,rhovs1,'vapor',visvs1,erx)       !saturation vapor
!                                                                       
!..............Store Vapor Properties as Gas #1.                                 
!                                                                       
               visg(1)=visvs 
               cong(1)=condvs 
               mw(1)=mwvap 
               mwgas=wmole_gas(i)
               wvapb=1.d0-cell%quala(i)
               wvref=wvapb
!
!..............Condensation model related (New condensation model for MARS)
!               
               IF(.FALSE.)THEN          ! IF the condensation model is used, an interface temperature calculated in the model should be used here.
                  timin=min(cell%tl(i),tsat-0.1d0)  !Interative Calculation is needed
                  Tintf=timin+0.1d0*(tsat-timin)   
                  !Tintf=cell%tg(i)              
                  CALL psatpd_cupid(Tintf,pvapi,dpvdt,1,erx)
!                                                                       
!.................Evaluate boundary layer reference mass fraction.                  
!                 ref : Knuth & Dershin, "Use of Reference States in Predicting Transport Rates in High-Speed Turbulent Flows with       
!                 Mass Transfers", IJHMT, Vol. 6, pp. 999-1018, 1963.                                                                           
                  mwmixb=mwvap*mwgas/(mwvap+(mwgas-mwvap)*(1.d0-cell%quala(i))) 
                  wvapi=(mwvap/mwgas)*pvapi/(cell%p(i)-(1.d0-(mwvap/mwgas))*pvapi) 
                  mwmixi=mwvap*mwgas/(mwvap+(mwgas-mwvap)*wvapi) 
                  wvref=(mwvap/(mwvap-mwgas))*log(mwmixb/mwmixi)/log((wvapb*mwmixb)/(wvapi*mwmixi))       !wvref=mass fraction of vapor
               ENDIF
!
               xmf(1)=wvref/mwvap 
               wgas=1.d0-wvref 
!                                                                       
!..............NC gas properties for each species 
!                                                                       
               trat=tsat**1.5d0/(tsat+treff) 
               delt=max(0.d0,tsat-250.d0) 
               dugdt=0.d0
!                                                                      
               DO k=1,n_ncg_sp
                  ncg=ncg_species(k)
                  visg(k+1)=visa(ncg)*trat
                  cong(k+1)=thcax(ncg)*tsat**thcbx(ncg)
                  dugndt=cvaox(ncg)+dcvax(ncg)*delt
                  dugdt=dugdt+qn_cell(i,k)*dugndt
                  mw(k+1)=wmole(ncg) 
                  xmf(k+1)=qn_cell(i,k)*wgas/mw(k+1)
               ENDDO 
!                                                                       
!..............Vapor+NC gas Mixture properties (mixture rule of Wilke)
!     
               conmix=0.d0
               vismix=0.d0
               DO k=1,n_ncg_sp+1
                  xphi=0.d0
                  DO j=1,n_ncg_sp+1
                     phiij=(1.d0+sqrt((visg(k)/visg(j))*sqrt(mw(j)/mw(k))))**2/sqrt(8.d0*(1.d0+(mw(k)/mw(j))))
                     xphi=xphi+xmf(j)*phiij
                  END DO
                  conmix=conmix+xmf(k)*cong(k)/xphi
                  vismix=vismix+xmf(k)*visg(k)/xphi
               END DO
               cell%lviscosg(i)=vismix
               cell%lcondg(i)=conmix
               cell%cpg(i)=wvref*cpg+wgas*(dugdt+ra_cell(i))

            ENDIF !i_ncg_vis
!
         ENDIF
!
!........Mixture properties
!
         cell%rhom(i)=cell%alphag(i)*cell%rhog(i) +(1.d0-cell%alphag(i))*cell%rhol(i)
         cell%rhomr(i)=1.d0/cell%rhom(i)
         cell%quals(i)=cell%alphag(i)*cell%rhog(i)*cell%rhomr(i)
         cell%rhod(i) =cell%rhol(i)
      ENDDO
!
!........transport properties
!
      IF(nfluid.eq.1)THEN
         DO i=1,ncell_fluid         
            cell%tl(i)=tlimit(cell%tl(i)) 
         ENDDO         
         IF(i_ncg_vis.eq.0)THEN
            CALL viscos_lw_cupid(cell%tg,cell%rhog,cell%ts,cell%lviscosg,'vap',ncell_fluid)
            CALL cond_lw_cupid(cell%tg,cell%rhog,cell%lcondg,ncell_fluid)
         ENDIF
         CALL viscos_lw_cupid(cell%tl,cell%rhol,cell%ts,cell%lviscosl,'liq',ncell_fluid)
         CALL cond_lw_cupid(cell%tl,cell%rhol(:),cell%lcondl(:),ncell_fluid)
      ELSE
         DO i=1,ncell_fluid         
            CALL thcond_cupid(i,FluidType,iones,iones,1,1,cell%tl(i),cell%p(i),cell%rhol(i),'liquid',cell%lcondl(i),erx)             !saturation vapor
            CALL viscos_cupid(i,FluidType,iones,iones,1,1,cell%ts(i),cell%rhol(i),'liquid',cell%lviscosl(i),erx)     !saturation vapor
            CALL thcond_cupid(i,FluidType,iones,iones,1,1,cell%tg(i),cell%p(i),cell%rhog(i),'vapor',cell%lcondg(i),erx)              !saturation vapor
            CALL viscos_cupid(i,FluidType,iones,iones,1,1,cell%ts(i),cell%rhog(i),'vapor',cell%lviscosg(i),erx)      !saturation vapor
         ENDDO
         
      ENDIF
!
      CALL surftn_cupid(ncell_fluid,FluidType,cell%ts,cell%sigma)
!
      IF(np.gt.1) CALL communicate_1d(cell%lviscosg, &
                                      cell%lviscosl)
! 
      END SUBROUTINE steamtable
