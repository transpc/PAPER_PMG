!
      SUBROUTINE convert_temp2erg_gas(pi,tgi,qualai,egi,rhogi,pps_oi,estm_oi, &
                                  taoi,cvaoi,uaoi,dcvai,rai,q_n)
!
!     Convert temparature into internal energy
!
      USE STM_TBL_cupid  , ONLY: st_tbl,              &
                                 nt,np,ns,ns2,ndxstd, &
                                 nfluid,              &
                                 wmole
      USE Zncg     , ONLY: n_ncg_sp,ncg_species
!
      IMPLICIT NONE
!
!.....Input
      REAL(8) pi,tgi,qualai
      REAL(8) taoi,cvaoi,uaoi,dcvai,rai
!.....Output
      REAL(8) egi,rhogi
      REAL(8) pps_oi,estm_oi
!.....Local variables
      INTEGER it,ncg,p
      LOGICAL erx,SSS
      REAL(8) psats,tsat,usubfs,usubgs,term,vsubfs,vsubgs,ua,pps,             &
              hsubfs,hsubgs,betafs,betags,kapafs,kapags,cpfs,cpgs,entfs,entgs
      REAL(8) press,vbars,ubars,hbars,betas,kapas,cps,quals,uerm
      REAL(8) q_n(n_ncg_sp),mole_weight_ncg,mole_weight_steam,unitmmole
!.....Local arrays
      REAL(8) :: s(36)
!      
      EQUIVALENCE(s( 1),tsat),    &
                 (s( 2),press),   &
                 (s( 3),vbars),   &
                 (s( 4),ubars),   &
                 (s( 5),hbars),   &
                 (s( 6),betas),   &
                 (s( 7),kapas),   &
                 (s( 8),cps),     &
                 (s( 9),quals),   &
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
      DATA SSS/.FALSE./   !ACTUAL   
      DATA mole_weight_steam/18.02d0/
! 
!.....Judge state of steam 
! 
!     Initialize s for sth2x3_cupid
      s(:)=0.d0
      s(2)=pi
      pps=1.0d0 !initial steam molecular fraction =1.d0
      pps_oi=pi !ST-pik
      IF(qualai.gt.0.99999d0)qualai=0.99999d0
!      IF(qualai.gt.0.9999999999d0)qualai=0.9999999999d0
      IF(qualai.gt.1.0e-5)THEN 
         mole_weight_ncg=0.0d0
         unitmmole=0.0d0               !20181031 LCW mixture molecular weight
         DO p=1,n_ncg_sp
            ncg=ncg_species(p)
            unitmmole=unitmmole+q_n(p)/wmole(ncg)
         ENDDO
         DO p=1,n_ncg_sp
            ncg=ncg_species(p)
            mole_weight_ncg=mole_weight_ncg+q_n(p)/wmole(ncg)/unitmmole*wmole(ncg)
         ENDDO
         pps=qualai/(1.0d0-qualai)*mole_weight_steam/mole_weight_ncg
         pps=1.0d0/(1.0d0+pps)
         s(2)=pi*pps  !steam partial pressure 
         pps_oi=s(2)                 
      ENDIF
      s(1)=tgi
      IF(nfluid.eq.1)then 
         CALL sth2x3_cupid(s,it,erx,                          &
                           st_tbl(ndxstd),                    &
                           st_tbl(ndxstd+nt),                 &
                           st_tbl(ndxstd+nt+np+13*ns+13*ns2))
      ELSEIF(nfluid.eq.2)then 
         CALL std2x3_cupid(st_tbl(ndxstd),s,it,erx) 
      ELSEIF(nfluid.eq.15)then 
         CALL nth2x3_cupid(st_tbl(ndxstd),s,it,erx)          
      ELSE 
         CALL strtp_cupid(st_tbl(ndxstd),s,it,erx) 
      ENDIF            
      IF(erx)then
         print *, '#### ERROR2: sth2x3_cupid called from convert_temp2erg'
         stop
      ENDIF
      !here s(1) is not tsat, but tgi.
!
!.....Calculate eg,rhog 
!        
!.....Super-saturated steam
      IF(it.eq.1)THEN 
         IF(nfluid.eq.1)then 
            CALL sth2x2_cupid(st_tbl(ndxstd),st_tbl(ndxstd+nt),s,erx)
         ELSEIF(nfluid.eq.2)then 
            CALL std2x2_cupid(st_tbl(ndxstd),s,erx)
         ELSEIF(nfluid.eq.15)then 
            CALL nth2x2_cupid(st_tbl(ndxstd),s,erx) 
         ELSE 
            CALL strpx_cupid(st_tbl(ndxstd),s,erx) 
         ENDIF           
!         
!........Supersaturated steam         
         IF(SSS)THEN 
             term=tgi-tsat
             uerm=term*(cpgs-s(2)*betags*vsubgs)
             usubgs=usubgs+uerm
             vsubgs=vsubgs*(1.0d0+betags*term)
             term=dmax1(tgi-taoi,0.0d0) 
             ua=cvaoi*tgi+0.5d0*dcvai*term*term+uaoi 
             egi=qualai*ua+(1.0d0-qualai)*usubgs
             rhogi=1.0d0/vsubgs+(pi-s(2))/(rai*tgi) 
!
!........Supersaturated steam -> saturated 
         ELSE  
             tgi=tsat              
             term=dmax1(tgi-taoi,0.0d0) 
             ua=cvaoi*tgi+0.5d0*dcvai*term*term+uaoi
             egi=qualai*ua+(1.0d0-qualai)*usubgs  
             rhogi=1.0d0/vsubgs+(pi-s(2))/(rai*tgi)
         ENDIF
!         
!.....Saturated steam 
      ELSEIF(it.eq.2)THEN 
         IF(nfluid.eq.1)then 
            CALL sth2x2_cupid(st_tbl(ndxstd),st_tbl(ndxstd+nt),s,erx)
         ELSEIF(nfluid.eq.2)then 
            CALL std2x2_cupid(st_tbl(ndxstd),s,erx)
         ELSEIF(nfluid.eq.15)then 
            CALL nth2x2_cupid(st_tbl(ndxstd),s,erx) 
         ELSE 
            CALL strpx_cupid(st_tbl(ndxstd),s,erx) 
         ENDIF           
         term=dmax1(tgi-taoi,0.0d0) 
         ua=cvaoi*tgi+0.5d0*dcvai*term*term+uaoi 
         egi=qualai*ua+(1.0d0-qualai)*usubgs     
         rhogi=1.0d0/vsubgs+(pi-s(2))/(rai*tgi)                   
!         
!.....Superheated steam 
      ELSEIF(it.eq.3.or.it.eq.4)THEN 
         term=dmax1(tgi-taoi,0.0d0) 
         ua=cvaoi*tgi+0.5d0*dcvai*term*term+uaoi 
         egi=qualai*ua+(1.0d0-qualai)*ubars
         rhogi=1.0d0/vbars+(pi-s(2))/(rai*tgi)
      ELSE
         WRITE(*,*)'Error in initial gas state!'
         WRITE(*,*)'p,tg=',pi,tgi
         STOP           
      ENDIF
!
      estm_oi=egi   
! 
!.....Judge state of liquid 
! 
      END SUBROUTINE convert_temp2erg_gas
