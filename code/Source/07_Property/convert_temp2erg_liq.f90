!
      SUBROUTINE convert_temp2erg_liq(pi,tli,eli)
!
!     Convert temparature into internal energy
!
      USE STM_TBL_cupid  , ONLY: st_tbl,              &
                                 nt,np,ns,ns2,ndxstd, &
                                 nfluid
!
      IMPLICIT NONE
!
!.....Input
      REAL(8) pi,tli
!.....Output
      REAL(8) eli,rholi
!.....Local variables
      INTEGER it
      LOGICAL erx,SHL
      REAL(8) psats,tsat,usubfs,usubgs,term,vsubfs,vsubgs,             &
              hsubfs,hsubgs,betafs,betags,kapafs,kapags,cpfs,cpgs,entfs,entgs
      REAL(8) press,vbars,ubars,hbars,betas,kapas,cps,quals,uerm
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
      DATA SHL/.FALSE./ 
! 
!.....Judge state of steam 
! 
! 
!.....Judge state of liquid 
! 
!     Initialize s for sth2x3_cupid
      s(:)=0.d0
      s(1)=tli
      s(2)=pi
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
         print *, '#### ERROR1: sth2x3_cupid called from convert_temp2erg'
         stop
      ENDIF
!      
!.....Calculate el,rhol         
!
!.....Subcooled liquid
      IF(it.eq.1)THEN 
         eli=ubars
         rholi=1.0d0/vbars
!         
!.....Saturated liquid
      ELSEIF(it.eq.2)THEN 
         s(2)=pi         
         IF(nfluid.eq.1)then 
            CALL sth2x2_cupid(st_tbl(ndxstd),st_tbl(ndxstd+nt),s,erx)
         ELSEIF(nfluid.eq.2)then 
            CALL std2x2_cupid(st_tbl(ndxstd),s,erx)
         ELSEIF(nfluid.eq.15)then 
            CALL nth2x2_cupid(st_tbl(ndxstd),s,erx) 
         ELSE 
            CALL strpx_cupid(st_tbl(ndxstd),s,erx) 
         ENDIF          
         eli=usubfs !
         rholi=1.0d0/vsubfs
!
!.....Superheated liquid
      ELSEIF(it.eq.3.or.it.eq.4)THEN 
         s(2)=pi         
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
!........Superheated liquid       
         IF(SHL)THEN 
             term=tli-tsat
             uerm=term*cpfs
             eli=usubfs+uerm
             vsubfs=vsubfs*(1.0d0+betafs*term)
             rholi=1.0d0/vsubfs 
!
!........Superheated liquid -> saturated       
         ELSE 
             tli=tsat
             eli=usubfs 
             rholi=1.0d0/vsubfs             
         ENDIF
      ELSE
         WRITE(*,*)'Error in initial liquid state!'
         WRITE(*,*)'p,tl=',pi,tli         
         STOP   
      ENDIF       
! 
      END SUBROUTINE convert_temp2erg_liq  