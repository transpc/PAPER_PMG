!
      SUBROUTINE udfn_prop_liquid
!
!     Computes various thermodynamic properties.
!
      USE STM_TBL_cupid , ONLY: st_tbl,              &
                                nt,np,ns,ns2,ndxstd, &
                                nfluid
      USE Zsg           , ONLY: n_group,n_1d,p_1d,t_1d,en_1d,rho_1d, &
                                vis_1d,cond_1d,cp_1d,drde_1d
!
      IMPLICIT NONE
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
      INTEGER i,iq,j
!
      LOGICAL s3
      LOGICAL erx
!
      REAL(8) :: ts
      REAL(8) :: dte,dv,rdp,rdu,tf,vb,vsat,term
      REAL(8) :: tt,pres,vbar,ubar,hbar,beta,kapa,cp,psat,vsubf,vsubg, &
                 usubf,hsubf,hsubg,betaf,kapaf,cpf
      REAL(8) :: tsat,psats,vsubfs,vsubgs,usubfs,hsubfs,hsubgs, &
                 betafs,kapafs,cpfs
      REAL(8) :: prop(36),s(36)
!
      EQUIVALENCE(prop( 1),tt),     &
                 (prop( 2),pres),   &
                 (prop( 3),vbar),   &
                 (prop( 4),ubar),   &
                 (prop( 5),hbar),   &
                 (prop( 6),beta),   &
                 (prop( 7),kapa),   &
                 (prop( 8),cp),     &
                 (prop(10),psat),   &
                 (prop(11),vsubf),  &
                 (prop(12),vsubg),  &
                 (prop(13),usubf),  &
                 (prop(15),hsubf),  &
                 (prop(16),hsubg),  &
                 (prop(17),betaf),  &
                 (prop(19),kapaf),  &
                 (prop(21),cpf)
!
      EQUIVALENCE(s( 1),tsat),    &
                 (s(10),psats),   &
                 (s(11),vsubfs),  &
                 (s(12),vsubgs),  &
                 (s(13),usubfs),  &
                 (s(15),hsubfs),  &
                 (s(16),hsubgs),  &
                 (s(17),betafs),  &
                 (s(19),kapafs),  &
                 (s(21),cpfs)
!                                                                       
!.....Set up for loops over components and volumes.                        
!
      DO i=1,n_group
         DO j=1,n_1d(i)
!
!...........save btflag is input option for thermal condition ex) 3=P,T
!
            pres=p_1d(i,j)
            s(2)=p_1d(i,j)
            s(9)=0.0d0 
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
            ts=tsat 
!
            ubar=en_1d(i,j) 
            prop(10)=tsat 
!
!...........Make (P,uf) call, already have sat properties 
!
            IF(nfluid.eq.1)then 
               s3=.true.
               CALL sth2x6_cupid(s3,prop,iq,erx,                    &
                                 st_tbl(ndxstd),                    &
                                 st_tbl(ndxstd+nt),                 &
                                 st_tbl(ndxstd+nt+np+13*ns+13*ns2))
            ELSEIF(nfluid.eq.2)then 
               CALL std2xf_cupid(st_tbl(ndxstd),prop,iq,erx)
            ELSEIF(nfluid.eq.15)then 
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
               vb=vsubf*betaf 
               cpf=cpfs
               vsat=vsubf 
               term=(ubar-usubf)/(cpf-p_1d(i,j)*vb) 
               tf=tt+term 
               vsubf=vsubf+vb*term 
               betaf=vb/vsubf 
               kapaf=vsat*kapaf/vsubf 
!         
               IF(vsubf.gt.0.0d0.and.term.le.50.0d0) GOTO 130 
            ENDIF
!
!...........Subcooled liquid state.
!
            tf=tt 
            cpf=cp 
            betaf=beta 
            kapaf=kapa 
            vsubf=vbar 
!
!...........Liquid related Derivatives.
!
  130       term=-1.0d0/(vsubf*vsubf)          ! -1/v^2
            dv=vsubf*betaf                     ! v*beta
            rdu=1.0d0/(cpf-dv*p_1d(i,j))       ! 1/(Cp-v*beta*P)
!
!...........dRho/dU=-v*beta/[v^2*(Cp-v*beta*P)]
!
            drde_1d(i,j)=term*dv*rdu 
!
            dv=cpf*vsubf*kapaf-tf*(vsubf*betaf)**2    !Cp*v*K-T(v*beta)^2
            dte=p_1d(i,j)*vsubf*kapaf-tf*vsubf*betaf  !P*v*K-T*v*beta
            rdp=-rdu 
!
            rho_1d(i,j)=1.0d0/vsubf 
            t_1d(i,j)=tf
            cp_1d(i,j)=cpf !hkcho 20081020
!
            IF(t_1d(i,j).gt.ts) STOP '#### primary coolant is superheated ####'
!
!...........transport properties
!
            CALL udfn_sg_viscos_lw_cupid(t_1d(i,j),rho_1d(i,j),vis_1d(i,j))
            CALL usfn_sg_cond_lw_cupid(t_1d(i,j),rho_1d(i,j),cond_1d(i,j))
!
         ENDDO
      ENDDO
! 
      END SUBROUTINE udfn_prop_liquid
