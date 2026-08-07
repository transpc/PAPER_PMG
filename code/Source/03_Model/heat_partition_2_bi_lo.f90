!
      SUBROUTINE heat_partition_2_bi_lo
!
!     This routine obtained partitioned heat at specific wall boundary cells
!     whose property is solid-fluid interface, constant temperature, constant 
!     heatflux(nbcon=-2,-3&-4,-5), using bi-section method. This routine do 
!     bi-section to twall until (q1-q2(twall))/q1 < 1.d-6 is satisfied.
!     q1=assinged heat, q2=partitioned heat, twall=wall temperature, err=q1-q2
!
      USE VOL_DATA        , ONLY: cell
      USE Wall_DATA       , ONLY: face
      USE Zmpi            , ONLY: jperm  
      USE Zzone           , ONLY: ncell_fluid
      USE Zcore           , ONLY: myrank        
      USE Ztimecon        , ONLY: time
      USE Zparam          , ONLY: nin_max,pi
      USE Znum_cell       , ONLY: i_neigh
      USE Zb_condition    , ONLY: twall      
      USE Zconst1         , ONLY: iheatpart,iat
      USE Zcoord3         , ONLY: volp
      USE Zface           , ONLY: q1cell,qqcell,qecell,qclcell,qcgcell,ndensitycell
      USE Zheat_partition , ONLY: q1,qe,qq,qcg,qcl,d_depart,ndensity,bfreq,a_two,dlo
      USE Ziat            , ONLY: dbubble_init,iat_nucl
      USE Zqvol           , ONLY: gamma_wall,qwall_solid,dry_weight
      USE Zbc_index       , ONLY: nbcon,icell_type,iface_wall,iface_wall1
      USE Znormal         , ONLY: sa_walll
      USE Zuserdefined    , ONLY: udfl_psbt_cfx_model
      USE Zio_unit        , ONLY: unit_log
      
!
      IMPLICIT NONE
!
!.....Input
      INTEGER i,j,j0,nfcondition
!.....Local variables
      INTEGER count,dtoption,iter,iter_limit
      INTEGER fr_lift,fr_nden,fr_depfreq,fr_depdia
      REAL(8) d_depart_init,hi_gas,hi_liq
      REAL(8) deltats,deltatl,deltatg,deltats_old,deltatl_old,deltatg_old
      REAL(8) deltarho,deltahlg !give
      REAL(8) ra_coeff,heat_partition_2_bi_lo_err !take
      REAL(8) err1,err2,err3
      REAL(8) dtwall,tw1,tw2,tw3
      REAL(8) erra,errb,errc,twa,twb,twc
      DATA iter_limit,ra_coeff/10000,1.0d0/
      DATA fr_lift,fr_nden,fr_depfreq,fr_depdia/0,0,0,0/
!
      DO i=1,ncell_fluid
         IF(icell_type(i).ne.1) cycle
         j0=i_neigh(i)-1
         j=iface_wall(i)
         nfcondition=nbcon(j+j0)
         IF(nfcondition.eq.-1)THEN
            face%wall_fluxl_diff(i)=0.0d0
            face%wall_fluxg_diff(i)=0.0d0
            face%wall_fluxd_diff(i)=0.0d0
            face%ddepartw(i)=0.0d0
            face%ratio_evap(i)=0.0d0
            GOTO 100
!           RETURN
         ENDIF
         j=iface_wall1(i)
         IF(j.eq.0) RETURN
         nfcondition=nbcon(j+j0)
!........Initialize papameters
!    
         d_depart_init=dbubble_init
         iter=1    
         count=1  
         dtwall=0.0d0 
         dtoption=1
         err1=0.0d0
         err2=0.0d0
         err3=0.0d0
!
!........Set the sub-models in HPM
!      
         fr_lift=mod(iheatpart,10000)/1000        ! Lift off model;                  0=off,  1=On
         fr_nden=mod(iheatpart,1000)/100          ! Nucleation site density model:   0=Cole, 1=Lemmert and Chwala, 2=ocamustafaogullary, 3=Hibiki, 4=Modified 
         fr_depfreq=mod(iheatpart,100)/10         ! Departure frequency model:       0=Cole, 1=Situ
         fr_depdia=mod(iheatpart,10)              ! Depart diameter model;           1: Cole and Rosenhow (1968) 2: Fritz (1935) 3: Tolubinsky model 4: Unal  5: Kocamusta  6:Constant      
!
!........Do bi-section according to boundary properties
!
         IF(nfcondition.eq.-2)THEN
!
!........AT fluid-solid interface (nbcon=-2)
!
            !face%twall_partition(i)     !! old step°ª »ç¿ë
            deltaTs=face%twall_partition(i)-cell%ts(i)
            deltaTl=face%twall_partition(i)-cell%tl(i)
            deltaTg=face%twall_partition(i)-cell%tg(i)     
            deltarho=DMAX1(0.0d0,cell%rhol(i)-cell%rhog(i))
            deltahlg=DMAX1(0.0d0,cell%hgsat(i)-cell%hl(i))       
         
            q1=heat_partition_2_bi_lo_err(deltaTs,deltaTl,deltaTg,deltarho,deltahlg &
                ,d_depart_init,fr_lift,fr_nden,fr_depfreq,fr_depdia,i,nfcondition,ra_coeff,face%twall_partition(i))
!            
         ELSEIF(nfcondition.ge.-nin_max)THEN     
!
!........AT constant temperature (nbcon=-3,-4)
!
            face%twall_partition(i)=twall(-nfcondition)      
            IF(face%twall_partition(i).lt.cell%ts(i))THEN         
               q1cell(i)=0.0d0
               qqcell(i)=0.0d0
               qecell(i)=0.0d0
               qclcell(i)=0.0d0
               qcgcell(i)=0.0d0
               ndensitycell(i)=0.0d0
               cell%Ddepart(i)=0.0d0
               cell%Dlift(i)=0.0d0   
               cell%twall(i)=cell%tl(i)    
               face%twall_partition(i)=cell%tl(i)
               face%wall_fluxl_diff(i)=0.0d0
               face%wall_fluxg_diff(i)=0.0d0
               face%wall_fluxd_diff(i)=0.0d0
               face%ddepartw(i)=0.0d0
               face%ratio_evap(i)=0.0d0
               face%bfreq(i)=0.0d0            
!              gamma_wall(i)=0.0d0   
               GOTO 100
!              RETURN      
            ENDIF
            deltaTs=face%twall_partition(i)-cell%ts(i)
            deltaTl=face%twall_partition(i)-cell%tl(i)
            deltaTg=face%twall_partition(i)-cell%tg(i)     
            deltarho=DMAX1(0.0d0,cell%rhol(i)-cell%rhog(i))
            deltahlg=DMAX1(0.0d0,cell%hgsat(i)-cell%hl(i))       
         
            q1=heat_partition_2_bi_lo_err(deltaTs,deltaTl,deltaTg,deltarho,deltahlg &
              ,d_depart_init,fr_lift,fr_nden,fr_depfreq,fr_depdia,i,nfcondition,ra_coeff,face%twall_partition(i))
              
         
!            IF(myrank.eq.0)THEN
!               PRINT*,'Heat partition using bi-section is not available in constant temperature boundary!'
!               PAUSE
!               STOP
!            ENDIF

         ELSE
!
!........At constant heat flux (nbcon=-5) 
!
            q1=qwall_solid(-nfcondition)         
            deltarho=DMAX1(0.0d0,cell%rhol(i)-cell%rhog(i))
            deltahlg=DMAX1(0.0d0,cell%hgsat(i)-cell%hl(i))
!
!...........Set all parameters Zero if wall heat flux is zero
!         
            IF(q1.le.0.0d0)THEN         
               q1cell(i)=0.0d0
               qqcell(i)=0.0d0
               qecell(i)=0.0d0
               qclcell(i)=0.0d0
               qcgcell(i)=0.0d0
               ndensitycell(i)=0.0d0
               cell%Ddepart(i)=0.0d0
               cell%Dlift(i)=0.0d0   
               cell%twall(i)=cell%tl(i)    
               face%twall_partition(i)=cell%tl(i)
               face%wall_fluxl_diff(i)=0.0d0
               face%wall_fluxg_diff(i)=0.0d0
               face%wall_fluxd_diff(i)=0.0d0
               face%ddepartw(i)=0.0d0
               face%ratio_evap(i)=0.0d0
               face%bfreq(i)=0.0d0            
!              gamma_wall(i)=0.0d0   
               GOTO 100
!              RETURN      
            ENDIF
!
!........Wetted area fraction (1-dry area fraction) calculation
!     
            IF(udfl_psbt_cfx_model)THEN
               dry_weight(i)=0.0d0
            ELSE
               dry_weight(i)=DMAX1(0.0d0,DMIN1(1.0d0,(cell%alphag(i)-0.95d0)/(0.99d0-0.95d0)))
            ENDIF 
!
!..........Bi-section method until (q1-q2(twall))/q1 < 1.d-6, err=q1-q2
!                  
            DO
            IF(ITER.eq.1)THEN
!
!..............Prepare bi-section by obtaining two end values (tw1,tw2) of 'err(tw1)*err(tw2)<0.0'
!                  

               tw1=face%twall_partition(i)
               IF(time.lt.1.0d-4)tw1=cell%ts(i)+0.5d0
          2    deltaTs=DMAX1(0.0d0,tw1-cell%ts(i))
               deltaTl=(tw1-cell%tl(i))
               deltaTg=(tw1-cell%tg(i))
!
!..............First calculation of err=q1-q2=assigned heat-calculated heat              
!
               err1=heat_partition_2_bi_lo_err(deltaTs,deltaTl,deltaTg,deltarho,deltahlg &
                                    ,d_depart_init,fr_lift,fr_nden,fr_depfreq,fr_depdia,i,nfcondition,ra_coeff,tw1) 
               IF(DABS(err1)/DMAX1(1.0d0,q1).le.1.0d-6)THEN
                  tw3=tw1
                  EXIT
               ENDIF   
               deltaTs_old=tw1-cell%ts(i)
               deltaTl_old=deltaTl
               deltaTg_old=deltaTg
!
!..............Calculate dtwall which is added to twall                
!
          1    IF(err1.ge.0)THEN
                  dtwall=count*0.5d0
               ELSE
                  dtwall=-1.0d0*count*0.5d0
               ENDIF
!
!..............Update parameters                            
!
               tw2=tw1+dtwall  
               IF(tw2.lt.273.15)THEN
                  dtwall=273.15-tw1
                  tw2=tw1+dtwall
               ENDIF 
               deltaTs=DMAX1(0.0d0,deltaTs_old+dtwall)
               deltaTl=deltaTl_old+dtwall
               deltaTg=deltaTg_old+dtwall
!
!..............Second calculation of err=q1-q2(twall)             
!
               err2=heat_partition_2_bi_lo_err(deltaTs,deltaTl,deltaTg,deltarho,deltahlg &
                                       ,d_depart_init,fr_lift,fr_nden,fr_depfreq,fr_depdia,i,nfcondition,ra_coeff,tw2)
               IF(DABS(err2)/DMAX1(1.0d0,q1).le.1.0d-6)THEN
                  tw3=tw2
                  EXIT
               ENDIF                                       
!                   
!..............Iteration until err1(tw1) and err2(tw2) have the opposite sign                   
!
               IF(err1*err2 .gt. 0.) THEN
                   count=count+1
                   IF(count.gt.2000)THEN
                      WRITE(*,*)'Error in the first step in Bi-section iteration of HPM'
                      WRITE(*,*)'time, cell, tw1, tw2, err1, err2, dry_weight'
                      WRITE(*,*) time,jperm(i),tw1,tw2,err1,err2,dry_weight(i)
                      WRITE(unit_log,*)'Error in the first step in Bi-section iteration of HPM'
                      WRITE(unit_log,*)'time, cell, tw1, tw2, err1, err2, dry_weight'
                      WRITE(unit_log,*) time,jperm(i),tw1,tw2,err1,err2,dry_weight(i)                    
                      dtwall=0.0d0
                      tw2=tw1
                      STOP
                      !GOTO 3
                   ENDIF
                   GOTO 1
               ENDIF
!
!..............Third calculation of err using tw3 between tw1 and tw2
!
             3 erra=err1
               errb=err2    
               twa=tw1
               twb=tw2          
               tw3=tw1+(tw2-tw1)/(err2-err1)*(0.0d0-err1)
               IF(isnan(tw3))tw3=tw1+(tw2-tw1)/2.0d0
               deltaTs=DMAX1(0.0d0,tw3-cell%ts(i))
               deltaTl=(tw3-cell%tl(i))
               deltaTg=(tw3-cell%tg(i))
               err3=heat_partition_2_bi_lo_err(deltaTs,deltaTl,deltaTg,deltarho,deltahlg &
                                       ,d_depart_init,fr_lift,fr_nden,fr_depfreq,fr_depdia,i,nfcondition,ra_coeff,tw3)
               twc=tw3
               errc=err3                                               
            ELSE
!
!..............Do bi-section with tw1,tw2,tw3 of err1*err2<0, from iter.ge.2.
!            
               IF(DABS(err3)/DMAX1(1.0d0,q1).le.1.0d-6) EXIT
               IF(DABS(dtwall/tw3).lt.1.d-4 .and. iter.gt.2) EXIT
               dtwall=tw3  
               IF(err3*err2 .lt. 0.) THEN !kill point1; point3 -> point1; 
                 err1=err3
                 tw1=tw3
                 tw3=tw1+(tw2-tw1)/(err2-err1)*(0.d0-err1) 
                 IF(isnan(tw3))tw3=tw1+(tw2-tw1)/2.0d0                 
               ELSEIF(err3*err1 .lt. 0.) THEN !kill point2; point3->point2;
                 err2=err3
                 tw2=tw3
                 tw3=tw1+(tw2-tw1)/(err2-err1)*(0.d0-err1)
                 IF(isnan(tw3))tw3=tw1+(tw2-tw1)/2.0d0                 
               ENDIF 
               dtwall=dtwall-tw3
               deltaTs=DMAX1(0.0d0,tw3-cell%ts(i))
               deltaTl=(tw3-cell%tl(i))
               deltaTg=(tw3-cell%tg(i))
               err3=heat_partition_2_bi_lo_err(deltaTs,deltaTl,deltaTg,deltarho,deltahlg &
                                       ,d_depart_init,fr_lift,fr_nden,fr_depfreq,fr_depdia,i,nfcondition,ra_coeff,tw3)
            ENDIF  
            iter=iter+1
!
!...........Stop bi-section when the iteration exceeds tolerance.
!            
            IF(iter.gt.iter_limit)THEN
               IF(myrank.eq.0)THEN
                  WRITE(*,*)'Convergence failure in heat_partition_bi'
                  WRITE(*,*)'time,cell,twa,tw3,twb,erra,err3,errb,voidg'
                  WRITE(*,*) time,jperm(i),twa,tw3,twb,erra,err3,errb,cell%alphag(i)
                  WRITE(unit_log,*)'Convergence failure in heat_partition_bi'
                  WRITE(unit_log,*)'time,cell,twa,tw3,twb,erra,err3,errb,voidg'
                  WRITE(unit_log,*) time,jperm(i),twa,tw3,twb,erra,err3,errb,cell%alphag(i)                  
               ENDIF 
               
               tw3=face%twall_partition(i)
               IF(dry_weight(i).eq.1.0d0)THEN
                  qq=0.0d0
                  qe=0.0d0
                  qcl=0.0d0
                  qcg=qcgcell(i)  
                  IF(time.lt.1.0d-4)qcg=q1   
               ELSEIF(dry_weight(i).eq.0.0d0)THEN
                  qq=qqcell(i)
                  qe=qecell(i)
                  qcl=qclcell(i)
                  qcg=0.0d0  
                  IF(time.lt.1.0d-4)qcl=q1           
               ELSE
                  qq=qqcell(i)/(1.0d0-dry_weight(i))
                  qe=qecell(i)/(1.0d0-dry_weight(i))
                  qcl=qclcell(i)/(1.0d0-dry_weight(i))
                  qcg=qcgcell(i)/dry_weight(i) 
                  IF(time.lt.1.0d-4)THEN
                     qcl=q1/(1.0d0-dry_weight(i))
                     qcg=q1/dry_weight(i)               
                  ENDIF
               ENDIF 
               EXIT
            ENDIF 
            ENDDO
!          
         ENDIF
!
!........Calculate heat assinged to cell
!
         q1cell(i)=q1
         qqcell(i)=qq*(1.0d0-dry_weight(i))
         qecell(i)=qe*(1.0d0-dry_weight(i))
         qclcell(i)=qcl*(1.0d0-dry_weight(i))
         qcgcell(i)=qcg*dry_weight(i)
         ndensitycell(i)=ndensity
         IF(d_depart.ge.0) cell%Ddepart(i)=d_depart
         cell%Dlift(i)=dlo   
         cell%twall(i)=tw3
!  
!........Calculate diffusive heat flux which will be used in energy diffusion
!    
         IF(iter.lt.iter_limit.and.nfcondition.ne.-3.and.nfcondition.ne.-4) face%twall_partition(i)=tw3

         face%wall_fluxl_diff(i)=(1.0d0-dry_weight(i))*(qq+qcl)*sa_walll(i)
         face%wall_fluxg_diff(i)=dry_weight(i)*qcg*sa_walll(i)      
         face%wall_fluxd_diff(i)=0.0d0
         face%ddepartw(i)=d_depart
         face%bfreq(i)=bfreq
!
         IF(q1.gt.0.0d0)THEN
            face%ratio_evap(i)=qe*(1.0d0-dry_weight(i))/q1
         ELSE
            face%ratio_evap(i)=0.0d0
         ENDIF
!
         IF((1.0d0-dry_weight(i))*qe.gt.0.0d0)THEN
            hi_gas=cell%hgsat(i)
            hi_liq=cell%hl(i)
            gamma_wall(i)=gamma_wall(i)+(1.0d0-dry_weight(i))*qe/(hi_gas-hi_liq)*sa_walll(i)/volp(i)
!
!...........IAT source by nuculate boiling
!
            IF(iat.gt.0)IAT_nucl(i)=(1.0d0-dry_weight(i))*pi*d_depart**2*ndensity*bfreq*a_two*sa_walll(i)/volp(i)         
            IF(fr_lift.eq.1)IAT_nucl(i)=pi*dlo**2*ndensity*bfreq*a_two*sa_walll(i)/Volp(i)*ra_coeff                 
         ENDIF
100      CONTINUE
      ENDDO
!
      END SUBROUTINE heat_partition_2_bi_lo
