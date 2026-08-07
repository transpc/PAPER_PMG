!
      SUBROUTINE rv_gap_conductance
!
      USE Vol_Data       , ONLY: cell
      USE Zconst2        , ONLY: dt
      USE STM_TBL_cupid  , ONLY: wmole,thcax,thcbx
      USE Zrv_gap_cond   , ONLY: block,nr_gapi,nr_gapo,hte,GapIntW,irupt,iplas,cltave,dtdt,               &
                                 hte_clad_o,StrPlas,GapDisp_fission,GapWdth,CladExR,ts_old,tg_old,tl_old, &
                                 cond_gap,h_gap,width_gap,block_gap,Gap_P0,GapRough,                      &
                                 ngas,indgas,molgas
      USE Zrv_hts_2d     , ONLY: nr_2d,t_fuel,ri_2d,ri_2d_opt,ri_2d_input 
      USE Zrv_ncell      , ONLY: ncell_fuel_rod,cupid_cell_hts2d,nrod_fuel_rod
!
      IMPLICIT NONE
!
      INTEGER i,j,inde,kk,ii,indy,k,jj,ndivs,No_r
      REAL(8) sig,utf,eplas,emissc,emissf 
      REAL(8) gprinc1,gprinc2,xnu,ey,term1,term2,pgas,pfluid
      REAL(8) hte_clad,term3,ucl,gapwid,rclad,rfuel,xkg,hgap,fe,hrad,htot,tfuel,tclad
      REAL(8) sqrtmw(7)
      REAL(8) xkg_i,xkg_k,wmole_i,wmole_k,phij,phij_up,phij_down
      LOGICAL initial_read,initial_irupt,initial_hte
      LOGICAL,SAVE:: initial_temp
!                                                                       
!.....Sqrt of molecular weight: he, h2, n2, kr, xe, o2, ar, 
!                  
      DATA sqrtmw/2.00065d0, 1.41982d0, 5.29277d0, 9.15423d0, 11.45862d0, 5.65680d0, 6.32044d0/                  
!
!.....NCG property: Helium, Hydrogen, Nitrogen, Krypton, Xenon, Air, Argon
!     wmolea=molecular mass of non-condensible gas. = sqrtmw^2 (that is, sqrtmw=sqrt(wmolea))
!     thcax=A constants
!     thcbx=B constatns
!
!.....Stefan-Boltzmann constant (w/m**2-k**4)                              
!
      DATA sig/5.6697d-8/     
!
!.....Logical variale
!    
      DATA initial_read,initial_irupt,initial_hte,initial_temp/.TRUE.,.TRUE.,.TRUE.,.TRUE./
!
!.....initialize: need to be input variables
!
      CALL rv_read_gapcond(initial_read)
!
!.....Initialize old temperatures for pgas
!
      IF(initial_temp) THEN
         ALLOCATE(ts_old(ncell_fuel_rod),tg_old(ncell_fuel_rod),tl_old(ncell_fuel_rod)) 
      ENDIF
!
!.....Check the rod is ruptured or not with blockage ratio
!
      IF(initial_irupt) THEN
         block=0.0
         irupt=0
         StrPlas=0.d0
         initial_irupt=.FALSE.
         width_gap(:)=1000.d0
         CladExR=0.d0         
      ENDIF
!      IF(block.gt.0.d0) THEN
!         irupt=1
!      ENDIF
!
!.....Calculate average fuel temperature
!
      IF(initial_hte) THEN
         ALLOCATE(hte(nr_2d),GapIntW(nr_2d))
         hte(:)=0.d0
         initial_hte=.FALSE.
      ENDIF
!
!      
      DO kk=1,ncell_fuel_rod                             !ncell_fuel_rod=nrod_2d*nz_2d*nz_fine (total rod number divied by z direction)
!
!........geometrical data 
         IF(ri_2d_opt.eq.1)THEN
            No_r=nrod_fuel_rod(k)
            ri_2d(:)=ri_2d_input(No_r,:)            
         ENDIF
!         
         hte_clad_o=hte(nr_gapo)  !averaged clad temperature for old time step  !If cald has more than two division, this should be modified
         DO i=1,nr_2d-1   
            hte(i)=0.5d0*(t_fuel(kk,i)+t_fuel(kk,i+1))
            GapIntW(i)=ri_2d(i+1)-ri_2d(i)
         ENDDO   
!
         tfuel=t_fuel(kk,nr_gapi)
         tclad=t_fuel(kk,nr_gapo)      
!
!........Not ruptured Case                                                                       
!
         irupt=0
!         IF(block_gap(kk).gt.0.d0.or.width_gap(kk).le.1.d-5) THEN
         IF(block_gap(kk).gt.0.d0) THEN
            irupt=1
         ENDIF  
!
         IF(irupt.eq.0)then 
!
!...........UTF: radial displacement of fuel by thermal expansion
!     
            inde=nr_gapi-1
            utf=0.d0
            DO i=1,inde
               utf=utf+(1.0d-5*hte(i)-3.0d-3+4.0d-2*DEXP(-5.0d3/hte(i)))*GapIntW(i)   !UTF
            ENDDO
!
!...........Computation of cladding deformation.                                 
!
            inde=nr_gapo
            gprinc1=ri_2d(inde)                               !ri=inner radius of clad
            gprinc2=ri_2d(inde+1)                             !ro=outer radius of clad            
! 
            xnu=gprinc1**2                                    !ri**2 
            ey=gprinc2**2                                     !ro**2 
!            
            ii=cupid_cell_hts2d(kk)                        
!!        
!!           save the inital temperature values     
!            IF(initial_temp) THEN
!               ts_old(kk)=cell%ts(ii)
!               tg_old(kk)=cell%Tg(ii)
!            ENDIF             
!            pfluid=cell%p(ii)                                 !coolant pressure
!            pgas=gap_P0*MAX(cell%ts(ii),cell%Tg(ii))/MAX(ts_old(kk),tg_old(kk))  !gap pressure
!        
!           save the inital temperature values     
            IF(initial_temp) THEN
               ts_old(kk)=cell%ts(ii)
               tl_old(kk)=cell%tl(ii)
!               ts_old(57)=cell%ts(19)
!               tl_old(57)=cell%tl(19)
            ENDIF             
            pfluid=cell%p(ii)                                 !coolant pressure
            pgas=gap_P0*MAX(cell%ts(ii),cell%tl(ii))/MAX(ts_old(kk),tl_old(kk))  !gap pressure
!            pgas=gap_P0*MAX(cell%ts(19),cell%tl(19))/MAX(ts_old(57),tl_old(57))  !gap pressure
            IF(irupt.eq.1) pgas=pfluid
!            
            term1=(pgas*gprinc1-pfluid*gprinc2)/(gprinc2-gprinc1)   !sig_n
            term2=(pgas*xnu-pfluid*ey)/(ey-xnu)                     !sig_z                    
!
!...........Young's modulus (ey) and poisson ratio (xnu) for cladding.     
!
            hte_clad=hte(nr_gapo) !clad avg. temperature
            CALL rv_celmdr(hte_clad,ey,xnu)  
!
!...........Cladding thermal expansion (eps_TC)
!
            CALL rv_cthxpr(hte_clad,term3)   !term3=eps_TC                     
!     
!...........Plastic deformation: Ep, Blockage ratio
!       
            eplas=0.0d0 
            block=0.0d0
            cltave=hte(nr_gapo)      !avg. clad temperature
            dtdt=(hte(nr_gapo)-hte_clad_o)/dt
            IF(iplas.eq.1) CALL rv_cplexp(cltave,dtdt,term1,irupt,eplas,block) !eplas=Ep, block=blockage ratio
!
!...........Plastic strain into change of clad radius.                           
!
            eplas=MAX(eplas,StrPlas)
            StrPlas=eplas
            ucl=(term3+eplas+(term1-xnu*term2)/ey)*(gprinc1+gprinc2)*0.5d0   !(eps_TC+Ep+ue/rcm)*rcm=uTC+ucc+ue=uC 
!
!...........Average gap width.                                                   
!
            inde=nr_gapi
            ey=4.0d-3*GapIntW(inde)
            gapwid=MAX(GapIntW(inde)-utf+ucl+GapDisp_fission,ey) !tg=to-uF+uC (to=GapIntW, uF=uTF+us+ur, uC=ucl, )
!
!...........Calculate clad and fuel radii.                                       
!
            rclad=GapIntW(nr_gapo)+ucl  
!
!...........Do not let the gap width or clad radius diminish after plastic strain begins.
!
            rfuel=rclad-gapwid 
!
!...........Save the gap width and clad outer radius for editing in majout.      
!
            GapWdth=gapwid 
!
!...........Keep the sign on cladex if it is ever set.                           
!
            CladExR=sign(gprinc2+ucl,CladExR)
!
!...........Set clad radius negative if it has burst.                            
!
           IF(block.gt.0.0d0.and.CladExR.gt.0.0d0) CladExR=-CladExR
!
!........Ruptured Case                                                                       
!
         ELSE
            inde=nr_gapo
            ey=4.0d-3*GapIntW(inde) 
            gapwid=GapWdth
            rclad=DABS(CladExR)
            rfuel=rclad-gapwid
         ENDIF
!
!........Gap thermal conductivity model: xkg
!
         IF(ngas.eq.1) THEN
            k=indgas(1)
            xkg=thcax(k)*(hte(nr_gapi))**thcbx(k)
         ELSE
            xkg=0.d0
            DO i=1,ngas
               k=indgas(i)
               xkg_i=thcax(k)*(hte(nr_gapi))**thcbx(k)
               wmole_i=wmole(k)
               phij=0.d0
               DO j=1,ngas
                  IF(j.eq.i) CYCLE
                  jj=indgas(j)
                  xkg_k=thcax(kk)*(hte(nr_gapi))**thcbx(jj)
                  wmole_k=wmole(jj)
                  phij_up=(1.d0+DSQRT(xkg_i/xkg_k)*(wmole_i/wmole_k)**0.25)
                  phij_up=phij_up*phij_up
                  phij_down=DSQRT(8.d0*(1.d0+wmole_i/wmole_k))
                  phij=phij+phij_up/phij_down*molgas(jj)
               ENDDO
               xkg=xkg+xkg_i*molgas(k)/(molgas(k)+phij)
            ENDDO
         ENDIF
!
!........Now have geometry figured, so calculate conductance.                 
!        Temperature jump distance.                                           
!
         indy=nr_gapi  !For gas temperature in the Gap
         term1=0.425d0-2.3d-4*MIN(hte(indy),1000.0d0)    !aHe
         term2=0.740d0-2.5d-4*MIN(hte(indy),1000.0d0)    !aXe
         term3=0.0d0 
         DO i=1,ngas 
            k=indgas(i) 
            term3=term3+molgas(i)*(term1+7.85562d-3*(sqrtmw(k)-sqrtmw(1))  &                !MARS code.  USE THIS!
                 *(sqrtmw(k)+sqrtmw(1))*(term2-term1))/sqrtmw(k)                 
!           term3=term3+molgas(k)*(term1+(wmole(k)-wmole(1))/(wmole(5)-wmole(1))  &        !MARS report: fi*ai*Mi^-0.5 for (g1+g2)
!                *(term2-term1))/sqrtmw(k)                                    
         END DO 
!         
         xnu=0.024688d0*xkg*DSQRT(hte(indy))/(pgas*term3)+3.2d0*GapRough   !(g1+g2)+3.2(RF+RC)=0.024688*kg*T^0.5 / (Pg*fi*ai*Mi^-0.5) +3.2(RF+RC)
!
!........Calculation of effective gap conductace, using the deformation information.      
!
         term1=0.0d0 
         term3=ey*0.5d0 
         ey=2.0d0*gapwid 
         ndivs=8
         DO i=1,ndivs 
!        
!           tn=tg+[-1+(2n-1)/N]t0         
            term2=MIN(MAX(term3,gapwid+(-1.0d0+(DFLOAT(i)*2.0d0-1.0d0)*0.125d0)*GapIntW(nr_gapi)),ey)   
!
!           1/(tn+(g1+g2)+3.2(RF+RC))             
            term1=term1+1.0d0/(term2+xnu)   
         END DO 
!
!        hg=kg/N * 1/(tn+(g1+g2)+3.2(RF+RC))  
!         hgap=term1*0.125d0*xkg 
         hgap=term1*1.d0/DFLOAT(ndivs)*xkg 
!
!........Set fuel and clad surface emissivities to reasonable value.          
!
!        F=1 / [ 1/esp_f + Rf/Rc*(1/(esp_c-1)) ]
         emissf=0.6d0 
         emissc=0.6d0 
!         fe=1.0d0/(1.0d0/emissf+(rfuel/rclad)*(1.0d0/(emissc-1.0d0)))   !MARS report
         fe=1.0d0/(1.0d0/emissf+(rfuel/rclad)*(1.0d0/emissc-1.0d0))      !MARS code: bug . BUT! USE THIS!
!
!        hr=sig*F*(Tf^2+Tc^2)(Tf+Tc)         
         hrad=sig*fe*(tfuel**2+tclad**2)*(tfuel+tclad) 
!
!        htot         
         htot=hgap+hrad 
!
!        Effective conductivity.                                              
!        k * initial width of the gap.                                        
!
         h_gap(kk)=htot
!         h_gap(kk)=hgap
         xkg=htot*GapIntW(nr_gapi)
         cond_gap(kk)=xkg
!
         width_gap(kk)=GapWdth 
         block_gap(kk)=block        

         IF(irupt.ne.0) print*,'irupt=',irupt,kk
      ENDDO
      initial_temp=.false.
!           
      END SUBROUTINE rv_gap_conductance



