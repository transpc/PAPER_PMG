!
      SUBROUTINE udfn_porous_user(ncell,x,vol,porosity,nmaterial,nzone,sl,hydraulicd,icore)
!
!.....This routine change the cell value of somaGrid.
!
      USE Zcore         , ONLY: np,myrank
      USE Zzone         , ONLY: ncell_fluid_all,ncell_cond_all
      USE Zparam        , ONLY: nn,ndim,pi,mesh_openfoam
      USE Znum_cell     , ONLY: i_neigh_tmp,j_neigh_tmp, & 
                                perm_tmp1
      USE Zconst1       , ONLY: vv_prob
      USE Zconst2       , ONLY: hydraulicd_init
      USE Zio_unit      , ONLY: unit_log
      USE unitManager   , ONLY: createUnit      
      USE Zrv_subchan   , ONLY: subchannel_type_tmp    
      USE KSMR          , ONLY: vol_sg       
      ! fwl - ldh
!
      IMPLICIT NONE
!     
!     input
      INTEGER ncell
      REAL(8) x(nn,ndim),vol(nn)
!     output
      INTEGER nmaterial(nn),nzone(nn),icore(nn)
      REAL(8) porosity(nn)
      REAL(8) sl(nn,ndim),hydraulicd(nn) 
!     local variables
      INTEGER i,j,k,na,runit,n_core
      INTEGER err
      REAL(8) r
!     local arrray
      INTEGER :: itmp(2)
      REAL(8) permeability(nn,ndim)
!     apr1400_lbloca
      INTEGER num_ch,nz0_2d,ch_opt,io
      INTEGER,ALLOCATABLE::nf_input(:)
      INTEGER,ALLOCATABLE::num_core(:)
!     halden650_5
      REAL(8) :: r1,r2,r3,r4,depth,w1,w2,w3,w4,x1,x2,x3,x4
      
      ! fwl - ldh

      IF(myrank.eq.0) THEN
!
!.....NuScale
!      
      IF(vv_prob.eq.'Nuscale-RVV')THEN
         OPEN(222,file='rv_core.in')
         READ(222,*) n_core
         ALLOCATE(num_core(n_core))
         DO i=1,n_core
            READ(222,*) num_core(i) !global core index
         ENDDO
         CLOSE(222)
         icore(:)=0
         DO i=1,nn
            DO j=1,n_core 
               IF(i.eq.j) icore(i)=1  !icore=1 --> RV model, icore=0 --> simple model
            ENDDO
         ENDDO
      ENDIF       
!  
!.....KSMR
!      
      IF(vv_prob.eq.'KSMR')THEN
            hydraulicd_init=1            
         WRITE(*,*)'>>>>> sl, porosity, permeability, hd is calculating for KSMR'
         DO i=1,nn 
            porosity(i)=1.d0
            permeability(i,:)=1.d0
            sl(i,:)=1.d0 
            hydraulicd(i)=0.012729d0
            IF(subchannel_type_tmp(i).eq.1) THEN !center
               porosity(i)=0.56753d0
               sl(i,1)=76.499d0
               sl(i,2)=76.499d0
               permeability(i,1)=1.0d0
               permeability(i,2)=1.0d0
               permeability(i,3)=0.56753d0
               hydraulicd(i)=0.012729d0
            ELSEIF(subchannel_type_tmp(i).eq.2) THEN !side
               porosity(i)=0.78376d0
               sl(i,1)=76.499d0
               sl(i,2)=76.499d0
               permeability(i,1)=1.0d0
               permeability(i,2)=1.0d0
               permeability(i,3)=0.78376d0
               hydraulicd(i)=0.018923d0
            ELSEIF(subchannel_type_tmp(i).eq.3) THEN !center               
               porosity(i)=0.89188d0
               sl(i,1)=76.499d0
               sl(i,2)=76.499d0
               permeability(i,1)=1.0d0
               permeability(i,2)=1.0d0
               permeability(i,3)=0.89188d0
               hydraulicd(i)=0.018055d0
            ENDIF       
         ENDDO            
      ENDIF
!      
      IF(vv_prob.eq.'KSMR-SG-porous')THEN
!        total volume of SG
         vol_sg=0.d0
         DO i=1,nn 
            IF(i.ge.8396.and.i.le.10595) THEN !SG 
               vol_sg=vol_sg+vol(i) 
            ENDIF
         ENDDO    
!         
!        porosity of SG         
         DO i=1,nn 
            porosity(i)=1.d0
            permeability(i,:)=1.d0
            sl(i,:)=1.d0 
            hydraulicd(i)=1.d0
            IF(i.ge.8396.and.i.le.10595) THEN !SG 
               porosity(i)=30.822d0/vol_sg
               permeability(i,:)=1.d0
               sl(i,:)=1.d0 
               hydraulicd(i)=4.d0*(9.65289d0-4.16575d0)/555.434d0 !4Af/P=0.03  (Af=Atot-Asolid)
            ENDIF
         ENDDO            
      ENDIF             
!
!.....sgp_separator
!      
      IF(vv_prob.eq.'sgp_separator')THEN
         OPEN(222,file='porosity.in')
         DO i=1,nn
            READ(222,*) porosity(i)
         ENDDO
         CLOSE(222)
      ENDIF     
!     
!.....apr1400_lbloca
!
      IF(vv_prob.eq.'apr1400_lbloca'   .or.&
         vv_prob.eq.'opr1000_mc_rv_lbloca' )THEN
         
         runit=createUnit('rv_core.in')
         runit=4
         OPEN(runit,file='rv_core.in',status='old',iostat=err)
         porosity(:)=1.0d0
         perm_tmp1(:)=1.0d0
         IF(err.eq.0)THEN
             READ(runit,*)na
             IF(na.ne.ncell)THEN
                IF(myrank.eq.0)WRITE(*,*)'Error! na and ncell in rv_core.in!',na,ncell
                CLOSE(runit)
                PAUSE
                STOP          
             ENDIF    
             IF(myrank.eq.0)WRITE(*,"(11x,a)")'Reading icore, porosity and permeability in udfn_porous_user...'
             DO i=1,ncell
                READ(runit,*)icore(i),porosity(i),(perm_tmp1(j),j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1) !1.core,2.reflector,3.barrel-in,4.lower plenum,5.upper plenum
                IF(0)THEN !no solid or ncell_cond_all=0 due to nmaterial==0    
                   IF(icore(i).ne.1.and.porosity(i).lt.1.0d0)THEN !apr1400_lbloca_debug_porous
                      nzone(i)=2
                      nmaterial(i)=-4
                   ENDIF
                ENDIF   
             ENDDO
             CLOSE(runit)
         ELSE
            icore(:)=0
            IF(myrank.eq.0)WRITE(*,"(11x,a)")'Use ht_str_2d.in instead of rv_core.in!'
            runit=createUnit('ht_str_2d.in')
            runit=52
            OPEN(runit,file='ht_str_2d.in',status='old',iostat=io)
            IF(io.gt.0) RETURN
            READ(runit,*) num_ch,nz0_2d,ch_opt
            IF(nz0_2d.eq.0)RETURN
            IF(myrank.eq.0)WRITE(*,"(11x,a)")'Reading ht_str_2d.in...'
            IF(myrank.eq.0)WRITE(unit_log,"(11x,a)")'Reading ht_str_2d.in...'      
            IF(ch_opt.eq.1)THEN
               DO i=1,num_ch
                  READ(runit,*) itmp(1)
               ENDDO
            ELSE
               READ(runit,*) itmp(1)
            ENDIF     
            ALLOCATE(nf_input(nz0_2d))
            DO i=1,num_ch
               READ(runit,*) (nf_input(j),j=1,nz0_2d)
               DO j=1,nz0_2d
                  icore(nf_input(j))=1
                  porosity(nf_input(j))=0.54d0 !make it diverge
               ENDDO 
            ENDDO
            DEALLOCATE(nf_input)
            CLOSE(runit)
         ENDIF
!
         IF(mesh_openfoam.eq.1)THEN
           IF(myrank.eq.0)WRITE(*,"(11x,a)")'Skip setting porosity to D.C. outer wall.'
           ELSE
           IF(myrank.eq.0)WRITE(*,"(11x,a)")'Skip setting porosity to D.C. outer wall.'
         ENDIF          
!
         RETURN
      ENDIF
!     
!.....fs_31701_3D
!        
      IF(vv_prob.eq.'fs_31203_3D'    .or.&
         vv_prob.eq.'fs_31302_3D'    .or.&            
         vv_prob.eq.'fs_31701_3D'    .or.&            
         vv_prob.eq.'fs_31805_3D'    ) THEN
         IF(myrank.eq.0)WRITE(*,"(11x,a)")'Set porosity and permeability of fs_31701_3D to 0.54!'
         IF(myrank.eq.0)WRITE(unit_log,"(11x,a)")'Set porosity and permeability of fs_31701_3D to 0.54!'
         DO i=1,nn
            porosity(i)=0.54d0
!            DO j=1,num_neigh(i)           
!               IF(ABS(xn(j,i,ndim)).gt.0.9d0)THEN
!                  perm(j,i)=0.54d0
!               ELSE   
!                  perm(j,i)=1.0d0
!               ENDIF   
!            ENDDO   
         ENDDO
      ENDIF               
!     
!.....icarus2002
!        
      IF(vv_prob.eq.'icarus2002') THEN          
          DO i=1,nn
            IF(nzone(i).le.3) THEN
               porosity(i)=0.687d0  
               hydraulicd(i)=0.009886d0                
            ELSE
               porosity(i)=1.0d0  
            ENDIF
            permeability(i,:)=1.d0
         ENDDO          
      ENDIF      
!     
!.....CUPID-MARS: ATLAS
!        
      IF(vv_prob.eq.'atlas_mc'.or.vv_prob.eq.'atlas')THEN
         DO i=1,nn
            r=DSQRT(x(i,1)**2+x(i,2)**2)
            IF(r.le.0.168d0)then
               porosity(i)=0.787d0
            ELSE
               porosity(i)=1.0d0   
            ENDIF
         ENDDO   
      ENDIF   
!
!.....PAFS-POOL
!
      IF(vv_prob.eq.'PAFS-POOL') THEN
          DO i=1,ncell
             DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
                k=j_neigh_tmp(j)
                IF(i.eq.1700.and.(k.eq.1811.or.k.eq.1806)) THEN
                   perm_tmp1(j)=1.d0
                ELSEIF(i.eq.1702.and.(k.eq.1812.or.k.eq.1807)) THEN
                   perm_tmp1(j)=1.d0
                ELSEIF(i.eq.1705.and.(k.eq.1813.or.k.eq.1808)) THEN
                   perm_tmp1(j)=1.d0            
                ELSEIF(i.eq.1709.and.(k.eq.1814.or.k.eq.1809)) THEN
                   perm_tmp1(j)=1.d0            
                ELSEIF(i.eq.1713.and.(k.eq.1815.or.k.eq.1810)) THEN
                   perm_tmp1(j)=1.d0                           
                ELSEIF(i.eq.1683.and.(k.eq.1801.or.k.eq.1796)) THEN
                   perm_tmp1(j)=1.d0
                ELSEIF(i.eq.1685.and.(k.eq.1802.or.k.eq.1797)) THEN
                   perm_tmp1(j)=1.d0                       
                ELSEIF(i.eq.1688.and.(k.eq.1803.or.k.eq.1798)) THEN
                   perm_tmp1(j)=1.d0                       
                ELSEIF(i.eq.1692.and.(k.eq.1804.or.k.eq.1799)) THEN
                   perm_tmp1(j)=1.d0                       
                ELSEIF(i.eq.1696.and.(k.eq.1805.or.k.eq.1800)) THEN
                   perm_tmp1(j)=1.d0                       
                ELSEIF(i.eq.1811.and.k.eq.1700) THEN
                   perm_tmp1(j)=1.d0
                ELSEIF(i.eq.1812.and.k.eq.1702) THEN
                   perm_tmp1(j)=1.d0         
                ELSEIF(i.eq.1813.and.k.eq.1705) THEN
                   perm_tmp1(j)=1.d0         
                ELSEIF(i.eq.1814.and.k.eq.1709) THEN
                   perm_tmp1(j)=1.d0         
                ELSEIF(i.eq.1815.and.k.eq.1713) THEN
                   perm_tmp1(j)=1.d0         
                ELSEIF(i.eq.1806.and.k.eq.1700) THEN
                   perm_tmp1(j)=1.d0
                ELSEIF(i.eq.1807.and.k.eq.1702) THEN
                   perm_tmp1(j)=1.d0         
                ELSEIF(i.eq.1808.and.k.eq.1705) THEN
                   perm_tmp1(j)=1.d0         
                ELSEIF(i.eq.1809.and.k.eq.1709) THEN
                   perm_tmp1(j)=1.d0         
                ELSEIF(i.eq.1810.and.k.eq.1713) THEN
                   perm_tmp1(j)=1.d0         
                ELSEIF(i.eq.1801.and.k.eq.1683) THEN
                   perm_tmp1(j)=1.d0
                ELSEIF(i.eq.1802.and.k.eq.1685) THEN
                   perm_tmp1(j)=1.d0         
                ELSEIF(i.eq.1803.and.k.eq.1688) THEN
                   perm_tmp1(j)=1.d0         
                ELSEIF(i.eq.1804.and.k.eq.1692) THEN
                   perm_tmp1(j)=1.d0         
                ELSEIF(i.eq.1805.and.k.eq.1696) THEN
                   perm_tmp1(j)=1.d0                        
                ELSEIF(i.eq.1796.and.k.eq.1683) THEN
                   perm_tmp1(j)=1.d0
                ELSEIF(i.eq.1797.and.k.eq.1685) THEN
                   perm_tmp1(j)=1.d0         
                ELSEIF(i.eq.1798.and.k.eq.1688) THEN
                   perm_tmp1(j)=1.d0         
                ELSEIF(i.eq.1799.and.k.eq.1692) THEN
                   perm_tmp1(j)=1.d0         
                ELSEIF(i.eq.1800.and.k.eq.1696) THEN
                   perm_tmp1(j)=1.d0                        
                ENDIF
             ENDDO
          ENDDO 
      ENDIF           
!     
!.....rbht1196_3d 
!        
      IF(vv_prob.eq.'rbht1196_3d') THEN
         DO i=1,nn
            porosity(i)=0.583302d0
            permeability(i,:)=1.d0
            !permeability(i,3)=0.583302d0
         ENDDO
      ENDIF          
!               
!.....atlas-mslb
!
      IF(vv_prob.eq.'atlas_mc_porous'.or.vv_prob.eq.'pwr_mc_poro'.or.vv_prob.eq.'apr1400_mc_poro'.or.vv_prob.eq.'opr1000_mc_poro')THEN
!
         CALL atlas_udfn_porous_user(vol,porosity,nmaterial,nzone)
            IF(np.gt.1) THEN
               itmp(1)=ncell_fluid_all
               itmp(2)=ncell_cond_all
            ENDIF
!         
      ENDIF
!
      ENDIF ! myrank
      IF(vv_prob.eq.'atlas_mc_porous'.or.vv_prob.eq.'pwr_mc_poro'.or.vv_prob.eq.'apr1400_mc_poro'.or.vv_prob.eq.'opr1000_mc_poro')THEN
!
!     atlas_udfn_porous_user changes nmaterial rank 0 only 
!
         IF(np.gt.1) THEN
            CALL broadcast_i(itmp,2)
            ncell_fluid_all=itmp(1)
            ncell_cond_all=itmp(2)
         ENDIF
      ENDIF
!             
!
!.....OPR1000 Full vessel
!
      IF(vv_prob.eq.'OPR1000_fullvessel_1x1'           .or. &
         vv_prob.eq.'opr1000_mc_rv'                    .or. &
         vv_prob.eq.'opr1000_rv'                       .or. &
         vv_prob.eq.'opr1000_rv_lbloca'                       .or. &
         vv_prob.eq.'OPR1000_fullcore_modmesh02_rv'    .or. &
         vv_prob.eq.'OPR1000_fullcore_modmesh02_rv_vessel'      )then
         IF(myrank.eq.0) THEN
            porosity=1.0d0
            permeability=1.0d0
            perm_tmp1=1.0d0
         ENDIF ! myrank
         CALL udfn_opr1000_zone_porosity(x,nzone,porosity,sl)
        
      ENDIF                  
!
!.....OPR1000 Single assembly 
!
      IF(vv_prob.eq.'OPR1000_single_assem')then
         IF(myrank.eq.0) THEN
            porosity=1.0d0
            permeability=1.0d0
            perm_tmp1=1.0d0
         ENDIF ! myrank
         CALL udfn_opr1000_snglassm_porosity(nzone,porosity,sl)
      ENDIF
!     
!
!.....APR1400 Full vessel
!
      IF(vv_prob.eq.'apr1400_mc_rv' .or. &
         vv_prob.eq.'apr1400_rv')then
         IF(myrank.eq.0) THEN
            porosity=1.0d0
            permeability=1.0d0
            perm_tmp1=1.0d0
         ENDIF
         CALL udfn_apr1400_zone_porosity(x,nzone,porosity)
      ENDIF
!
!.....Hydraulic Diameter
!     
      IF(myrank.eq.0) THEN
         IF(vv_prob.eq.'OPR1000_fullvessel_1x1'               .or. &
            vv_prob.eq.'opr1000_mc_rv'                        .or. &
            vv_prob.eq.'opr1000_rv'                           .or. &
            vv_prob.eq.'opr1000_rv_lbloca'                    .or. &
            vv_prob.eq.'OPR1000_fullcore_modmesh02_rv'        .or. &
            vv_prob.eq.'OPR1000_fullcore_modmesh02_rv_vessel' .or. &
            vv_prob.eq.'OPR1000_single_assem'                 .or. &
            vv_prob.eq.'apr1400_mc_rv' .or. &
            vv_prob.eq.'apr1400_rv'                                 ) then
            hydraulicd_init=1            
            CALL udfn_porous_hyd1(ncell,nzone,hydraulicd)
         ENDIF
      ENDIF ! myrank
!
!.....halden650_5       
!
      IF(myrank.eq.0) THEN
      IF(vv_prob.eq.'halden650_5')THEN
         WRITE(*,"(11x,a)")'--vol_tmp is changed in udfn_porous_user!'
         x1=0.014955
         x2=0.014955+0.0031
         x3=0.014955+0.0031+0.02466
         x4=0.014955+0.0031+0.024660+0.003
         r1=0.01
         r2=0.0131
         r3=0.017
         r4=0.02
         w1=0.014955
         w2=0.0031
         w3=0.02466
         w4=0.003
         depth=0.014955
         DO i=1,nn !pik-halden-debug, not exact, should consider annulas situation.
            IF(x(i,1).lt.x1)THEN
            ELSEIF(x(i,1).lt.x2)THEN !heater
               vol(i)=vol(i)/depth/w2*3.1415d0*(r2*r2-r1*r1)*0.5d0
            ELSEIF(x(i,1).lt.x3)THEN !fluid
               vol(i)=vol(i)/depth/w3*3.1415d0*(r3*r3-r2*r2)      
            ELSEIF(x(i,1).lt.x4)THEN !flask
               vol(i)=vol(i)/depth/w4*3.1415d0*(r4*r4-r3*r3)*0.5d0   
            ENDIF  
         ENDDO
      ENDIF
      ENDIF ! myrank
!      
      END SUBROUTINE udfn_porous_user
!
!-----------------------------------------------------------------------------------------------------------------------------------------      
            
      SUBROUTINE udfn_opr1000_zone_porosity(x,nzone,porosity,sl)
!
!.....This routine change the cell value of somaGrid.
!
      USE Zmpi          , ONLY: maxmt_cell
      USE Zcore         , ONLY: np,myrank       
      USE Zparam        , ONLY: nn,ndim,pi
      USE Znum_cell     , ONLY: i_neigh_tmp,j_neigh_tmp, &
                                perm_tmp1,sv_tmp1
      USE Zconst1       , ONLY: vv_prob
      USE Zmodel        , ONLY: vfwl_k,fsar
      !rod-scale
      USE Zrv_subchan   , ONLY: subchannel_type_tmp
!      
      IMPLICIT NONE
!     
!     input
      REAL(8) x(nn,ndim)
!     output      
      INTEGER nzone(nn)
      REAL(8) porosity(nn)
      REAL(8) sl(nn,ndim)
      
!     local variables
      INTEGER i,j,k,ix,iy
      INTEGER,ALLOCATABLE:: assem_ix(:),assem_iy(:)
!     local arrray
      REAL(8) radius
      REAL(8) y0,wid_fa,widx0,widy0,widx1,widy1
      INTEGER i1
      !rod-scale
      REAL(8) xn(maxmt_cell,ndim)     
      REAL(8) sa 
!
!.....OPR1000 Full vessel
!
      ALLOCATE(vfwl_k(20),fsar(16,16)) !See udfn_mom_wall
      FSAR(:,:)=0.1d0
      vfwl_k(:)=0.0d0      
      IF(myrank.eq.0) THEN
         ALLOCATE(assem_ix(nn),assem_iy(nn))
         assem_ix(:)=16
         assem_iy(:)=16
         wid_fa=0.208736
         y0=wid_fa*7.d0+wid_fa*0.5d0
         DO i=1,nn
            IF(x(i,3).gt.0.0d0.and.x(i,3).lt.4.53d0) THEN
               DO i1=1,15
                  widx0=-y0+wid_fa*dble(i1-1)
                  widx1=-y0+wid_fa*dble(i1  )
                  widy0= y0-wid_fa*dble(i1-1)
                  widy1= y0-wid_fa*dble(i1  )
                  IF(x(i,1).ge.widx0 .and. x(i,1).le.widx1)then
                     assem_ix(i)=i1
                  ENDIF
                  IF(x(i,2).le.widy0 .and. x(i,2).ge.widy1)then
                     assem_iy(i)=i1
                  ENDIF
               ENDDO
            ENDIF   
         ENDDO
!
!.........Input for k factor
!
         OPEN(1234, file='fsar_vfwl.in')
          DO iy=1,15
             READ(1234,*) (FSAR(ix,iy),ix=1,15)
          ENDDO
          DO ix=1,20
             READ(1234,*) vfwl_k(ix)
          ENDDO
          CLOSE(1234)
      ENDIF ! myrank
      IF(np.gt.1) THEN
         CALL broadcast_r(vfwl_k,20)
         CALL broadcast_r(FSAR,256)
      ENDIF
!
      IF(myrank.eq.0) THEN
!.........nzone
!          
          DO i=1,nn
             radius=DSQRT(x(i,1)**2.0d0+x(i,2)**2.0d0)
             !
             !Upper head            
             IF(x(i,3).gt.10.5d0) THEN
                nzone(i)=1
                porosity(i)=0.854d0
             !
             !Upper guide structure assembly
             ELSEIF(x(i,3).gt.6.35d0.and.x(i,3).lt.10.5d0.and.radius.lt.1.76d0) THEN
                nzone(i)=2
                porosity(i)=0.167d0
             !
             !Upper plenum connected hot leg                
             ELSEIF(x(i,3).gt.5.35.and.x(i,3).lt.6.35) THEN
                IF(radius.lt.1.35d0) THEN
                   nzone(i)=3 
                   porosity(i)=0.167d0
                ELSE
                   nzone(i)=4 
                   porosity(i)=0.95d0
                ENDIF
             !
             !Hot leg    
             ELSEIF(x(i,3).gt.5.35d0.and.x(i,3).lt.6.35d0.and.radius.gt.1.76d0.and.x(i,2).gt.-0.55d0.and.x(i,2).lt.0.55d0) THEN
                nzone(i)=4
                porosity(i)=0.95d0
             !
             !Upper plenum below hot leg
             ELSEIF(x(i,3).gt.4.53d0.and.x(i,3).lt.5.3d0.and.radius.lt.1.76d0) THEN
                nzone(i)=5
                porosity(i)=0.167d0
             ELSEIF(x(i,3).gt.0.3d0.and.x(i,3).lt.4.53d0.and.radius.lt.1.76d0) THEN  !1cell:0.3, old mesh:0.5
                !
                !Reflector(Core)                 
                IF(fsar(assem_ix(i),assem_iy(i)).eq.0.0d0.or.assem_ix(i).eq.16.or.assem_iy(i).eq.16) THEN
                   nzone(i)=7
                   porosity(i)=0.95d0
                ELSE            
                   !
                   !Core
                   nzone(i)=6
                   porosity(i)=0.57d0
                ENDIF

             ELSEIF(x(i,3).gt.0.0d0.and.x(i,3).lt.0.3d0.and.radius.lt.1.76d0) THEN   !1cell:0.3, old mesh:0.5
                !
                !Reflector(Core inlet)
                IF(fsar(assem_ix(i),assem_iy(i)).eq.0.0d0.or.assem_ix(i).eq.16.or.assem_iy(i).eq.16) THEN
                   nzone(i)=9
                   porosity(i)=0.95d0
                ELSE
                   !
                   !Core inlet                    
                   nzone(i)=8
                   porosity(i)=0.57d0
                ENDIF
             ELSEIF(x(i,3).lt.0.0d0) THEN
                IF(x(i,3).gt.-1.2) THEN
                   !
                   !Lower support plate                    
                   IF(radius.lt.1.35d0) THEN
                      nzone(i)=10
                      porosity(i)=0.889d0
                   ! 
                   !Flow skirt                      
                   ELSEIF(radius.ge.1.35d0) THEN
                      nzone(i)=11
                      porosity(i)=0.889d0
                   ENDIF  
                ELSE
                   !
                   !Lower head
                   nzone(i)=12
                   porosity(i)=0.889d0
                ENDIF
            !
            !Downcomer
            ELSEIF(x(i,3).gt.0.0d0.and.x(i,3).lt.5.3d0.and.radius.gt.1.76d0) THEN
               nzone(i)=13
               porosity(i)=0.95d0
            !
            !Downcomer connected to cold leg
            ELSEIF(x(i,3).gt.5.3d0.and.x(i,3).lt.6.35d0.and.radius.gt.1.76d0.and.radius.lt.2.05d0) THEN                  
               IF(x(i,2).gt.0.55d0.or.x(i,2).lt.-0.55d0) THEN  !hot leg region                
                  nzone(i)=14
               porosity(i)=0.95d0
               ENDIF      
            !
            !Cold leg
            ELSEIF(x(i,3).gt.5.45d0.and.x(i,3).lt.6.2d0.and.radius.gt.2.05d0) THEN
               nzone(i)=15
               porosity(i)=0.95d0
            !   
            !Vessel annulus above cold leg                  
            ELSEIF(x(i,3).gt.6.35d0.and.x(i,3).lt.10.5d0.and.radius.gt.1.76d0) THEN
               nzone(i)=16
               porosity(i)=0.95d0

            ENDIF
            IF(nzone(i).ne.6)porosity(i)=1.0d0
         ENDDO
!                              
         DEALLOCATE(assem_ix,assem_iy)
!
!        Fuel Assemly ONLY
         IF(vv_prob.eq.'OPR1000_fullcore_modmesh02_rv'       .or. &
            vv_prob.eq.'OPR1000_fullcore_modmesh02_rv_vessel')then
!
!............xn definition
!
            xn=0.0d0
            IF(ndim.eq.2) THEN
               DO i=1,nn
                  DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
                     sa=sv_tmp1(j,1)**2+sv_tmp1(j,2)**2
                     sa=1.d0/DSQRT(sa)
                     xn(j,1)=sv_tmp1(j,1)*sa
                     xn(j,2)=sv_tmp1(j,2)*sa
                  ENDDO
               ENDDO
            ELSE
               DO i=1,nn
                  DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
                     sa=sv_tmp1(j,1)**2+sv_tmp1(j,2)**2+sv_tmp1(j,3)**2
                     sa=1.d0/DSQRT(sa)
                     xn(j,1)=sv_tmp1(j,1)*sa
                     xn(j,2)=sv_tmp1(j,2)*sa
                     xn(j,3)=sv_tmp1(j,3)*sa
                  ENDDO
               ENDDO
            ENDIF         
!
!           rod-scale   
!
            sl=1.d0
            perm_tmp1=1.d0
            DO i=1,nn
               IF(subchannel_type_tmp(i).eq.1)THEN
                  porosity(i)=0.570863299d0
                  sl(i,1)=1/0.012852d0
                  sl(i,2)=1/0.012852d0
               ELSEIF(subchannel_type_tmp(i).eq.2)THEN
                  porosity(i)=0.654432d0
                  DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
                     k=j_neigh_tmp(j)
                     IF(k.ne.0)then
                        IF(subchannel_type_tmp(k).eq.2)THEN
                           IF(abs(xn(j,1).gt.0.5d0))THEN
                              sl(i,1)=1/(((0.00798d0+0.012852d0)/2+(0.00798d0+0.003352d0)/2)/2)
                              sl(i,2)=1/0.012852d0
                           ELSE
                              sl(i,1)=1/0.012852d0
                              sl(i,2)=1/(((0.00798d0+0.012852d0)/2+(0.00798d0+0.003352d0)/2)/2)
                           ENDIF
                        ENDIF
                     ENDIF   
                  ENDDO
               ELSEIF(subchannel_type_tmp(i).eq.3)THEN
                  porosity(i)=0.721726841d0
                  sl(i,1)=1/(((0.00798d0+0.012852d0)/2+(0.00798d0+0.003352d0)/2)/2)
                  sl(i,2)=1/(((0.00798d0+0.012852d0)/2+(0.00798d0+0.003352d0)/2)/2)
               ELSEIF(subchannel_type_tmp(i).eq.4)THEN
                  porosity(i)=0.05d0
                  sl(i,1)=1/0.012852d0
                  sl(i,2)=1/0.012852d0
               ELSEIF(subchannel_type_tmp(i).eq.5)THEN
                  porosity(i)=0.3196209766d0
                  sl(i,1)=1/0.012852d0
                  sl(i,2)=1/0.012852d0
               ELSEIF(subchannel_type_tmp(i).eq.6)THEN
                  porosity(i)=0.5855619953d0  !KSB 180123 modi
                  sl(i,1)=1/0.012852d0
                  sl(i,2)=1/0.012852d0
               ENDIF
               DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
                  IF(ABS(xn(j,3)).gt.0.5d0) perm_tmp1(j)=porosity(i)
               ENDDO
            ENDDO
         ENDIF   

         IF(vv_prob.eq.'OPR1000_fullcore_modmesh02_rv')then
            nzone=1
         ENDIF      

      ENDIF ! myrank

      RETURN
      ENDSUBROUTINE udfn_opr1000_zone_porosity      
!
!-----------------------------------------------------------------------------------------------------------------------------------------      
            
      SUBROUTINE udfn_apr1400_zone_porosity(x,nzone,porosity)
!
!.....This routine change the cell value of somaGrid.
!
      USE Zcore         , ONLY: np,myrank       
      USE Zparam        , ONLY: nn,ndim,pi
      USE Zmodel        , ONLY: vfwl_k,fsar
!      
      IMPLICIT NONE
!     
!     input
      REAL(8) x(nn,ndim)
!     output      
      INTEGER nzone(nn)
      REAL(8) porosity(nn)
      
!     local variables
      INTEGER i,ix,iy
      INTEGER,ALLOCATABLE:: assem_ix(:),assem_iy(:)
!     local arrray
      REAL(8) radius
      REAL(8) y0,wid_fa,widx0,widy0,widx1,widy1
      INTEGER i1
!
!.....OPR1000 Full vessel
!
      ALLOCATE(vfwl_k(20),fsar(18,18)) !See udfn_mom_wall
      FSAR(:,:)=0.1d0
      vfwl_k(:)=0.0d0      
      IF(myrank.eq.0) THEN
         ALLOCATE(assem_ix(nn),assem_iy(nn))
         assem_ix(:)=18
         assem_iy(:)=18
         wid_fa=0.208736
         y0=wid_fa*8.d0+wid_fa*0.5d0
         DO i=1,nn
            IF(x(i,3).gt.0.0d0.and.x(i,3).lt.4.53d0) THEN
               DO i1=1,17
                  widx0=-y0+wid_fa*dble(i1-1)
                  widx1=-y0+wid_fa*dble(i1  )
                  widy0= y0-wid_fa*dble(i1-1)
                  widy1= y0-wid_fa*dble(i1  )
                  IF(x(i,1).ge.widx0 .and. x(i,1).le.widx1)then
                     assem_ix(i)=i1
                  ENDIF
                  IF(x(i,2).le.widy0 .and. x(i,2).ge.widy1)then
                     assem_iy(i)=i1
                  ENDIF
               ENDDO
            ENDIF   
         ENDDO
!
!.........Input for k factor
!
         OPEN(1234, file='fsar_vfwl.in')
          DO iy=1,17
             READ(1234,*) (FSAR(ix,iy),ix=1,17)
          ENDDO
          DO ix=1,20
             READ(1234,*) vfwl_k(ix)
          ENDDO
          CLOSE(1234)
      ENDIF ! myrank
      IF(np.gt.1) THEN
         CALL broadcast_r(vfwl_k,20)
         CALL broadcast_r(FSAR,324)
      ENDIF
!
      IF(myrank.eq.0) THEN
!.........nzone
!          
          DO i=1,nn
             radius=DSQRT(x(i,1)**2.0d0+x(i,2)**2.0d0)
             !
             !Upper head            
             IF(x(i,3).gt.10.5d0) THEN
                nzone(i)=1
                porosity(i)=0.854d0
             !
             !Upper guide structure assembly
             ELSEIF(x(i,3).gt.6.35d0.and.x(i,3).lt.10.5d0.and.radius.lt.1.99d0) THEN
                nzone(i)=2
                porosity(i)=0.167d0
             !
             !Upper plenum connected hot leg                
             ELSEIF(x(i,3).gt.5.3.and.x(i,3).lt.6.35) THEN
                IF(radius.lt.1.57d0) THEN
                   nzone(i)=3 
                   porosity(i)=0.167d0
                ELSE
                   nzone(i)=4 
                   porosity(i)=0.95d0
                ENDIF
             !
             !Hot leg    
             ELSEIF(x(i,3).gt.5.3d0.and.x(i,3).lt.6.35d0.and.radius.gt.1.76d0) THEN
                nzone(i)=4
                porosity(i)=0.95d0
             !
             !Upper plenum below hot leg
             ELSEIF(x(i,3).gt.4.41d0.and.x(i,3).lt.5.3d0.and.radius.lt.1.99d0) THEN
                nzone(i)=5
                porosity(i)=0.167d0
             ELSEIF(x(i,3).gt.0.3d0.and.x(i,3).lt.4.41d0.and.radius.lt.1.99d0) THEN  !1cell:0.3, old mesh:0.5
                !
                !Reflector(Core)                 
                IF(fsar(assem_ix(i),assem_iy(i)).eq.0.0d0.or.assem_ix(i).eq.18.or.assem_iy(i).eq.18) THEN
                   nzone(i)=7
                   porosity(i)=0.95d0
                ELSE            
                   !
                   !Core
                   nzone(i)=6
                   porosity(i)=0.57d0
                ENDIF

             ELSEIF(x(i,3).gt.0.0d0.and.x(i,3).lt.0.3d0.and.radius.lt.1.99d0) THEN   !1cell:0.3, old mesh:0.5
                !
                !Reflector(Core inlet)
                IF(fsar(assem_ix(i),assem_iy(i)).eq.0.0d0.or.assem_ix(i).eq.18.or.assem_iy(i).eq.18) THEN
                   nzone(i)=9
                   porosity(i)=0.95d0
                ELSE
                   !
                   !Core inlet                    
                   nzone(i)=8
                   porosity(i)=0.57d0
                ENDIF
             ELSEIF(x(i,3).lt.0.0d0) THEN
                IF(x(i,3).gt.-1.2) THEN
                   !
                   !Lower support plate                    
                   IF(radius.lt.1.57d0) THEN
                      nzone(i)=10
                      porosity(i)=0.889d0
                   ! 
                   !Flow skirt                      
                   ELSEIF(radius.ge.1.57d0) THEN
                      nzone(i)=11
                      porosity(i)=0.889d0
                   ENDIF  
                ELSE
                   !
                   !Lower head
                   nzone(i)=12
                   porosity(i)=0.889d0
                ENDIF
            !
            !Downcomer
            ELSEIF(x(i,3).gt.0.0d0.and.x(i,3).lt.5.3d0.and.radius.gt.1.99d0) THEN
               nzone(i)=13
               porosity(i)=0.95d0
            !
            !Downcomer connected to cold leg
            !ELSEIF(x(i,3).gt.5.3d0.and.x(i,3).lt.6.35d0.and.radius.gt.1.76d0.and.radius.lt.2.05d0) THEN                  
            !   IF(x(i,2).gt.0.55d0.or.x(i,2).lt.-0.55d0) THEN  !hot leg region                
            !      nzone(i)=14
            !   porosity(i)=0.95d0
            !   ENDIF      
            !
            !Cold leg
            !ELSEIF(x(i,3).gt.5.45d0.and.x(i,3).lt.6.2d0.and.radius.gt.2.05d0) THEN
            !   nzone(i)=15
            !   porosity(i)=0.95d0
            !   
            !Vessel annulus above cold leg                  
            ELSEIF(x(i,3).gt.6.35d0.and.x(i,3).lt.10.5d0.and.radius.gt.1.99d0) THEN
               nzone(i)=16
               porosity(i)=0.95d0

            ENDIF
            IF(nzone(i).ne.6)porosity(i)=1.0d0
         ENDDO
!                              
         DEALLOCATE(assem_ix,assem_iy)
!
      ENDIF ! myrank

      RETURN
      ENDSUBROUTINE udfn_apr1400_zone_porosity      
!
!======================================================================
!======================================================================
!      
      SUBROUTINE udfn_opr1000_snglassm_porosity(nzone,porosity,sl)
      
      USE Zmpi          , ONLY: maxmt_cell
      USE Zcore         , ONLY: myrank       
      USE Zparam        , ONLY: nn,ndim,pi
      USE Znum_cell     , ONLY: i_neigh_tmp,j_neigh_tmp, &
                                perm_tmp1,sv_tmp1
      USE Zconst1       , ONLY: vv_prob
      !rod-scale
      USE Zrv_subchan   , ONLY: subchannel_type_tmp
!      
      IMPLICIT NONE
!     
!     input
!     output      
      INTEGER nzone(nn)
      REAL(8) porosity(nn)
      REAL(8) sl(nn,ndim)
      
!     local variables
      INTEGER i,j,k
      !rod-scale
      REAL(8) xn(maxmt_cell,ndim)     
      REAL(8) sa 
!      
      IF(myrank.eq.0) THEN
!
!........Zone 
         nzone=1         
!
!        Single assembly 
         IF(vv_prob.eq.'OPR1000_single_assem')then
!
!...........xn definition
!
            xn=0.0d0
            IF(ndim.eq.2) THEN
               DO i=1,nn
                  DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
                     sa=sv_tmp1(j,1)**2+sv_tmp1(j,2)**2
                     sa=1.d0/DSQRT(sa)
                     xn(j,1)=sv_tmp1(j,1)*sa
                     xn(j,2)=sv_tmp1(j,2)*sa
                  ENDDO
               ENDDO
            ELSE
               DO i=1,nn
                  DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
                     sa=sv_tmp1(j,1)**2+sv_tmp1(j,2)**2+sv_tmp1(j,3)**2
                     sa=1.d0/DSQRT(sa)
                     xn(j,1)=sv_tmp1(j,1)*sa
                     xn(j,2)=sv_tmp1(j,2)*sa
                     xn(j,3)=sv_tmp1(j,3)*sa
                  ENDDO
               ENDDO
            ENDIF         
!
!...........rod-scale (subchannel resolution)
!
            sl=1.d0
            perm_tmp1=1.d0
            DO i=1,nn
               IF(subchannel_type_tmp(i).eq.1)THEN
                  porosity(i)=0.570863299d0
                  sl(i,1)=1/0.012852d0
                  sl(i,2)=1/0.012852d0
               ELSEIF(subchannel_type_tmp(i).eq.2)THEN
                  porosity(i)=0.654432d0
                  DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
                     k=j_neigh_tmp(j)
                     IF(k.ne.0)then
                        IF(subchannel_type_tmp(k).eq.2)THEN
                           IF(abs(xn(j,1).gt.0.5d0))THEN
                              sl(i,1)=1/(((0.00798d0+0.012852d0)/2+(0.00798d0+0.003352d0)/2)/2)
                              sl(i,2)=1/0.012852d0
                           ELSE
                              sl(i,1)=1/0.012852d0
                              sl(i,2)=1/(((0.00798d0+0.012852d0)/2+(0.00798d0+0.003352d0)/2)/2)
                           ENDIF
                        ENDIF
                     ENDIF   
                  ENDDO
               ELSEIF(subchannel_type_tmp(i).eq.3)THEN
                  porosity(i)=0.721726841d0
                  sl(i,1)=1/(((0.00798d0+0.012852d0)/2+(0.00798d0+0.003352d0)/2)/2)
                  sl(i,2)=1/(((0.00798d0+0.012852d0)/2+(0.00798d0+0.003352d0)/2)/2)
               ELSEIF(subchannel_type_tmp(i).eq.4)THEN
                  porosity(i)=0.05d0
                  sl(i,1)=1/0.012852d0
                  sl(i,2)=1/0.012852d0
               ELSEIF(subchannel_type_tmp(i).eq.5)THEN
                  porosity(i)=0.3196209766d0
                  sl(i,1)=1/0.012852d0
                  sl(i,2)=1/0.012852d0
               ELSEIF(subchannel_type_tmp(i).eq.6)THEN
                  porosity(i)=0.5855619953d0  !KSB 180123 modi
                  sl(i,1)=1/0.012852d0
                  sl(i,2)=1/0.012852d0
               ENDIF
               DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
                  IF(ABS(xn(j,3)).gt.0.5d0) perm_tmp1(j)=porosity(i)
               ENDDO
            ENDDO
         ENDIF   

      ENDIF ! myrank
!
      END SUBROUTINE  udfn_opr1000_snglassm_porosity
              
