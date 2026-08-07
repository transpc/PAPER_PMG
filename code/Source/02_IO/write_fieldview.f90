!
      SUBROUTINE write_fieldview
!
!     SAVE field view DATA for Paraview: 3 vectors, 48 scalars = 51 variables
!
      USE VOL_DATA         , ONLY: cell 
      USE Solid_DATA       , ONLY: solid                 
      USE Wall_DATA        , ONLY: face      
      USE Zmpi             , ONLY: jperm,celem
      USE Zzone            , ONLY: ncell_fluid,ncell_fluid_all,ncell_cond,ncell_cond_all,icore,nzone
      USE Zrv_ncell        , ONLY: ncell_fluid_core,cupid_cell_channel,ncell_fluid_core_all,qvol_mas
      USE Zparam           , ONLY: nn,ndim
      USE Zcore            , ONLY: myrank
      USE Znum_cell        , ONLY: n_fluid,i_neigh_tmp,j_nbcon_tmp
      USE Zcoord3          , ONLY: porosity
      USE Zface            , ONLY: q1cell,qqcell,qecell,qclcell,qcgcell,ndensitycell   
      USE Zgrad_ls_c3d     , ONLY: lsindex
      USE Ziat             , ONLY: ia
      USE Zndforce         , ONLY: d_bfc,f_wl,dis_closewall,clift       
      USE Zncg             , ONLY: qn_cell 
      USE Zpress           , ONLY: p,dpdx_o
      USE Zqvol            , ONLY: gamma_wall,t_bulk,dry_weight,H_il,H_ig,H_gf,qporous_gas,qporous_liq,qvol_liq,gamma
      USE Ztimecon         , ONLY: time
      USE Zturb            , ONLY: turb_ke,turb_dp,turb_keg,turb_dpg,utau,yplus,yplusg,tauw
      USE Zmodel           , ONLY: qconden,molefr,qrad,resist
      USE Zvector          , ONLY: vg_n,vl_n,vd_n
      USE Zwall_HTC        , ONLY: mode,mul_o,twall_rv
      USE Zporous          , ONLY: chn_type
      USE Zrv_ncell        , ONLY: asm_nx,asm_ny,asm_nz,asm_ni,asm_ni2,chn_nx,chn_ny
      USE Zrv_ncell        , ONLY: dnbr_cupid1
      USE Zconst2          , ONLY: hydraulicd,gfactor
      USE viewData_common  , ONLY: nvector,nscalar,nframe,viwname,viwUnit,crit_zero, &
                                   cupid_rv_jperm, &
                                   alphal_all,alphag_all, &
                                   viewField
      USE Zio_unit        , ONLY: unit_log
      USE Zbc_index       , ONLY: ngrad
      USE KSMR            , ONLY: zone_comp
!
      IMPLICIT NONE
! 
!.....Local variables
      INTEGER :: i,ivar
      INTEGER :: ii,j,k,n_solid_porous
      INTEGER :: nncup,na
      REAL(8) :: xi_he,xi_steam,mol_he,mol_steam,ni_he,ni_steam
!.....Local arrays
      REAL(8) :: tmp1(ncell_fluid)
!.....Local allocatable arrays
      INTEGER,DIMENSION(:),ALLOCATABLE :: n_fluid_tmp,htmode_all  
      INTEGER,DIMENSION(:),ALLOCATABLE :: idat_rv,idat_rv_all,cupid_cell_channel_all,cupid_jperm
      INTEGER,DIMENSION(:),ALLOCATABLE :: itmp1_all
      REAL(8),DIMENSION(:),ALLOCATABLE :: tmp1_all,tmp2_all
      REAL(8),DIMENSION(:),ALLOCATABLE :: mul_o_all
      REAL(8),DIMENSION(:),ALLOCATABLE :: tsolid,tsolid_all
      REAL(8),DIMENSION(:),ALLOCATABLE :: dat_rv_all,dat_rv
!
      nncup=nn
      na=ncell_fluid_all
      !IF(cupid_mars)nncup=ncell_old(1) 
      nvector=0
      nscalar=0          
!        
!.....OPEN file and WRITE header        
!
      IF (nframe == 0) THEN 
         IF(myrank.eq.0) THEN
            viwname ='somaPlot.viw'
            OPEN (unit=viwUnit, file=trim(viwname),status='replace', form='unformatted')
            WRITE(viwUnit) nncup, ncell_fluid
            WRITE(viwUnit) viewField%nVectors
            WRITE(viwUnit) (viewField%vectorVar(i),i=1,viewField%nVectors)
            WRITE(viwUnit) viewField%nScalars
            WRITE(viwUnit) (viewField%scalarVar(i),i=1,viewField%nScalars)
            nframe = nframe + 1
         ENDIF   
      ENDIF
!      
      IF(myrank.eq.0) WRITE(viwUnit) nframe,time
!
!.....Allocate temp array
!
      IF(ncell_fluid_core_all.gt.0) THEN
         ALLOCATE(cupid_rv_jperm(ncell_fluid_core))
         DO i=1,ncell_fluid_core
            ii=cupid_cell_channel(i)
            cupid_rv_jperm(i)=jperm(ii)
         ENDDO   
      ELSE
         ALLOCATE(cupid_rv_jperm(1))
         cupid_rv_jperm(:)=0
      ENDIF
!
!.....Save alphag&alphal to check the zero fraction
!
      crit_zero=1.0d-4    ! criterion to check whether the phase fraction is zero or not
      IF(myrank.eq.0) THEN
         ALLOCATE(alphag_all(nncup))
         ALLOCATE(alphal_all(nncup))
      ELSE
         ALLOCATE(alphag_all(1))
         ALLOCATE(alphal_all(1))
      ENDIF
      CALL gatherv_r(cell%alphag,ncell_fluid,alphag_all,na,0)        
      CALL gatherv_r(cell%alphal,ncell_fluid,alphal_all,na,0)       
!
      IF(myrank.eq.0) THEN
         DO i=na+1,nncup
            alphag_all(i)=0.d0
            alphal_all(i)=0.d0
         ENDDO
      ENDIF
!
      DO ivar = 1, viewField%nVectors
!
!----------------------------------------------------------------------
!-----------------------Basic vector variables-------------------------
!----------------------------------------------------------------------               
!     
!........1. Gas-phase velocity
!
         IF(TRIM(ADJUSTL(viewField%vectorVar(ivar))) == "vg" ) THEN
            CALL wr_ndim_z(vg_n,alphag_all) 
!
!........2. Liquid-phase velocity
!
         ELSEIF(TRIM(ADJUSTL(viewField%vectorVar(ivar))) == "vl" ) THEN
            CALL wr_ndim_z(vl_n,alphal_all) 
!
!........3. Drop-phase velocity
!
         ELSEIF(TRIM(ADJUSTL(viewField%vectorVar(ivar))) == "vd" ) THEN
            CALL wr_ndim(vd_n) 
!
!........4. Gradient of pressure at cell center
!
         ELSEIF(TRIM(ADJUSTL(viewField%vectorVar(ivar))) == "dpdx" ) THEN
            CALL wr_ndim_z(dpdx_o,alphal_all)             
!!!!
!!!!........4. wall shear
!!!!
!!!         ELSEIF(TRIM(ADJUSTL(viewField%vectorVar(ivar))) == "wall_shear" ) THEN
!!!            CALL wr_ndim(wall_shear) 
!     
!........5. wall lubrication force
!
         ELSEIF(TRIM(ADJUSTL(viewField%vectorVar(ivar))) == "wall_lub" ) THEN
            CALL wr_ndim(f_wl) 
         ENDIF             
!
      END DO
!
      DO ivar = 1, viewField%nScalars
!      
!----------------------------------------------------------------------
!-----------------------Basic scalar variables-------------------------
!----------------------------------------------------------------------               
!         
!........4. Gas-phase volume fraction
!
         IF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "alphag" ) THEN
            nscalar=nscalar+1
            IF(myrank.eq.0) WRITE(viwUnit) (alphag_all(i),i=1,nncup)
!
!.......5. Liquid-phase volume fraction
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "alphal" ) THEN
            nscalar=nscalar+1
            IF(myrank.eq.0) WRITE(viwUnit) (alphal_all(i),i=1,nncup)   
!
!........6. Gas-phase volume fraction
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "alphad" ) THEN
            CALL wr_1d(cell%alphad) 
!
!........7. Gas-phase density
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "rhog" ) THEN
            CALL wr_1d(cell%rhog) 
!
!........8. Liquid-phase density
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "rhol" ) THEN
            CALL wr_1d(cell%rhol) 
!
!........9. Gas-phase density
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "rhoa" ) THEN
            CALL wr_1d(cell%rhoa) 
!
!........10. Pressure
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "p" ) THEN
            CALL wr_1d_scale(p,1.d-6)
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "ppa" ) THEN
            CALL wr_1d(p) 
!
!........11. Gas-phase temperature
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "tg" ) THEN
            CALL wr_1d_z(cell%tg,alphag_all) 
!
!........12. Liquid-phase temperature
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "tl" ) THEN
            CALL wr_1d_z(cell%tl,alphal_all) 
!
!........13. Gas-phase energy
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "eg" ) THEN
            CALL wr_1d(cell%eg) 
!
!........14. Liquid-phase energy
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "el" ) THEN
            CALL wr_1d(cell%el) 
!        
!........15. Saturation temperature
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "ts" ) THEN
            CALL wr_1d(cell%ts) 
!        
!........16. Boron concentration
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "cboron" ) THEN
            CALL wr_1d_scale(cell%cboron,1.d6)
!
!........17. Drop-phase energy
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "ed" ) THEN
            CALL wr_1d(cell%ed) 
!
!........18. Saturation gas enthalpy
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "hg" ) THEN
            CALL wr_1d(cell%hg) 
!
!........19. Saturation liquid enthalpy
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "hl" ) THEN
            CALL wr_1d(cell%hl) 
!
!----------------------------------------------------------------------
!-----------------------Interaface-Transfer-related variables----------
!----------------------------------------------------------------------
!        
!........20. Bubble diameter
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "dbubble" ) THEN
            CALL wr_1d(cell%D1) 
!        
!........21. Interfacial area concentration
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "iac" ) THEN
            CALL wr_1d(ia) 
!        
!........22. Subcooling 
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "subcooling" ) THEN
            DO i=1,ncell_fluid
               tmp1(i)=cell%ts(i)-cell%tl(i)
            ENDDO
            CALL wr_1d(tmp1) 
!        
!........23. Quality 
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "quality" ) THEN
            CALL wr_1d(cell%quala) 
!        
!........24. Interfacial drag 
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "vfgl" ) THEN
            CALL wr_1d(cell%vfgl) 
!        
!........25. Interfacial heat transfer coefficient, H_il 
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "hil" ) THEN
            CALL wr_1d(H_il) 
!        
!........26. Interfacial heat transfer coefficient, H_ig 
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "hig" ) THEN
            CALL wr_1d(H_ig) 
!        
!........27. Interfacial heat transfer coefficient, H_gf
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "hgf" ) THEN
            CALL wr_1d(H_gf) 
!
!----------------------------------------------------------------------
!-----------------------Turbulence-related variables-------------------
!----------------------------------------------------------------------
!
!........28. Turbulent kinetic energy
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "ke" ) THEN
            CALL wr_1d(turb_ke) 
!
!........28. Turbulent kinetic energy
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "keg" ) THEN
            CALL wr_1d(turb_keg) 
!
!........29. Turbulent dissipation rate
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "dp" ) THEN
            CALL wr_1d(turb_dp) 
!
!........29. Turbulent dissipation rate
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "dpg" ) THEN
            CALL wr_1d(turb_dpg) 
!
!........30. Gas-phase turbulent viscosity
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "mutg" ) THEN
            CALL wr_1d(cell%tviscosg) 
!
!........31. Liquid-phase turbulent viscosity
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "mutl" ) THEN
            CALL wr_1d(cell%tviscosl) 
!
!........32. Drop-phase turbulent viscosity
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "mutd" ) THEN
            CALL wr_1d(cell%tviscosd) 
!
!........33. Gas-phase kinematic turbulent viscosity
!
        ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "nutg" ) THEN
            CALL wr_1d_div(cell%tviscosg,cell%rhog) 
!
!........34. Liquid-phase  kinomatic turbulent viscosity
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "nutl" ) THEN
            CALL wr_1d_div(cell%tviscosl,cell%rhol) 
!
!........35. Gas-phase viscosity
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "mug" ) THEN
            CALL wr_1d(cell%lviscosg) 
!
!........36. Liquid-phase viscosity
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "mul" ) THEN
            CALL wr_1d(cell%lviscosl) 
!
!........37. Gas-phase kinematic viscosity
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "nug" ) THEN
            CALL wr_1d_div(cell%lviscosg,cell%rhog) 
!
!........38. Liquid-phase kinematic viscosity
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "nul" ) THEN
            CALL wr_1d_div(cell%lviscosl,cell%rhol) 
!        
!........39. Friction velocity
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "utau" ) THEN
            CALL wr_1d(utau) 
!        
!........40. Yplus
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "yplus" ) THEN
            CALL wr_1d(yplus) 
!        
!........40. Yplusg
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "yplusg" ) THEN
            CALL wr_1d(yplusg) 
!         
!----------------------------------------------------------------------
!-----------------------Heat-partition-related variables---------------
!----------------------------------------------------------------------
!        
!........41. Dry area weight
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "dry_weight" ) THEN
            CALL wr_1d(dry_weight) 
!        
!........42. Wall temperature
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "twall" ) THEN
            CALL wr_1d(face%twall_partition) 
!
!........43. Nucleation site density
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "ndensity" ) THEN
            CALL wr_1d(ndensitycell) 
!
!........44. Bubble departure diameter
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "ddepart" ) THEN
            CALL wr_1d(cell%Ddepart) 
!
!........45. Quenching heat flux
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "qq" ) THEN
            CALL wr_1d(qqcell) 
!
!........46. Evaporation heat flux
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "qe" ) THEN
            CALL wr_1d(qecell) 
!        
!........47. Liquid convectionm heat flux
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "qcl" ) THEN
            CALL wr_1d(qclcell) 
!        
!........48. Gas convection heat flux
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "qcg" ) THEN
            CALL wr_1d(qcgcell) 
!
!........49. Total heat flux
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "qin" ) THEN
            CALL wr_1d(q1cell) 
!        
!........50. Wall Evaporation rate (kg/s)
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "gamma_wall" ) THEN
            CALL wr_1d(gamma_wall) 
!        
!........51. Estimated bulk temperature
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "tbulk" ) THEN
            CALL wr_1d(t_bulk) 
!        
!........52. Flow regime
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "regime" ) THEN
            CALL wr_1d_i(cell%regime) 
!                          
!         
!----------------------------------------------------------------------
!-----------------------Addition---------------------------------------
!----------------------------------------------------------------------
!                   
!........53. Boundary conditions
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "bc" ) THEN
            nscalar=nscalar+1
            IF(myrank.eq.0) THEN
               ALLOCATE(itmp1_all(nncup))
               itmp1_all(:)=0
               DO i=1,na
                  DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
                     k=j_nbcon_tmp(j)
                     IF(k.ne.0) THEN
                        IF(ABS(k).ge.ABS(itmp1_all(i)).and.k.ne.-1) itmp1_all(i)=k
                     ENDIF
                  ENDDO
               ENDDO
               WRITE(viwUnit) (REAL(itmp1_all(i)),i=1,nncup)      
               DEALLOCATE(itmp1_all)
            ENDIF
!        
!........54. Distance from wall
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "walldis" ) THEN
            CALL wr_1d(d_bfc) 
!        
!........55. Distance from wall group1
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "walldis1" ) THEN
            CALL wr_1d(dis_closewall(1,1))
!        
!........56. Distance from wall wall group2
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "walldis2" ) THEN
            CALL wr_1d(dis_closewall(1,2)) 
!........57. Distance from wall wall group3
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "walldis3" ) THEN
            CALL wr_1d(dis_closewall(1,3)) 
! 
!........58. Wall drag of liquid phase
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "vfwl" ) THEN
            CALL wr_1d(cell%vfwl) 
!        
!........59. Color of subdomain
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "celem" ) THEN
            nscalar=nscalar+1
            IF(myrank.eq.0) WRITE(viwUnit) (real(celem(i)),i=1,nncup)
!
!........60. Ratio of viscosity
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "r_mul" ) THEN
            CALL wr_1d_div(cell%tviscosl,cell%lviscosl) 
!
!........61. Ratio of viscosity
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "r_mug" ) THEN
            CALL wr_1d_div(cell%tviscosg,cell%lviscosg) 
!
!........62. betal
!               
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "betal" ) THEN
            CALL wr_1d(cell%betal) 
!
!........63.betag           
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "betag" ) THEN
            CALL wr_1d(cell%betag) 
!
!........64.buoyancy coefficient
!         
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "ced33" ) THEN
            CALL wr_1d(cell%ced33) 
!
!........65.sigma
!               
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "sigma" ) THEN
            CALL wr_1d(cell%sigma) 
!
!........66.lcondl
!               
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "lcondl" ) THEN
            CALL wr_1d(cell%lcondl) 
!
!........67.lcondg
!               
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "lcondg" ) THEN
            CALL wr_1d(cell%lcondg) 
!
!........68.clift
!               
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "clift" ) THEN
            CALL wr_1d(clift) 
!
!........69.solidity
!               
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "solidity" ) THEN
            DO i=1,ncell_fluid
               tmp1(i)=1.d0-porosity(i)
            ENDDO
            CALL wr_1d(tmp1) 
!
!........70.tsolid
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "tsolid" ) THEN     
            nscalar=nscalar+1
            ALLOCATE(cupid_jperm(ncell_cond))
            IF(myrank.eq.0) THEN
               ALLOCATE(tsolid(ncell_cond_all),n_fluid_tmp(ncell_cond_all))
            ELSE
               ALLOCATE(tsolid(1),n_fluid_tmp(1))
            ENDIF
            DO i=1,ncell_cond
               IF(n_fluid(i).eq.0) THEN
                  cupid_jperm(i)=n_fluid(i)
               ELSE
                  cupid_jperm(i)=jperm(n_fluid(i))
               ENDIF
            ENDDO 
            CALL gatherv_i(cupid_jperm,ncell_cond,n_fluid_tmp,ncell_cond_all,1)
            CALL gatherv_r(solid%tsol ,ncell_cond,tsolid     ,ncell_cond_all,1)
!
            IF(myrank.eq.0) THEN
               ALLOCATE(tsolid_all(nncup))
               tsolid_all(:)=0.d0
                n_solid_porous=0
                DO i=1,ncell_cond_all
                   ii=n_fluid_tmp(i)
                   IF(ii.ge.1.and.ii.le.ncell_fluid_all)THEN
                       tsolid_all(ii)=tsolid(i)
                       n_solid_porous=i                   
                   ELSE 
                       tsolid_all(ncell_fluid_all+i-n_solid_porous)=tsolid(i)
                   ENDIF 
                ENDDO          
                WRITE(viwUnit) (tsolid_all(i),i=1,nncup)
               DEALLOCATE(tsolid_all)
            ENDIF    
            DEALLOCATE(tsolid,n_fluid_tmp,cupid_jperm)
!
!........71.tliquid&tsolid
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "tliqsol" ) THEN     
            nscalar=nscalar+1
            ALLOCATE(cupid_jperm(ncell_cond))
            IF(myrank.eq.0) THEN
               ALLOCATE(tsolid(ncell_cond_all),n_fluid_tmp(ncell_cond_all))
               ALLOCATE(tmp1_all(na))
            ELSE
               ALLOCATE(tsolid(1),n_fluid_tmp(1))
               ALLOCATE(tmp1_all(1))
            ENDIF
            DO i=1,ncell_cond
               IF(n_fluid(i).eq.0) THEN
                  cupid_jperm(i)=n_fluid(i)
               ELSE
                  cupid_jperm(i)=jperm(n_fluid(i))
               ENDIF
            ENDDO 
            CALL gatherv_i(cupid_jperm,ncell_cond,n_fluid_tmp,ncell_cond_all,1)
            CALL gatherv_r(solid%tsol ,ncell_cond,tsolid     ,ncell_cond_all,1)
!           
            CALL gatherv_r(cell%tl,ncell_fluid,tmp1_all,na,0)       
!
            IF(myrank.eq.0)THEN
               ALLOCATE(tsolid_all(nncup))
               tsolid_all(:)=0.d0
               n_solid_porous=0
               DO i=1,ncell_cond_all
                  ii=n_fluid_tmp(i)
                  IF(ii.ge.1.and.ii.le.ncell_fluid_all)THEN
!bug                     tsolid_all(ii)=tsolid(ii)    ! tsolid(i)
                      n_solid_porous=i                   
                  ELSE 
                      tsolid_all(ncell_fluid_all+i-n_solid_porous)=tsolid(i)
                  ENDIF 
               ENDDO        
               tsolid_all(1:na)=tmp1_all(1:na)
               WRITE(viwUnit) (tsolid_all(i),i=1,nncup)
               DEALLOCATE(tsolid_all)
            ENDIF    
            DEALLOCATE(tsolid,n_fluid_tmp,cupid_jperm)
            DEALLOCATE(tmp1_all)
!
!
!........69-1.tliquid +twall + tsolid
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "tall" ) THEN     
            nscalar=nscalar+1
            ALLOCATE(cupid_jperm(ncell_cond))
            IF(myrank.eq.0) THEN
               ALLOCATE(tsolid(ncell_cond_all),n_fluid_tmp(ncell_cond_all))
               ALLOCATE(tmp1_all(na),tmp2_all(na))
            ELSE
               ALLOCATE(tsolid(1),n_fluid_tmp(1))
               ALLOCATE(tmp1_all(1),tmp2_all(1))
            ENDIF
            DO i=1,ncell_cond
               IF(n_fluid(i).eq.0) THEN
                  cupid_jperm(i)=n_fluid(i)
               ELSE
                  cupid_jperm(i)=jperm(n_fluid(i))
               ENDIF
            ENDDO 
            CALL gatherv_i(cupid_jperm,ncell_cond,n_fluid_tmp,ncell_cond_all,1)
            CALL gatherv_r(solid%tsol ,ncell_cond,tsolid     ,ncell_cond_all,1)
!            
            CALL gatherv_r(cell%tl             ,ncell_fluid,tmp1_all,na,0)
            CALL gatherv_r(face%twall_partition,ncell_fluid,tmp2_all,na,0)
!
            IF(myrank.eq.0)THEN
               DO i=1,na
                  IF(tmp2_all(i).ne.0.and.ISNAN(tmp2_all(i)).eq..FALSE.) tmp1_all(i)=tmp2_all(i)
               ENDDO            
               ALLOCATE(tsolid_all(nncup))
               tsolid_all(:)=0.d0
               DO i=1,ncell_cond_all
                  ii=n_fluid_tmp(i)
                  IF(ii.ge.1.and.ii.le.ncell_fluid_all)THEN
!bug                    tsolid_all(ncell_fluid_all)=tsolid(ii) !tsolid(i
                  ELSE 
                      tsolid_all(ncell_fluid_all+i)=tsolid(i)
                  ENDIF 
               ENDDO        
               tsolid_all(1:na)=tmp1_all(1:na)
               WRITE(viwUnit) (tsolid_all(i),i=1,nncup)
               DEALLOCATE(tsolid_all)
            ENDIF    
            DEALLOCATE(tsolid,n_fluid_tmp)
            DEALLOCATE(tmp1_all,tmp2_all)
!
!........73. cell index
!               
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "cellind" ) THEN
            nscalar=nscalar+1
            IF(myrank.eq.0) WRITE(viwUnit) (i+0.0d0,i=1,nncup)     
!        
!........74. Shear stress
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "tauw" ) THEN
            CALL wr_1d(tauw) 
!        
!........75. Non-dimensional velocity 
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "uplus" ) THEN
            nscalar=nscalar+1
            IF(myrank.eq.0) THEN
               ALLOCATE(tmp1_all(na),tmp2_all(na))
            ELSE
               ALLOCATE(tmp1_all(1),tmp2_all(1))
            ENDIF
            CALL gatherv_r(vl_n(1,ndim),ncell_fluid,tmp1_all,na,0)
            CALL gatherv_r(utau        ,ncell_fluid,tmp2_all,na,0)
!
            IF(myrank.eq.0) THEN
!bug            DO i=1,nncup
!bug               IF(DABS(utau(i)).ge.1.0d-10)tmp1_all(i)=tmp1_all(i)/tmp2_all(i)
               DO i=1,ncell_fluid_all
                  IF(ABS(tmp2_all(i)).ge.1.0d-10) tmp1_all(i)=tmp1_all(i)/tmp2_all(i)
               ENDDO
               DO i=na+1,nn
                  tmp1_all(i)=0.d0
               ENDDO
               WRITE(viwUnit) (tmp1_all(i),i=1,nncup)
            ENDIF
            DEALLOCATE(tmp1_all,tmp2_all)
!
!........76. volume heat from solid to gas in porous media
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "qporous_g" ) THEN
            CALL wr_1d(qporous_gas) 
!
!........77. volume heat from solid to gas in porous media
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "qporous_l" ) THEN
            CALL wr_1d(qporous_liq) 
!
!........78.gamma multiplier
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "mul_o" ) THEN
            nscalar=nscalar+1
            ALLOCATE(dat_rv_all(ncell_fluid_core_all),cupid_cell_channel_all(ncell_fluid_core_all))
            IF(ncell_fluid_core.gt.0)THEN
               ALLOCATE(dat_rv(ncell_fluid_core),cupid_jperm(ncell_fluid_core))
            ELSE
               ALLOCATE(dat_rv(1),cupid_jperm(1))
            ENDIF   
            DO i=1,ncell_fluid_core
               dat_rv(i)=mul_o(i)
               cupid_jperm(i)=jperm(cupid_cell_channel(i))
            ENDDO   
            CALL gatherv_i(cupid_jperm,ncell_fluid_core,cupid_cell_channel_all,ncell_fluid_core_all,2)
            CALL gatherv_r(dat_rv     ,ncell_fluid_core,dat_rv_all            ,ncell_fluid_core_all,2)
!
            IF(myrank.eq.0) THEN
               ALLOCATE(mul_o_all(nncup))
               mul_o_all(:)=0.d0
               DO i=1,ncell_fluid_core_all  
                  ii=cupid_cell_channel_all(i)
                  mul_o_all(ii)=dat_rv_all(i)
               ENDDO                      
               WRITE(viwUnit) (mul_o_all(i),i=1,nncup)     
               DEALLOCATE(mul_o_all)
            ENDIF
               DEALLOCATE(dat_rv,dat_rv_all)
               DEALLOCATE(cupid_jperm,cupid_cell_channel_all)
!
!........79.Wall heat transfer mode
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "htmode" ) THEN
            nscalar=nscalar+1  
            ALLOCATE(cupid_cell_channel_all(ncell_fluid_core_all),idat_rv_all(ncell_fluid_core_all))
            IF(ncell_fluid_core.gt.0)THEN
               ALLOCATE(cupid_jperm(ncell_fluid_core),idat_rv(ncell_fluid_core))
            ELSE
               ALLOCATE(cupid_jperm(1),idat_rv(1))
            ENDIF   
            DO i=1,ncell_fluid_core
               idat_rv(i)=mode(i)
               cupid_jperm(i)=jperm(cupid_cell_channel(i))
            ENDDO   
            CALL gatherv_i(cupid_jperm,ncell_fluid_core,cupid_cell_channel_all,ncell_fluid_core_all,2)
            CALL gatherv_i(idat_rv    ,ncell_fluid_core,idat_rv_all           ,ncell_fluid_core_all,2)
            IF(myrank.eq.0) THEN
               ALLOCATE(htmode_all(nncup))
               htmode_all(:)=0
               DO i=1,ncell_fluid_core_all  
                  ii=cupid_cell_channel_all(i)
                  htmode_all(ii)=idat_rv_all(i)
               ENDDO     
               WRITE(viwUnit) (REAL(htmode_all(i)),i=1,nncup)
               DEALLOCATE(htmode_all)
            ENDIF 
            DEALLOCATE(idat_rv,idat_rv_all)
            DEALLOCATE(cupid_jperm,cupid_cell_channel_all)      
!
!........79-2.Wall heat transfer mode
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "lsindex" ) THEN
            CALL wr_1d_i(lsindex)               
! 
!........80. Wall drag of gas phase
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "vfwg" ) THEN
            CALL wr_1d(cell%vfwg) 
!        
!........81. Wall temperature in RV_model1
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "twall_rv1" ) THEN
            CALL wr_1d_rv(twall_rv(1,1)) 
!        
!........82. Wall temperature in RV_model2
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "twall_rv2" ) THEN
            CALL wr_1d_rv(twall_rv(1,2)) 
!        
!........83. Wall temperature in RV_model3
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "twall_rv3" ) THEN
            CALL wr_1d_rv(twall_rv(1,3)) 
!        
!........84. Wall temperature in RV_model4
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "twall_rv4" ) THEN
            CALL wr_1d_rv(twall_rv(1,4)) 
!        
!........85. Wall temperature in RV_model5
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "twall_rv5" ) THEN
            CALL wr_1d_rv(twall_rv(1,5)) 
!        
!........86. Non-condensable gas fractions 1
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "qn_1" ) THEN
            CALL wr_1d(qn_cell(1,1)) 
!        
!........87. Non-condensable gas fractions 2
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "qn_2" ) THEN
            CALL wr_1d(qn_cell(1,2)) 
!        
!........88.icore
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "icore" ) THEN
            CALL wr_1d_i(icore) 
!        
!........89.nzone
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "zone" ) THEN
            CALL wr_1d_i(nzone) 
!        
!........90.zone_comp (KSMR PSH edit)
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "zone_comp" ) THEN
            IF(nframe>1) THEN
                CALL wr_1d_i(zone_comp)         !(PSH)
            ELSE
                CALL wr_1d_i(nzone)
            ENDIF 
!        
!........##. wall condensation heat flux
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "qconden" ) THEN
            CALL wr_1d(qconden) 
!        
!........##. wall radiation heat flux
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "qrad" ) THEN
            CALL wr_1d(qrad)                 
!        
!........##. cpg
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "cpg" ) THEN
            CALL wr_1d(cell%cpg) 
!        
!........##. pps
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "pps" ) THEN
            CALL wr_1d(cell%pps) 
!        
!........##. qvol_mas
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "qvol_mas" ) THEN
            CALL wr_1d(qvol_mas) 
!        
!........##. qvol_liq
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "qvol_liq" ) THEN
            CALL wr_1d(qvol_liq) 
!        
!........##. gfactor
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "gfactor" ) THEN
            CALL wr_1d(gfactor) 
!        
!........##. hd
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "hd" ) THEN
            CALL wr_1d(hydraulicd) 
!
!........##. molar fraction
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "molefr" ) THEN
            mol_he=4.003d0
            mol_steam=18.02d0
            DO i=1,ncell_fluid
               xi_he=cell%quala(i)
               xi_steam=1.0d0-cell%quala(i)
               ni_he=xi_he/mol_he
               ni_steam=xi_steam/mol_steam
               molefr(i)=ni_he/(ni_he+ni_steam)
            ENDDO
            CALL wr_1d(molefr) 
!
!........##. mass diffusivity
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "mdiff" ) THEN
            CALL wr_1d(cell%mdiff) 
!        
!........##. subchannel type
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "sc_type" ) THEN
            CALL wr_1d_i(chn_type) 
!        
!........##. nx_asm
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "asm_nx" ) THEN
            nscalar=nscalar+1
            IF(myrank.eq.0) WRITE(viwUnit) (real(asm_nx(i)),i=1,nncup)
!        
!........##. asm_ny
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "asm_ny" ) THEN
            nscalar=nscalar+1
            IF(myrank.eq.0) WRITE(viwUnit) (real(asm_ny(i)),i=1,nncup)
!        
!........##. asm_nz
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "asm_nz" ) THEN
            nscalar=nscalar+1
            IF(myrank.eq.0) WRITE(viwUnit) (real(asm_nz(i)),i=1,nncup)
!        
!........##. asm_ni
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "asm_ni" ) THEN
            nscalar=nscalar+1
            IF(myrank.eq.0) WRITE(viwUnit) (real(asm_ni(i)),i=1,nncup)
!        
!........##. asm_ni2
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "asm_ni2" ) THEN
            nscalar=nscalar+1
            IF(myrank.eq.0) WRITE(viwUnit) (real(asm_ni2(i)),i=1,nncup)
!        
!........##. chn_nx
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "chn_nx" ) THEN
            nscalar=nscalar+1
            IF(myrank.eq.0) WRITE(viwUnit) (real(chn_nx(i)),i=1,nncup)
!        
!........##. chn_ny
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "chn_ny" ) THEN
            nscalar=nscalar+1
            IF(myrank.eq.0) WRITE(viwUnit) (real(chn_ny(i)),i=1,nncup)
!        
!........##. DNBR
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "dnbr_cupid1" ) THEN
            CALL wr_1d(dnbr_cupid1) 
!        
!........##. gamma
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "gamma" ) THEN
            CALL wr_1d(gamma)
!        
!........lsindex
!            
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "lsindex" ) THEN
            CALL wr_1d_i(lsindex)             
!        
!........ngrad
!            
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "ngrad" ) THEN
            CALL wr_1d_i(ngrad)  
!            
!........exit quality
!               
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "quals" ) THEN
            CALL wr_1d(cell%quals)     
            
!........resist
!
         ELSEIF(TRIM(ADJUSTL(viewField%scalarVar(ivar))) == "resist" ) THEN
            CALL wr_1d(resist)   
         ENDIF
!
      ENDDO
!      
!
!.....Deallocate temp array
!
      DEALLOCATE(cupid_rv_jperm)
!
      DEALLOCATE(alphag_all)
      DEALLOCATE(alphal_all)
!
      IF(nvector.ne.viewField%nVectors.or.nscalar.ne.viewField%nScalars)THEN
         IF(myrank.eq.0)THEN
            WRITE(*,*)'          Check variable names using Variable Info and write_fieldview.f90 !!!'
            WRITE(*,*)nvector,viewField%nVectors,nscalar,viewField%nScalars
            WRITE(unit_log,*)'          Check variable names using Variable Info and write_fieldview.f90 !!!'
            WRITE(unit_log,*)nvector,viewField%nVectors,nscalar,viewField%nScalars
         ENDIF
         CALL finalize_mpi
         STOP
      ENDIF
!   

      END SUBROUTINE write_fieldview
!
!DIR$ ATTRIBUTES NOINLINE :: wr_ndim_z
      SUBROUTINE wr_ndim_z(x,alpha_all)
!
      USE Zmpi             , ONLY: ncell_fp
      USE Zzone            , ONLY: ncell_fluid,ncell_fluid_all
      USE Zcore            , ONLY: myrank
      USE Zparam           , ONLY: nn,ndim
      USE viewData_common  , ONLY: nvector,crit_zero,viwUnit
!
      IMPLICIT NONE
! 
!.....Input
      REAL(8) x(ncell_fp,ndim)
      REAL(8) alpha_all(ncell_fluid_all)
!.....Local variables
      INTEGER i,ix
      INTEGER nncup,na
!.....Local arrays
      REAL(8),DIMENSION(:),ALLOCATABLE :: tmp
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: tmp0
!
      nncup=nn
      na=ncell_fluid_all
      !IF(cupid_mars)nncup=ncell_old(1) 
!
      nvector=nvector+1
      IF(myrank.eq.0) THEN
         ALLOCATE(tmp0(na,ndim))
      ELSE
         ALLOCATE(tmp0(1,1))
      ENDIF
      CALL gatherv_r_2d(x,ncell_fp,tmp0,ncell_fluid,na,0)
!
      IF(myrank.eq.0) THEN                   
         ALLOCATE(tmp(nncup))
         DO i=na+1,nncup
            tmp(i)=0.d0
         ENDDO
!........Make zero if the phase fraction is almost zero
!            
         DO ix=1,ndim
            DO i=1, na
               IF(alpha_all(i).le.crit_zero) THEN
                  tmp(i)=0.0d0
               ELSE
                  tmp(i)=tmp0(i,ix)
               ENDIF
            ENDDO
            WRITE(viwUnit) (tmp(i),i=1,nncup)
         ENDDO
         DEALLOCATE(tmp)
      ENDIF   
      DEALLOCATE(tmp0)
!
      END SUBROUTINE wr_ndim_z
!
!DIR$ ATTRIBUTES NOINLINE :: wr_ndim
      SUBROUTINE wr_ndim(x)
!
      USE Zmpi             , ONLY: ncell_fp
      USE Zzone            , ONLY: ncell_fluid,ncell_fluid_all
      USE Zcore            , ONLY: myrank
      USE Zparam           , ONLY: nn,ndim
      USE viewData_common  , ONLY: nvector,viwUnit
!
      IMPLICIT NONE
! 
!.....Input
      REAL(8) x(ncell_fp,ndim)
!.....Local variables
      INTEGER i,ix
      INTEGER nncup,na
!.....Local arrays
      REAL(8),DIMENSION(:),ALLOCATABLE :: tmp
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: tmp0
!
      nncup=nn
      na=ncell_fluid_all
      !IF(cupid_mars)nncup=ncell_old(1) 
!
      nvector=nvector+1
      IF(myrank.eq.0) THEN
         ALLOCATE(tmp0(na,ndim))
      ELSE
         ALLOCATE(tmp0(1,1))
      ENDIF
      CALL gatherv_r_2d(x,ncell_fp,tmp0,ncell_fluid,na,0)
!
      IF(myrank.eq.0) THEN                   
         ALLOCATE(tmp(nncup))
         DO i=na+1,nncup
            tmp(i)=0.d0
         ENDDO
         DO ix=1,ndim
            DO i=1, na
               tmp(i)=tmp0(i,ix)
            ENDDO
            WRITE(viwUnit) (tmp(i),i=1,nncup)
         ENDDO
         DEALLOCATE(tmp)
      ENDIF   
!
      END SUBROUTINE wr_ndim
!
!DIR$ ATTRIBUTES NOINLINE :: wr_1d_z
      SUBROUTINE wr_1d_z(x,alpha_all)
!
      USE Zzone            , ONLY: ncell_fluid,ncell_fluid_all
      USE Zcore            , ONLY: myrank
      USE Zparam           , ONLY: nn
      USE viewData_common  , ONLY: nscalar,viwUnit,crit_zero
!
      IMPLICIT NONE
! 
!.....Input
      REAL(8) x(ncell_fluid)
      REAL(8) alpha_all(ncell_fluid_all)
!.....Local variables
      INTEGER i
      INTEGER nncup,na
!.....Local arrays
      REAL(8),DIMENSION(:),ALLOCATABLE :: tmp
!
      nncup=nn
      na=ncell_fluid_all
      !IF(cupid_mars)nncup=ncell_old(1) 
!
      nscalar=nscalar+1
      IF(myrank.eq.0) THEN
         ALLOCATE(tmp(nncup))
      ELSE
         ALLOCATE(tmp(1))
      ENDIF
      CALL gatherv_r(x,ncell_fluid,tmp,na,0)
!
      IF(myrank.eq.0) THEN
!
!........Make zero if the phase fraction is almost zero
!            
         DO i=1, na
            IF(alpha_all(i).le.crit_zero) tmp(i)=0.0d0
         ENDDO
         DO i=na+1,nncup
            tmp(i)=0.d0
         ENDDO
         WRITE(viwUnit) (tmp(i),i=1,nncup)
      ENDIF                              
      DEALLOCATE(tmp)
!
      END SUBROUTINE wr_1d_z
!
!DIR$ ATTRIBUTES NOINLINE :: wr_1d
      SUBROUTINE wr_1d(x)
!
      USE Zzone        , ONLY: ncell_fluid,ncell_fluid_all
      USE Zcore        , ONLY: myrank
      USE Zparam       , ONLY: nn
      USE viewData_common  , ONLY: nscalar,viwUnit
!
      IMPLICIT NONE
! 
!.....Input
      REAL(8) x(ncell_fluid)
!.....Local variables
      INTEGER i
      INTEGER nncup,na
!.....Local arrays
      REAL(8),DIMENSION(:),ALLOCATABLE :: tmp
!
      nncup=nn
      na=ncell_fluid_all
      !IF(cupid_mars)nncup=ncell_old(1) 
!
      nscalar=nscalar+1
      IF(myrank.eq.0) THEN
         ALLOCATE(tmp(nncup))
      ELSE
         ALLOCATE(tmp(1))
      ENDIF
      CALL gatherv_r(x,ncell_fluid,tmp,na,0)
!
      IF(myrank.eq.0) THEN
         DO i=na+1,nncup
            tmp(i)=0.d0
         ENDDO
         WRITE(viwUnit) (tmp(i),i=1,nncup)
      ENDIF
      DEALLOCATE(tmp)
!
      END SUBROUTINE wr_1d
!
!DIR$ ATTRIBUTES NOINLINE :: wr_1d_scale
      SUBROUTINE wr_1d_scale(x,scale)
!
      USE Zzone        , ONLY: ncell_fluid,ncell_fluid_all
      USE Zcore        , ONLY: myrank
      USE Zparam       , ONLY: nn
      USE viewData_common  , ONLY: nscalar,viwUnit
!
      IMPLICIT NONE
!
!.....Input
      REAL(8) scale
      REAL(8) x(ncell_fluid)
!.....Local variables
      INTEGER i
      INTEGER nncup,na
!.....Local arrays
      REAL(8),DIMENSION(:),ALLOCATABLE :: tmp
!
      nncup=nn
      na=ncell_fluid_all
      !IF(cupid_mars)nncup=ncell_old(1)
!
      nscalar=nscalar+1
      IF(myrank.eq.0) THEN
         ALLOCATE(tmp(nncup))
      ELSE
         ALLOCATE(tmp(1))
      ENDIF
      CALL gatherv_r(x,ncell_fluid,tmp,na,0)
!
      IF(myrank.eq.0) THEN
         DO i=na+1,nncup
            tmp(i)=0.d0
         ENDDO
         WRITE(viwUnit) (tmp(i)*scale,i=1,nncup)
      ENDIF
      DEALLOCATE(tmp)
!
      END SUBROUTINE wr_1d_scale
!
!DIR$ ATTRIBUTES NOINLINE :: wr_1d_div
      SUBROUTINE wr_1d_div(x1,x2)
!
      USE Zzone            , ONLY: ncell_fluid,ncell_fluid_all
      USE Zcore            , ONLY: myrank
      USE Zparam           , ONLY: nn
      USE viewData_common  , ONLY: nscalar,viwUnit
!
      IMPLICIT NONE
! 
!.....Input
      REAL(8) x1(ncell_fluid),x2(ncell_fluid)
!.....Local variables
      INTEGER i
      INTEGER nncup,na
!.....Local arrays
      REAL(8),DIMENSION(:),ALLOCATABLE :: tmp1,tmp2
!
      nncup=nn
      na=ncell_fluid_all
      !IF(cupid_mars)nncup=ncell_old(1) 
!
      nscalar=nscalar+1
      IF(myrank.eq.0) THEN
         ALLOCATE(tmp1(nncup),tmp2(na))
      ELSE
         ALLOCATE(tmp1(1),tmp2(1))
      ENDIF
      CALL gatherv_r(x1,ncell_fluid,tmp1,na,0)
      CALL gatherv_r(x2,ncell_fluid,tmp2,na,0)
!
      IF(myrank.eq.0) THEN
         DO i=1,na
            IF(tmp2(i).ne.0.d0) tmp1(i)=tmp1(i)/tmp2(i)
         ENDDO
         DO i=na+1,nncup
            tmp1(i)=0.d0
         ENDDO
         WRITE(viwUnit) (tmp1(i),i=1,nncup)
      ENDIF
      DEALLOCATE(tmp1,tmp2)
!
      END SUBROUTINE wr_1d_div
!
!DIR$ ATTRIBUTES NOINLINE :: wr_1d_i
      SUBROUTINE wr_1d_i(x)
!
      USE Zzone            , ONLY: ncell_fluid,ncell_fluid_all
      USE Zcore            , ONLY: myrank
      USE Zparam           , ONLY: nn
      USE viewData_common  , ONLY: nscalar,viwUnit
!
      IMPLICIT NONE
! 
!.....Input
      INTEGER x(ncell_fluid)
!.....Local variables
      INTEGER i
      INTEGER nncup,na
!.....Local arrays
      INTEGER,DIMENSION(:),ALLOCATABLE :: itmp
!
      nncup=nn
      na=ncell_fluid_all
      !IF(cupid_mars)nncup=ncell_old(1) 
!
      nscalar=nscalar+1
      IF(myrank.eq.0) THEN
         ALLOCATE(itmp(nncup))
      ELSE
         ALLOCATE(itmp(1))
      ENDIF
      CALL gatherv_i(x,ncell_fluid,itmp,na,0)
!
      IF(myrank.eq.0) THEN
         DO i=na+1,nncup
            itmp(i)=0
         ENDDO
         WRITE(viwUnit) (REAL(itmp(i)),i=1,nncup)
      ENDIF
      DEALLOCATE(itmp)
!
      END SUBROUTINE wr_1d_i
!
      SUBROUTINE wr_1d_rv(x) 
!
      USE Zrv_ncell        , ONLY: ncell_fluid_core
      USE Zrv_ncell        , ONLY: ncell_fluid_core,ncell_fluid_core_all
      USE Zcore            , ONLY: myrank
      USE Zparam           , ONLY: nn
      USE viewData_common  , ONLY: nscalar,viwUnit,cupid_rv_jperm
!
      IMPLICIT NONE
!
!.....Input
      REAL(8) x(ncell_fluid_core)
!.....Local variables
      INTEGER i,ii
      INTEGER nncup,na
!.....Local arrays
      INTEGER,DIMENSION(:),ALLOCATABLE :: itmp
      REAL(8),DIMENSION(:),ALLOCATABLE :: tmp,tmp1
!
      nncup=nn
      na=ncell_fluid_core_all
      !IF(cupid_mars)nncup=ncell_old(1) 
!
      nscalar=nscalar+1
      IF(myrank.eq.0) THEN
         ALLOCATE(itmp(na),tmp(na))
      ELSE
         ALLOCATE(itmp(1),tmp(1))
      ENDIF
      CALL gatherv_i(cupid_rv_jperm,ncell_fluid_core,itmp,na,2)
      CALL gatherv_r(x             ,ncell_fluid_core,tmp ,na,2)
!
      IF(myrank.eq.0) THEN
         ALLOCATE(tmp1(nncup))
         tmp1(:)=0.d0
         DO i=1,na
            ii=itmp(i)
            tmp1(ii)=tmp(i)
         ENDDO                
         WRITE(viwUnit) (tmp1(i),i=1,nncup)     
         DEALLOCATE(tmp1)
      ENDIF
      DEALLOCATE(itmp,tmp)
!
      END SUBROUTINE wr_1d_rv
