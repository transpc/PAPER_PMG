!
      SUBROUTINE turbulence_mod
!
!     This routine select the turbulence model and wall model.
!
      USE VOL_DATA     , ONLY: cell            
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: ndim
      USE Ztimecon     , ONLY: itim      
      USE Znum_cell    , ONLY: i_neigh
      USE Zbc_index    , ONLY: nbcon,iface_wall
      USE Zb_condition , ONLY: vb_liq,vin_liq
      USE Zconst1      , ONLY: iturb,restart,turb_phase
      USE Zface        , ONLY: Twall_Model,                         &
                               Kepsilon,Kepsilon_real,Kepsilon_RNG, &
                               Laminar,Zequation,SST,Free_slip,     &
                               liq_only,gas_only,both,LES_WALE
      USE Zndforce     , ONLY: dvdxl
      USE Zturb        , ONLY: turb_ke,turb_dp,utau,yplus,walln,wcd_gas,wcd_liq,wvis_gas,wvis_liq, &
                               ustar_ke,ustar_keg,w_real_ke,w_real_keg,cmu_Real,cmug_Real
      USE Zturbzeq     , ONLY: tlengs
      USE Zvector      , ONLY: vg_n,vl_n
      USE Znormal      , ONLY: wall_cell_l,num_wallcells_l,xn_wallcell_l,sa_wallcell_l
      USE Zio_unit     , ONLY: unit_log
!
      IMPLICIT NONE
!
!.....Local variable
      INTEGER i,j,n,ii,j0
      INTEGER Turbulence_Model      
      LOGICAL, SAVE :: INITIAL=.TRUE.
      REAL(8) wvisl,wcdl,wvisg,wcdg
      REAL(8) tPrg,tPrl
!.....Local arrays
      REAL(8) vgn_i(3),vln_i(3)
      REAL(8) utaul,yplusl,cv(3)
!
      n=ncell_fp
      IF(iturb.eq.Kepsilon_real)THEN
         IF(INITIAL.eq..TRUE.)THEN
            ALLOCATE(ustar_ke(n),ustar_keg(n),w_real_ke(n),w_real_keg(n),cmu_Real(n),cmug_Real(n))      
            DO i=1,n 
               ustar_ke(i)  =0.0d0
               ustar_keg(i) =0.0d0
               w_real_ke(i) =0.0d0
               w_real_keg(i)=0.0d0
               cmu_Real(i)  =0.09d0
               cmug_Real(i) =0.09d0
            ENDDO
            INITIAL=.FALSE.
         ENDIF 
      ENDIF 
!
      IF(itim.eq.1 .and. restart.eq.0) THEN
         cell%eviscosl(:)=cell%lviscosl(:)
         cell%eviscosg(:)=cell%lviscosg(:)
         cell%eviscosd(:)=cell%lviscosl(:)
         cell%tviscosl(:)=0.0d0
         cell%tviscosg(:)=0.0d0
         cell%tviscosd(:)=0.0d0
         cell%condl(:)   =cell%lcondl(:)         
         cell%condg(:)   =cell%lcondg(:)
      ENDIF        
!
!.....Select turbulence model and wall model
!
      Turbulence_Model=iturb
      Twall_Model=Turbulence_Model
!   
!.....1) Laminar
!      
      IF(Turbulence_Model.eq.Laminar) THEN
!
         DO i=1,ncell_fluid
            cell%eviscosg(i)=cell%lviscosg(i)
            cell%eviscosl(i)=cell%lviscosl(i)
            cell%eviscosd(i)=cell%lviscosl(i)
            cell%tviscosg(i)=0.d0
            cell%tviscosl(i)=0.d0
            cell%tviscosd(i)=0.d0
            cell%condg(i)   =cell%lcondg(i)
            cell%condl(i)   =cell%lcondl(i)
         ENDDO
!   
!.....2) Zero equation model
!         
      ELSEIF(Turbulence_Model.eq.Zequation) THEN
!
!........Wall Function for Zero Equation Model
!
         DO ii=1, num_wallcells_l
            i=wall_cell_l(ii)
            j0=i_neigh(i)-1 
            j=iface_wall(i)
            vgn_i(1)=vg_n(i,1)
            vgn_i(2)=vg_n(i,2)
            vln_i(1)=vl_n(i,1)
            vln_i(2)=vl_n(i,2)
            cv(1)   =xn_wallcell_l(ii,1)
            cv(2)   =xn_wallcell_l(ii,2)
            IF(ndim.eq.3) THEN
               vgn_i(3)=vg_n(i,3)
               vln_i(3)=vl_n(i,3)
               cv(3)   =xn_wallcell_l(ii,3)
            ENDIF
            utaul=utau(i)
            CALL Wall_Zeq_i(i,nbcon(j+j0),vgn_i,vln_i,                &
                            cell%lviscosl(i),cell%tviscosl(i),        &
                            cell%lcondl(i),cell%cpl(i),cell%rhol(i),  &
                            cell%lviscosg(i),cell%tviscosg(i),        &
                            cell%lcondg(i),cell%cpg(i),cell%rhog(i),  &
                            cv,sa_wallcell_l(ii),walln(i),            &
                            wvisl,wcdl,wvisg,wcdg,                    &
                            utaul,yplusl)
!               
             wvis_liq(i)=wvisl
             wcd_liq(i)=wcdl
             wvis_gas(i)=wvisg
             wcd_gas(i)=wcdg
             utau(i)=utaul
             yplus(i)=yplusl
!               
             IF(.true.) THEN
                wvis_gas(i)=cell%lviscosg(i)/cell%lviscosl(i)*wvisl
                wcd_gas(i)=cell%lcondg(i)/cell%lcondl(i)*wcdl
             ENDIF
         ENDDO
!
!........Viscosity and kinetic energy and dissipation for zero equation model
!
         CALL grad_vel(2,vl_n,dvdxl,vb_liq,vin_liq)
         CALL Turb_Zeq
!
         DO i=1,ncell_fluid
            cell%eviscosg(i)=cell%lviscosg(i)+cell%tviscosg(i)
            cell%eviscosl(i)=cell%lviscosl(i)+cell%tviscosl(i)
! bug cell%lviscosd never computed replaced by cell%lviscosl
!           cell%eviscosd(i)=cell%lviscosd(i)
            cell%eviscosd(i)=cell%lviscosl(i)
            tPrg=0.9d0
            tPrl=0.9d0
            cell%condg(i)=cell%lcondg(i)+cell%tviscosg(i)*cell%cpg(i)/tPrg
            cell%condl(i)=cell%lcondl(i)+cell%tviscosl(i)*cell%cpl(i)/tPrl
         ENDDO
!
         IF(ndim.eq.2) THEN         
            DO i=1,ncell_fluid
               turb_ke(i)=tlengs*DSQRT(dvdxl(i,1,ndim)**2)
               turb_dp(i)=2.0d0*cell%lviscosl(i)/cell%rhol(i)*(dvdxl(i,1,ndim)**2)
            ENDDO
         ENDIF
!         
         IF(ndim.eq.3) THEN
            DO i=1,ncell_fluid
               turb_ke(i)=tlengs*DSQRT(dvdxl(i,1,ndim)**2+dvdxl(i,2,ndim)**2)
               turb_dp(i)=2.0d0*cell%lviscosl(i)/cell%rhol(i)*(dvdxl(i,1,ndim)**2+dvdxl(i,2,ndim)**2)
            ENDDO
         ENDIF
!    
!.....3)k-epsilon turbulence model and wall funciton
!     * Flag
!       Kepsilon_l = Flag for liquid ONLY
!       Kepsilon_g = Flag for gas    ONLY
!       Kepsilon_lg= Flag for BOTH liquid & gas
!
      ELSE IF(Turbulence_Model.eq.Kepsilon.or.Turbulence_Model.eq.Kepsilon_real.or.Turbulence_Model.eq.Kepsilon_RNG) THEN
!      
!........k-e for liquid ONLY
         IF(turb_phase.eq.liq_only) THEN
            CALL Turb_ke_liq
            CALL Turb_ke_vis_liq
!            
!........k-e for gas ONLY
         ELSE IF(turb_phase.eq.gas_only) THEN
            CALL Turb_ke_gas
            CALL Turb_ke_vis_gas
!            
!........k-e for both liquid and gas
         ELSE IF(turb_phase.eq.both) THEN
            CALL Turb_ke_liq
            CALL Turb_ke_gas
            CALL Turb_ke_vis_liq
            CALL Turb_ke_vis_gas
            CALL Turb_ke_vis_lg
         ENDIF
!
!........Calculate thermal conductivity
!         
         CALL Turb_ke_cond
!         
!......4)SST k-w model         
      ELSE IF(Turbulence_Model.eq.SST) THEN
!      
!........k-w for liquid ONLY
         IF(turb_phase.eq.liq_only) THEN
            CALL Turb_SST_liq
            CALL Turb_SST_vis_liq
!
!........k-w for gas ONLY
         ELSE IF(turb_phase.eq.gas_only) THEN
            CALL Turb_SST_gas
            CALL Turb_SST_vis_gas         
!
!........k-w for both liquid and gas
         ELSE IF(turb_phase.eq.both) THEN
            CALL Turb_SST_liq
            CALL Turb_SST_gas
            CALL Turb_SST_vis_liq
            CALL Turb_SST_vis_gas
            CALL Turb_ke_vis_lg        ! Identical function to k-e model         
         ENDIF
!
!........Calculate thermal conductivity
!        
         CALL Turb_ke_cond
!   
!.....7) LES model
!      
      ELSEIF(Turbulence_Model.eq.LES_WALE)THEN
         CALL Turb_LES_WALE(turb_phase)
         CALL Turb_LES_cond(turb_phase)
!
!.....5) Free slip model
!         
      ELSEIF(Turbulence_Model.eq.Free_slip) THEN
!
!.....6) Not proper model          
!
      ELSE
         WRITE(*,*) '          ITURB should be (-2,-1,0,1,2,3,4) !'
         WRITE(unit_log,*)'          ITURB should be (-2,-1,0,1,2,3,4) !'
!         
      ENDIF
!           
      END SUBROUTINE turbulence_mod
