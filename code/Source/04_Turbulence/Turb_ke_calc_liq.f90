!
      SUBROUTINE turb_ke_calc_liq
!
!     This routine discretizes and solves the linearized
!     turbulence kinetic energy equations explicitly.
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: np
      USE Zparam       , ONLY: ke_small,ced1,ced2,cmu,cappa,ke_cff,dp_cff,RNG_cff,ndim,prt,  &
                                ced1_RNG,ced2_RNG,RNG_cff,ced2_Real,ke_cff_Real,dp_cff_Real, &
                               RNG_cffr,ke_cff_Realr,dp_cff_Realr,ke_cffr,dp_cffr
      USE Zare         , ONLY: ar_liq
      USE Zb_condition , ONLY: turb_keb,turb_dpb,alphab_liq,rhob_liq
      USE Zbc_index    , ONLY: icell_type
      USE Zconst1      , ONLY: turboil,buoyancy_turb,lowreynolds,iturb
      USE Zconst2      , ONLY: dt,grav,hydraulicd
      USE Zcoord3      , ONLY: volr,porosity
      USE Zface        , ONLY: qecell,Kepsilon_RNG,Kepsilon_real
      USE Zimplicit    , ONLY: imp_ke_diff,imp_ke_conv,eps_imp_ke,eps_imp_dp,max_iter_ke,max_iter_dp
      USE Zndforce     , ONLY: d_bfc      
      USE Zturb        , ONLY: turb_ke,turb_ke_o,turb_dp,turb_dp_o,pro_ke,diff_ke,diff_dp,yplus,strn_ke,cmu_real
      USE Zturb        , ONLY: s_macroturb_source
      USE Zvector      , ONLY: vrel_o,vl_o
      USE Zvec_major   , ONLY: flux_l_nf
      USE Zio_unit     , ONLY: unit_log
!
      IMPLICIT NONE
!
!     local varaibles
      INTEGER i,iml
      INTEGER idim      
      INTEGER kill
!
      REAL(8) sk_kata,sk_void,se_void,cdk
      REAL(8) hfg,gravity
      REAL(8) am
      REAL(8) ut,un,gravm,beta                                  !for buoyancy
      REAL(8) ced3,a,b,c,d                                      !for buoyancy
      REAL(8) f1,f2,fmu,Ret,a_square,DlowReynolds,ElowReynolds  !for low-Reynolds-number k-e model  
      REAL(8) ced1_mod,ced2_mod,etha                     !for RNG, Realizable k-e model
      REAL(8) hd,rel,fric,velo,velo2,velo3,velo4,yplim,pro_porous(ncell_fluid),kinf,coeff_eps   !LSJ 161122 porous
      REAL(8) KK,bb,cd                                                                           !LSJ 161122 porous       
!     local arrays 
      REAL(8) :: diff_k(ncell_fluid),diff_e(ncell_fluid),conv_k(ncell_fluid),conv_e(ncell_fluid)
      REAL(8) :: bm(ncell_fluid),cm(ncell_fluid)
      REAL(8) :: dtldx(ncell_fluid,ndim) 
!
      ar_liq(:)=cell%alphal(:)*cell%rhol(:) !This should be defined before k-e calculation
      conv_k(:)=0.d0
      conv_e(:)=0.d0
      diff_k(:)=0.d0
      diff_e(:)=0.d0      
!
!.....Initialize constants
!
      CDk=0.44d0
      gravity=9.8d0
! 
!.....Standard k-e turbulence model
!         
      f1=1.d0
      f2=1.d0
      fmu=1.d0
      DlowReynolds=0.d0
      ElowReynolds=0.d0         
!
!.....Diffusivity of kinetic energy and dissipation
!
      kill=0
      IF(iturb.ne.Kepsilon_RNG.and.iturb.ne.Kepsilon_real)THEN
         DO i=1,ncell_fluid
            diff_ke(i)=(cell%lviscosl(i)+cell%tviscosl(i)*ke_cffr)*cell%alphal(i)
            diff_dp(i)=(cell%lviscosl(i)+cell%tviscosl(i)*dp_cffr)*cell%alphal(i)
         ENDDO
      ELSEIF(iturb.eq.Kepsilon_RNG.and.iturb.ne.Kepsilon_real)THEN   !RNG_ke model
         DO i=1,ncell_fluid
            diff_ke(i)=(cell%lviscosl(i)+cell%tviscosl(i)*RNG_cffr)*cell%alphal(i)
            diff_dp(i)=(cell%lviscosl(i)+cell%tviscosl(i)*RNG_cffr)*cell%alphal(i)
         ENDDO     
      ELSEIF(iturb.ne.Kepsilon_RNG.and.iturb.eq.Kepsilon_real)THEN   !Realizable_ke model
         DO i=1,ncell_fluid
            diff_ke(i)=(cell%lviscosl(i)+cell%tviscosl(i)*ke_cff_Realr)*cell%alphal(i)
            diff_dp(i)=(cell%lviscosl(i)+cell%tviscosl(i)*dp_cff_Realr)*cell%alphal(i)
         ENDDO    
      ELSE
         kill=1
         WRITE(*,*)  'Input of Turbulence model option is incorrect!' 
         WRITE(unit_log,*) 'Input of Turbulence model option is incorrect!'        
      ENDIF
      IF(np.gt.1) THEN
         CALL allreducei_i1(kill)
         IF(kill.gt.0) THEN
            CALL barrier_mpi
            CALL finalize_mpi
            STOP
         ENDIF
      ELSE
         IF(kill.gt.0) THEN
            STOP
         ENDIF
      ENDIF
!
!.....Communicate transport variables for parallel computing
!
      IF(np.gt.1) CALL communicate_1d(turb_ke_o, &
                                      turb_dp_o, &
                                      diff_ke,   &
                                      diff_dp,   &
                                      ar_liq)
!
!.....Convection and diffusion of k-e
!
      IF(imp_ke_diff.eq.0) THEN  
         CALL turb_ke_diffusion(turb_ke_o,turb_dp_o,turb_keb,turb_dpb,diff_ke,diff_dp,diff_k,diff_e)
      ENDIF
      IF(imp_ke_conv.eq.0) THEN  
         CALL turb_ke_convection(turb_ke_o,turb_dp_o,ar_liq,turb_keb,turb_dpb,alphab_liq,rhob_liq,flux_l_nf,conv_k,conv_e)      
      ENDIF
!
!.....For buoyancy effect
!
      IF(buoyancy_turb.eq.1) THEN
         CALL grad_scalar(cell%tl_o,dtldx,ncell_fluid)
      ENDIF
!
!.....Kinetic energy, k
!
      bm(:)=0.0d0    
      cm(:)=0.0d0 
      DO i=1,ncell_fluid
!
!........Coeff. for Low Reynolds model (Chien): DlowReynolds
!
         IF(lowreynolds.eq.1) DlowReynolds=-2.d0*cell%lviscosl(i)*turb_ke_o(i)/d_bfc(i)/d_bfc(i)
!      
         sk_void=0.0d0
         sk_kata=0.0d0
!
!........Source term by bubble motion (Lahey,2005)
!
         IF(turboil.eq.1) THEN 
            sk_void=0.25d0*(1.0d0+CDk**(4.d0/3.d0))*cell%alphag(i)*vrel_o(i)*vrel_o(i)*vrel_o(i)/cell%d1(i)
            sk_void=cell%alphal(i)*sk_void
!
!...........Kinetic energy source term by wall boiling (Kataoka and Serizawa,1997)
!
            IF(icell_type(i).eq.1) THEN
               hfg=cell%hgsat(i)-cell%hlsat(i)
               sk_kata=8.0d0/9.d0*gravity*(cell%rhol(i)-cell%rhog(i))/cell%rhol(i)*qecell(i)/(cell%rhog(i)*hfg)
            ENDIF
         ENDIF
!
!........Source term by convection and diffusion
!
         IF(imp_ke_diff.eq.0) THEN      
            diff_k(i)=diff_k(i)*volr(i)
            diff_e(i)=diff_e(i)*volr(i)
         ENDIF            
         IF(imp_ke_conv.eq.0) THEN      
            conv_k(i)=conv_k(i)*volr(i)
            conv_e(i)=conv_e(i)*volr(i)            
         ENDIF   
!
!........bm for the Kinetic energy, k
!
!choi         bm(i)=-conv_k(i)+diff_k(i)+cell%alphal(i)*pro_ke(i)-ar_liq(i)*turb_dp_o(i)+sk_void+sk_kata+DlowReynolds
         IF(imp_ke_diff.eq.0.and.imp_ke_conv.eq.0) THEN !explicit source
            bm(i)=-conv_k(i)+diff_k(i)+cell%alphal(i)*pro_ke(i)+sk_void+sk_kata+DlowReynolds -ar_liq(i)*turb_dp_o(i)
         ELSE       !implicit source 
            bm(i)=-conv_k(i)+diff_k(i)+cell%alphal(i)*pro_ke(i)+sk_void+sk_kata+DlowReynolds !-ar_liq(i)*turb_dp_o(i)
            cm(i)=-ar_liq(i)*turb_dp_o(i)/turb_ke_o(i)
         ENDIF
!
!........Buoyancy Effect for k
!  
         IF(buoyancy_turb.eq.1) THEN
            DO idim=1,ndim
               beta=cell%betal(i)
               bm(i)=bm(i)+cell%alphal(i)*beta*grav(idim)*cell%tviscosl(i)/prt*dtldx(i,idim)
            ENDDO
         ENDIF
!         
!........udfl_turb_source for k (porous media): macroscopic turbulence source (LSJ)
!
         IF(s_macroturb_source.eq.'chandesris') THEN
            hd=hydraulicd(i)
            velo=MAXVAL(DABS(vl_o(i,:)))
            rel=DMAX1(1.0d0,(cell%rhol(i)*velo*hd/cell%lviscosl(i)))
            fric=0.25d0*rel**(-0.25d0)
            velo3=velo*velo*velo
            yplim=7.d0 !user-defined
!            yplim=30.d0 !user-defined
            pro_porous(i)=2.d0*fric*velo3/hd*(1.d0-yplim*DSQRT(fric/2.d0)) !production of k (Pk) for porous media
!            
            bm(i)=bm(i)+cell%alphal(i)*cell%rhol(i)*pro_porous(i)
         ELSEIF(s_macroturb_source.eq.'nakayama') THEN             
            hd=hydraulicd(i)
            KK=porosity(i)*hd*hd/32.d0 
            velo=MAXVAL(DABS(vl_o(i,:)))               
            velo3=velo*velo*velo
            bb=0.3164d0/(2.d0*hd*porosity(i)*porosity(i)*(hd*DMAX1(1.e-8,velo)/cell%lviscosl(i)*cell%rhol(i))**0.25d0) 
            pro_porous(i)=porosity(i)*porosity(i)*bb*velo3
!            
            bm(i)=bm(i)+cell%alphal(i)*cell%rhol(i)*pro_porous(i)
         ENDIF   
!         
         IF(imp_ke_diff.eq.0.and.imp_ke_conv.eq.0) THEN
            am=ar_liq(i)/dt-cm(i)          
            IF(cell%alphal(i).gt.1.0d-8) THEN
               turb_ke(i)=turb_ke_o(i)+bm(i)/am
            ELSE
               turb_ke(i)=turb_ke_o(i)
            ENDIF
            turb_ke(i)=DMAX1(turb_ke(i),1.d-15)            
         ENDIF   
      ENDDO  
!      
      IF(imp_ke_diff.eq.1.or.imp_ke_conv.eq.1) THEN
!
!........1:liquid-phase, 2:gas-phase        
!
         iml=1
         CALL imp_diffusion_ke(1,bm,cm,turb_ke_o,turb_keb,diff_ke,turb_ke,ar_liq,flux_l_nf,alphab_liq,rhob_liq,eps_imp_ke,max_iter_ke,iml)
         turb_ke(:)=DMAX1(turb_ke(:),1.d-15)            
      ENDIF   
!
!.....Dissipation, e
!
      bm(:)=0.0d0   
      cm(:)=0.0d0    
      DO i=1,ncell_fluid
         sk_void=0.0d0
         se_void=0.0d0
!
!........Coeff. for Low Reynolds model (Chien): f1,f2,ElowReynolds
!
         IF(lowreynolds.eq.1) THEN
            f1=1.d0
            Ret=turb_ke_o(i)*turb_ke_o(i)*cell%rhol(i)/cell%lviscosl(i)/turb_dp_o(i)
            a_square=Ret*Ret/6.d0/6.d0
            f2=1.d0-2.d0/9.d0*DEXP(-a_square)
            ElowReynolds=-2.d0*cell%lviscosl(i)*turb_dp_o(i)/d_bfc(i)/d_bfc(i)*DEXP(-0.5d0*yplus(i))
         ENDIF   
!                  
         IF(turboil.eq.1) THEN 
            sk_void=0.25d0*(1.0d0+CDk**(4.d0/3.d0))*cell%alphag(i)*vrel_o(i)*vrel_o(i)*vrel_o(i)/cell%d1(i)
            sk_void=cell%alphal(i)*sk_void
            IF(turb_ke(i).gt.ke_small) se_void=ced2*turb_dp_o(i)/turb_ke(i)*sk_void
         ENDIF   
!
         IF(iturb.ne.Kepsilon_RNG.and.iturb.ne.Kepsilon_real)THEN       ! Standard k-e
            ced1_mod=ced1
            ced2_mod=ced2
         ELSEIF(iturb.eq.Kepsilon_RNG.and.iturb.ne.Kepsilon_real)THEN   ! RNG k-e
            ced1_mod=ced1_RNG
            IF(turb_dp_o(i).ne.0.0d0)THEN
               etha=strn_ke(i)*turb_ke(i)/turb_dp_o(i)
               ced2_mod=ced2_RNG+cmu*etha**3*(1.0d0-etha/4.38d0)/(1.0d0+0.012d0*etha**3)
            ELSE
               ced2_mod=ced2_RNG
            ENDIF            
         ELSEIF(iturb.ne.Kepsilon_RNG.and.iturb.eq.Kepsilon_real)THEN   ! Realizable k-e
            IF(turb_dp_o(i).ne.0.0d0)THEN
               etha=strn_ke(i)*turb_ke(i)/turb_dp_o(i)
               ced1_mod=DMAX1(0.43d0,etha/(etha+5.0d0))
            ELSE
               ced1_mod=0.43d0
            ENDIF           
            ced2_mod=ced2_Real
         ELSE
            WRITE(*,*) 'Input of Turbulence model option is incorrect!' 
            WRITE(unit_log,*) 'Input of Turbulence model option is incorrect!'
         ENDIF
!         
         IF(iturb.ne.Kepsilon_real)THEN
            IF(imp_ke_diff.eq.0.and.imp_ke_conv.eq.0) THEN !explicit source
               bm(i)=cell%alphal(i)/turb_ke(i)*(f1*ced1_mod*pro_ke(i)-f2*ced2_mod*cell%rhol(i)*turb_dp_o(i))
               bm(i)=-conv_e(i)+diff_e(i)+se_void +bm(i)*turb_dp_o(i) +ElowReynolds
            ELSE      !implicit source
               bm(i)=-conv_e(i)+diff_e(i)+se_void +ElowReynolds  !+bm(i)*turb_dp_o(i)
               cm(i)=cell%alphal(i)/turb_ke(i)*(f1*ced1_mod*pro_ke(i)-f2*ced2_mod*cell%rhol(i)*turb_dp_o(i))
            ENDIF
         ELSE
            IF(imp_ke_diff.eq.0.and.imp_ke_conv.eq.0) THEN !explicit source
               bm(i)=cell%alphal(i)*(f1*cell%rhol(i)*ced1_mod*strn_ke(i)-f2*ced2_mod*cell%rhol(i)*turb_dp_o(i)/(turb_ke(i)+DSQRT(cell%lviscosl(i)/cell%rhol(i)*turb_dp_o(i))))
               bm(i)=-conv_e(i)+diff_e(i)+se_void +bm(i)*turb_dp_o(i) +ElowReynolds
            ELSE      !implicit source
               bm(i)=-conv_e(i)+diff_e(i)+se_void +ElowReynolds  !+bm(i)*turb_dp_o(i)
               cm(i)=cell%alphal(i)*(f1*cell%rhol(i)*ced1_mod*strn_ke(i)-f2*ced2_mod*cell%rhol(i)*turb_dp_o(i)/(turb_ke(i)+DSQRT(cell%lviscosl(i)/cell%rhol(i)*turb_dp_o(i))))
            ENDIF         
         ENDIF
!
!........udfl_turb_source for epsilon (porous media): macroscopic source
!       
         IF(s_macroturb_source.eq.'chandesris') THEN         
            hd=hydraulicd(i)
            velo=MAXVAL(DABS(vl_o(i,:)))
            velo2=velo*velo
            rel=DMAX1(1.0d0,(cell%rhol(i)*velo*hd/cell%lviscosl(i)))
            coeff_eps=0.0368d0 !user defined           
!            coeff_eps=0.05d0 !user defined           
            kinf=coeff_eps*velo2*rel**(-1.d0/6.d0)
!              
            bm(i)=bm(i)+cell%alphal(i)*cell%rhol(i)*ced2*pro_porous(i)*pro_porous(i)/kinf
         ELSEIF(s_macroturb_source.eq.'nakayama') THEN   
            hd=hydraulicd(i)
            velo=MAXVAL(DABS(vl_o(i,:)))
            velo4=velo*velo*velo*velo
            KK=DMAX1(1.e-8,porosity(i)*hd*hd/32.d0) 
            bb=0.3164d0/(2.d0*hd*porosity(i)*porosity(i)*(hd*DMAX1(1.e-8,velo)/cell%lviscosl(i)*cell%rhol(i))**0.25d0)                
            cd=0.09d0 !user-defined
            pro_porous(i)=(ced1*ced1/ced2)*porosity(i)*porosity(i)*bb*DSQRT(cd*porosity(i)/2.d0/KK)*velo4
!            pro_porous(i)=ced2*porosity(i)*porosity(i)*bb*DSQRT(cd*porosity(i)/2.d0/KK)*velo4
!
            bm(i)=bm(i)+cell%alphal(i)*cell%rhol(i)*pro_porous(i)
         ENDIF                   
!
!........Buyoancy Effect for e
!         
         IF(buoyancy_turb.eq.1) THEN
            ut=0.d0
            un=0.d0
            gravm=0.d0
            DO idim=1,ndim
               ut=ut+grav(idim)*vl_o(i,idim)
               un=un+vl_o(i,idim)*vl_o(i,idim)
               gravm=gravm+grav(idim)*grav(idim)
            ENDDO   
            gravm=gravm**0.5
            ut=ut/gravm
            un=(DMAX1(1d-9,un-ut*ut))**0.5
            a=ut/un
            b=DMAX1(0.d0,DMIN1(20.d0,DABS(a)))
            c=DEXP(b)
            d=DEXP(-b)
            ced3=(c-d)/(c+d)
            IF(lowreynolds.eq.1)ced3=1
!
            DO idim=1,ndim
               beta=cell%betal(i)
               cm(i)=cm(i)+cell%alphal(i)/turb_ke(i)*ced1*ced3*beta*grav(idim)*cell%tviscosl(i)/prt*dtldx(i,idim)
            ENDDO
        ENDIF   
!
        IF(imp_ke_diff.eq.0.and.imp_ke_conv.eq.0) THEN 
            am=ar_liq(i)/dt-cm(i)        
            IF(cell%alphal(i).gt.1.0d-8) THEN
              turb_dp(i)=ar_liq(i)/dt/am*turb_dp_o(i)+bm(i)/am
            ELSE
              turb_dp(i)=turb_dp_o(i)
            ENDIF  
            turb_dp(i)=DMAX1(turb_dp(i),1.d-15)      
        ENDIF           
      ENDDO 
!      
      IF(imp_ke_diff.eq.1.or.imp_ke_conv.eq.1) THEN
         iml=2
         CALL imp_diffusion_ke(1,bm,cm,turb_dp_o,turb_dpb,diff_dp,turb_dp,ar_liq,flux_l_nf,alphab_liq,rhob_liq,eps_imp_dp,max_iter_dp,iml)    
         turb_dp(:)=DMAX1(turb_dp(:),1.d-15)      
      ENDIF   
!
!.....Dissipation at the boundary cells next to the wall
!
      DO i=1,ncell_fluid
         IF(icell_type(i).eq.1)THEN
            IF(Kepsilon_real.ne.1)THEN
               turb_dp(i)=(cmu**0.75*turb_ke(i)**1.5d0)/(cappa*d_bfc(i))
            ELSE
               turb_dp(i)=(cmu_Real(i)**0.75*turb_ke(i)**1.5d0)/(cappa*d_bfc(i))
            ENDIF
         ENDIF
      ENDDO
!
      RETURN
      END SUBROUTINE turb_ke_calc_liq
