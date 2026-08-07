!
      SUBROUTINE turb_ke_calc_gas
!
!     This routine discretizes and solves the linearized
!     turbulence kinetic energy equations explicitly.
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: np
      USE Zparam       , ONLY: ke_small,ced1,ced2,cmu,cappa,ndim,prtr,             &
                               ced1_RNG,ced2_RNG,ced2_Real,                       &
                               RNG_cffr,ke_cff_Realr,dp_cff_Realr,ke_cffr,dp_cffr
      USE Zare         , ONLY: ar_gas
      USE Zb_condition , ONLY: turb_kegb,turb_dpgb,alphab_gas,rhob_gas,vb_gas,vin_gas
      USE Zbc_index    , ONLY: icell_type
      USE Zconst1      , ONLY: buoyancy_turb,lowreynolds,iturb
      USE Zface        , ONLY: Kepsilon_RNG,Kepsilon_real
      USE Zconst2      , ONLY: dtr,grav,hydraulicd
      USE Zcoord3      , ONLY: volr,porosity
      USE Zimplicit    , ONLY: imp_ke_diff,imp_ke_conv,eps_imp_ke,eps_imp_dp,max_iter_ke,max_iter_dp
      USE Zndforce     , ONLY: d_bfc          
      USE Zturb        , ONLY: turb_keg,turb_keg_o,turb_dpg,turb_dpg_o,pro_keg,diff_ke,diff_dp,yplusg,strn_keg,cmug_real
      USE Zturb        , ONLY: s_macroturb_source
      USE Zvector      , ONLY: vg_o
      USE Zvec_major   , ONLY: flux_g_nf
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER i,iml
      INTEGER ix
      INTEGER kill
      REAL(8) t1,t2
      REAL(8) am
      REAL(8) ut,un,gravm
      REAL(8) ced3,a,b,c,d                                      !for buoyancy
      REAL(8) f1,f2,fmu,Ret,a_square,DlowReynolds,ElowReynolds  !for low Reynolds k-e model    
      REAL(8) ced1_mod,ced2_mod,etha                     !for RNG, Realizable k-e model
      REAL(8) hd,reg,fric,velo,velo2,velo3,velo4,yplim,kinf,coeff_eps                           !LSJ 161122 porous
      REAL(8) KK,bb,cd                                                                          !LSJ 161122 porous       
      REAL(8) Eii,Eij,gradE                                     !for lowreynolds=2
      REAL(8) Gb
!.....Local arrays 
      REAL(8) :: diff_k(ncell_fluid),diff_e(ncell_fluid),conv_k(ncell_fluid),conv_e(ncell_fluid)
      REAL(8) :: bm(ncell_fluid),cm(ncell_fluid)
      REAL(8) :: drgdx(ncell_fluid,ndim)
      REAL(8) :: pro_porous(ncell_fluid)                                                        !LSJ 161122 porous
      REAL(8) :: dvdx(ncell_fp,ndim,ndim)
! 
      ar_gas(:)=cell%alphag(:)*cell%rhog(:) !This should be defined before k-e calculation    
      DO i=1,ncell_fluid
         conv_k(i)=0.d0
         conv_e(i)=0.d0
         diff_k(i)=0.d0
         diff_e(i)=0.d0      
      ENDDO     
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
            diff_ke(i)=(cell%lviscosg(i)+cell%tviscosg(i)*ke_cffr)*cell%alphag(i)
            diff_dp(i)=(cell%lviscosg(i)+cell%tviscosg(i)*dp_cffr)*cell%alphag(i)
         ENDDO
      ELSEIF(iturb.eq.Kepsilon_RNG.and.iturb.ne.Kepsilon_real)THEN     !RNG_ke model
         DO i=1,ncell_fluid
            diff_ke(i)=(cell%lviscosg(i)+cell%tviscosg(i)*RNG_cffr)*cell%alphag(i)
            diff_dp(i)=(cell%lviscosg(i)+cell%tviscosg(i)*RNG_cffr)*cell%alphag(i)
         ENDDO
      ELSEIF(iturb.ne.Kepsilon_RNG.and.iturb.eq.Kepsilon_real)THEN     !Realizable_ke model
         DO i=1,ncell_fluid
            diff_ke(i)=(cell%lviscosg(i)+cell%tviscosg(i)*ke_cff_Realr)*cell%alphag(i)
            diff_dp(i)=(cell%lviscosg(i)+cell%tviscosg(i)*dp_cff_Realr)*cell%alphag(i)
         ENDDO
      ELSE
         kill=1
         WRITE(*,*)  'Input of Turbulence model option is incorrect!'
         WRITE(97,*) 'Input of Turbulence model option is incorrect!'        
         STOP
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
      IF(np.gt.1) CALL communicate_1d(turb_keg_o, &
                                      turb_dpg_o, &
                                      diff_ke,    &
                                      diff_dp,    &
                                      ar_gas)
!
!.....Convection and diffusion of k-e
!
      IF(imp_ke_diff.eq.0) THEN      
         CALL turb_ke_diffusion(turb_keg_o,turb_dpg_o,turb_kegb,turb_dpgb,diff_ke,diff_dp,diff_k,diff_e)
!
!........Source term by convection and diffusion
!        
         DO i=1,ncell_fluid   
            diff_k(i)=diff_k(i)*volr(i)
            diff_e(i)=diff_e(i)*volr(i)
         ENDDO
      ENDIF      
      IF(imp_ke_conv.eq.0) THEN      
         CALL turb_ke_convection(turb_keg_o,turb_dpg_o,ar_gas,turb_kegb,turb_dpgb,alphab_gas,rhob_gas,flux_g_nf,conv_k,conv_e)      
!
!........Source term by convection and diffusion
!        
         DO i=1,ncell_fluid   
            conv_k(i)=conv_k(i)*volr(i)
            conv_e(i)=conv_e(i)*volr(i)
         ENDDO
      ENDIF      
!
!.....For buoyancy effect
!
      IF(buoyancy_turb.eq.1) THEN
         CALL grad_scalar(cell%rhog,drgdx,ncell_fluid)
      ENDIF
!      
!.....Kinetic energy, k
!         
!
!.....Coeff. for Low Reynolds model (Chien): DlowReynolds
!
      IF(lowreynolds.ge.1) THEN
         IF(imp_ke_diff.eq.0.and.imp_ke_conv.eq.0) THEN !explicit source
            DlowReynolds=-2.d0*cell%lviscosg(i)*turb_keg_o(i)/(d_bfc(i)**2)
            DO i=1,ncell_fluid   
               !bm(i)=-conv_k(i)+diff_k(i)+cell%alphag(i)*pro_keg(i)+DlowReynolds-ar_gas(i)*turb_dpg_o(i)
               bm(i)=-conv_k(i)+diff_k(i)+cell%alphag(i)*pro_keg(i)+DlowReynolds -ar_gas(i)*turb_dpg_o(i)
               cm(i)=0.d0
            ENDDO
         ELSE
            DO i=1,ncell_fluid   
               !bm(i)=-conv_k(i)+diff_k(i)+cell%alphag(i)*pro_keg(i)+DlowReynolds-ar_gas(i)*turb_dpg_o(i)
               bm(i)=cell%alphag(i)*pro_keg(i)+DlowReynolds !-ar_gas(i)*turb_dpg_o(i)
               cm(i)=-ar_gas(i)*turb_dpg_o(i)/turb_keg_o(i)
            ENDDO
         ENDIF
      ELSE
         IF(imp_ke_diff.eq.0.and.imp_ke_conv.eq.0) THEN !explicit source
            DO i=1,ncell_fluid   
               !bm(i)=-conv_k(i)+diff_k(i)+cell%alphag(i)*pro_keg(i)-ar_gas(i)*turb_dpg_o(i)
               bm(i)=-conv_k(i)+diff_k(i)+cell%alphag(i)*pro_keg(i) -ar_gas(i)*turb_dpg_o(i)
               cm(i)=0.d0
            ENDDO
         ELSE
            DO i=1,ncell_fluid   
               !bm(i)=-conv_k(i)+diff_k(i)+cell%alphag(i)*pro_keg(i)-ar_gas(i)*turb_dpg_o(i)
               bm(i)=cell%alphag(i)*pro_keg(i) !-ar_gas(i)*turb_dpg_o(i)
               cm(i)=-ar_gas(i)*turb_dpg_o(i)/turb_keg_o(i)
            ENDDO
         ENDIF
      ENDIF
!
!........Buoyancy Effect for k
!  
      IF(buoyancy_turb.eq.1) THEN
         IF(ndim.eq.2) THEN
            DO i=1,ncell_fluid
               Gb=grav(1)*drgdx(i,1)+grav(2)*drgdx(i,2)
               Gb=-cell%tviscosg(i)*Gb/cell%rhog(i)*prtr
               bm(i)=bm(i)+cell%alphag(i)*Gb
            ENDDO
         ELSE
            DO i=1,ncell_fluid
               Gb=grav(1)*drgdx(i,1)+grav(2)*drgdx(i,2)+grav(3)*drgdx(i,3)
               Gb=-cell%tviscosg(i)*Gb/cell%rhog(i)*prtr
               bm(i)=bm(i)+cell%alphag(i)*Gb
            ENDDO
         ENDIF 
      ENDIF 
!         
!........udfl_turb_source for k (porous media): macroscopic turbulence source (LSJ)
!
      IF(s_macroturb_source.eq.'chandesris') THEN
         DO i=1,ncell_fluid
            hd=hydraulicd(i)
            velo=MAXVAL(ABS(vg_o(i,:)))
            reg=MAX(1.0d0,(cell%rhog(i)*velo*hd/cell%lviscosg(i)))
            fric=0.25d0*reg**(-0.25d0)
            velo3=velo*velo*velo
            yplim=7.d0 !user-defined
!            yplim=30.d0 !user-defined
            pro_porous(i)=2.d0*fric*velo3/hd*(1.d0-yplim*SQRT(fric/2.d0)) !production of k (Pk) for porous media
!
            bm(i)=bm(i)+cell%alphag(i)*cell%rhog(i)*pro_porous(i)
         ENDDO
!            
      ELSEIF(s_macroturb_source.eq.'nakayama') THEN             
         DO i=1,ncell_fluid
            hd=hydraulicd(i)
            KK=porosity(i)*hd*hd/32.d0 
            velo=MAXVAL(ABS(vg_o(i,:)))               
            velo3=velo*velo*velo
            t1=sqrt(hd*MAX(1.d-8,velo)/cell%lviscosg(i)*cell%rhog(i))
            bb=0.3164d0/(2.d0*hd*porosity(i)*porosity(i)*t1*sqrt(t1)) 
            pro_porous(i)=porosity(i)*porosity(i)*bb*velo3
!
            bm(i)=bm(i)+cell%alphag(i)*cell%rhog(i)*pro_porous(i)
         ENDDO
      ENDIF              
!         
      IF(imp_ke_diff.eq.0.and.imp_ke_conv.eq.0) THEN  
         DO i=1,ncell_fluid
            am=ar_gas(i)*dtr
            IF(cell%alphag(i).gt.1.0d-8) THEN
               turb_keg(i)=turb_keg_o(i)+bm(i)/am            
            ELSE
               turb_keg(i)=turb_keg_o(i)
            ENDIF
            turb_keg(i)=MAX(turb_keg(i),0.d0)            
         ENDDO
      ENDIF
!
      IF(imp_ke_diff.eq.1.or.imp_ke_conv.eq.1) THEN       
!
!........1:liquid-phase, 2:gas-phase
!
         iml=1
         CALL imp_diffusion_ke(2,bm,cm,turb_keg_o,turb_kegb,diff_ke,turb_keg,ar_gas,flux_g_nf,alphab_gas,rhob_gas,eps_imp_ke,max_iter_ke,iml)        
         turb_keg(:)=MAX(turb_keg,1.d-15)           
      ENDIF            
!
!.....Dissipation, e
!
      IF(lowreynolds.eq.2) THEN
         CALL grad_vel(1,vg_o,dvdx,vb_gas,vin_gas)
      ENDIF
!
      DO i=1,ncell_fluid
!
!........Coeff. for Low Reynolds model (Chien): f1,f2,ElowReynolds
!
         IF(lowreynolds.eq.1) THEN
            f1=1.d0
            IF(turb_dpg_o(i).ge.ke_small) THEN 
               Ret=turb_keg(i)*turb_keg(i)*cell%rhog(i)/cell%lviscosg(i)/turb_dpg_o(i)
               a_square=(Ret/6.d0)**2
               f2=1.d0-2.d0/9.d0*EXP(-a_square)
            ENDIF   
            ElowReynolds=-2.d0*cell%lviscosg(i)*turb_dpg_o(i)/(d_bfc(i)**2)*EXP(-0.5d0*yplusg(i))
         ELSEIF(lowreynolds.eq.2)THEN
            f1=1.d0
            IF(turb_dpg_o(i).ge.ke_small)THEN 
               Ret=turb_keg(i)*turb_keg(i)*cell%rhog(i)/cell%lviscosg(i)/turb_dpg_o(i)
               a_square=Ret*Ret
               f2=1.0d0-0.3d0*EXP(-a_square)
            ENDIF
!
            Eii=0.0d0
            Eij=0.0d0
            DO ix=1,ndim
               Eii=Eii+dvdx(i,ix,ix)*dvdx(i,ix,ix)
            ENDDO
            Eij=dvdx(i,1,1)*dvdx(i,2,1)+dvdx(i,2,2)*dvdx(i,1,2)
            IF(ndim.eq.3)THEN
               Eij=Eij+dvdx(i,1,1)*dvdx(i,3,1)+dvdx(i,2,2)*dvdx(i,3,2)+dvdx(i,3,3)*(dvdx(i,1,3)+dvdx(i,2,3))
            ENDIF
            gradE=Eii+Eij
            ElowReynolds=2.0d0*cell%lviscosg(i)*cell%tviscosg(i)*gradE*gradE/cell%rhog(i)
         ENDIF   
!
         IF(iturb.eq.Kepsilon_RNG) THEN    ! RNG k-e
            ced1_mod=ced1_RNG
            IF(turb_dpg_o(i).ne.0.0d0)THEN
               etha=strn_keg(i)*turb_keg(i)/turb_dpg_o(i)
               ced2_mod=ced2_RNG+cmu*etha**3*(1.0d0-etha/4.38d0)/(1.0d0+0.012d0*etha**3)
            ELSE
               ced2_mod=ced2_RNG
            ENDIF            
         ELSEIF(iturb.eq.Kepsilon_real)THEN    ! Realizable k-e
            IF(turb_dpg_o(i).ne.0.0d0)THEN
               etha=strn_keg(i)*turb_keg(i)/turb_dpg_o(i)
               ced1_mod=MAX(0.43d0,etha/(etha+5.0d0))
            ELSE
               ced1_mod=0.43d0
            ENDIF           
            ced2_mod=ced2_Real
         ELSE
            ced1_mod=ced1
            ced2_mod=ced2
         ENDIF
!         
         IF(iturb.ne.Kepsilon_real)THEN
            IF(imp_ke_diff.eq.0.and.imp_ke_conv.eq.0) THEN !explicit source
               bm(i)= cell%alphag(i)/turb_keg(i)                                      &
                     *(f1*ced1_mod*pro_keg(i)-f2*ced2_mod*cell%rhog(i)*turb_dpg_o(i))
               bm(i)=-conv_e(i)+diff_e(i)+bm(i)*turb_dpg_o(i) +ElowReynolds
               cm(i)=0.d0
            ELSE      !implicit source
               bm(i)=-conv_e(i)+diff_e(i)+ElowReynolds  !+bm(i)*turb_dp_o(i)
               cm(i)= cell%alphag(i)/turb_keg(i)                                      &
                     *(f1*ced1_mod*pro_keg(i)-f2*ced2_mod*cell%rhog(i)*turb_dpg_o(i))
            ENDIF
         ELSE
            IF(imp_ke_diff.eq.0.and.imp_ke_conv.eq.0) THEN !explicit source
               bm(i)= cell%alphag(i)                                                      &
                     *( f1*cell%rhog(i)*ced1_mod*strn_keg(i)                              &
                       -f2*ced2_mod*cell%rhog(i)*turb_dpg_o(i)                            &
                       /(turb_keg(i)+SQRT(cell%lviscosg(i)/cell%rhog(i)*turb_dpg_o(i))) )
               bm(i)=-conv_e(i)+diff_e(i)+bm(i)*turb_dpg_o(i) +ElowReynolds
               cm(i)=0.d0
            ELSE      !implicit source
               bm(i)=-conv_e(i)+diff_e(i)+ElowReynolds  !+bm(i)*turb_dp_o(i)
               cm(i)= cell%alphag(i)                                                     &
                     *( f1*cell%rhog(i)*ced1_mod*strn_keg(i)                             &
                       -f2*ced2_mod*cell%rhog(i)*turb_dpg_o(i)                           &
                       /(turb_keg(i)+SQRT(cell%lviscosg(i)/cell%rhog(i)*turb_dpg_o(i))))
            ENDIF         
         ENDIF
      ENDDO 
!
!........udfl_turb_source for epsilon (porous media)
!       
      IF(s_macroturb_source.eq.'chandesris') THEN         
         DO i=1,ncell_fluid
            hd=hydraulicd(i)
            velo=MAXVAL(ABS(vg_o(i,:)))
            velo2=velo*velo
            reg=MAX(1.0d0,(cell%rhog(i)*velo*hd/cell%lviscosg(i)))
            coeff_eps=0.0368d0 !user defined           
!            coeff_eps=0.05d0 !user defined           
            kinf=coeff_eps*velo2*reg**(-1.d0/6.d0)
!              
            bm(i)=bm(i)+cell%alphag(i)*cell%rhog(i)*ced2*pro_porous(i)*pro_porous(i)/kinf
         ENDDO 
      ELSEIF(s_macroturb_source.eq.'nakayama') THEN   
         DO i=1,ncell_fluid
            hd=hydraulicd(i)
            velo=MAXVAL(ABS(vg_o(i,:)))
            velo4=velo*velo*velo*velo
            KK=MAX(1.d-8,porosity(i)*hd*hd/32.d0) 
            t1=sqrt(hd*MAX(1.d-8,velo)/cell%lviscosg(i)*cell%rhog(i))
            bb=0.3164d0/(2.d0*hd*porosity(i)*porosity(i)*t1*sqrt(t1)) 
            cd=0.09d0 !user-defined
            pro_porous(i)=(ced1*ced1/ced2)*porosity(i)*porosity(i)*bb*SQRT(cd*porosity(i)/2.d0/KK)*velo4
!            pro_porous(i)=ced2*porosity(i)*porosity(i)*bb*SQRT(cd*porosity(i)/2.d0/KK)*velo4
!
            bm(i)=bm(i)+cell%alphag(i)*cell%rhog(i)*pro_porous(i)
         ENDDO 
      ENDIF          
!
!........Buyoancy Effect for e
!            
      IF(buoyancy_turb.eq.1) THEN
         IF(ndim.eq.2) THEN
            DO i=1,ncell_fluid
               ut=grav(1)*vg_o(i,1)+grav(2)*vg_o(i,2)
               un=vg_o(i,1)*vg_o(i,1)+vg_o(i,2)*vg_o(i,2)
               gravm=grav(1)*grav(1)+grav(2)*grav(2)
               gravm=sqrt(gravm)
               ut=ut/gravm
               un=sqrt(MAX(1d-9,un-ut*ut))
               a=ut/un
               b=MAX(0.d0,MIN(20.d0,ABS(a)))
               c=EXP(b)
!              d=EXP(-b)
               d=1.d0/c
               ced3=(c-d)/(c+d)
               cell%ced33(i)=ced3
!
               Gb=grav(1)*drgdx(i,1)+grav(2)*drgdx(i,2)
               Gb=-cell%tviscosg(i)*Gb/cell%rhog(i)*prtr
               IF(imp_ke_diff.eq.0.and.imp_ke_conv.eq.0) THEN !explicit source
                  bm(i)=bm(i)+cell%alphag(i)*f1*ced1*ced3*Gb*turb_dpg_o(i)/turb_keg(i)
               ELSE                                           !implicit source
                  cm(i)=cm(i)+cell%alphag(i)*f1*ced1*ced3*Gb/turb_keg(i)
               ENDIF
            ENDDO 
         ELSE
            DO i=1,ncell_fluid
               ut=grav(1)*vg_o(i,1)+grav(2)*vg_o(i,2)+grav(3)*vg_o(i,3)
               un=vg_o(i,1)*vg_o(i,1)+vg_o(i,2)*vg_o(i,2)+vg_o(i,3)*vg_o(i,3)
               gravm=grav(1)*grav(1)+grav(2)*grav(2)+grav(3)*grav(3)
               gravm=sqrt(gravm)
               ut=ut/gravm
               un=sqrt(MAX(1d-9,un-ut*ut))
               a=ut/un
               b=MAX(0.d0,MIN(20.d0,ABS(a)))
               c=EXP(b)
!              d=EXP(-b)
               d=1.d0/c
               ced3=(c-d)/(c+d)
               cell%ced33(i)=ced3
!
               Gb=grav(1)*drgdx(i,1)+grav(2)*drgdx(i,2)+grav(3)*drgdx(i,3)
               Gb=-cell%tviscosg(i)*Gb/cell%rhog(i)*prtr
               IF(imp_ke_diff.eq.0.and.imp_ke_conv.eq.0) THEN !explicit source
                  bm(i)=bm(i)+cell%alphag(i)*f1*ced1*ced3*Gb*turb_dpg_o(i)/turb_keg(i)
               ELSE                                           !implicit source
                  cm(i)=cm(i)+cell%alphag(i)*f1*ced1*ced3*Gb/turb_keg(i)
               ENDIF
            ENDDO 
         ENDIF
      ENDIF   
!
      IF(imp_ke_diff.eq.0.and.imp_ke_conv.eq.0) THEN 
         DO i=1,ncell_fluid
            am=ar_gas(i)*dtr-cm(i)                
            IF(cell%alphag(i).gt.1.0d-8) THEN
               turb_dpg(i)=ar_gas(i)*dtr/am*turb_dpg_o(i)+bm(i)/am
            ELSE
               turb_dpg(i)=turb_dpg_o(i)
            ENDIF
            turb_dpg(i)=MAX(turb_dpg(i),0.d0)            
         ENDDO 
      ENDIF
!
!
!.....1:liquid-phase, 2:gas-phase        
!
      IF(imp_ke_diff.eq.1.or.imp_ke_conv.eq.1) THEN 
         iml=2
         CALL imp_diffusion_ke(2,bm,cm,turb_dpg_o,turb_dpgb,diff_dp,turb_dpg,ar_gas,flux_g_nf,alphab_gas,rhob_gas,eps_imp_dp,max_iter_dp,iml)
         turb_dpg(:)=MAX(turb_dpg(:),1.d-20)          
      ENDIF  
!
!.....Dissipation at the boundary cells next to the wall
!
      IF(iturb.ne.Kepsilon_real)THEN
         t1=cmu**0.75
         DO i=1,ncell_fluid
            IF(icell_type(i).eq.1) THEN
               t2=turb_keg(i)*sqrt(turb_keg(i))
               turb_dpg(i)=(t1*t2)/(cappa*d_bfc(i))
            ENDIF
         ENDDO
      ELSE
         DO i=1,ncell_fluid
            IF(icell_type(i).eq.1) THEN
               t1=cmug_Real(i)**0.75
               t2=turb_keg(i)*sqrt(turb_keg(i))
               turb_dpg(i)=(t1*t2)/(cappa*d_bfc(i))
            ENDIF
         ENDDO
      ENDIF
!
      END SUBROUTINE turb_ke_calc_gas
