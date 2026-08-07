!
      SUBROUTINE turb_SST_calc_gas
!
!     This routine discretizes and solves the linearized
!     turbulence kinetic energy equations explicitly.
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: np
      USE Zparam       , ONLY: ke_small,ced1,ced2,cmu,cappa,ke_cff,dp_cff,ndim,prt
      USE Zare         , ONLY: ar_gas
      USE Zb_condition , ONLY: turb_kegb,turb_dpgb,alphab_gas,rhob_gas
      USE Zbc_index    , ONLY: icell_type
      USE Zconst1      , ONLY: buoyancy_turb
      USE Zconst2      , ONLY: dt,grav
      USE Zcoord3      , ONLY: volr
      USE Zimplicit    , ONLY: imp_ke_diff,imp_ke_conv,eps_imp_ke,eps_imp_dp,max_iter_ke,max_iter_dp
      USE Zndforce     , ONLY: d_bfc          
      USE Zturb        , ONLY: turb_keg,turb_keg_o,turb_dpg,turb_dpg_o,pro_keg,diff_ke,diff_dp,f_b2,strn_keg
!
      USE Zvec_major    , ONLY: flux_g_nf
!
      IMPLICIT NONE
!
      INTEGER i,ix,iml
!      
      REAL(8) am
!     REAL(8) ced3,a,b,c,d                                      !for buoyancy
      REAL(8) f1,f2,fmu,DlowReynolds,ElowReynolds  !for low Reynolds k-e model    
      REAL(8) cd_kw,arg1,arg1_min,arg1_max1,arg1_max2,arg2                                                 !CYJ k-w
      REAL(8) sigma_k1,sigma_k2,sigma_k3,sigma_w1,sigma_w2,sigma_w3,beta_d                                 !CYJ k-w
      REAL(8) alpha1,alpha2,alpha3,beta1,beta2,beta3                                                       !CYJ k-w
      REAL(8) Gb
!     local arrays
      REAL(8) :: conv_k(ncell_fluid),conv_e(ncell_fluid)
      REAL(8) :: diff_k(ncell_fluid),diff_e(ncell_fluid)
      REAL(8) :: bm(ncell_fluid),cm(ncell_fluid)   
      REAL(8) :: kw_div(ncell_fluid),f_b1(ncell_fluid),dkedx(ncell_fluid,ndim),ddpdx(ncell_fluid,ndim)              !CYJ k-w
      REAL(8) :: drgdx(ncell_fluid,ndim)
!
      ar_gas(:)=cell%alphag(:)*cell%rhog(:) !This should be defined before k-e calculation    
      conv_k(:)=0.d0
      conv_e(:)=0.d0
      diff_k(:)=0.d0
      diff_e(:)=0.d0      
! 
!.....Standard k-e turbulence model
!     
      f1=1.d0
      f2=1.d0
      fmu=1.d0
      DlowReynolds=0.d0
      ElowReynolds=0.d0
!
!.....Additional constants - CYJ k-w
!
      sigma_k1=0.85d0
      sigma_k2=1.d0
      sigma_w1=0.5d0
      sigma_w2=0.856d0
      beta_d=0.09d0
      alpha1=5.d0/9.d0
      alpha2=0.44d0
      beta1=0.075d0
      beta2=0.0828d0       
!
!.....Calculation of Blending function (F1,F2) - CYJ k-w
!
      kw_div(:)=0.d0
      CALL grad_scalar(turb_keg_o,dkedx,ncell_fluid)
      CALL grad_scalar(turb_dpg_o,ddpdx,ncell_fluid)
!      
      DO i=1,ncell_fluid  

         DO ix=1,ndim
            kw_div(i)=kw_div(i)+2.d0*cell%rhog(i)*sigma_w2/turb_dpg_o(i)*dkedx(i,ix)*ddpdx(i,ix)
         ENDDO  
         !IF(turb_dp_o(i).eq.0)kw_div(i)=0.d0
         !kw_div(i)=DMAX1(kw_div(i),1.d-10)
         cd_kw=DMAX1(kw_div(i),1.d-10)
         arg1_min=4.d0*cell%rhog(i)*sigma_w2*turb_keg_o(i)/(cd_kw*d_bfc(i)*d_bfc(i))
         arg1_max1=SQRT(turb_keg_o(i))/(beta_d*turb_dpg_o(i)*d_bfc(i))
         arg1_max2=500.d0*cell%lviscosg(i)/cell%rhog(i)/(d_bfc(i)*d_bfc(i)*turb_dpg_o(i))                  !!!CYJ: cell%lviscosl(i) vs. cell%tviscosl(i)
         arg1=DMIN1(DMAX1(arg1_max1,arg1_max2),arg1_min)
         f_b1(i)=DTANH(arg1**4)   
!         
         arg2=DMAX1(2.d0*arg1_max1,arg1_max2)
         f_b2(i)=DTANH(arg2*arg2)        
      ENDDO      
!
!.....Diffusivity of kinetic energy and dissipation, f1,f2,fmu definition for Low Reynolds k-e model
!
      DO i=1,ncell_fluid
!         sigma_k3=1.d0/(f_b1(i)/sigma_k1+(1.d0-f_b1(i))/sigma_k2)           ! Fluent 12.0
!         sigma_w3=1.d0/(f_b1(i)/sigma_w1+(1.d0-f_b1(i))/sigma_w2)
         sigma_k3=f_b1(i)*sigma_k1+(1.d0-f_b1(i))*sigma_k2                    ! Flunet 14.0
         sigma_w3=f_b1(i)*sigma_w1+(1.d0-f_b1(i))*sigma_w2
         diff_ke(i)=(cell%lviscosg(i)+cell%tviscosg(i)*sigma_k3)*cell%alphag(i)
         diff_dp(i)=(cell%lviscosg(i)+cell%tviscosg(i)*sigma_w3)*cell%alphag(i)
      ENDDO
!
!.....Communicate transport variables for parallel computing
!
      IF(np.gt.1) CALL communicate_1d(turb_keg_o, &
                                      turb_dpg_o, &
                                      diff_ke,    &
                                      diff_dp,    &
                                      ar_gas)
!
!.....Convection and diffusion of k-e & k-w
!
      IF(imp_ke_diff.eq.0) THEN      
         CALL turb_ke_diffusion(turb_keg_o,turb_dpg_o,turb_kegb,turb_dpgb,diff_ke,diff_dp,diff_k,diff_e)
      ENDIF      
      IF(imp_ke_conv.eq.0) THEN      
         CALL turb_ke_convection(turb_keg_o,turb_dpg_o,ar_gas,turb_kegb,turb_dpgb,alphab_gas,rhob_gas,flux_g_nf,conv_k,conv_e)      
      ENDIF      
!
!.....For buoyancy effect
!
      IF(buoyancy_turb.eq.1)THEN
         CALL grad_scalar(cell%rhog,drgdx,ncell_fluid)
      ENDIF
!      
!.....Kinetic energy, k
!         
      bm(:)=0.d0
      cm(:)=0.d0  
      DO i=1,ncell_fluid   
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
         !bm(i)=-conv_k(i)+diff_k(i)+cell%alphag(i)*pro_keg(i)-ar_gas(i)*turb_dpg_o(i)+DlowReynolds
!
!........bm for the Kinetic energy, k
!
         pro_keg(i)=DMIN1(pro_keg(i),10.d0*beta_d*cell%rhog(i)*turb_keg_o(i)*turb_dpg_o(i))
!        
         IF(imp_ke_diff.eq.0.and.imp_ke_conv.eq.0) THEN
            bm(i)=-conv_k(i)+diff_k(i)+cell%alphag(i)*pro_keg(i)-beta_d*ar_gas(i)*turb_keg_o(i)*turb_dpg_o(i)            !+sk_void+sk_kata+DlowReynolds
         ELSE
            bm(i)=-conv_k(i)+diff_k(i)+cell%alphag(i)*pro_keg(i) !+sk_void+sk_kata+DlowReynolds
            cm(i)=-beta_d*ar_gas(i)*turb_dpg_o(i) 
         ENDIF            
!
!........Buoyancy Effect for k
!
         IF(buoyancy_turb.eq.1)THEN
            IF(ndim.eq.2)THEN
               Gb=grav(1)*drgdx(i,1)+grav(2)*drgdx(i,2)
            ELSEIF(ndim.eq.3)THEN
               Gb=grav(1)*drgdx(i,1)+grav(2)*drgdx(i,2)+grav(3)*drgdx(i,3)
            ENDIF
            Gb=-cell%tviscosg(i)*Gb/cell%rhog(i)/prt
            bm(i)=bm(i)+cell%alphag(i)*Gb
         ENDIF
!         
         IF(imp_ke_diff.eq.0.and.imp_ke_conv.eq.0) THEN  
            am=ar_gas(i)/dt
            IF(cell%alphag(i).gt.1.d-8) THEN
               turb_keg(i)=turb_keg_o(i)+bm(i)/am            
            ELSE
               turb_keg(i)=turb_keg_o(i)
            ENDIF
            turb_keg(i)=DMAX1(turb_keg(i),0.d0)            
         ENDIF
!
      ENDDO
!
      IF(imp_ke_diff.eq.1.or.imp_ke_conv.eq.1) THEN       
!
!........1:liquid-phase, 2:gas-phase
!
         iml=1
         CALL imp_diffusion_ke(2,bm,cm,turb_keg_o,turb_kegb,diff_ke,turb_keg,ar_gas,flux_g_nf,alphab_gas,rhob_gas,eps_imp_ke,max_iter_ke,iml)          
         turb_keg(:)=DMAX1(turb_keg,1.d-15)           
      ENDIF            
!
!.....Dissipation, w
!
      bm(:)=0.d0   
      cm(:)=0.d0
      DO i=1,ncell_fluid
!
!         alpha3=alpha3_inf/alpha_star*(1.d0/9.d0+re_t/r_w)/(1.d0+re_t/r_w)            !Fluent 12.0
!         
         alpha3=f_b1(i)*alpha1+(1.d0-f_b1(i))*alpha2                                   !Fluent 14.0
         beta3=f_b1(i)*beta1+(1.d0-f_b1(i))*beta2
!         
         IF(cell%tviscosg(i).eq.0)cell%tviscosg(i)=cell%lviscosg(i)
         IF(imp_ke_diff.eq.0.and.imp_ke_conv.eq.0) THEN
            bm(i)=(1.d0-f_b1(i))*cell%alphag(i)*kw_div(i) +  &
                   cell%alphag(i)*alpha3*cell%rhog(i)*strn_keg(i)*strn_keg(i) -  &
                   beta3*ar_gas(i)*turb_dpg_o(i)*turb_dpg_o(i)   
         ELSE
            bm(i)=(1.d0-f_b1(i))*cell%alphag(i)*kw_div(i) +  &
                   cell%alphag(i)*alpha3*cell%rhog(i)*strn_keg(i)*strn_keg(i)
            cm(i)=-beta3*ar_gas(i)*turb_dpg_o(i)
         ENDIF
!         
         bm(i)=-conv_e(i)+diff_e(i)+bm(i) 
!
         IF(imp_ke_diff.eq.0.and.imp_ke_conv.eq.0) THEN 
            am=ar_gas(i)/dt-cm(i)                
            IF(cell%alphag(i).gt.1.d-8) THEN
               turb_dpg(i)=ar_gas(i)/dt/am*turb_dpg_o(i)+bm(i)/am
            ELSE
               turb_dpg(i)=turb_dpg_o(i)
            ENDIF
            turb_dpg(i)=DMAX1(turb_dpg(i),0.d0)            
         ENDIF
!
      ENDDO 
!
      IF(imp_ke_diff.eq.1.or.imp_ke_conv.eq.1) THEN 
!
!........1:liquid-phase, 2:gas-phase        
!
         iml=2
         CALL imp_diffusion_ke(2,bm,cm,turb_dpg_o,turb_dpgb,diff_dp,turb_dpg,ar_gas,flux_g_nf,alphab_gas,rhob_gas,eps_imp_dp,max_iter_dp,iml)
         turb_dpg(:)=DMAX1(turb_dpg(:),1.d-20)          
      ENDIF  
!
!...Dissipation at the boundary cells next to the wall
!
      DO i=1,ncell_fluid
         IF(icell_type(i).eq.1) THEN
            turb_keg(i)=0.d0
            turb_dpg(i)=10.d0*6.d0*cell%lviscosg(i)/cell%rhog(i)/(beta1*d_bfc(i)*d_bfc(i))
         ENDIF
      ENDDO
!
      END SUBROUTINE turb_SST_calc_gas
!
