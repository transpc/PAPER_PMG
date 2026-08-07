!
      SUBROUTINE turb_SST_calc_liq
!
!     This routine discretizes and solves the linearized
!     turbulence kinetic energy equations explicitly.
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: np
      USE Zimplicit    , ONLY: imp_ke_diff,imp_ke_conv,eps_imp_ke,eps_imp_dp,max_iter_ke,max_iter_dp
      USE Zparam       , ONLY: ke_small,ced1,ced2,cmu,cappa,ke_cff,dp_cff,ndim,prt
      USE Zare         , ONLY: ar_liq
      USE Zb_condition , ONLY: turb_keb,turb_dpb,alphab_liq,rhob_liq
      USE Zbc_index    , ONLY: icell_type
      USE Zconst2      , ONLY: dt
      USE Zcoord3      , ONLY: volr
      USE Zndforce     , ONLY: d_bfc      
      USE Zturb        , ONLY: turb_ke,turb_ke_o,turb_dp,turb_dp_o,pro_ke,diff_ke,diff_dp,f_b2,strn_ke
!
      USE Zvec_major    , ONLY: flux_l_nf
!
      IMPLICIT NONE
!
      INTEGER i,ix,iml
!
      REAL(8) sk_void,se_void,cdk
      REAL(8) gravity
      REAL(8) am
!     REAL(8) ut,un,gravm,beta                                  !for buoyancy
!     REAL(8) ced3,a,b,c,d                                      !for buoyancy
      REAL(8) f1,f2,fmu,DlowReynolds,ElowReynolds  !for low-Reynolds-number k-e model  
!     REAL(8) hd,rel,fric,velo,velo2,velo3,velo4,yplim,pro_porous(ncell_fluid),kinf,coeff_eps !LSJ 161122 porous
!     REAL(8) KK,bb,cd                                                                        !LSJ 161122 porous  
      !REAL(8) kw_div(ncell_fp),f_b1(ncell_fp),dkedx(ndim,ncell_fp),ddpdx(ndim,ncell_fp)                    !CYJ k-w
      REAL(8) cd_kw,arg1,arg1_min,arg1_max1,arg1_max2,arg2                                                 !CYJ k-w
      REAL(8) sigma_k1,sigma_k2,sigma_k3,sigma_w1,sigma_w2,sigma_w3,beta_d                                 !CYJ k-w
      REAL(8) alpha1,alpha2,alpha3,beta1,beta2,beta3                                                       !CYJ k-w
!     REAL(8) alpha1_inf,alpha2_inf,alpha3_inf,alpha_star,re_t,r_k,r_w                                     !CYJ k-w
!     local arrays
      REAL(8) :: diff_k(ncell_fluid),diff_e(ncell_fluid),conv_k(ncell_fluid),conv_e(ncell_fluid)
!     REAL(8) :: temp(ncell_fp),dtldx(ncell_fp,ndim)         
      REAL(8) :: bm(ncell_fluid),cm(ncell_fluid)
      REAL(8) :: kw_div(ncell_fluid),f_b1(ncell_fluid),dkedx(ncell_fluid,ndim),ddpdx(ncell_fluid,ndim)
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
!.....Additional constants - CYJ k-w
!
      sigma_k1=0.85d0
      sigma_k2=1.0d0
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
      CALL grad_scalar(turb_ke_o,dkedx,ncell_fluid)
      CALL grad_scalar(turb_dp_o,ddpdx,ncell_fluid)      
!      
      DO i=1,ncell_fluid  
         DO ix=1,ndim
            kw_div(i)=kw_div(i)+2.d0*cell%rhol(i)*sigma_w2/turb_dp_o(i)*dkedx(i,ix)*ddpdx(i,ix)
         ENDDO  
         cd_kw=DMAX1(kw_div(i),1.d-10)
         arg1_min=4.d0*cell%rhol(i)*sigma_w2*turb_ke_o(i)/(cd_kw*d_bfc(i)*d_bfc(i))
         arg1_max1=SQRT(turb_ke_o(i))/(beta_d*turb_dp_o(i)*d_bfc(i))
         arg1_max2=500.d0*cell%lviscosl(i)/cell%rhol(i)/(d_bfc(i)*d_bfc(i)*turb_dp_o(i))                  !!!CYJ: cell%lviscosl(i) vs. cell%tviscosl(i)
         arg1=DMIN1(DMAX1(arg1_max1,arg1_max2),arg1_min)
         f_b1(i)=DTANH(arg1**4)   
!         
         arg2=DMAX1(2.d0*arg1_max1,arg1_max2)
         f_b2(i)=DTANH(arg2*arg2)        
      ENDDO
!
!.....Diffusivity of kinetic energy and dissipation
!
      DO i=1,ncell_fluid
!         sigma_k3=1.0d0/(f_b1(i)/sigma_k1+(1.d0-f_b1(i))/sigma_k2)           ! Fluent 12.0
!         sigma_w3=1.0d0/(f_b1(i)/sigma_w1+(1.d0-f_b1(i))/sigma_w2)
         sigma_k3=f_b1(i)*sigma_k1+(1.d0-f_b1(i))*sigma_k2                    ! Flunet 14.0
         sigma_w3=f_b1(i)*sigma_w1+(1.d0-f_b1(i))*sigma_w2
         diff_ke(i)=(cell%lviscosl(i)+cell%tviscosl(i)*sigma_k3)*cell%alphal(i)
         diff_dp(i)=(cell%lviscosl(i)+cell%tviscosl(i)*sigma_w3)*cell%alphal(i)
      ENDDO
!
!.....Communicate transport variables for parallel computing
!
      IF(np.gt.1) CALL communicate_1d(turb_ke_o, &
                                      turb_dp_o, &
                                      diff_ke,   &
                                      diff_dp,   &
                                      ar_liq)
!
!.....Convection and diffusion of k-e & k-w
!
      IF(imp_ke_diff.eq.0) THEN  
         CALL turb_ke_diffusion(turb_ke_o,turb_dp_o,turb_keb,turb_dpb,diff_ke,diff_dp,diff_k,diff_e)
      ENDIF
!      
      IF(imp_ke_conv.eq.0) THEN  
         CALL turb_ke_convection(turb_ke_o,turb_dp_o,ar_liq,turb_keb,turb_dpb,alphab_liq,rhob_liq,flux_l_nf,conv_k,conv_e)   
      ENDIF     
!
!.....For buoyancy effect
!
!      temp(:)=cell%tl_o(:)
!      CALL grad_scalar(temp,dtldx)
!
!.....Kinetic energy, k
!
      bm(:)=0.0d0    
      cm(:)=0.0d0 
      DO i=1,ncell_fluid
!!
!!........Source term by bubble motion (Lahey,2005)
!!
!         IF(turboil.eq.1) THEN 
!            sk_void=0.25d0*(1.0d0+CDk**(4.d0/3.d0))*cell%alphag(i)*vrel_o(i)*vrel_o(i)*vrel_o(i)/cell%d1(i)
!            sk_void=cell%alphal(i)*sk_void
!!
!!...........Kinetic energy source term by wall boiling (Kataoka and Serizawa,1997)
!!
!            IF(icell_type(i).eq.1) THEN
!               hfg=cell%hgsat(i)-cell%hlsat(i)
!               sk_kata=8.0d0/9.d0*gravity*(cell%rhol(i)-cell%rhog(i))/cell%rhol(i)*qecell(i)/(cell%rhog(i)*hfg)
!            ENDIF
!         ENDIF
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
         pro_ke(i)=DMIN1(pro_ke(i),10.d0*beta_d*cell%rhol(i)*turb_ke_o(i)*turb_dp_o(i))
!        
         IF(imp_ke_diff.eq.0.and.imp_ke_conv.eq.0) THEN
            bm(i)=-conv_k(i)+diff_k(i)+cell%alphal(i)*pro_ke(i)-beta_d*ar_liq(i)*turb_ke_o(i)*turb_dp_o(i)            !+sk_void+sk_kata+DlowReynolds
         ELSE
            bm(i)=-conv_k(i)+diff_k(i)+cell%alphal(i)*pro_ke(i) !+sk_void+sk_kata+DlowReynolds
            cm(i)=-beta_d*ar_liq(i)*turb_dp_o(i) 
         ENDIF                 
         
!!
!!........Buoyancy Effect for k
!!  
!         IF(buoyancy_turb.eq.1) THEN
!            DO idim=1,ndim
!               beta=cell%betal(i)
!               bm(i)=bm(i)+cell%alphal(i)*beta*grav(idim)*cell%tviscosl(i)/prt*dtldx(i,idim)
!            ENDDO
!         ENDIF
!         
         IF(imp_ke_diff.eq.0.and.imp_ke_conv.eq.0) THEN
            am=ar_liq(i)/dt          
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
!----------------------------------------------------------------------------------------------------------------------------------------------      
!
!.....Dissipation, w
!
      bm(:)=0.0d0   
      cm(:)=0.0d0    
      DO i=1,ncell_fluid
         sk_void=0.0d0
         se_void=0.0d0
!
!         alpha3=alpha3_inf/alpha_star*(1.d0/9.d0+re_t/r_w)/(1.d0+re_t/r_w)            !Fluent 12.0
!         
         alpha3=f_b1(i)*alpha1+(1.d0-f_b1(i))*alpha2                                   !Fluent 14.0
         beta3=f_b1(i)*beta1+(1.d0-f_b1(i))*beta2
!
         IF(cell%tviscosl(i).eq.0)cell%tviscosl(i)=cell%lviscosl(i)
         IF(imp_ke_diff.eq.0.and.imp_ke_conv.eq.0) THEN
            bm(i)=(1.d0-f_b1(i))*cell%alphal(i)*kw_div(i) +  &
                   cell%alphal(i)*alpha3*cell%rhol(i)*strn_ke(i)*strn_ke(i) -  &
                   beta3*ar_liq(i)*turb_dp_o(i)*turb_dp_o(i)   
         ELSE
            bm(i)=(1.d0-f_b1(i))*cell%alphal(i)*kw_div(i)+  &
                   cell%alphal(i)*alpha3*cell%rhol(i)*strn_ke(i)*strn_ke(i)
            cm(i)=-beta3*ar_liq(i)*turb_dp_o(i)    !+cell%alphal(i)*alpha3/turb_ke(i)*pro_ke(i) 
         ENDIF
!         
         bm(i)=-conv_e(i)+diff_e(i)+bm(i)        
!!
!!........Buyoancy Effect for w
!!         
!         IF(buoyancy_turb.eq.1) THEN
!            ut=0.d0
!            un=0.d0
!            gravm=0.d0
!            DO idim=1,ndim
!               ut=ut+grav(idim)*vl_o(idim,i)
!               un=un+vl_o(idim,i)*vl_o(idim,i)
!               gravm=gravm+grav(idim)*grav(idim)
!            ENDDO   
!            gravm=gravm**0.5
!            ut=ut/gravm
!            un=(un-ut*ut)**0.5
!            un=DMAX1(DABS(un),0.000001)*DSIGN(1.d0,un)
!            a=ut/un
!            b=DMAX1(0.d0,DMIN1(20.d0,DABS(a)))
!            c=DEXP(b)
!            d=DEXP(-b)
!            ced3=(c-d)/(c+d)
!!
!            beta=cell%betal(i)
!            DO idim=1,ndim
!               pkb=beta*grav(idim)*cell%tviscosl(i)/prt*dtldx(i,idim)
!               cm(i)=cm(i)+cell%alphal(i)/turb_ke(i)*ced1*ced3*   ((alpha1+1.d0)*DMAX1(pkb,0.d0)-pkb)
!            ENDDO
!        ENDIF   
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
      iml=2
      IF(imp_ke_diff.eq.1.or.imp_ke_conv.eq.1) THEN
         CALL imp_diffusion_ke(1,bm,cm,turb_dp_o,turb_dpb,diff_dp,turb_dp,ar_liq,flux_l_nf,alphab_liq,rhob_liq,eps_imp_dp,max_iter_dp,iml)
         turb_dp(:)=DMAX1(turb_dp(:),1.d-15)       !!!CYJ Min=0.0d0 -> 1d-15
      ENDIF       
!
!.....Dissipation at the boundary cells next to the wall
!
      DO i=1,ncell_fluid
         IF(icell_type(i).eq.1)THEN   !!! CYJ k-w
!              turb_dp(i)=6.d0*cell%lviscosl(i)/cell%rhol(i)/(beta1*d_bfc(i)*d_bfc(i))
!              dp_w=cmu**(-0.25)*utau(i)/(cappa*d_bfc(i))
!              turb_dp(i)=SQRT(turb_dp(i)*turb_dp(i)+dp_w*dp_w)
            turb_ke(i)=0.0d0
            turb_dp(i)=10.0d0*6.d0*cell%lviscosl(i)/cell%rhol(i)/(beta1*d_bfc(i)*d_bfc(i))
         ENDIF
!!         IF(icell_type(i).eq.1) turb_dp(i)=SQRT(turb_ke(i)/(cmu**(-0.25)*cappa*d_bfc(i))
!         IF(icell_type(i).eq.1)turb_dp(i)=6.d0*cell%lviscosl(i)/cell%rhol(i)/(beta1*d_bfc(i)*d_bfc(i))
!!          IF(icell_type(i).eq.1)turb_dp(i)=turb_ke(i)**0.5d0/(0.41d0*0.09d0**0.25d0*d_bfc(i))
      ENDDO
!
      END SUBROUTINE turb_SST_calc_liq
!
