      SUBROUTINE rv_int_fric_model(n,vfgl_i,vfgd_i,ia_i)
!
!.....interfacial drag coefficients for n cells.
!
      USE VOL_DATA       , ONLY: cell
      USE Zzone          , ONLY: ncell_fluid     
      USE Zcore          , ONLY: np,myrank
      USE Zparam         , ONLY: pi,ndim
      USE Zvector        , ONLY: vg_o,vl_o,vrel_o
      USE Zconst2        , ONLY: hydraulicd
      USE Zcoord1        , ONLY: xloc
      USE Zcoord2        , ONLY: cell_leng      
      USE Znum_cell      , ONLY: i_neigh,neigh
      USE Zmodel         , ONLY: drift_c0,drift_c1
      USE STM_TBL_cupid  , ONLY: pcrit
      USE Zrv_int_friction , ONLY: iam,vr_drift,iregime,drift_model,drift_model,drag_model,            &
                                   ia_bubbly,ia_slug_tb,ia_slug_sb,ia_churn,ia_annular_ann,            &
                                   ia_annular_drp,ia_MPR,ia_invann_ann,ia_invann_sb,ia_invchn,         &
                                   ia_invslg_ann,ia_invslg_drp,ia_mist,ia_MPO,ia_VST,vfgl_liquid,      &
                                   vfgl_gas,vfgl_drift,vfgl_bubbly,vfgl_slug,vfgl_churn,               &
                                   vfgl_annular,vfgl_MPR,vfgl_invann,vfgl_invchn,vfgl_invslg,          &
                                   vfgl_mist,vfgl_MPO,vfgl_VST,alpha_am,alpha_min,alpha_max,           &
                                   alphag_lower,alphag_upper,alphabs,alphacd,alphade,alphasa,          &
                                   alphags,alphaam,alphamin,alphamax,dbub_min,dbub_max,                &
                                   ddrp_min,ddrp_max,wf_dry,wf_VST,pres,xnn,fdir,dh,dht,ddrp,lcell,    &
                                   drop_max,rhog,rhol,rhogoverl,rholoverg,rhol_g,mug,mul,sigma,        &
                                   red_min,alphag,alphal,alphag,vr,vrd,vl,vg,gravity,vm,vtb,           &
                                   upward,sideward,wf_liquid,wf_gas,wf_bubbly,wf_slug,wf_drift,        &
                                   wf_churn,wf_annular,wf_MPR,wf_invann,wf_invchn,wf_invslg,           &
                                   wf_mist,wf_MPO,jg,jl,liquid,gas,bubbly,slug,drift,churn,annular,    &
                                   MPR,invann,invchn,invslg,mist,MPO,VST,dcon,web,xa,wfactor,mum,      &
                                   alphadrp,red,cd,dbub,sigma_gr,sigmagr,vgjs,alphagupper,alphaglower, &
                                   alphagsb,fmixlevel,ia_VST_sb,ia_VST_surf,wf_preCHF,wf_postDRY,      &
                                   wf_unVST,wf_sum,ia_preCHF,ia_postDRY,ia_unVST,vfgl_preCHF ,         &
                                   vfgl_postDRY,vfgl_unVST,mdrag_preCHF,mdrag_postDRY,mdrag_unVST,     &
                                   mdrag_whole,c0,c1,                                                  &
                                   i   
      USE Zrv_model      , ONLY: rv_fric_i      
      USE Zwall_HTC      , ONLY: reflod
      USE Zporous        , ONLY: s_subchannel_mixing
      USE Zrv_int_friction , ONLY: s_rv_int_fric
      USE Zio_unit         , ONLY: unit_log
!      
      IMPLICIT NONE
!
!.....Input
      INTEGER,INTENT(IN)  :: n           ! cell number.
!.....Output
      REAL(8),INTENT(OUT) :: vfgl_i(n)   ! interfacial drag coefficient for gas to liquid.
      REAL(8),INTENT(OUT) :: vfgd_i(n)   ! interfacial drag coefficient for gas to droplet.
      REAL(8),INTENT(OUT) :: ia_i(n)     ! interfacial area.
!.....Local variables
      INTEGER :: j,k,j0,mdrag
      INTEGER :: iacrit,face_alpha_option
      LOGICAL,SAVE :: initial=.TRUE.
      LOGICAL,SAVE :: mars_coding_only              
!     REAL(8) :: rv_int_fric_weight,alphawid
      REAL(8) :: rv_int_fric_weight
      REAL(8) :: vfg,vfg2,preduc,drmin    !MPR
      REAL(8) :: alpdrp,voidx,surfa,rey,fic      
!
      IF(s_subchannel_mixing.eq.'evvd')s_rv_int_fric='drag_coeff'
!
      IF(rv_fric_i.eq.1) THEN
         mars_coding_only=.FALSE.         !MARS manual model
      ELSE
         mars_coding_only=.TRUE.          !MARS code model
      ENDIF       
!
      face_alpha_option=0                 !rv face-upwind void fraction option (0=not use, 1=use)
!
      iam=0    !0=calculate ia here, 1=obtain ia from other subroutine
      iacrit=1 !0=set alphag criterion, 1=obtain alphag criterion from other subroutine
!.....pcrit is constant read in stread.f90 why change?
!     pcrit=22.4d6
!
      IF(initial)THEN
         initial=.FALSE.
         vfgl_i(:)=0.d0
         vfgd_i(:)=0.d0
         IF(myrank.eq.0)THEN
            WRITE(*,"(11x,a,1x,a)")'--Interfacial drag model in RV is ',s_rv_int_fric
!            READ(*,*)mdrag
         ENDIF            
         IF(np.gt.1) CALL broadcast_i1(mdrag)         
         ALLOCATE(vr_drift(ndim))
         ALLOCATE(iregime(14,n))
         IF(s_subchannel_mixing.eq.'evvd') THEN
            IF(myrank.eq.0) PRINT*,'Drag coefficient model is selected for EVVD model'
         ENDIF
      ENDIF
!
!........Drag model
!
      IF(s_rv_int_fric.eq.'constant')THEN
         vfgl_i(:)=1.d4
         vfgd_i(:)=1.d4
         RETURN
      ELSEIF(s_rv_int_fric.eq.'drag_coeff')THEN
         drift_model=0
         drag_model=1
      ELSEIF(s_rv_int_fric.eq.'drift_flux')THEN
         drift_model=1
         drag_model=1
      ELSE
         IF(myrank.eq.0)WRITE(*,"(11x,a)")'Check s_rv_int_fric= ',s_rv_int_fric
         IF(myrank.eq.0)WRITE(unit_log,"(11x,a)")'Check s_rv_int_fric= ',s_rv_int_fric
         STOP
      ENDIF
!
!.....face-upwind void fraction for interfacial friction model
!
      IF(face_alpha_option.eq.1) CALL face_alpha      
! 
      DO i=1,n
!
!........Initialize
!
         ia_bubbly     =0.d0
         ia_slug_tb    =0.d0
         ia_slug_sb    =0.d0
         ia_churn      =0.d0 
         ia_annular_ann=0.d0
         ia_annular_drp=0.d0
         ia_MPR        =0.d0         
         ia_invann_ann =0.d0
         ia_invann_sb  =0.d0
         ia_invchn     =0.d0 
         ia_invslg_ann =0.d0
         ia_invslg_drp =0.d0
         ia_mist       =0.d0
         ia_MPO        =0.d0
         ia_VST        =0.d0
         vfgl_liquid   =0.d0 
         vfgl_gas      =0.d0
         vfgl_drift    =0.d0  
         vfgl_bubbly   =0.d0
         vfgl_slug     =0.d0
         vfgl_churn    =0.d0 
         vfgl_annular  =0.d0
         vfgl_MPR      =0.d0         
         vfgl_invann   =0.d0
         vfgl_invchn   =0.d0
         vfgl_invslg   =0.d0 
         vfgl_mist     =0.d0
         vfgl_MPO      =0.d0
         vfgl_VST      =0.d0       
!
!........temperary assignment 
!
         IF(iacrit.eq.0)then  
             cell%alpha_bs(i)=0.3d0 !bubbly-slug  invann-invchn
             cell%alpha_cd(i)=0.4d0 !             invchn-invslg
             cell%alpha_de(i)=0.6d0 !slug-churn,  discontinuity occurs when drift flux model is used. 
             cell%alpha_sa(i)=0.8d0 !slug-annular invslg-mist
             cell%alpha_gs(i)=0.3d0 
         ENDIF
!         alpha_am=0.99d0 !next
         alpha_am=0.9999d0 !next
!         
         alpha_min=1.d-8
         alpha_max=1.d0
         
         IF(.true.)THEN !next
            j0=i_neigh(i)-1
            DO j=i_neigh(i),i_neigh(i+1)-1
               k=neigh(j)
               IF(k.le.0.or.k.gt.ncell_fluid)CYCLE
               IF(xloc(i,3).gt.xloc(k,3))then
                  alphag_lower=cell%alphag(k)
               ELSE
                  alphag_upper=cell%alphag(k)
               ENDIF
            ENDDO   
         ELSE
            alphag_upper=0.9d0 !for VST
            alphag_lower=0.2d0 !for VST     
         ENDIF      
!      
         cell%ddrp(i)=0.0001d0
!         cell%length(i)=0.1d0 
!         cell%length(i)=0.5d0*(cell_leng(i,1)+cell_leng(i,2))  !mod. by LSJ
         cell%length(i)=cell_leng(i,ndim)  ! cell height
!         cell%fdir(i)=1     !Commment-out by LSJ since cell%fdir(:) is defined in rv_flow_regime.f90.
!
!........common variables
!
!........void fraction criteria
!
         alphabs =cell%alpha_bs(i) !bubbly-slug,  iannular-ichurn
         alphacd =cell%alpha_cd(i) !              ichurn-islug
         alphade =cell%alpha_de(i) !slug-churn 
         alphasa =cell%alpha_sa(i) !churn-annular islug-mist
         alphags =cell%alpha_gs(i) !next, small bubble fraction? !alpha_gs(i)=MIN( alphag, alpha_bs(i) ) * beta
         alphaam =alpha_am
         alphamin=alpha_min
         alphamax=alpha_max
!      
         dbub_min=0.001d0   !next
         dbub_max=0.01d0    !next
         ddrp_min=0.00001d0 !next-critical 
         ddrp_max=0.01d0    !next  
         IF(mars_coding_only) ddrp_max=0.5d0 
!
!........weighting factor for CHF and VST
!
         wf_dry=cell%wf_dry(i)
         wf_VST=cell%wf_vst(i) 
         pres  =cell%p(i) 
         xnn   =cell%quala(i)     
!      
!........flow direction for drift flux model, next      
!
         fdir=cell%fdir(i)      
!      
!........geometrical factor
!
         dh      =hydraulicd(i) 
         dht     =dh !physical cell diamter, next
         ddrp    =cell%ddrp(i)
         lcell   =cell%length(i)
         drop_max=dh       ! default
!      
!........properties
!
         rhog     =cell%rhog(i)
         rhol     =cell%rhol(i)
         rhogoverl=rhog/rhol
         rholoverg=1.d0/rhogoverl
         rhol_g   =rhol-rhog      
         rhol_g   =MAX(1.d-5,rhol_g)
!         mug     =cell%eviscosg(i)
!         mul     =cell%eviscosl(i)
         mug      =cell%lviscosg(i)     !LSJ modification
         mul      =cell%lviscosl(i)     !LSJ modification  
         sigma    =cell%sigma(i)
         red_min  =1.d0 ! yjm
!
!........T/H factor  
!
         IF(face_alpha_option.eq.1) THEN
            alphag=MIN(1.d0,MAX(cell%alphagf(i),1.d-8))
            alphal=MIN(1.d0,MAX(cell%alphalf(i),1.d-8))          
         ELSE
            alphag=cell%alphag(i)
            alphag=MIN(1.d0,MAX(alphag,1.d-8))
            alphal=MIN(1.d0,MAX(1.d0-alphag,1.d-8))
            alphag=1.d0-alphal           
         ENDIF   
!         
         vr=MAX(1.d-8,vrel_o(i))
         vrd=0.d0 !drift
         vl=vl_o(i,1)*vl_o(i,1)+vl_o(i,2)*vl_o(i,2)
         IF(ndim.eq.3)vl=vl+vl_o(i,3)*vl_o(i,3)
         vl=SQRT(vl)
         vg=vg_o(i,1)*vg_o(i,1)+vg_o(i,2)*vg_o(i,2)
         IF(ndim.eq.3)vg=vg+vg_o(i,3)*vg_o(i,3)
         vg=SQRT(vg)
         gravity=9.81d0  !next
         vm=(vl*alphal*rhol+vg*alphag*rhog)/(alphal*rhol+alphag*rhog)
         vtb=0.35d0*SQRT(gravity*dh*rhol_g/rhol)
!
         jg=alphag*vg !next
         jl=alphal*vl !next
!
!........flow direction with cell%fdir(i), upward, sideward, downward, countcurrent, next
!
         IF(fdir.eq.1)THEN
            upward=1
            sideward=0
         ELSEIF(fdir.eq.2)THEN
            upward=0
            sideward=1      
         ELSE
            upward=0
            sideward=0         
         ENDIF   
!   
!........weighting factor for each regime using void fraction
!
!        alphawid=0.01d0
         wf_liquid =rv_int_fric_weight(alphag,  -1.d0,  0.d0)
         wf_gas    =rv_int_fric_weight(alphag,   1.d0,  2.d0)
!         
         wf_bubbly =rv_int_fric_weight(alphag,   0.d0, alphabs)
         wf_slug   =rv_int_fric_weight(alphag,alphabs, alphade)
         wf_drift  =rv_int_fric_weight(alphag,  0.d0,alphade)
         wf_churn  =rv_int_fric_weight(alphag,alphade,alphasa)
         wf_annular=rv_int_fric_weight(alphag,alphasa,alphaam)
         wf_MPR    =rv_int_fric_weight(alphag,alphaam,  1.d0)
!         
         wf_invann =rv_int_fric_weight(alphag,  0.d0,alphabs)
         wf_invchn =rv_int_fric_weight(alphag,alphabs,alphacd)
         wf_invslg =rv_int_fric_weight(alphag,alphacd,alphasa)
         wf_mist   =rv_int_fric_weight(alphag,alphasa,alphaam)
         wf_MPO    =rv_int_fric_weight(alphag,alphaam,  1.d0)
!
!........weighting factor is 0 at inner and outer edge in churn and invchn
!      
         IF(wf_churn.gt.0.d0.and.wf_churn.lt.0.5d0)THEN
            wf_churn=0.d0
            IF(wf_slug   .gt.0.5d0 .and. wf_slug   .lt.1.d0) wf_slug   =1.d0
            IF(wf_drift  .gt.0.5d0 .and. wf_drift  .lt.1.d0) wf_drift  =1.d0
            IF(wf_annular.gt.0.5d0 .and. wf_annular.lt.1.d0) wf_annular=1.d0
         ENDIF   
         IF(wf_churn.ge.0.5d0 .and. wf_churn.lt.1.d0)THEN
            wf_churn=1.d0
            IF(wf_slug   .gt.0.d0 .and. wf_slug   .le.0.5d0) wf_slug   =0.d0
            IF(wf_drift  .gt.0.d0 .and. wf_drift  .le.0.5d0) wf_drift  =0.d0
            IF(wf_annular.gt.0.d0 .and. wf_annular.le.0.5d0) wf_annular=0.d0
         ENDIF                 
         IF(wf_invchn.gt.0.d0.and.wf_invchn.lt.0.5d0)THEN
            wf_invchn=0.d0
            IF(wf_invslg.gt.0.5d0 .and. wf_invslg.lt.1.d0) wf_invslg=1.d0
            IF(wf_invann.gt.0.5d0 .and. wf_invann.lt.1.d0) wf_invann=1.d0
         ENDIF 
         IF(wf_invchn.ge.0.5d0.and.wf_invchn.lt.1.d0)THEN
            wf_invchn=1.d0
            IF(wf_invslg.gt.0.d0 .and. wf_invslg.le.0.5d0) wf_invslg=0.d0
            IF(wf_invann.gt.0.d0 .and. wf_invann.le.0.5d0) wf_invann=0.d0
         ENDIF       
!
!........adjust weighting factor according to drag model
!
         IF(drift_model.eq.1)THEN
            wf_bubbly=0.d0
            wf_slug  =0.d0
         ELSE
            wf_drift =0.d0
         ENDIF      
!      
!........adjust weighting factor according to dTsat (CHF curve)
!
         IF(wf_dry.eq.0.d0)THEN !Pre-CHF
            wf_invann=0.d0
            wf_invchn=0.d0
            wf_invslg=0.d0
            wf_mist  =0.d0
            wf_MPO   =0.d0 
         ELSEIF(wf_dry.eq.1.d0)THEN !Post-DRY
            wf_bubbly =0.d0
            wf_slug   =0.d0
            wf_drift  =0.d0
            wf_churn  =0.d0
            wf_annular=0.d0
            wf_MPR    =0.d0    
         ELSE !Transition between Pre-CHF and Post-DRY
         ENDIF  
!  
!........adjust weighting factor according to stratification
!
!         IF(vm.le.0.5d0*vtb)THEN
!            wf_VST=1.d0
!         ELSEIF(vm.ge.vtb)THEN
!            wf_VST=0.d0
!         ELSE
!            !wf_VST=(0.d0-1.d0)/(vtb-0.5d0*vtb)*(vm-vtb)+0.d0
!            wf_VST=2.d0/vtb*(vtb-vm)
!         ENDIF        
!   
         IF(wf_VST.eq.1.d0)THEN !completely stratified
            wf_bubbly=0.d0
            wf_slug   =0.d0
            wf_drift  =0.d0
            wf_churn  =0.d0
            wf_annular=0.d0
            wf_MPR    =0.d0  
            wf_invann =0.d0
            wf_invchn =0.d0
            wf_invslg =0.d0
            wf_mist   =0.d0
            wf_MPO    =0.d0          
         ENDIF
! 
!........Select regime for friction to be calculated, with weighting factor
!
!        intialize single-phase regimes
         liquid=0
         gas=0
!         
!        intialize normal 2-phase regimes     
         bubbly =0
         slug   =0
         drift  =0
         churn  =0
         annular=0
         MPR    =0
!
!        intialize inverted 2-phase regimes      
         invann=0
         invchn=0
         invslg=0
         mist  =0
         MPO   =0
!
!        vertically stratified regime
         VST=0      
!         
!        turn on each drag model with a weighting factor      
         IF(wf_liquid .gt.0.d0 .and. wf_liquid .le.1.d0) liquid =1
         IF(wf_gas    .gt.0.d0 .and. wf_gas    .le.1.d0) gas    =1
         IF(wf_bubbly .gt.0.d0 .and. wf_bubbly .le.1.d0) bubbly =1
         IF(wf_slug   .gt.0.d0 .and. wf_slug   .le.1.d0) slug   =1
         IF(wf_drift  .gt.0.d0 .and. wf_drift  .le.1.d0) drift  =1
         IF(wf_churn  .gt.0.d0 .and. wf_churn  .le.1.d0) churn  =1
         IF(wf_annular.gt.0.d0 .and. wf_annular.le.1.d0) annular=1
         IF(wf_MPR    .gt.0.d0 .and. wf_MPR    .le.1.d0) MPR    =1
         IF(wf_invann .gt.0.d0 .and. wf_invann .le.1.d0) invann =1
         IF(wf_invchn .gt.0.d0 .and. wf_invchn .le.1.d0) invchn =1
         IF(wf_invslg .gt.0.d0 .and. wf_invslg .le.1.d0) invslg =1
         IF(wf_mist   .gt.0.d0 .and. wf_mist   .le.1.d0) mist   =1
         IF(wf_MPO    .gt.0.d0 .and. wf_mist   .le.1.d0) MPO    =1
         IF(wf_VST    .gt.0.d0 .and. wf_VST    .le.1.d0) VST    =1
!
         iregime( 1,i)=liquid
         iregime( 2,i)=gas
         iregime( 3,i)=bubbly
         iregime( 4,i)=slug
         iregime( 5,i)=drift
         iregime( 6,i)=churn
         iregime( 7,i)=annular
         iregime( 8,i)=MPR
         iregime( 9,i)=invann
         iregime(10,i)=invchn
         iregime(11,i)=invslg
         iregime(12,i)=mist
         iregime(13,i)=MPO
         iregime(14,i)=VST
!
         IF(bubbly.eq.1.or.slug.eq.1.or.invann.eq.1 &                     !bubbles
             .or.MPR.eq.1.or.annular.eq.1 &                  !preCHF, droplets 
             .or.MPO.eq.1.or.mist.eq.1.or.invslg.eq.1)THEN   !postCHF,droplets
            IF(reflod)THEN
                dcon(1)=5.d-3
!                dcon(2)=1.5d-3
                dcon(2)=2.5d-4  ! MARS source interphaseDrag Line 1962
                dcon(3)=2.d-4
                web(1)=5.d0
                web(2)=1.5d0
!                web(3)=1.5d0
!
!                use TRACE 5.0 We number of reflood We (based on dmax) =12.0
!                sauter mean diam = 1/3 * dmax
!                We = 3/1*12.0 = 4.0
!
                web(3)=4.d0  ! yjm ori: web(3)=4.d0
            ELSE
                dcon(1)=5.d-3
                dcon(2)=2.5d-3
                dcon(3)=2.d-4
                web(1)=5.d0
                web(2)=1.5d0
                web(3)=6.d0          
            ENDIF    
            xa(:)=MAX(web(:)*sigma,1.d-10)
         ENDIF        
!
!........interpolation drag coefficient model for churn between slug and annular
!
         IF(churn.eq.1)THEN
            CALL bubblyslug_annular_for_churn(1)
            wfactor=(alphasa-alphag)/(alphasa-alphade)
            wfactor=MIN(1.d0,MAX(0.d0,wfactor))  
            IF(drift_model.eq.1)THEN    
               vfgl_churn=vfgl_drift*wfactor+vfgl_annular*(1.d0-wfactor)
!               mdrag_annular=vfgl_annular*vr
!               mdrag_churn=mdrag_drift*wfactor+mdrag_annular*(1.d0-wfactor)
!               vfgl_churn=mdrag_churn/vrd !next
            ELSE
               vfgl_churn=vfgl_slug*wfactor+vfgl_annular*(1.d0-wfactor)
               ia_churn=(ia_slug_tb+ia_slug_sb)*wfactor+(ia_annular_ann+ia_annular_drp)*(1.d0-wfactor)
            ENDIF   
            vfgl_slug   =0.d0
            vfgl_annular=0.d0
            vfgl_drift  =0.d0
         ENDIF 
!
!........interpolation drag coefficient model for invchurn between slug and annular
!
         IF(invchn.eq.1)THEN
            CALL invann_invslg_for_invchn(1)
            wfactor=(alphacd-alphag)/(alphacd-alphabs)
            wfactor=MIN(1.d0,MAX(0.d0,wfactor))          
            IF(.true.)THEN
               vfgl_invchn=vfgl_invann**wfactor*vfgl_invslg**(1.d0-wfactor)
               ia_invchn=(ia_invann_ann+ia_invann_sb)*wfactor+(ia_invslg_ann+ia_invslg_drp)*(1.d0-wfactor)
            ELSE
               vfgl_invchn=vfgl_invann*wfactor+vfgl_invslg*(1.d0-wfactor)
               ia_invchn=(ia_invann_ann+ia_invann_sb)*wfactor+(ia_invslg_ann+ia_invslg_drp)*(1.d0-wfactor)
            ENDIF            
            vfgl_invann=0.d0
            vfgl_invslg=0.d0
         ENDIF                                                                       
!
!........drag coefficients for pure gas and pure liquid
!
         IF(liquid.eq.1)THEN
!            vfgl_i(i)=1.d4
            vfgl_i(i)=1.d1
!            
            IF(mars_coding_only) vfgl_i(i)=1.d1   !MARS code             
            vfgl_liquid=vfgl_i(i)
         ENDIF
!         
         IF(gas.eq.1)THEN
!            vfgl_i(i)=1.d4
            vfgl_i(i)=1.d1
!            
            IF(mars_coding_only) vfgl_i(i)=1.d1   !MARS code             
            vfgl_gas=vfgl_i(i)
         ENDIF
!
!........drag coefficient model for bubbly,slug,annular, 
!........and drift flux model for bubbly&slug
!      
         CALL bubblyslug_annular_for_churn(0)
!
!........drag coefficient model for MPO
!
         IF(MPR.eq.1)THEN
            IF(mars_coding_only) THEN  
!               velo_g=DABS(vg_o(i,ndim))
!               velo_l=DABS(vl_o(i,ndim))
!               
               alpdrp=MAX(alphal,1.d-4)
!               vfg=(velo_g-velo_l)
               vfg=MAX(ABS(vr),0.001d0)  ! yjm MARS
!
!              BubbleDropDrag: fic            
               vfg2=MAX(vfg*vfg,xa(2)/(rhog*MIN(dcon(2)*alpdrp**0.333333333d0,dh)))
!
!               rdiam=MIN(rhog*vfg2/xa(2),1.d0/84.d-6)
               ddrp=xa(2)/(rhog*vfg2)
               preduc=cell%p(i)/pcrit
               IF(preduc.lt.0.25d0)THEN
                  IF(preduc.lt.0.025d0)THEN
                     drmin=dcon(2)
                  ELSE
                     drmin=dcon(2)+(dcon(3)-dcon(2))*4.44444d0*(preduc-0.025d0)
                  ENDIF   
               ELSE
                  drmin=dcon(3)
               ENDIF
               ddrp=MAX(ddrp,drmin)
               ddrp=MIN(dh,ddrp)
               ddrp=MAX(MIN(ddrp,drop_max),8.6d-5)
               cell%ddrp(i)=ddrp  
!
               voidx=MAX(1.d-5,1.d0-alpdrp)
!               
!               surfa=10.8d0*alpdrp*rdiam
               surfa=3.6d0*alpdrp/ddrp
!               
               rey=xa(2)*voidx**3/(mug*SQRT(vfg2*voidx)) 
!               t1=SQRT(rey)
!               t2=SQRT(t1)
!               t=t1*t2
!               fic=rhog*surfa*(1.d0+0.1d0*rey**0.75d0)/rey
!               fic=rhog*surfa*(1.d0+0.1d0*t)/rey
               fic=rhog*surfa*MAX(0.05625d0,3.d0*(1.d0+0.1d0*rey**0.75d0)/rey)
!         
               vfgl_MPR=fic*vr            
            ELSE           
!               IF(xnn.ne.0.d0.and.alphag.eq.1.d0)THEN
!                  alphadrp=MAX(1.d-3,alphal)
!               ELSE
!                  alphadrp=MAX(1.d-4,alphal)
!               ENDIF
               alphadrp=MAX(1.d-4,alphal)                   !LSJ modification    
!            
               vfg=vr
!
               vfg2=MAX(vfg*vfg,xa(2)/rhog*MIN(dcon(2)*alphal**0.333333333d0,dht))
               ddrp=xa(2)/(rhog*vfg2)
               ddrp=MAX(ddrp,84.d-6)            
               ddrp=MIN(ddrp,dh)  
!            
!              mum=mug/(alphag**2*SQRT(alphag)) !d0
               mum=mug/(alphag**2.5d0) !d0
!               red=(1.d0-alphadrp)**2.5d0*rhog*ddrp*vr/mum
!               red=(1.d0-alphadrp)**2.5*rhog*ddrp*vr/mum
!              red=(1.d0-alphadrp)**2*SQRT(1.d0-alphadrp)*rhog*ddrp*vr/mug      !LSJ modification
               red=(1.d0-alphadrp)**2.5d0*rhog*ddrp*vr/mug      !LSJ modification
               red=MAX(red_min,red)
!               cd=24.d0/red*(1.d0+0.1d0*red**0.75d0)
!              t1=SQRT(red)
!              t2=SQRT(t1)
!              t=t1*t2
               cd=24.d0/red*(1.d0+0.1d0*red**0.75)
!              cd=24.d0/red*(1.d0+0.1d0*t)
               ia_MPR=3.6d0*alphadrp/ddrp
               IF(iam.eq.1)ia_MPR=cell%ia_MPR(i)
               IF(iam.eq.0)cell%ia_MPR(i)=ia_MPR
               vfgl_MPR=1.d0/8.d0*rhog*cd*ia_MPR*vr
!               vfgl_MPR=0.125d0*rhog*cd*ia_MPR*vr
!               vfgl_MPR=rhog*cd*ia_MPR*vr/8.d0
!               vfgl_MPR=vr*ia_MPR/8.d0
            ENDIF 
!
         ENDIF  
!      
!........drag coefficient model for inverted annular and inverted slug
!      
         CALL invann_invslg_for_invchn(0)
!      
!........drag coefficient model for mist
!      
         IF(mist.eq.1)THEN
!
!           For the mist flow regime, the MARS code model is not used. (only MARS manual model is used.)
!           If MARS code model is applied, RBHT result will be much deviated from the real situation.
!
            alphadrp=MAX(1.d-4,alphal)
!            
            vfg=vr
            vfg=MAX(DABS(vfg),0.001d0)
            vfg2=vfg*vfg !next
!            ddrp=xa(2)/(rhog*vfg2)
            ddrp=xa(3)/(rhog*vfg2)
            preduc=cell%p(i)/pcrit
            IF(preduc.lt.0.25d0)THEN
               IF(preduc.lt.0.025d0)THEN
                  drmin=dcon(2)
               ELSE
                  drmin=dcon(2)+(dcon(3)-dcon(2))*4.44444d0*(preduc-0.025d0)
               ENDIF   
            ELSE
               drmin=dcon(3)
            ENDIF
!            ddrp=MAX(ddrp,drmin)
!            ddrp=MIN(dh,dcon(2),ddrp)
            ddrp=MIN(dh,ddrp)         ! yjm  ori: on
            ddrp=MAX(MIN(ddrp,drop_max),8.4d-5)
            cell%ddrp(i)=ddrp 
!                       
            mum=mug
!            red=(1.d0-alphadrp)**2*SQRT(1.d0-alphadrp)*rhog*ddrp*vr/mum
!            red=(1.d0-alphadrp)**2.5d0*rhog*ddrp*vr/mum
            red=(1.d0-alphadrp)**2.5d0*rhog*ddrp*vfg/mum
            red=MAX(red_min,red)
!            t1=SQRT(red)
!            t2=SQRT(t1)
!            t=t1*t2
!            cd=MIN(24.d0/red*(1.d0+0.1d0*red**0.75d0),0.05626*8.d0) !next
            cd=MAX(24.d0/red*(1.d0+0.1d0*red**0.75d0),0.05626*8.d0) !next ! yjm  DMIN -> DMAX
!            cd=MIN(24.d0/red*(1.d0+0.1d0*t),0.05626*8.d0) !next
            ia_mist=3.6d0*alphadrp/ddrp
            IF(iam.eq.1)ia_mist=cell%ia_mist(i)
            IF(iam.eq.0)cell%ia_mist(i)=ia_mist
            vfgl_mist=1.d0/8.d0*rhog*cd*ia_mist*vr
!            
         ENDIF    
!
!........drag coefficient model for vertically stratified regime
!
         IF(VST.eq.1)THEN
!        
            dbub=MAX(dbub_min,5.d0*sigma/(rhol*vr*vr)) !We=10.0
            dbub=MIN(dbub_max,dbub)
            sigma_gr=sigma/gravity/rhol_g
            sigmagr=sigma*gravity*rhol_g
            vgjs=MIN(1.53d0,0.345d0*SQRT((dh/sigma_gr)*(sigmagr/rhol)))
            alphagupper=alphag_upper
            alphaglower=alphag_lower
            alphagsb=jg/(1.2d0*(jg+jl)+vgjs)
            alphagsb=MIN(alphag,MIN(alphagupper,MAX(alphagsb,alphaglower)))
            alphagsb=MIN(alphabs,alphagsb)
            !fmixlevel=(alphagupper-alphag)/(alphagupper-alphagsb)
            fmixlevel=(alphagupper-alphag)/MAX(1.d-5,alphagupper-alphagsb)
            fmixlevel=MIN(1.d0,MAX(0.d0,fmixlevel))
            alphagsb=alphagsb*fmixlevel
            ia_VST_sb=3.6d0*alphagsb/dbub
            ia_VST_surf=1.d0/lcell
            ia_VST=ia_VST_sb+ia_VST_surf
            IF(iam.eq.1)ia_VST=cell%ia_VST(i)
            IF(iam.eq.0)cell%ia_VST(i)=ia_VST         
            vfgl_VST=ia_VST*1000.d0*vr
!            
         ENDIF         
!
!........drag coefficient model for MPO
!
         IF(MPO.eq.1)THEN
!      
            vfg=vr
            alphadrp=MAX(1.d-4,alphal)
!         
            vfg=MAX(DABS(vfg),0.001d0)
            vfg2=vfg*vfg !next
            ddrp=xa(3)/(rhog*vfg2)
            preduc=cell%p(i)/pcrit
            IF(preduc.lt.0.25d0)THEN
               IF(preduc.lt.0.025d0)THEN
                  drmin=dcon(2)
               ELSE
                  drmin=dcon(2)+(dcon(3)-dcon(2))*4.44444d0*(preduc-0.025d0)
               ENDIF   
            ELSE
               drmin=dcon(3)
            ENDIF
            ddrp=MAX(ddrp,drmin)
!            ddrp=MIN(dh,dcon(2),ddrp)
            ddrp=MIN(dh,ddrp)
            ddrp=MAX(MIN(ddrp,drop_max),8.6d-5)
            cell%ddrp(i)=ddrp
!           
            mum=mug
!            red=(1.d0-alphadrp)**2*SQRT(1.d0-alphadrp)*rhog*ddrp*vr/mum
!            red=(1.d0-alphadrp)**2.5d0*rhog*ddrp*vr/mum
            red=(1.d0-alphadrp)**2.5d0*rhog*ddrp*vfg/mum
            red=MAX(red_min,red)
!            t1=SQRT(red)
!            t2=SQRT(t1)
!            t=t1*t2
!            cd=MIN(24.d0/red*(1.d0+0.1d0*red**0.75d0),0.05626d0*8.d0) !next
            cd=MAX(24.d0/red*(1.d0+0.1d0*red**0.75d0),0.05626d0*8.d0) !next  ! yjm MIN -> MAX
!            cd=MIN(24.d0/red*(1.d0+0.1d0*t),0.05626d0*8.d0) !next
            ia_MPO=3.6d0*alphadrp/ddrp
            IF(iam.eq.1)ia_MPO=cell%ia_MPO(i)
            IF(iam.eq.0)cell%ia_MPO(i)=ia_MPO
            vfgl_MPO=1.d0/8.d0*rhog*cd*ia_MPO*vr
!            
         ENDIF   
!
!.......sum of weighting factor
!      
         wf_preCHF =wf_liquid+wf_drift+ wf_bubbly+wf_slug  +wf_churn+wf_annular+wf_MPR+wf_gas
         wf_postDRY=wf_liquid+wf_invann+wf_invchn+wf_invslg+wf_mist            +wf_MPO+wf_gas
         wf_unVST  =wf_preCHF*(1.d0-wf_dry)+wf_postDRY*wf_dry
         wf_sum    =wf_unVST *(1.d0-wf_VST)+wf_VST       
!   
!........Obtain interfacial area
!
         ia_preCHF= wf_bubbly*ia_bubbly                        &
                   +wf_slug*(ia_slug_tb+ia_slug_sb)            &
                   +wf_churn*ia_churn                          &
                   +wf_annular*(ia_annular_ann+ia_annular_drp) &
                   +wf_MPR*ia_MPR         
         ia_postDRY= wf_invann*(ia_invann_ann+ia_invann_sb)  &
                    +wf_invchn*ia_invchn                     &
                    +wf_invslg*(ia_invslg_ann+ia_invslg_drp) &
                    +wf_mist*ia_mist                         &
                    +wf_MPO*ia_MPO
         ia_unVST=ia_preCHF*(1.d0-wf_dry)+ia_postDRY*wf_dry
         ia_i(i) =ia_unVST *(1.d0-wf_VST)+ia_VST*wf_VST
!   
!........Obtain drag coefficient using drag coefficient

         vfgl_preCHF= wf_slug*vfgl_slug       &
                     +wf_churn*vfgl_churn     &
                     +wf_annular*vfgl_annular &
                     +wf_MPR*vfgl_MPR         &
                     +wf_liquid*vfgl_liquid   &
                     +wf_gas*vfgl_gas         &
                     +wf_drift*vfgl_drift     &
                     +wf_bubbly*vfgl_bubbly 
         vfgl_postDRY= wf_liquid*vfgl_liquid &
                      +wf_gas*vfgl_gas       &
                      +wf_invann*vfgl_invann &
                      +wf_invchn*vfgl_invchn &
                      +wf_invslg*vfgl_invslg &
                      +wf_mist*vfgl_mist     &
                      +wf_MPO*vfgl_MPO
         vfgl_unVST=vfgl_preCHF*(1.d0-wf_dry)+vfgl_postDRY*wf_dry
         vfgl_i(i)=vfgl_unVST*(1.d0-wf_VST)+vfgl_VST*wf_VST
!
!........Obtain drag coefficient using drag force
!
         mdrag_preCHF= wf_liquid*vfgl_liquid*vr   &
                      +wf_gas*vfgl_gas*vr         &
                      +wf_drift*vfgl_drift*vrd    &
                      +wf_bubbly*vfgl_bubbly*vr   &
                      +wf_slug*vfgl_slug*vr       &
                      +wf_churn*vfgl_churn*vr     &
                      +wf_annular*vfgl_annular*vr &
                      +wf_MPR*vfgl_MPR*vr         
         mdrag_postDRY= wf_liquid*vfgl_liquid*vr &
                       +wf_gas*vfgl_gas*vr       &
                       +wf_invann*vfgl_invann*vr &
                       +wf_invchn*vfgl_invchn*vr &
                       +wf_invslg*vfgl_invslg*vr &
                       +wf_mist*vfgl_mist*vr     &
                       +wf_MPO*vfgl_MPO*vr
         mdrag_unVST=mdrag_preCHF*(1.d0-wf_dry)+mdrag_postDRY*wf_dry
         mdrag_whole=mdrag_unVST*(1.d0-wf_VST)+vfgl_VST*vr*wf_VST
!      
!........c1,c0
!
         IF(wf_drift.ge.1.d0)THEN !drift-flux fomulation of momentum eq.
            drift_c1(i)=c1
            drift_c0(i)=c0
!            IF(vrd.ne.0.d0) THEN
!               vfgl_i(i)=mdrag_whole/vrd         
!            ELSE
!               vfgl_i(i)=mdrag_whole        !temporary pik
!            ENDIF   
         ELSE !drag-cofficient fomulation of momentum eq.
            drift_c1(i)=1.d0
            drift_c0(i)=1.d0 
!            vfgl_i(i)=mdrag_whole/vr         
         ENDIF  
!
!........droplet & gas
!   
         vfgd_i(i)=1.d4  !next
!       
!........limit the magnitude of drag coefficient
!
!         vfgl_i(i)=MAX(100.d0,MIN(vfgl_i(i),1.d6))
         vfgl_i(i)=MAX(1.d-2,MIN(vfgl_i(i),1.d6)) 
!      
      ENDDO
!   
      ENDSUBROUTINE rv_int_fric_model
