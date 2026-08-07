!
      MODULE Zrv_int_friction
      IMPLICIT NONE
      SAVE
!.....will be removed
      INTEGER,ALLOCATABLE :: iregime(:,:)      
!.....control
      INTEGER :: i,iam
      INTEGER :: drift_model,drag_model !input
      INTEGER :: drag,drift,bestion,epri,fore,asali
      INTEGER :: liquid,gas,bubbly,slug,churn,annular,MPR
      INTEGER :: invann,invchn,mist,MPO,VST,invslg
      INTEGER :: fdir,upward,sideward
!
      REAL(8) :: alpha_am !parameter
      REAL(8) :: alphabs,alphacd,alphade,alphasa,alphaam,alphags
      REAL(8) :: alphag,alphal,alphab,alphagb,alphadrp
      REAL(8) :: alpha_min,alpha_max !parameter
      REAL(8) :: alphamin,alphamax
      REAL(8) :: alphaglower,alphagupper,alphag_lower,alphag_upper,fmixlevel !VST
      REAL(8) :: ia_bubbly,ia_slug_tb,ia_slug_sb,ia_churn,ia_annular_drp,ia_annular_ann,ia_MPR
      REAL(8) :: ia_invann_ann,ia_invann_sb,ia_invchn,ia_invslg_drp,ia_invslg_ann,ia_mist,ia_mpo
      REAL(8) :: ia_VST_sb,ia_VST_surf,ia_VST
      REAL(8) :: ia_PRECHF,ia_POSTDRY,ia_unVST
      REAL(8) :: lcell
      REAL(8),ALLOCATABLE :: alphaff(:),alpha_b(:)
!                
      REAL(8) :: pres
      REAL(8) :: ddrp_min,ddrp_max,dbub_max,dbub_min,red_min,reb_min !parameter      
      REAL(8) :: ddrp,red,dbub,reb,dgas,cd
      REAL(8) :: mum,mug,mul,rhog,rhol,rhol_g,sigma
      REAL(8) :: vtb,vm,vg,vl,rhogoverl,rholoverg,gravity,phij
      REAL(8) :: vrd,drop_max
      REAL(8),ALLOCATABLE :: vr_drift(:)
!.....weighting factor      
      REAL(8) :: wf_dry,wf_VST,wf_unVST,wf_preCHF,wf_postDRY
      REAL(8) :: wf_liquid,wf_gas,wf_drift,wf_bubbly,wf_slug,wf_churn,wf_annular,wf_MPR
      REAL(8) :: wf_invann,wf_invchn,wf_invslg,wf_mist,wf_MPO
      REAL(8) :: wf_epri,wf_bestion
      REAL(8) :: wfactor,wf_sum
!.....drag_coefficient
      REAL(8) :: vfgl_liquid,vfgl_gas,vfgl_bubbly,vfgl_slug,vfgl_churn,vfgl_annular,vfgl_MPR
      REAL(8) :: vfgl_slug_tb,vfgl_slug_sb,vfgl_annular_ann,vfgl_annular_drp
      REAL(8) :: vfgl_invann,vfgl_invchn,vfgl_invslg,vfgl_mist,vfgl_MPO
      REAL(8) :: vfgl_invann_ann,vfgl_invann_sb,vfgl_invslg_ann,vfgl_invslg_drp
      REAL(8) :: vfgl_VST,vfgl_unVST
      REAL(8) :: vfgl_drift,vfgl_sum
      REAL(8) :: vfgl_preCHF,vfgl_postDRY
      REAL(8) :: mdrag_drift,mdrag_annular,mdrag_churn
      REAL(8) :: mdrag_preCHF,mdrag_postDRY,mdrag_unVST,mdrag_whole
!.....drift  
      REAL(8) :: vr,dh,dht !common         
      REAL(8) :: cp,L,re,reg,rel,jg,jl,k0,r,a1,b1,b2,d1,d2     !epri 
      REAL(8) :: vgj0,fdebug,c2,c3,c4,c5,c6,c7,c8,c9,c10             !epri
      REAL(8) :: c0,c1,vgj,ci                                        !epri&bestion
!.....drag
      REAL(8) :: fi,f1,f2,f3,vgmul,cann,uc,kug_crit,vlabs,cf,alphalf,alphald !annular
      REAL(8) :: alphab_1 !invslg
      REAL(8) :: inertens,dstar,dstar1,dstar2,delta,deltastar !invann
      REAL(8) :: sigmagr,sigma_gr,vgjs,alphagsb !VST
!.....droplet diameter
      REAL(8) :: we,xnn         
!.....reflod
      REAL(8) :: xa(3),dcon(3),web(3)             
!.....input
      CHARACTER(30) s_rv_int_fric
!      
      ENDMODULE Zrv_int_friction
