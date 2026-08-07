      MODULE Zconst1
! 
      IMPLICIT NONE
      SAVE
!
      CHARACTER*70 vv_prob
      INTEGER noutput,restart,save_option,iat
      INTEGER iprofbc,mtopol,mdrag,mHTC,mboron,iheatpart,parallel
      INTEGER mdiffoff,lsquareoff,cplmars,cplmaster,turbubble,turboil,topolsurface
      INTEGER mdiffscheme
      INTEGER buoyancy_turb
      INTEGER iturb,turb_phase,lowreynolds
      INTEGER i_macroturb_source,i_turb_zero,i_wall_fric,i_turb_disp,i_wall_lub, &
              i_bubble_diameter,i_drop_diameter,i_rv_int_fric,i_critical_flow, &
              i_subchannel_fric,i_subchannel_fric_axial,i_subchannel_fric_cross,i_subchannel_mixing,i_2p_multiplier
      INTEGER nd_face,fric_face,p_work_face
      INTEGER wconden
      INTEGER iVisRatio
      INTEGER rv_htmodel_forCFD        !use  RV heat transfer model with solid grid (not with RV heat structures)      
      REAL(8) vis_ratio
      REAL(8) Nlift,Nwlf,Ntdf
!
      LOGICAL lwconden_alphal0,lrestart_changed_nbcon,lrestart_overwrite      
!
      END MODULE Zconst1
