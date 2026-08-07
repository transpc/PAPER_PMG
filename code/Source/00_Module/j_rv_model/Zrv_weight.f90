MODULE Zrv_weight
!
   IMPLICIT NONE
   SAVE
!
   INTEGER liquid,gas
   INTEGER bubbly,slug,churn,annular,mpr
   INTEGER invann,invchn,invslg,mist,mpo
   INTEGER vst
!
   REAL(8) wf_liquid,wf_gas                                   !single-phase
   REAL(8) wf_bubbly,wf_slug,wf_churn,wf_annular,wf_MPR       !pre-CHF
   REAL(8) wf_invann,wf_invchn,wf_invslg,wf_mist,wf_MPO       !post-dry 
   REAL(8) wf_VST
   REAL(8) wf_sum
!
   REAL(8) pfnrgj
!
END MODULE Zrv_weight
