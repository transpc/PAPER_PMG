MODULE Zrv_gap_cond
!
   IMPLICIT NONE
   SAVE
!
   INTEGER nr_gapi,nr_gapo  
   INTEGER irupt,iplas
   INTEGER ngas
   INTEGER indgas(7)  
!
   REAL(8) block,cltave,dtdt,hte_clad_o,StrPlas,GapDisp_fission,GapWdth,CladExR
   REAL(8) gap_P0,GapRough
   REAL(8) molgas(7)   
   REAL(8) rhocp_gap
   REAL(8), ALLOCATABLE:: hte(:),GapIntW(:),ts_old(:),tg_old(:),tl_old(:),cond_gap(:),h_gap(:),width_gap(:),block_gap(:)
   
!
!   LOGICAL GapConductance   
!  
END MODULE Zrv_gap_cond