      Module Zsbloca
! 
      IMPLICIT NONE
      SAVE
!      
      INTEGER ncell_si(2)
      REAL(8) pzr_level,      &
               rg_break, rl_break, eg_break, el_break, cg_break, cl_break, break_flow,   &
               si_flow, break_flow_eng, si_flow_eng,        &
               q_liq,ge_err, & 
               break_flow_eng_g, break_flow_eng_l,     & 
               eng_gg,pw_g,pw_l,pw_l_int       
! bug created in sbloca_out_user uninitialized
      REAL(8) :: eng_gg_int=0.d0
      REAL(8) :: ge_err_int=0.d0
      REAL(8) :: pw_g_int=0.d0
      REAL(8) :: break_flow_int=0.d0
      REAL(8) :: si_flow_int=0.d0
      REAL(8) :: break_flow_eng_int=0.d0
      REAL(8) :: si_flow_eng_int=0.d0
      REAL(8) :: q_liq_int=0.d0
      REAL(8) :: break_flow_eng_g_int=0.d0
      REAL(8) :: break_flow_eng_l_int=0.d0
! bug created in sbloca_out_user uninitialized called only if time.gt.300
      REAL(8) :: p_break=0.d0
      REAL(8) :: t_break=0.d0
      REAL(8) :: q_break=0.d0
      REAL(8) :: rhom_break=0.d0
      REAL(8) :: h_break=0.d0
      REAL(8) :: a_break=0.d0
!
      ENDMODULE Zsbloca
