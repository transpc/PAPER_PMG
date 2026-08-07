      MODULE Zrv_model
! 
      IMPLICIT NONE
      SAVE
!
      INTEGER::rv_model,rv_fw_reg,rv_ht_str,rv_ht_w,rv_fric_i,rv_ht_i,rv_fric_w
      INTEGER::rv_choke,rv_gapcond
      INTEGER::ia_option       
!
      LOGICAL:: free_model
      LOGICAL:: lfric_swap,lhtc_swap,lfricw_swap
!
!...mcp model, valve model
!      
      INTEGER::rv_mcp,rv_valve       
!
      END MODULE Zrv_model
