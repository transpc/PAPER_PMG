      MODULE Zdel_scalar
!      
      IMPLICIT NONE
      SAVE
!
      INTEGER smac3_pres_eng,prn_div_eng
      INTEGER stmtbl_repeat_for_nc,max_ihtc_opt,relax_interface_dtemp
      INTEGER limit_eng_src_opt,limit_iht_opt,suspend_iht_opt,suspend_erg_opt
      REAL(8) ag_min_hig,al_min_hil,qu_min_hgf,dsrc
      REAL(8) max_ihtc_opt_coeff 
!      
      REAL(8),Allocatable::del_eg(:),del_el(:),del_x(:),del_ag(:),del_ad(:),del_rhog(:),del_rhol(:)
!
      INTEGER,Allocatable::err_stm(:)
!
      END MODULE Zdel_scalar
      