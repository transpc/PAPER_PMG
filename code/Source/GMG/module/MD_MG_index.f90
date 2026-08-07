      MODULE MD_MG_index
! 
      IMPLICIT NONE
      SAVE
!
      INTEGER(4) nlevel,ncycle,mxnbne,maxit_1,ip_nmax,nmax,iter_mg,n_GC, isth,AR_hi,mxnbne_mg
      REAL(8):: relax
      REAL(8) crit_1,crit_bcg_mg
      INTEGER(4) iter_max,nlevel_N,n1_min,n2_min,ioplv,ip_lev,isol_mg, id_GS_sym
      CHARACTER :: report_text*100
      INTEGER(4) itergs(20)
      INTEGER(4),DIMENSION(:),ALLOCATABLE :: isend_m,irecv_m
! NEW for polynormial
      integer (4) icheb(5)
      real(8) rcheb(5)
! re-calculate: stiffness_MG
      INTEGER icase_MG
      INTEGER :: ihybrid
! hibrid solver
      INTEGER(4) isol_start, i_precond
      INTEGER(4) nsol_start, nsol_start_mg
!
    END MODULE MD_MG_index
