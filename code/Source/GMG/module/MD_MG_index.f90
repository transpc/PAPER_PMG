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
! G2 셋업 분배 이중 모드 (LOOP C011): 0=파일(MG_tmp), 1=MPI 통신
!   stg_*: rank0 이 writer(subdomain_infor_mg)에서 채우고 reader(read_mesh_MPI)가
!          BCAST 로 분배하는 스테이징 메타 (PMG_infor 파일 내용과 1:1)
      INTEGER(4) :: isetup_comm
      INTEGER(4),DIMENSION(:,:),ALLOCATABLE :: stg_iintf,stg_inodegl,stg_inbdc,stg_ialvP
      INTEGER(4),DIMENSION(:),ALLOCATABLE :: stg_inmax,stg_nnzc0,stg_nnzi,stg_nnzr
!
    END MODULE MD_MG_index
