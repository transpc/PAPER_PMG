      MODULE MD_MG_index
! 
      IMPLICIT NONE
      SAVE
!
      INTEGER(4) nlevel,ncycle,mxnbne,maxit_1,ip_nmax,nmax,iter_mg,n_GC,mxnbne_mg
      REAL(8) crit_1,crit_bcg_mg
      INTEGER(4) iter_max,nlevel_N,n1_min,n2_min,ioplv,ip_lev
      CHARACTER :: report_text*100
      INTEGER(4) itergs(20)
! icase_MG: pressure_solve 2차(non-orth) 호출의 SOLVE_GMG 인자 — 고정 2 (행렬 불변)
      INTEGER(4) :: icase_MG
      INTEGER(4),DIMENSION(:),ALLOCATABLE :: isend_m,irecv_m
! NEW for polynormial
      integer (4) icheb(5)
      real(8) rcheb(5)
! l1-보정 스무더 (해결책 A, Baker et al. 2011): 코어스 레벨 GS 대각에
!   랭크 밖(고스트 열) 연결 |a_ij| 합을 더해 파티션 무관 수렴 보장.
!   0 = 기존(순수 대각), 1 = l1 보정. np=1 에서는 고스트가 없어 두 모드 동일.
      INTEGER(4) :: il1_gs
! G2 셋업 분배 이중 모드 (LOOP C011): 0=파일(MG_tmp), 1=MPI 통신
!   stg_*: rank0 이 writer(subdomain_infor_mg)에서 채우고 reader(read_mesh_MPI)가
!          BCAST 로 분배하는 스테이징 메타 (PMG_infor 파일 내용과 1:1)
      INTEGER(4) :: isetup_comm
      INTEGER(4),DIMENSION(:,:),ALLOCATABLE :: stg_iintf,stg_inodegl,stg_inbdc,stg_ialvP
      INTEGER(4),DIMENSION(:),ALLOCATABLE :: stg_inmax,stg_nnzc0,stg_nnzi,stg_nnzr
!   finest fan-out (part###.out) 스테이징 (C011-3): prc 순 연접 정수/실수 스트림 + prc별 길이
      INTEGER(4),DIMENSION(:),ALLOCATABLE :: stg_fibuf,stg_ficnt,stg_frcnt
      REAL(8),DIMENSION(:),ALLOCATABLE :: stg_frbuf
!   coarse fan-out (part_MG###.out) 스테이징 (C011-4): 레벨-major 루프라 카운트
!   선계산이 불가 — prc별 성장형(배증) 스트림에 pack, reader 가 연접→SCATTERV
      TYPE :: stg_buf_t
         INTEGER(4),DIMENSION(:),ALLOCATABLE :: ib
         REAL(8),DIMENSION(:),ALLOCATABLE :: rb
         INTEGER(4) :: ni = 0
         INTEGER(4) :: nr = 0
      END TYPE stg_buf_t
      TYPE(stg_buf_t),DIMENSION(:),ALLOCATABLE :: stg_mg
!
    CONTAINS
!
      SUBROUTINE stg_mg_init(nprc)
      INTEGER(4), INTENT(IN) :: nprc
      INTEGER(4) :: p
      IF(ALLOCATED(stg_mg)) DEALLOCATE(stg_mg)   ! ioplv 재진입(GOTO 500) 대비
      ALLOCATE(stg_mg(nprc))
      DO p = 1, nprc
         ALLOCATE(stg_mg(p)%ib(4096), stg_mg(p)%rb(4096))
         stg_mg(p)%ni = 0
         stg_mg(p)%nr = 0
      ENDDO
      END SUBROUTINE stg_mg_init
!
      SUBROUTINE stg_pushi(p, iv)
      INTEGER(4), INTENT(IN) :: p, iv
      INTEGER(4),DIMENSION(:),ALLOCATABLE :: t
      IF(stg_mg(p)%ni .EQ. SIZE(stg_mg(p)%ib)) THEN
         ALLOCATE(t(2*stg_mg(p)%ni))
         t(1:stg_mg(p)%ni) = stg_mg(p)%ib
         CALL MOVE_ALLOC(t, stg_mg(p)%ib)
      ENDIF
      stg_mg(p)%ni = stg_mg(p)%ni + 1
      stg_mg(p)%ib(stg_mg(p)%ni) = iv
      END SUBROUTINE stg_pushi
!
      SUBROUTINE stg_pushr(p, xv)
      INTEGER(4), INTENT(IN) :: p
      REAL(8), INTENT(IN) :: xv
      REAL(8),DIMENSION(:),ALLOCATABLE :: t
      IF(stg_mg(p)%nr .EQ. SIZE(stg_mg(p)%rb)) THEN
         ALLOCATE(t(2*stg_mg(p)%nr))
         t(1:stg_mg(p)%nr) = stg_mg(p)%rb
         CALL MOVE_ALLOC(t, stg_mg(p)%rb)
      ENDIF
      stg_mg(p)%nr = stg_mg(p)%nr + 1
      stg_mg(p)%rb(stg_mg(p)%nr) = xv
      END SUBROUTINE stg_pushr
!     list-directed ASCII 왕복과 bitwise 동일한 라운딩 (C011-1 실험 검증).
!     과도기 통신 모드가 파일 모드와 골든 bitwise 를 유지하기 위한 shim —
!     제거(정확값 전달)는 C011-5 에서 별도 수치 사이클로 수행
      FUNCTION rt_ascii(x) RESULT(y)
      REAL(8), INTENT(IN) :: x
      REAL(8) :: y
      CHARACTER(32) :: buf
      WRITE(buf,*) x
      READ(buf,*) y
      END FUNCTION rt_ascii
!
    END MODULE MD_MG_index
