!=======================================================================
!     dump_pmg — PMG 골든 회귀용 덤프 훅 (LOOP C008, PLAN §4-2)
!
!     환경변수 CUPID_PMG_DUMP="<itim>[,<itim>...]" 가 설정된 경우에만 활성.
!     미설정 시 첫 호출에서 비활성 판정 후 즉시 반환 (평상시 오버헤드 ~0).
!
!     기록 위치: 실행 디렉토리의 pmg_dump/ (rank 0 이 생성)
!       s<itim>_k<site>_r<rank>_c<call>.pre  : 솔브 입력
!           [ver, itim, isite, icall, n, nnz0c] + src(n) + au0c(nnz0c) + diag(n)
!       s<itim>_k<site>_r<rank>_c<call>.post : 수렴해
!           [ver, itim, isite, icall, n] + u(n)
!       setup_r<rank>.bin (덤프 최초 1회)     : 셋업 입력 (standalone 재생용)
!           [ver, ndim, nelem_mg, nf_max] + num_neigh_mg + neigh_mg
!           + xloc_tmp + celem
!     포맷: unformatted stream (-r8 정밀도 보존). its 는 fort.501(rank0) 참조.
!
!     유닛은 NEWUNIT= 사용 — 기존 코드의 고정 유닛 대역과 충돌하지 않음
!     (PLAN §5-1 "유닛 번호 위생" 참조).
!=======================================================================
      MODULE dump_pmg
!
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: dump_pmg_pre, dump_pmg_post, dump_pmg_setup_hook
!
      INTEGER, PARAMETER :: MAXSTEPS = 16       ! 덤프 대상 스텝 최대 수
      INTEGER, PARAMETER :: IVER = 1            ! 덤프 포맷 버전
      LOGICAL, SAVE :: checked = .FALSE.        ! env 파싱 완료 여부
      LOGICAL, SAVE :: active  = .FALSE.
      LOGICAL, SAVE :: setup_done = .FALSE.
      INTEGER, SAVE :: nsteps = 0
      INTEGER, SAVE :: steps(MAXSTEPS) = -1
      INTEGER, SAVE :: ncall(2) = 0             ! (사이트별) 스텝 내 호출 카운터
      INTEGER, SAVE :: itim_last = -1
!
      CONTAINS
!
!-----------------------------------------------------------------------
      LOGICAL FUNCTION dump_now(itim)
!     env 1회 파싱 + 현재 스텝이 덤프 대상인지 판정
      INTEGER, INTENT(IN) :: itim
      CHARACTER(256) :: buf
      INTEGER :: stat, lenb, i, i1, i2
!
      IF (.NOT. checked) THEN
         checked = .TRUE.
         CALL GET_ENVIRONMENT_VARIABLE('CUPID_PMG_DUMP', buf, lenb, stat)
         IF (stat .EQ. 0 .AND. lenb .GT. 0) THEN
            active = .TRUE.
            i1 = 1
            DO WHILE (i1 .LE. lenb .AND. nsteps .LT. MAXSTEPS)
               i2 = INDEX(buf(i1:lenb), ',')
               IF (i2 .EQ. 0) THEN
                  i2 = lenb + 1
               ELSE
                  i2 = i1 + i2 - 1
               END IF
               nsteps = nsteps + 1
               READ (buf(i1:i2-1), *) steps(nsteps)
               i1 = i2 + 1
            END DO
         END IF
      END IF
!
      dump_now = .FALSE.
      IF (.NOT. active) RETURN
      DO i = 1, nsteps
         IF (steps(i) .EQ. itim) dump_now = .TRUE.
      END DO
      END FUNCTION dump_now
!
!-----------------------------------------------------------------------
      SUBROUTINE dump_pmg_pre(isite, nnz0c, src, au0c, diag)
!     SOLVE_GMG 직전 호출 — 솔브 입력 기록 (+최초 1회 셋업 기록)
      USE Zcore,       ONLY: myrank
      USE Ztimecon,    ONLY: itim
      USE MD_MPI,      ONLY: nintf
!
      INTEGER, INTENT(IN) :: isite, nnz0c
      REAL(8), INTENT(IN) :: src(nintf), au0c(nnz0c), diag(nintf)
      INTEGER :: iu
      CHARACTER(64) :: fn
!
      IF (.NOT. dump_now(itim)) RETURN
!
      IF (itim .NE. itim_last) THEN            ! 새 스텝 → 호출 카운터 리셋
         ncall = 0
         itim_last = itim
      END IF
      ncall(isite) = ncall(isite) + 1
!
      WRITE (fn, '(A,I0,A,I0,A,I0,A,I0,A)') 'pmg_dump/s', itim, '_k',    &
             isite, '_r', myrank, '_c', ncall(isite), '.pre'
      OPEN (NEWUNIT=iu, FILE=fn, FORM='unformatted', ACCESS='stream',    &
            STATUS='replace')
      WRITE (iu) IVER, itim, isite, ncall(isite), nintf, nnz0c
      WRITE (iu) src(1:nintf)
      WRITE (iu) au0c(1:nnz0c)
      WRITE (iu) diag(1:nintf)
      CLOSE (iu)
      END SUBROUTINE dump_pmg_pre
!
!-----------------------------------------------------------------------
      SUBROUTINE dump_pmg_post(isite)
!     SOLVE_GMG 직후 호출 — 수렴해 u* 기록 (기준값)
      USE Zcore,       ONLY: myrank
      USE Ztimecon,    ONLY: itim
      USE MD_MPI,      ONLY: nintf
      USE MD_matrix,   ONLY: u
!
      INTEGER, INTENT(IN) :: isite
      INTEGER :: iu
      CHARACTER(64) :: fn
!
      IF (.NOT. dump_now(itim)) RETURN
!
      WRITE (fn, '(A,I0,A,I0,A,I0,A,I0,A)') 'pmg_dump/s', itim, '_k',    &
             isite, '_r', myrank, '_c', ncall(isite), '.post'
      OPEN (NEWUNIT=iu, FILE=fn, FORM='unformatted', ACCESS='stream',    &
            STATUS='replace')
      WRITE (iu) IVER, itim, isite, ncall(isite), nintf
      WRITE (iu) u(1:nintf)
      CLOSE (iu)
      END SUBROUTINE dump_pmg_post
!
!-----------------------------------------------------------------------
      SUBROUTINE dump_pmg_setup_hook
!     셋업 입력 기록 — pmg_standalone 이 프로덕션 체인을 재생할 때 필요한 전부.
!     반드시 subdomain_infor_MG **호출 전**에 불러야 함 — 그 루틴이 종료 시
!     num_neigh_mg/neigh_mg 를 DEALLOCATE 함 (6_subdomain_infor_mg.f90:644,
!     LOG C008-r1 에서 확인). 호출처: read_grid.f90 의 MG 블록.
!     env CUPID_PMG_DUMP 미설정 시 no-op.
      USE Zcore,       ONLY: myrank
      USE Zparam,      ONLY: ndim
      USE Zcoord1,     ONLY: xloc_tmp
      USE Zmpi,        ONLY: celem
      USE md_geometry, ONLY: nelem_mg, num_neigh_mg, neigh_mg
      USE MD_parameter, ONLY: nf_max
!
      INTEGER :: iu
      CHARACTER(64) :: fn
!
      IF (.NOT. dump_now(-999999)) THEN        ! env 파싱만 수행 (스텝 무관)
         IF (.NOT. active) RETURN
      END IF
      IF (setup_done) RETURN
      setup_done = .TRUE.
!
      IF (myrank .EQ. 0) CALL system('mkdir -p pmg_dump')
      WRITE (fn, '(A,I0,A)') 'pmg_dump/setup_r', myrank, '.bin'
      OPEN (NEWUNIT=iu, FILE=fn, FORM='unformatted', ACCESS='stream',    &
            STATUS='replace')
      WRITE (iu) IVER, ndim, nelem_mg, nf_max
      WRITE (iu) num_neigh_mg(1:nelem_mg)
      WRITE (iu) neigh_mg(1:nf_max, 1:nelem_mg)
      WRITE (iu) xloc_tmp(1:nelem_mg, 1:ndim)
      WRITE (iu) celem(1:nelem_mg)
      CLOSE (iu)
      END SUBROUTINE dump_pmg_setup_hook
!
      END MODULE dump_pmg
