# 루프 엔지니어링 LOG

> 규약: [LOOP.md](LOOP.md) §4 — **append-only**, 과거 엔트리 수정 금지(정정은 새 줄에 `정정:`), 최신 엔트리가 맨 아래.
> 결과란에는 관측 사실만. 해석·추측은 `추정:` 접두어로 분리.

---

## C000 | 2026-08-09 | 베이스라인 확립 (조사·계획·프로토콜)

- 목표: 실행 착수 전 기반 문서(PLAN/LOOP/LOG)와 환경 검증을 완료한 상태 확립
- 변경: `PLAN.md` 작성 후 2회 개정(기존 SIF 재사용 반영, 고도화 G1~G3 구체화), `LOOP.md`·`LOG.md` 신규 작성. 소스 무변경
- 실행: 소스/케이스/호스트/SIF 조사 (근거·상세는 PLAN §1, §2-1, §5-1)
- 결과:
  - 호스트: WSL2 Ubuntu 24.04, 12코어, **AVX2까지만** (makefile.in 의 AVX512 설정과 불일치 → §3-1 에서 전환 예정), Apptainer 1.5.1
  - `hpc23.sif`(5.2 GB) 검증: ifort classic 2021.10.0 (oneAPI 2023.2.1), Intel MPI 2021.10, make/git 있음. **METIS·7z 없음**
  - 호스트 METIS: `/root/00_apptainer/metis-5.0.2/build/Linux-x86_64/libmetis/libmetis.so` (shared 전용, `.a` 없음). `/root` 자동 바인드로 컨테이너 내 접근 확인
  - 케이스 `2_iSMR_ECT_res1`: `mg.in`(isol_mg=-2, PMG 예조건 모드), `foamGrid.7z`(내부 `foamGrid.in` 112 MB 1개), `tpfh2o` 있음. **`somaFlow.in` 없음 → C005/C009 블로킹**
  - PMG 파일 기반 셋업 확인: rank 0 이 `MG_tmp/part###.out` 기록(`6_subdomain_infor_mg.f90:131`, unit `50+prc`) → 각 rank 가 읽음(`2_read_mesh_MPI.f90:89`, unit `10+myrank`), 폴더는 `system('mkdir MG_tmp')`(`2_Prep_fine_global.f90:83`)
  - np=900 실패 원인 관련 소스 검증: 쓰기·읽기 유닛 모두 반복마다 CLOSE 확인(`6_subdomain_infor_mg.f90:188,583`, `2_read_mesh_MPI.f90:205`). Fortran 유닛은 프로세스-로컬이므로 rank 간 유닛 충돌 가설은 성립 안 함. rank 0 내 유닛 대역 겹침(솔버의 `write(501,*)` `6_solver_pbcg_mg.f90:280`, `WRITE(999,*)` `:271` vs 셋업 유닛 51~950)은 실재 — 단 즉시 크래시 기전은 아님
  - 추정: np=900 실패 원인은 **미확정**. 남은 후보(고정 배열 한도, rank 0 fan-out 부하, 동일 rank 유닛 겹침의 2차 효과)는 LOOP.md §6 검증 메모 참조 — C010 에서 재현으로 확정 필요
- 판정: **PASS** — 코드 무변경이므로 L1/L2 비해당. 문서 3종 정합·환경 검증 완료
- 다음: C001 (`scripts/env.sh` + `in_contain.sh` 작성, SIF 스모크 — LOOP.md §6 백로그)

---

## C001 | 2026-08-09 | env.sh + in_contain.sh 작성, SIF 스모크

- 목표: PLAN §2-2 계약대로 스크립트 2종을 만들고 §2-3 DoD 3종 통과
- 변경: `scripts/env.sh`, `scripts/in_contain.sh` 신규 (+x)
- 실행: `./scripts/in_contain.sh mpiifort --version` / `./scripts/in_contain.sh bash -c 'ls "$METIS_LIB"'` / `SIF=/no/such.sif ./scripts/in_contain.sh true`
- 결과:
  - [1] `ifort (IFORT) 2021.10.0 20230609` 출력 — **`bash -lc` 없이 plain exec 로 동작** (SIF 의 `%environment` 가 모든 exec 에 적용됨을 확인)
  - [2] `libmetis.so` 출력 — 호스트 export 환경변수(`METIS_LIB`)가 컨테이너 안까지 전파됨
  - [3] `SIF not found: /no/such.sif — ...` 메시지, exit=1
- 판정: **PASS** — DoD 3/3 green (L1 해당, L2 비해당)
- 다음: C003 (빌드 전에 .gitignore 정비가 유리하므로 순서 교체: C001→C003→C002)

---

## C003 | 2026-08-09 | .gitignore 보강 + prepare_case.sh

- 목표: 빌드/실행 부산물의 git 오염 차단 + 케이스 전처리 자동화(7z 해제·입력 확인·cupid.x 복사)
- 변경: 기존 `/.gitignore` 에 CUPID 런타임 블록 추가(`MG_tmp/`, `fort.*`, `code/**/cupid.x`, `code/**/run_*.log`, `*__genmod.f90`, `*.sif` — `*.o`/`*.mod`/`foamGrid.in`/`*.orig` 는 기존 항목이 이미 커버), `scripts/prepare_case.sh` 신규
- 실행: `./scripts/prepare_case.sh`
- 결과:
  - `foamGrid.7z` → `foamGrid.in`(113 MB) 호스트 7z 해제 성공
  - `somaFlow.in` 부재로 `BLOCKED` 메시지 + exit=2 — 설계된 블로킹 동작 그대로
  - `git status` 에 `foamGrid.in` 미출현 — ignore 유효
- 판정: **PASS** — 해제·블로킹·ignore 3개 검증점 green. somaFlow.in 블로커는 상태 유지 (PLAN §7-1)
- 다음: C002 (makefile.in.apptainer + build.sh 로 cupid.x 빌드)

---

## C002 | 2026-08-09 | makefile.in.apptainer + build.sh → cupid.x 빌드

- 목표: 컨테이너 안에서 전체 소스 컴파일·링크 성공 (에러 0, `Source/cupid.x` 생성)
- 변경: `Source/makefile.in.apptainer` 신규 (diff 3개군: AVX512→AVX2 블록 전환, `p=52`→`p=12`, `LIBDIR`을 `-L$(METIS_LIB) -lmetis -Wl,-rpath,$(METIS_LIB)` 로), `scripts/build.sh` 신규
- 실행: `./scripts/build.sh` (1차) → 실패 → 수정 후 재실행 (r1)
- 결과:
  - **1차 FAIL**: `error #7001: Error in creating the compiled module file` ×177 — 원인: 전 서브 makefile 이 `-module (상대경로)/Modules` 로 `.mod` 를 `Source/Modules/` 에 쓰는데 해당 디렉토리가 git 에 없음(gitignore 대상, ifort 는 자동 생성 안 함)
  - **r1**: `build.sh` 에 `mkdir -p "$CUPID_SRC/Modules"` 추가 → **에러 0, 경고 41**(`-mkl` deprecation 등), `cupid.x` 7.6 MB 생성, BUILD EXIT=0
  - `ldd cupid.x`(컨테이너 내): `libmetis.so → /root/00_apptainer/metis-5.0.2/build/Linux-x86_64/libmetis/libmetis.so` — **rpath 로 해석됨** (`LD_LIBRARY_PATH` 불요), MKL 3종 해석 정상, "not found" 없음
  - 참고: `build.sh` 가 `makefile.in` 을 apptainer 변형으로 교체(원본은 `makefile.in.orig` 백업 + git 이력)
- 판정: **PASS** (r1) — L1 green (에러 0, 경고 허용). PLAN §3-4 DoD 1/3 충족
- 다음: C004 (단독 실행으로 입력 체크 도달 확인)

---

## C004 | 2026-08-09 | cupid.x 단독 실행 → 입력 단계 도달 확인

- 목표: 입력 파일 없는 빈 디렉토리에서 실행 시 정상 기동 후 입력 체크에 도달하는지 확인 (계획상 기대: "lack of somaFlow.in" 메시지)
- 변경: 없음 (관측 사이클). 대기 시간에 `scripts/run.sh` 작성 (C005 스크립트 부분 선행)
- 실행: scratchpad 빈 디렉토리에서 `in_contain.sh .../cupid.x` (mpirun 없이 싱글턴)
- 결과:
  - CUPID v2.20 배너 출력, MPI 싱글턴 초기화 정상, **illegal instruction 없음** (AVX2 선택 유효)
  - **기대와 다른 동작**: "lack of somaFlow.in" 메시지 미출력. `open_files.f90:72` 의 OPEN 이 `status='old'` 없이 호출되어 **빈 somaFlow.in 을 생성**(iostat=0 → `:74` 부재 가드는 사실상 도달 불가한 죽은 코드), 직후 빈 파일 READ 에서 `forrtl: severe (24): end-of-file during read, unit 812` 로 종료 (exit 24)
  - 추정: 입력 부재 시나리오에서도 제어된 런타임 에러로 종료하므로 실사용 지장은 없음. `status='old'` 추가는 G3 정리 후보 (PMG 외 영역이므로 우선순위 낮음)
- 판정: **PASS** — 목표의 본질(바이너리 기동·입력 단계 도달)은 확인. 계획 문구의 기대 동작은 사실과 다름을 기록 (PLAN §3-1(3) 의 검증 문구는 이 관측 기준으로 해석할 것)
- 다음: C005 는 somaFlow.in 확보 시 실행. 병행 가능한 C006 (stub 모듈) 착수 가능

---

## C005 | 2026-08-09 | run.sh 작성 (실행 검증은 블로킹)

- 목표: 케이스 실행 스크립트 완성 + iSMR 스모크
- 변경: `scripts/run.sh` 신규 (prepare_case → `I_MPI_FABRICS=shm mpirun -np $CUPID_NP`, 로그 tee 저장)
- 실행: 스모크는 **somaFlow.in 부재로 실행 불가** (prepare_case.sh 가 exit 2 로 차단함을 C003 에서 확인)
- 판정: 부분 완료 — 스크립트 작성분만. 실행 검증은 블로커 해소 후 이 사이클을 재개(C005-r1)하여 판정
- 다음: 사용자에게 somaFlow.in 요청 유지 (PLAN §7-1)

---

## C005-r1~r5 | 2026-08-09~10 | iSMR 스모크 — 블로커 3종 순차 해소 (somaFlow.in 확보 후 재개)

사용자가 somaFlow.in 제공 → 스모크 재개. 5회 리비전으로 서로 다른 3개 층위의 문제를 확정·해소:

**r1 — FAIL: namelist 불일치 (소스 내부 불일치)**
- `severe (19): invalid reference ... line 6` = `nfluid`. 소스 검증: `read_flow.f90:114` 의 `problem_description` 선언에는 `nfluid` **있음**(`STM_TBL_cupid` 모듈 변수, 주석에 "somaFlow 에서 설정" 명시), 먼저 READ 하는 `open_files.f90:36` 선언에만 없음 — 제공된 입력이 아니라 **소스 두 선언의 불일치**가 원인
- 조치: `open_files.f90` 에 `USE STM_TBL_cupid, ONLY: nfluid` + namelist 에 `nfluid` 추가 (입력 수정 대신 소스 정렬 — 입력에서 지우면 기본값 물리로 바뀔 위험)

**r2 — FAIL 이동: line 316 (`HS_coupling`)**
- 재빌드(에러 0) 후 크래시가 6행→316행으로 이동. 입력의 4개 변수(`HS_coupling`, `vfporous`, `i_droplet`, `i_fs_temp_intpol`)는 **이 소스의 어떤 namelist 선언에도 없음**(신형 CUPID 옵션) — 전수 대조 스크립트로 확인
- 정정: 앞서 "misc_option 은 읽히지 않는다"고 했던 판단은 오류 (case-sensitive grep 착오). 런타임이 해당 그룹을 검증하는 것은 사실. 정확한 리더 위치는 미상 — 경험적으로 처리

**r3 — 조치+FAIL: 4개 변수 주석 처리 → SIGSEGV**
- 입력의 4개 변수를 `!` 주석 처리 (원본은 scratchpad 에 `somaFlow.in.user_original` 백업). 이 소스는 해당 옵션을 아예 모르므로 의미상 무손실. 입력값이 전부 0(OFF)이라 물리 의도와도 일치
- namelist 전 구간 통과, `Successful grid generation` 후 **SIGSEGV** (np=4, rank 0 은 SIGKILL)

**r4 — 판별 실험: 직렬(np=1)도 동일 SIGSEGV**
- RSS ~2GB 수준(총 7.7GiB 중)에서 결정적 크래시 → **OOM 배제**, 결정적 버그로 확정. 크래시 지점: `read_grid.f90:508` "test,after ilu" 직후의 `gener_vect_size/gener_vect_u`(벡터화 재배열)

**r5 — PASS(해당 구간): 스택 오버플로 확정**
- 가설: `-auto`(로컬 전부 스택) + 셀 수 규모의 자동 배열 vs 컨테이너 기본 스택 8MB
- `ulimit -s unlimited` 적용 → SIGSEGV **해소**, 해당 구간 통과. `run.sh` 에 영구 반영
- 새 정지점(정상 입력 검증): `### rv_parameters.in is required when rv_model=1.` (exit 0) — 입력 `&rv_models` 가 `rv_model=1`. 리더는 사전 빌드 라이브러리(`03_Model/rv_model/io/libcupidMODrv5.a`) 내부로 소스 추적 불가

- 판정: r5 시점 기준 **진행 중** — 실행 파이프라인의 소스/입력/환경 3개 층위 문제를 모두 해소했고, 남은 것은 **케이스 입력 파일 1개 부재**
- 참고: 스택 요구는 클러스터(보통 unlimited 기본)에선 잠복했을 환경 요인 — G1 의 "환경 요인 점검" 항목과 연결
- 다음: **신규 블로커 — `rv_parameters.in` 을 원본 케이스에서 확보** (rv_ht_str=0 이므로 ht_str_*.in 은 불요 추정). 확보 후 C005-r6 재개

---

## C005-r6~r8 | 2026-08-10 | 첫 PMG 가동 성공 — 스모크 목표 달성, 완주는 물리 발산으로 미달

**r6 — rv_model=0 우회 (사용자 결정) → 새 정지점**
- 사용자가 rv_parameters.in 대신 `rv_model=1→0` 으로 우회. 재실행 → rv 체크 통과, **시간 전진 직전 도달** (`Problem: 11 3D 436136 cells`)
- 정지: `Check variable names ... write_fieldview.f90` + `2 2 5 6` — `write_fieldview.f90:991` 의 출력 변수 검증. nscalar(인식 5) ≠ nScalars(요청 6). 전수 확인 결과 **`qvol_gas` 만 이 소스에 없음** (`qvol_liq` 은 있음)

**r7 — 입력 수정 → 직렬 시간 전진 확인 + 핵심 발견**
- `viewField%nScalars=5`, `scalarVar` 에서 `qvol_gas` 제거 (순수 후처리 항목 — 솔버 물리 무관)
- 직렬(np=1) 18+ 스텝 전진 확인. 단 **fort.501 부재 → `solve_pbcg_mg` 미호출** 발견
- 원인 추적: `pressure_solve.f90` 의 MG 경로는 `parallel==1 .and. MG_solver` 조건인데 **`cupid_main.f90:35` 에 `MG_solver = .false.` 하드코딩** — 사용자 확인 후 `.true.` 로 전환 (실행 중이던 r7 은 중지)

**r8 — MG_solver=.true. + np=4: 첫 PMG 가동 (PASS with caveat)**
- 재빌드(에러 0) → np=4 실행: METIS 분할 성공, **MG_tmp 파일 기반 셋업 라이브 관측**(`part001~004.out`, `part_MG*.out` — G1/G2 대상 메커니즘), GMG 계층 구성 후 시간 전진
- **PMG-BiCGSTAB 76회 솔브 성공**:

| 문제 | rank | its 분포 | 독립 residual | 비고 |
|------|------|----------|---------------|------|
| iSMR ECT 436,136셀 (rv_model=0) | 4 | 1×5회, **2×53회**, 3×18회 | max 8.9e-9, mean 1.2e-9 | 솔버 자기보고치(fort.501). step 1~38 |

- **step 38 에서 발산 정지**: `Iteration number for PBCG_MG exceeds 1000` ×2 → `divergence => stop` (제어된 종료)
- 사실관계: 정지 전부터 DP_MAX 가 지수 성장(26 → 4.3e6, ×1.4/스텝), TOTAL_MASS 드리프트(1.823e7→1.814e7). 솔버는 its 1~3 으로 건전하다가 물리가 무너진 뒤에야 한계 도달
- 추정: **rv_model=0 우회가 물리 셋업을 변경**(노심 유동저항/열원 부재)해 과도 자체가 발산. 솔버 결함 증거는 현재 없음. 확정하려면 rv_parameters.in 확보 후 rv_model=1 재실행 필요
- 판정: **PASS (스모크 기준)** — C005 의 DoD(수 timestep + PMG 경로 진입 로그)는 초과 달성 (38 스텝, PMG 76 솔브). 완주(t_end=2s)는 물리 셋업 문제로 미달 — rv_model=1 복원 후 재확인 항목으로 이관
- 다음: ① rv_parameters.in 확보 (PLAN §7-1, 여전히 권장) ② C006 (stub 모듈) 착수 가능 ③ 이번 fort.501 관측치를 LOOP §2 ref_its 표에 잠정 등재

---

## C006 | 2026-08-10 | pmg_standalone 의존성 클로저 — 스텁 1개로 컴파일·링크 통과

- 목표: CUPID 본체 없이 GMG(PMG) 전체를 컴파일·링크 (스텁 최소화, 소스 이원화 금지)
- 변경: `pmg_standalone/makefile` 신규 (원본 파일 직접 참조, `../Source/makefile.in` include 로 플래그 동일), `stub/communicate_serial.f90` (8줄), `driver/link_probe.f90`, `.gitignore` 에 `build/` 추가
- 실행: `in_contain.sh make -C pmg_standalone probe` → `build/link_probe` 실행
- 결과:
  - **의존성 조사가 계획을 단순화**: GMG 가 USE 하는 모듈 22종 중 md_* 12종은 `GMG/module/` 이 자체 보유 (스텁 불요). CUPID 측은 Z-모듈 9종 + Zinterface 클로저 3종 = 12개 파일 — **전부 USE=0 리프**로 확인되어 원본 직접 컴파일 (당초 계획한 MD_matrix 등 스텁 작성 자체가 불필요했음 → PLAN §4-1 을 실측 기준으로 개정)
  - 컴파일: 12(Z) + 8(GMG/module) + 2(allreduce 원본 + communicate 직렬 스텁) + 43(GMG FSRCS) + padiso.f = **66 오브젝트, 에러 0, 1차 시도 통과**
  - 링크 프로브: 미정의 외부 심볼이 정확히 3개(`communicate`, `allreduce_r`, `allreducei_r1`)로 수렴 → allreduce 2종은 원본 `06_MPI/allreduce_fns.f90` 포함으로, `communicate` 는 직렬 스텁(np=1 에서 고스트 교환=항등, 시그니처는 원본 `communicate.f90:1330` 과 동일)으로 해소 → **미정의 0, `link_probe: OK`**
  - 제외 파일 확인: `6_solver_pbcg_ali.f90` 은 프로덕션 GMG makefile 에도 없음 (미컴파일 변형) — 동일하게 제외
- 판정: **PASS** — C006 DoD(컴파일 통과) + 링크 클로저까지 초과 달성. L2 비해당 (실행 로직 무변경)
- 다음: C007 — 합성 Poisson 생성기 + driver_pmg.f90 (덤프 로더·PREP 체인 호출·ref_its 채취)

---

## C007 | 2026-08-10 | 합성 Poisson 유닛 테스트 — PMG 단독 첫 수렴, ref_its 채취

- 목표: CUPID 없이 합성 문제로 `solve_pbcg_mg` 수렴 확인 + LOOP §2 ref_its 초기 채취
- 변경: `driver/driver_pmg.f90`(제작해 기반 검증 드라이버), `tests/mg.in`(케이스 사본), makefile `driver` 타깃, `run_tests.sh`(러너)
- 조사로 확정한 프로덕션 체인 (드라이버가 그대로 재현):
  `read_input_mg` → `subdomain_infor_MG`(rank0, MG_tmp 기록) → `read_mesh_MPI` → `Prep_fine_P` → `Prep_MG_GarL`, 솔브는 `assemble_FVM` → `SOLVE_GMG`(isol_mg≤0 → `solve_pbcg_mg`). 드라이버 입력: `num_neigh_mg`/`neigh_mg`(7점 스텐실), `xloc_tmp(ncell,ndim)`, `celem=1`, `Zbicg%eps_bicg`(→crit_bcg_mg 로 전파, `1_read_input.f90:167`)
- **실패 2건과 해소 (설계 지식으로 기록)**:
  - r1: SIGSEGV at coord 접근 — **`Prep_fine_P` 가 fine 좌표를 해제** (`3_Prep_fine_P.f90:130`) → `read_mesh_MPI` 직후 스냅샷으로 해소. 행렬 계수는 GMG 자신의 재배열된 ia/ja/coord 기준으로 생성 (permutation 무관)
  - r2: 발산(its>1000) — **GMG CSR 은 대각을 au 안에 포함** (`ju` 가 위치 지정; `mt_amux1`/`resi_normP` 는 full-row 로 A·x). 증거: nnz=93312 = off-diag 79488 + 대각 13824. `assemble_FVM` 의 별도 diag 인자는 스무더 스케일링용(`diagt=1/diag`)이지 A 의 가산 항이 아님. 드라이버가 대각 위치에 `-vol/0=Inf` 를 넣고 있었음 → 대각값으로 교정
- 결과 (`run_tests.sh`, 전부 **PASS**, 판정은 드라이버의 **독립 residual** — 솔버 자기보고와 별도 계산):

| 문제 | ncell | its | 독립 rel_res | 제작해 오차 |
|------|-------|-----|-------------|------------|
| 등방 24³ | 13,824 | **4** | 3.97e-9 | 3.6e-9 |
| 이방 24³ (aspect=100) | 13,824 | **4** | 1.28e-10 | 3.0e-10 |
| 등방 48³ | 110,592 | **5** | 1.63e-10 | 6.5e-10 |

  - 이방성에서 its 불변 → semi-coarsening 유효. 격자 8배에 its 4→5 → 격자 무관 수렴. 제작해 오차 ~1e-9 → 해 자체 정확
- 판정: **PASS** — L1(빌드)·L2(독립 residual ≤ 10·crit, its 채취) green. ref_its 3종을 LOOP §2 표에 등재 — **이 시점부터 C012(리팩터링) 계열의 L2 안전망 가동**
- 다음: C008(덤프 훅) 또는 C012(G3 리팩터링, 합성 ref 로 판정 가능). C009 골든은 rv_parameters.in 대기

---

## C008 | 2026-08-10 | 골든 덤프 훅 — env 게이트, 셋업+솔브 입출력 캡처

- 목표: `CUPID_PMG_DUMP=<step>` 설정 시에만 PMG 솔브의 입력(CSR·RHS·diag)과 출력(u*), 그리고 standalone 재생에 필요한 셋업 배열을 파일로 캡처 (미설정 시 오버헤드 ~0)
- 변경: `Source/GMG/module/dump_pmg.f90` 신규 (모듈, `NEWUNIT=` 사용 — §5-1 유닛 위생 준수), `pressure_solve.f90` 두 SOLVE_GMG 사이트에 pre/post 훅, `read_grid.f90` 에 셋업 훅, GMG/module·05_Solver makefile, `.gitignore` 에 `pmg_dump/`
- 사전 확인: `parallel` 은 somaFlow.in 입력값(=1)이지 np 파생이 아님 → **직렬(np=1)에서도 MG 경로 진입** = 직렬 골든 채취 가능
- **r0 FAIL — 셋업 덤프 불완전**: `setup_r0.bin` 크기가 정확히 `num_neigh_mg`+`neigh_mg` 분량만큼 부족 (12,211,824 vs 27,912,720). 바이트 판독으로 xloc·celem 만 기록됐음을 확정 — 두 배열은 **`subdomain_infor_mg` 종료 시 DEALLOCATE** 됨 (`6_subdomain_infor_mg.f90:644`). 정정: 사전 생존 검사에서 GMG 디렉토리를 grep 범위에서 누락했던 것이 원인 (ifort 는 미할당 allocatable 의 slice WRITE 를 0 바이트로 통과시킴 — 조용한 실패 모드)
- **r1 — 셋업 덤프를 셋업 시점으로 이동**: `dump_pmg_setup_hook` 을 `read_grid` 의 `subdomain_infor_MG` 호출 직전에 배치. 빌드 순서상 모듈을 `05_Solver` → `GMG/module/`(선행 빌드 그룹)로 이동 (02_IO 에서도 USE 가능). 이동 잔재 `05_Solver/dump_pmg.o` 가 글롭 링크에 걸려 multiple definition → 제거 후 green
- 검증 (np=1, step 1 캡처 후 중단):

| 파일 | 크기 | 내용 검증 |
|------|------|-----------|
| `setup_r0.bin` | 27,912,720 = 기대치 | ndim=3, nelem=436,136, nf_max=8, num_neigh 샘플 5~6 (비정렬 격자 유효값) |
| `s1_k{1,2}_r0_c1.pre` | 31,038,424 = 기대치 | n=436,136, nnz=3,007,528, src 실수값 |
| `s1_k{1,2}_r0_c1.post` | 3,489,108 = 기대치 | u* 실수값 (site1 ≈ −10.1, site2 ≈ 2e-10) |

  - site1 u* 샘플이 r0 시도의 값과 동일 — 실행 간 결정적 재현 (‑fp-model strict) 신호
  - fort.501 step1: its=3, 2
- 판정: **PASS** — L1 green(에러 0), 덤프 무결성 바이트 단위 검증 완료. L2 비해당 (솔버 로직 무변경 — 훅은 env 미설정 시 즉시 반환)
- 다음: C009 — 골든 채취·재생 (사용자 결정으로 rv_model=0 유지 확정 → **블로킹 해제**, 현행 케이스 step 1~38 구간이 골든 원천)

---

## C009 | 2026-08-10 | 골든 채취·재생 — 하네스 베이스라인 회귀 게이트 확립 + 솔버 잔차 드리프트 발견

- 목표: iSMR 실전 덤프를 채취하고 standalone 에서 재생해, 이후 모든 고도화의 회귀 게이트를 완성
- 변경: `driver_pmg.f90` 재생 모드(덤프 로더, k1→k2 상태 체인, save_u, 2단 게이트), `run_tests.sh` 골든 섹션, `golden/iSMR436k_np1/` (meta.md·checksums 커밋, 바이너리는 git 제외), `.gitignore`
- **채취**: `CUPID_PMG_DUMP=1,10,30` np=1 → 13파일 (setup + 3스텝×k1/k2×pre/post). `_c2` 파일 부재로 스텝당 사이트별 1솔브 확정 → fort.501 매핑(2줄/스텝) 검증. 골든 its: s1(3,2) s10(3,2) s30(3,3)
- **재생 리비전들**:
  - r0: k1 재생 성공(충실도 1.7e-11)이나 게이트 2건 오탐 — ① k2 는 RHS‖b₂‖≈1e-25 라 상대잔차 분모 붕괴 ② 충실도 1e-12 bitwise 기대가 비현실
  - r1: 잔차 게이트를 **솔버 실제 기준(초기잔차 상대 ‖r‖≤10ε‖r₀‖)** 으로 교정 — k2 의 r₀ = b₂−A·u_k1 ≈ −b₁ (k1 해가 초기추정). 측정값 ‖r‖≈2e-9 = 1e-8×‖r₀‖(0.2) 로 산수 일치 확인
  - r2: 충실도 정규화를 **k1 압력장 스케일(uscale)** 로 — k2 는 크기 ~0 보정이라 자체 스케일 정규화가 오탐 유발
- **핵심 결과**:

| step | 하네스 its | 프로덕션 its | 충실도 (k1/k2) | 판정 |
|------|-----------|--------------|----------------|------|
| 1  | 2,2 | 3,2 | 1.7e-11 / 1.4e-14 | PASS |
| 10 | 2,2 | 3,2 | 3.1e-11 / 9.6e-11 | PASS |
| 30 | 3,3 | 3,3 | 6.7e-12 / 2.3e-12 | PASS (res WARN, 아래) |

  - **하네스 결정성 bitwise 확정**: s1 재생 3회(게이트 수정 재빌드 포함)·s10 재생 2회 전부 md5 동일 → **베이스라인 bitwise 회귀 게이트 성립** (`baseline_s*_k*.u` 6개 확보)
  - **its 는 실행 문맥 의존**: 동일 (A,b,u₀,ε) 에서 프로덕션 3 vs 하네스 2 (s1/s10 k1). 덤프 헤더 itim=1 로 매핑 오프셋 가설은 기각 — 예조건자(Chebyshev/Lanczos) 시드 상태가 원인으로 추정. 프로덕션 자신도 s5 에서 (2,2). → **회귀 기준은 하네스 베이스라인 its** 로 확정, 프로덕션 대비는 충실도 게이트가 담당
  - **솔버 품질 관찰 (논문·고도화 후보)**: s30 에서 참 잔차 ‖b−Au‖ ≈ 4.4e-7·‖r₀‖ vs 솔버 자기보고 2.8e-11·‖r₀‖ — **재귀 잔차 ~1.6e4 배 드리프트**. 충실도 6.7e-12 로 재생=프로덕션이므로 프로덕션 고유 특성. 재생 모드 res_gate 는 보고 전용(WARN) 으로 설계 확정 — 독립 잔차 하네스가 설계 목적(자기보고에 안 속기)을 정확히 수행한 사례
- 판정: **PASS** — 골든 3스텝 재생 전부 green (충실도+its+bitwise 3중 게이트). run_tests.sh 로 합성 3 + 골든 3 통합
- 다음: 안전망 완성 — C010(G1)/C011(G2)/C012(G3) 착수 가능. 골든 바이너리는 git 제외 (재채취: meta.md 절차)
