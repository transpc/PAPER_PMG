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

---

## C012-1 | 2026-08-10 | G3 리팩터링 1차 — 컴파일러 기반 미사용 코드 제거 (bitwise 게이트 첫 실전)

- 목표: 논문 핵심 경로 3파일의 미사용 변수·죽은 코드를 제거하되 run_tests 6/6 + bitwise green 유지
- 변경:
  - `-warn unused -syntax-only` 전수 스캔으로 인벤토리 작성 (오브젝트 무변경 검사) → GMG 전체에서 IMPLICIT NONE 공백은 `blaslapack_sub_n.f90`(1996줄) 1개뿐 — 별도 사이클로 이관
  - `6_solver_pbcg_mg.f90`: 미사용 import `ju, alu, crit` + 로컬 `i1` 제거
  - `1_read_input.f90`: 미사용 import `isol, ns, nd_max` 제거 (Znode USE 라인 자체 삭제)
  - `7_SOLVE_GMG.f90`: 20개 라인에서 미사용 36건 제거 (`status(mpi_status_size), tag` 미사용 MPI 상태 변수 3벌, 스무더 쌍둥이 루틴의 중복 미사용 등) + 죽은 주석 호출 1건
  - **`6_solver_pbcg_ali.f90` 삭제** (314줄 — 프로덕션 makefile 미등재 확인, 참조는 주석 1건뿐)
- 과정 기록: 자동 콤마 정리 스크립트가 연속행(`, &`)의 필수 콤마까지 제거해 9개 라인 파손 → 컴파일러 에러 목록 따라 수동 교정 → 3파일 에러 0 + 미사용 remark 0 확인
- 실행: `build.sh`(프로덕션, 에러 0) → `run_tests.sh` 전체
- 결과: **run_tests 6/6 PASS (exit 0)** — 합성 3/3 (its 4/4/5 불변, res_gate 수치 동일), 골든 3/3 (fid 5.3e-15/9.6e-11/2.3e-12, its 22/22/33) + **베이스라인 bitwise(cmp) 전부 통과** = 리팩터링 전과 해가 비트 단위 동일
- 판정: **PASS** — L1 green(빌드 에러 0), L2 green(its ±0), L3 상당(bitwise) green. 총 -317줄 순감 (죽은 파일 314 + 미사용 43건 정리 vs 서식 정리)
- 다음: C012-2 후보 — `blaslapack_sub_n.f90` IMPLICIT NONE 전환(1996줄), 나머지 GMG 파일 미사용 정리(06_solver_pcg_ilu 등), 죽은 주석 블록 정리. 유닛 번호 파라미터화는 G1 에서 (PLAN §5-4)

---

## C013 | 2026-08-11 | 환경 이주 (root→sjdo) — 빌드 복원 + 잠재 인터페이스 버그 발견·수정

- 목표: 새 환경(sjdo 계정, 신규 SIF)에서 빌드·합성 테스트 체계를 green 으로 복원
- 변경: `scripts/env.sh` (SIF=`~/00_apptainer/hpc2023_ubuntu_prc3.3.sif`, METIS_LIB=SIF 내장 `/usr/local/lib`), `Source/GMG/poly_smooth.f90` (아래 r1)
- 실행: `in_contain.sh mpiifort --version` → 2021.10.0 (구 SIF 와 동일) / `build.sh` → cupid.x 에러 0 / `run_tests.sh`
- 결과:
  - 구 환경과의 차이: 계정 root→sjdo, `/root/00_apptainer` 접근 불가, SIF 파일 교체(hpc23→hpc2023_ubuntu_prc3.3, ifort 버전 동일), METIS 는 SIF 내장 libmetis.so 로 대체(METIS_PartGraphKway 확인), CPU 12→20코어(AVX2까지 동일), git 제외 산출물(골든 바이너리·빌드물) 유실
  - **r1 (FAIL→원인 규명)**: standalone fresh 빌드가 `poly_smooth.f90` 에서 컴파일 에러 — `allreduce_r(a,b,n)` (3인자 배열) 를 2인자 스칼라로 호출하는 원본 잠재 버그 4건 (`np>1` 분기 한정 = np=1 하네스에선 죽은 코드). `-warn interface` 검사가 standalone 통합 build/ 에서만 발화 (프로덕션은 디렉토리별 컴파일이라 미발화, 구 머신에선 증분 빌드로 poly_smooth 재컴파일이 없어 미발화). 수정: 4개 호출부를 스칼라 전용 래퍼 `allreduce_r_s` (동일 파일군 내 기존 루틴, 솔버 본체와 동일 관용구) 로 교체
  - 합성 3종: its 4/4/5 (ref 동일), res_gate 수치 C007 과 유효숫자 일치 — 머신 교체 후에도 재현

| 문제 | rank | its (ref) | res_gate | 판정 |
|------|------|-----------|----------|------|
| iso24 | 1 | 4 (4) | 3.9707E-02 | PASS |
| aniso24 | 1 | 4 (4) | 1.2861E-03 | PASS |
| iso48 | 1 | 5 (5) | 1.6315E-03 | PASS |

- 판정: **PASS** — L1 green (프로덕션+standalone 빌드 에러 0), L2 green (합성 its ±0)
- 다음: C014 골든 재채취 (바이너리 유실). poly_smooth 의 np>1 경로는 C010 MPI 하네스에서 실증 필요 (프로덕션 np=4 가 이 버그로 어떻게 동작했는지 확인 대상)

---

## C014 | 2026-08-11 | 골든 재채취 — 머신 간 bitwise 재현 확인 + 재베이스라인

- 목표: 유실된 골든 바이너리를 meta.md 절차로 재채취하고 6/6 게이트를 green 으로 복원
- 변경: `golden/iSMR436k_np1/` (덤프 13개 + baseline 6개 재설치 — git 제외), `checksums.md5`·`fort501_production.txt`·`meta.md` 갱신
- 실행: `CUPID_NP=1 CUPID_PMG_DUMP=1,10,30 run.sh` ×2회 (결정성 검증) → replay ×3 → `run_tests.sh`
- 결과:
  - **C009 체크섬 대조**: setup+s1 4파일 md5 **일치** (다른 CPU·계정·재빌드 바이너리에서 bitwise 재현 — `-fp-model strict`+AVX2 고정의 효과 실증). s10/s30 8파일은 불일치
  - **fort.501 대조**: 76솔브 its 전부 동일, 자기보고 잔차만 step 2 부터 최하위 비트 드리프트
  - **결정성**: 신규 환경 2회 독립 실행 13/13 파일 bitwise 동일 → run-to-run 결정적. 추정: 드리프트는 SIF 교체에 따른 런타임 라이브러리(MKL/libimf) 미세 차이가 비-GMG 물리 구간(스텝 1 솔브 이후)에서 유입 — GMG 산술은 C012-1 bitwise 게이트로 불변 입증된 것과 정합
  - 재생 하네스 its: s1(2,2) s10(2,2) s30(3,3) — C009 기준과 동일. baseline .u 6개 재설치

| 문제 | rank | its (ref) | 충실도 | 판정 |
|------|------|-----------|--------|------|
| gold_s1 | 1 | 2,2 (2,2) | 5.3e-15 | PASS |
| gold_s10 | 1 | 2,2 (2,2) | 9.6e-11 | PASS |
| gold_s30 | 1 | 3,3 (3,3) | 2.3e-12 | PASS |

- 판정: **PASS** — run_tests 6/6 (합성 3 + 골든 3, bitwise 게이트 포함) exit 0. 골든 재베이스라인 사유: 환경 이주 (its 불변, 궤적 비트 드리프트만)
- 다음: C015 베이스라인 태그 + perf_log.md 기준표 → C010 착수

---

## C015 | 2026-08-11 | Phase 4 베이스라인 고정 — 태그 + 성능 기준표

- 목표: 골든 green 상태에 태그를 붙이고 perf_log.md 성능 기준표를 시작 (PLAN §5-4-1)
- 변경: `pmg_standalone/perf_log.md` 신설, git tag `baseline-phase4`
- 실행: 6개 케이스 각각 `/usr/bin/time` 벽시계 채취 (컨테이너 기동 포함 프로세스 시간)
- 결과: 합성 1.01/0.99/3.83 s, 골든 재생 12.35/11.83/12.10 s (환경: sjdo/WSL2 20코어, hpc2023 SIF). 드라이버 내부 솔버 단독 타이밍 훅은 부재 — 향후 도입 항목으로 perf_log.md 에 명시
- 판정: **PASS** — 인프라 사이클 DoD 충족 (기준표 존재 + 태그)
- 다음: C010 착수 — MPI 하네스 확장 (communicate 스텁 → 실물, np=2/4 골든), np=900 재현

---

## C010-1 | 2026-08-11 | MPI 하네스 확장 1단계 — communicate 스텁 → 실물 교체

- 목표: communicate_serial 스텁을 06_MPI/communicate.f90 실물로 교체하되 np=1 run_tests 6/6 + bitwise green 유지
- 변경: `pmg_standalone/makefile` (SSRC 교체), `06_MPI/communicate_allb.f90` 신설 (communicate.f90 에서 verbatim 분리 — CUPID 본체 모듈 VOL_DATA/Zpress/Zvector/Zare 의존이 이 초기화 루틴뿐이라 파일 단위 분리로 하네스 클로저 성립), `06_MPI/makefile` FSRC1 추가
- 실행: `build.sh` (프로덕션, 에러 0) → standalone clean 재빌드 → `run_tests.sh`
- 결과: 6/6 PASS + 베이스라인 bitwise 전부 통과 — 실물 communicate 는 np=1 에서 항등(niut=0), allb 분리는 동작 불변
- 관측: GMG 의 `communicate` 호출은 `06_solver_pcg_ilu.f90` (coarsest 경로) 7건뿐. 주 솔버 `6_solver_pbcg_mg` 는 자체 `communicate_s`, 스무더 계층은 `send_receive*` 사용 — np>1 확장 시 세 경로 모두 실물이 필요하며 이제 클로저에 전부 포함됨
- 판정: **PASS** — L1/L2/L3(bitwise) green
- 다음: C010-2 — 드라이버 MPI 모드 (합성 Poisson np=2 분할 실행, mpirun 하네스)

---

## C010-2 | 2026-08-11 | 드라이버 MPI 확장 — np 스케일 가동 + 수렴 붕괴 재현 (np≥6)

- 목표: 합성 Poisson 을 np>1 로 실행 가능하게 하고 (셋업 fan-out 경로 실전 가동), np 스케일 특성을 관측
- 변경: `driver/driver_pmg.f90` (np 가드 해제, k-슬랩 celem 분할, subdomain_infor_MG rank0 전용 + BARRIER, nnode 로컬 규격 조립, uex 전역 id 복원 방식 — np=1 bitwise 불변 설계, verify/res_norm 전역 리덕션), `run_tests.sh` (iso24_np2/np4 회귀 케이스 추가, ref_its 4/4)
- 실행: run_tests 8/8 + np 스캔 (24³/48³ × np 2~16)
- 결과:
  - np=1 회귀: 6종 전부 PASS, res_gate 수치 불변 + 베이스라인 bitwise green — MPI 확장이 np=1 산술 불변임을 게이트로 입증
  - **np>1 첫 가동**: MG_tmp 파일 fan-out(6_subdomain_infor_mg→2_read_mesh_MPI, G1 대상 기계장치)이 np=2~16 셋업에서 작동, sum(nintf)=ncell 정합

| 격자 | np=2 | np=4 | np=6 | np=8 | np=12 | np=16 |
|------|------|------|------|------|-------|-------|
| 24³ | 4 PASS | 4 PASS | 7 PASS | 773 PASS | **maxit 발산** | - |
| 48³ | 5 PASS | 5 PASS | **발산** | **발산** | **발산** | **발산** |

  - **np 스케일 수렴 붕괴 재현**: 등방 Poisson 인데 np≥6(48³)/np≥8(24³) 에서 its 급증→발산. 격자가 클수록 더 낮은 np 에서 붕괴 (48³ np6 슬랩 8층 발산 vs 24³ np6 4층 수렴 — 슬랩 두께 단독으로는 설명 안 됨)
  - 분리 실험: ioplv=0 + nlevel=4 수동으로도 48³ np6 발산 → 자동 레벨 선택(ioplv=1)이 원인 아님. nlevel=2 는 출력 없이 종료 (별도 실패 양상, 미조사)
  - 추정: 다중 도메인 coarse 계층 구성 혹은 병렬 스무더/코스닝의 np 결합 결함. np=900 "실행 불가"의 소규모 전조일 가능성 — 단, 슬랩 분할은 프로덕션 METIS 분할과 형상이 달라 분할 형상 민감성은 미분리
- 판정: **PASS** — 사이클 목표(np>1 가동+관측) 달성, run_tests 8/8. np≥6 붕괴는 신규 발견으로 백로그 등재
- 다음: C010-3 — ① 분할 형상 분리 (3D 블록 분할 옵션) ② 발산 시 레벨별 구조 덤프 (rank 별 coarse 셀 수, 빈 도메인 여부) ③ np=900 셋업 스모크 (붕괴 원인 확정 후)

---

## C010-3 | 2026-08-12 | 분할 형상 원인 분리 + np=900 크래시 재현·힙 오염 버그 수정

- 목표: np 수렴 붕괴의 분할 형상 의존 분리, np=900 실패 재현과 원인 규명
- 변경: `driver_pmg.f90` (5번째 인자 ipart: 0=k-슬랩, 1=3D 블록 분할 — np 정육면체 인수분해), `2_read_mesh_MPI.f90` (iar1/iai1 할당 크기 수정 — 아래 ③)
- 실행: 블록/슬랩 np 스캔, np=900 스모크 (24³ 블록), `-g -traceback` + `-check bounds,uninit` 빌드 재현, run_tests
- 결과:
  - ① **C010-2 발산의 원인 = 분할 형상**: 3D 블록 분할이면 24³·48³ × np 6/8/12/64 전부 its 4~5 PASS (슬랩 발산 구성 전부 해소). np 구조 결함 가설 기각. 단 48³·96³ np=128 은 블록도 발산 — 셋업 진단(fort.999)은 건강 (coarsest min=max 균등, zero cells 없음, 전역 체인 7200→1008→147→24) → solve 단계 문제, **미해결 경계로 등재**
  - ② **np=900 크래시 재현 성공** (24³ 블록, 15셀/rank): MG_tmp 1801파일 fan-out 은 완료, 판독·셋업에서 비결정적 크래시 — release 실행마다 양상 상이 (SIGSEGV 15개 / EOF(unit 11, part002.out) 3개 등) = **힙 오염 증상**. 프로덕션 np=900 "실행 불가"의 비결정성과 정합하는 메커니즘
  - ③ **원인 확정 1건 (bounds check)**: `2_read_mesh_MPI.f90:313` — `iar1(nnode+1)` 할당인데 `:532` 레벨 루프가 레벨별 `nnode1`(=ialv 차분) 행을 기록. 극소 도메인에서 coarse 레벨 폭 > fine nnode → 정확히 상한+1 오버런 22개 rank (56>55, 46>45, …). **수정**: 레벨 최대 폭 `MAXVAL(ialv 차분)` 로 할당. 게이트 8/8 + bitwise green (산술 불변 입증)
  - ④ 다음 계층 (수정 후 재실행): `3_Prep_fine_P.f90:186` — 어떤 로컬 셀의 `num_neigh=7 > nf_max=6`. 7점 스텐실에서 불가능한 값 = **writer(6_subdomain_infor_mg) 산출물 자체 이상**. 15셀/rank 극단 병리일 수 있어 프로덕션 규모(484셀/rank) 관련성은 격자 확대(48³+ np=900) 후 판단
  - ⑤ 과정 사고: tests/mg.in 에 실험 잔재(nlevel 6/glomax 4/ioplv 0)가 남아 게이트 전면 RED (its 변동+bitwise FAIL) → **잔재 원복만으로 8/8 green, 수치 기준 일치**. 게이트가 구성 오염을 정확히 포착한 사례. 교훈: mg.in 실험은 반드시 즉시 원복 검증까지 한 사이클로
- 판정: **PASS(부분)** — 형상 분리 완료, 실버그 1건 수정·검증. np=128 발산과 ④ 는 미해결
- 다음: C010-4 — ① 48³ np=900 (923셀/rank, 프로덕션 근접) 재현 ② writer 의 num_neigh>6 생성 경위 (극소 도메인 한정인지) ③ np=128 solve 발산 (nlv_glomax 상향 실험 포함)

---

## C010-4 | 2026-08-12 | np=900 완주 green + 발산 원인 = Chebyshev 병렬 고유값 추정으로 확정

- 목표: 프로덕션 근접 규모(48³, 123셀/rank)에서 np=900 재현, 잔여 발산·num_neigh>6 경위 규명
- 변경: `driver_pmg.f90` — 합성 조립의 열 공간을 nnode→**nnodegl**(전역 coarse 확장 포함)로 수정 (mycoord/uex 크기·루프). 대규모 np 에서 ja 가 nnodegl 까지 참조 → 이전엔 범위 밖 읽기로 행렬 오염 (하네스 자체 버그, np≤4 에선 nnodegl=nnode 라 잠복)
- 실행: 48³ np=900 (bounds/release 각 1회), 발산 구성 재스캔, 스무더 판별(GAS), run_tests
- 결과:
  - **np=900 완주 PASS, its=7** — bounds·release 빌드 모두 (결정적). C010-3 의 iar1 수정 + 본 드라이버 수정 후 셋업(1800파일 fan-out+판독+9레벨 계층)·솔브 전부 green. G1 의 "np=900 실행 가능" 이 합성 문제 기준으로 달성
  - 48³/np900 bounds 실행에서 read_mesh_MPI 클린 통과 → C010-3 ④(num_neigh=7, EOF)는 **15셀/rank 극단 병리**로 국한. 123셀/rank 에선 미발생
  - **발산 원인 확정 (스무더 판별)**: 발산 구성(48³ np6 슬랩, 48³ np128 블록)이 `smothing='GAS'` 로는 its=5 즉시 수렴 → **'POL'(Chebyshev) 스무더의 Lanczos eig_max 추정이 특정 분할(이방 도메인)에서 실제 λmax 를 하회 → 고주파 증폭 발산**. 소스 주석도 병렬 시 추정 차이를 인지 (`poly_smooth.f90` "in parallel ... may difference"). np 비단조성(np128 발산 vs np900 수렴)은 블록 인수분해 이방성(4×4×8 vs 9×10×10) 차이와 정합
  - run_tests 8/8 green. iso24_np4 res_gate 1.51e-2→2.37e-2 변화는 이전 조립이 범위 밖 garbage 를 읽던 것의 교정 (its ±0, np=1 계열은 bitwise 불변)
- 판정: **PASS** — G1 1차 목표(재현→원인→수정→np=900 green) 달성. 프로덕션 케이스 np=900 최종 검증은 클러스터 확보 시
- 다음: C010-5 후보 — ① 프로덕션 np=900 실증 (클러스터) ② Chebyshev eig 추정 강건화 (안전 계수/추정 실패 감지 → 논문 고도화 항목) ③ 15셀/rank 극단의 reader 정렬 문제 (우선순위 낮음) ④ np≥1000 파일명 인코딩 (PLAN 기지)

---

## C016 | 2026-08-12 | 상류 somaFlow.in 교체(ECT1) 수신 — 호환성 검증 + 신입력 스모크 재베이스라인

- 목표: 상류 커밋 d5e2177(somaFlow.in 전면 교체, blob d145c81→29ee238)이 기존 회귀 체계에 미치는 영향 확정 — 신입력 스모크로 새 ref_its 채취, max_bicg 10000→1000 축소에 따른 조기정지 여부 판정
- 변경: 소스·하네스 변경 없음 (입력만 상류에서 교체됨). 입력 요지: 15 MPa 고압 정지 → 101.3 kPa ECT1 실증 구성 (비등 열원 q0_liq=0.918e6/0.448e6, 벽열유속 1.44e6, iturb=2, mdrag=1, iheatpart=4, mHTC/mtopol=2, vfporous=1, dt=1e-7, cfl_ratio=0.05), max_bicg 10000→1000, lev/lev_type/levmpi/levmpi_type·nfluid·&smr_models 삭제
- 실행: 사전 소스 대조 + ifort namelist 관용성 테스트(컨테이너), `CUPID_NP=4 scripts/run.sh` (step 48 에서 수동 중단, 표본 97솔브 — 구 fort.501 은 fort.501.d145c81.bak 백업)
- 결과:
  - **호환성 green**: 전 네임리스트 정상 읽힘·기동 확인 (`Problem: ECT1`, PMG 경로 진입). lev 4종 삭제는 기본값 0(`read_flow.f90:756`)과 동일 = 기능 불변. nfluid 기본값 1 = 구값. 신규 키(imp_boron_trans, eps_imp_boron, max_iter_boron, rv_mcp, rv_valve)는 소스 네임리스트에 존재. `&smr_models` 는 소스가 읽지 않아 삭제 무해. vfporous/i_droplet/i_fs_temp_intpol 은 현 소스 misc_option 에 선언되어 있어 **C005-r3식 입력 수술 불필요**. `iprn=1.d0`(정수 변수에 실수 리터럴)은 ifort 관용 통과 (iostat=0, iprn=1 — 단독 재현 테스트로 확인)
  - **스모크 (np=4, step 1~48, 97솔브)**: 크래시·발산 없음. PMG its 분포 **1~7 (최빈 2: 31회)**, 전 솔브 수렴 잔차 ≤ 9.7e-9 (eps_bicg=1e-8 충족), 스텝당 k1+k2 2솔브 패턴 유지
  - **max_bicg=1000 조기정지 없음** — 관측 최대 its 7, 여유 ~143배
  - 구입력 대비: its 상한 3→7 (물리 모델 대거 활성화와 정합). 구입력의 step 38 발산은 신입력에서 미재현 (48스텝 정상 진행, DP_MAX ~9.2e3 상승 추세는 초기 과도 관찰만 — 완주 여부는 미확인, 우리가 중단)
  - 골든 재생 회귀(run_tests.sh)는 디스크 덤프 재생이라 **영향 없음**. 단 골든 재채취·프로덕션 대비 절차는 구입력에 결박 → meta.md 에 입력 버전 고정(provenance) 명시
- 판정: **PASS** — L1 green (기동·48스텝 무크래시) + 신입력 ref_its 채취 완료. 유의점 2건(max_bicg 조기정지 / 구 스모크 기준 무효화) 모두 판정·문서화
- 다음: ① 신입력 완주 확인 (t_end 20s — 장시간, 필요 시 별도 사이클) ② 논문 케이스 구성 확정 후 신입력 기준 골든 재채취 여부 결정 ③ C010-5 잔여 (프로덕션 np=900 실증)

---

## C017 | 2026-08-12 | ECT1 골든 재채취 + 안정성 러닝(~2s) — 사용자 결정: 20s 완주 대신 2s

- 목표: (사용자 지시) 20s 완주는 과함 — ~2s 까지만 진행해 안정성을 확인하고 ECT1 신입력 기준 골든을 재채취, run_tests 병행 게이트로 등재
- 변경: `golden/ECT1_436k_np1/` 신설 (meta.md·checksums.md5·fort501_production.txt 커밋, 덤프 17개+베이스라인 8개는 git 제외), `driver_pmg.f90` replay 4번째 인자 `fid_gate` (기본 1e-9 — 기존 동작 불변), `run_tests.sh` 골든 섹션을 `run_golden` 함수로 일반화 + ECT1 세트 4스텝 등재, LOOP §2 표 갱신
- 실행: 채취 np=1 ×3회 (s1/10/30 ×2 결정성 + s150 추가 채취), replay ×(4스텝×2회), 안정성 러닝 np=4 (t=0→2.23s), `run_tests.sh`
- 결과:
  - **채취 결정성**: 2회 독립 채취 13/13 md5 bitwise + fort.501 60솔브 동일. s150 런(3번째)도 공유 구간(setup+s1/10/30) 8/8 md5 동일 — 세 런이 bitwise 계보로 연결
  - **프로덕션 its (np=1)**: s1(3,1) s10(2,1) s30(4,3) s150(11,9) — 정착 상태가 초기 과도보다 어려움
  - **리플레이 (하네스 기준)**: its s1(1,1) s10(1,1) s30(3,3) s150(9,10), 전부 VERDICT PASS + 베이스라인 bitwise 성립(2회차 cmp). **충실도 관측 최대 2.4e-7 (s10 k1)** — its 문맥 의존(C009 기전)의 ECT1 발현이 구골든(≤1e-10)보다 큼 → 세트별 fid_gate 도입: 구골든 1e-9 유지, ECT1 1e-6 (관측 ×4 여유). 정밀 회귀는 bitwise 게이트 담당이므로 감지력 손실 없음
  - **안정성 러닝 (np=4, t 0→2.23s, 157스텝, 벽시계 ~9분)**: 크래시·발산 없음. DP_MAX 1.7e4 (peak, t≈0.1s) → ~1 (t>1.5s) **정착**, TOTAL_MASS 불변(1.814e7). dt 1e-7→5.6e-2 성장(상한 0.1 미도달). its 초기 1~7 → dt 포화 구간 8~11 (최빈 9), max 11 ≪ max_bicg=1000. 전 320솔브 잔차 ≤9.7e-9
  - **시간 결정 (사용자 위임)**: **2s 채택** — t≈0.7s 이후 DP_MAX 정착으로 과도 통과가 확인되고 벽시계 ~9분이라 계산 시간 지배 없음. (참고: dt 포화 추세면 20s 완주도 wall ~15~25분 추정 — 필요 시 저비용)
  - **run_tests 12/12 green**: 합성 5 + 구골든 3 (bitwise 불변 — 드라이버 fid_gate 추가가 산술 무영향 입증) + ECT1 골든 4
- 판정: **PASS** — 골든 이원화 완성 (구입력 3스텝 + ECT1 4스텝 병행 게이트), 안정성·시간 결정 완료
- 다음: ① perf_log.md 에 ect1 재생 케이스 시간 행 추가 (선택) ② 논문 케이스 확정 시 구골든 은퇴 여부 결정 ③ C010-5 잔여 (프로덕션 np=900 실증)

---

## C018 | 2026-08-14 | 문서 현행화 (PLAN 잔여 stale 4지점) + perf_log ECT1 행

- 목표: C016 검토에서 확인된 PLAN.md 의 stale 서술을 현행화하고 C017 잔여(perf 행)를 마감
- 변경: `PLAN.md` — §1-3 somaFlow.in 표(블로커→이력), §7 두 해소 항목(step 38 발산 구입력 한정 / C005-r3 입력 수술 무효), §4-4 골든 입력 결박·세트별 fid_gate 조항, §8 Phase 4 체크리스트(C010/C015/C016/C017 반영). `perf_log.md` — ect1 재생 4행
- 실행: ect1 재생 4케이스 `/usr/bin/time` 채취
- 결과: ect1_s1/s10/s30/s150 = 11.08/13.89/11.51/16.23 s (구골든 12s 급과 동급 — s150 은 its 9,10 이라 +30% 수준). its 전부 베이스라인과 일치
- 판정: **PASS** — 문서 사이클 (L1/L2 비해당), 잔여 stale 서술 소거
- 다음: C011 (G2 — MG_tmp 파일 왕복 → MPI 통신 대체) 착수

---

## C011-1 | 2026-08-14 | G2 설계 — 파일 왕복 인벤토리 + ASCII 라운딩 실험 → 이중 모드 전략 확정

- 목표: MG_tmp 파일 경유 분배의 전 데이터 흐름을 인벤토리하고, 통신 대체의 bitwise 게이트 전략을 실험으로 확정
- 변경: 소스 무변경 (조사·설계 사이클). PLAN §5-2 검증 전략 갱신
- 실행: writer/reader 정독 (탐색 에이전트 + 검증), ifort list-directed 왕복 실험 (컨테이너, 10만 표본 3스케일)
- 결과 (인벤토리 핵심 — 파일:라인 상세는 조사 기록):
  - 파일 3종: `part###.out`(finest, prc당 1개) / `part_MG###.out`(coarse 전 레벨 + A_GC) / `PMG_infor`(메타 공유 1개). 전부 **list-directed ASCII**, 일회성 fan-out (셋업 후 재사용 없음 — 전 소스 grep 로 확정)
  - writer(rank0 전용, `read_grid.f90:253`): 전역 배열이 이미 `(np,·)` 선행차원으로 슬라이스됨 → Scatterv 사상 자연스러움. coarse fan-out 은 **np개 파일 동시 OPEN**(`6_subdomain_infor_mg.f90:229-232`) — np=900 이면 fd 900개 (ulimit 1024 근접, G1 잔여 위험). `ioplv=1` 재진입(GOTO 500)으로 fan-out 최대 2회
  - reader(전 rank, `read_grid.f90:569`): allocate 크기가 파일 헤더 값 의존 → **2-phase(메타→allocate→페이로드) 프로토콜 강제**. PMG_infor 는 rank마다 앞 블록 더미 스킵 — **O(np²) 총 I/O**
  - 유일한 rank 간 중복 데이터 = `A_GC`(전역 coarsest 행렬, 전 rank 동일) → Bcast 대상. 나머지는 disjoint
  - writer/reader 사이 명시적 BARRIER 없음 — 사이 collective(`read_grid.f90:346` 등)가 우연히 동기화 (하네스는 명시 BARRIER 로 이미 문서화)
  - **ASCII 라운딩 실험 (설계 결정 근거)**: REAL(8) 10만 표본(1e-3/1e0/1e6 스케일) list-directed 왕복 → **67,066/100,000 bitwise 불일치** (16자리 출력 vs 왕복 보존에 17자리 필요). 내부 WRITE/READ 문자열 왕복도 정확히 동일 라운딩(67,066 일치) → **in-memory 라운딩 shim 으로 파일 모드와 bitwise 동일한 통신 모드 구성 가능**
- 확정 전략 (증분 로드맵, 각 단계 bitwise 게이트):
  1. C011-2: `isetup_comm` 스위치 (mg.in `&MG_MPI`, 기본 0=파일 — 기존 동작 불변) + **PMG_infor 통신화** (정수뿐 → 라운딩 무관, O(np²) 스킵 소거)
  2. C011-3: `part###.out`(finest) 2-phase Scatterv — REAL(8)(coord)은 문자열 왕복 shim 으로 파일 모드와 bitwise 유지
  3. C011-4: `part_MG###.out` (레벨 루프 + A_GC 는 Bcast) — Xintp/Xrest/coord1 shim
  4. C011-5: **shim 제거 (정확값 전달)** — 의도적 수치 개선으로 분리, its ±0 + 골든 재베이스라인 (PLAN §4-4 셋째 행 절차)
  5. C011-6: 파일 경로·MG_tmp 삭제 (통신 모드 default 승격)
- 판정: **PASS** — 설계 산출물(로드맵 5단계 + 게이트 전략) 확정, 실험 근거 확보
- 다음: C011-2 착수 — 선행 확인: PMG_infor 대상 배열(iintf/inodegl/inbdc/inmax/ialv_P/nnz*)의 모듈 상주 여부 (rank0 writer→reader 구간 생존 조건)

---

## C011-2 | 2026-08-14 | G2 증분 ① — isetup_comm 스위치 + PMG_infor 통신화 (BCAST)

- 목표: 이중 모드 스위치를 도입하고 PMG_infor(정수 메타)를 첫 통신화 — 파일 모드 기본 불변 + 통신 모드 bitwise 동일
- 변경:
  - `GMG/module/MD_MG_index.f90`: `isetup_comm`(0=파일/1=통신) + 스테이징 배열 8종 (`stg_*` — PMG_infor 내용과 1:1)
  - `GMG/1_read_input.f90`: `&MG_MPI` 네임리스트에 `isetup_comm` 추가, READ 전 기본값 0 (기존 mg.in 무수정 호환)
  - `GMG/6_subdomain_infor_mg.f90`(writer): PMG_infor 기록 분기 — 통신 모드는 파일 대신 스테이징 적재 (GOTO 500 재진입 대비 재할당 가드). 선행 확인: 대상 메타는 writer 지역 배열이라 모듈 승격이 필수였음 (inmax 만 MD_MG_coord 상주)
  - `GMG/2_read_mesh_MPI.f90`(reader): 통신 모드 분기 — 비루트 스테이징 할당 → `MPI_BCAST` 8건(전부 INTEGER, ndom>1 가드) → 자기 열(myrank+1) 추출 → 스테이징 해제. 파일 모드의 O(np) 더미 스킵·공유 파일 동시 OPEN 이 통신 모드에서 소거
- 실행: `build.sh`(에러 0) + standalone make → ① 파일 모드(기본) run_tests ② tests/mg.in `isetup_comm=1` 임시 설정 후 합성 np=1/2/4 + 골든 재생 4종 ③ mg.in 원복 후 run_tests 재확인
- 결과:
  - **파일 모드 12/12 green** (bitwise 게이트 포함) — 이중 모드 도입이 기존 경로 무변경임을 입증
  - **통신 모드 실증**: 합성 iso24 np=1/2/4 its 4/4/4, res_gate 유효숫자 일치. **MG_tmp 에 PMG_infor 미생성** (part 파일만 존재) = 통신 경로 실행 물증
  - **통신 모드 bitwise**: ect1_s1/s30/s150 + gold_s30 재생 전부 VERDICT PASS + 베이스라인 **cmp bitwise 동일** — 정수 메타 통신화는 수치 무영향 (설계 예측과 일치)
  - mg.in 원복 검증: run_tests 12/12 green (C010-3 ⑤ 교훈 절차 준수)
- 판정: **PASS** — L1(빌드)·L2(its ±0)·L3 상당(bitwise) green, 증분 ① 완료
- 다음: C011-3 — `part###.out`(finest) 2-phase Scatterv, REAL(8) coord 는 문자열 왕복 shim (C011-1 로드맵)

---

## C011-3 | 2026-08-14 | G2 증분 ② — finest Scatterv + 게이트 red 가 실버그 2건을 적발 (y0 힙 가비지 = C009 "시드" 정체)

- 목표: `part###.out`(finest)을 2-phase Scatterv 로 통신화, 파일 모드 대비 bitwise 유지
- 변경 (transport): `MD_MG_index` — finest 스테이징(stg_fibuf/frbuf/ficnt/frcnt) + `rt_ascii`(문자열 왕복 shim, CONTAINS 함수). `6_subdomain_infor_mg` — 카운트 선패스 + prc 루프 pack 분기 (coord 는 rt_ascii). `2_read_mesh_MPI` — 카운트 SCATTER → 페이로드 SCATTERV(INT/REAL 2스트림) → READ 순서 그대로 커서 unpack + **커서 정합 자기 검증**(kci/kcr==수신 길이 아니면 STOP)
- **게이트 red → 원인 규명 (r1)**:
  - 통신 모드 res_gate 가 실행 방식에 따라 이동 (np1 direct 3.9707e-2 = 파일과 일치 vs np1 mpirun 3.8457e-2), 골든 재생은 **k1 bitwise 일치·k2 만 불일치** ×4종 — 데이터 오류로 보기엔 패턴이 이상함
  - **전송 충실도 직접 입증**: 임시 계측으로 finest 전 배열(헤더·이웃·coord·SR×3) 체크섬 대조 → np=4 파일/통신 **완전 동일**. shim 자체도 3값/줄 왕복 실험 0/90,000 불일치 → transport 무결 확정, 원인은 솔버 측 레이아웃 의존
  - `-init=snan,arrays -fpe0 -traceback` 디버그 빌드로 추적 → **실버그 ①**: `6_solver_pbcg_mg.f90` 의 예조건 상태 벡터 `y0(n)` 이 할당 직후(첫 외부 반복) `u = y0` 로 읽힘 — **힙 가비지가 MG 예조건 내부 솔브의 초기 추정**. k2 솔브는 힙 재사용으로 직전 k1 의 y0 잔재를 물려받는 비공식 warm start 였음 = **C009 "예조건자 시드 실행 문맥 의존"의 정체** (C014 의 SIF 교체 후 잔차 최하위 비트 드리프트와도 정합). **실버그 ②**: `poly_smooth` Lanczos 워크스페이스 v/w/v_old 고스트 구간 미초기화 (동일 계열, 방어적 수정)
  - 수정: y0/z0·v/w/v_old 할당 직후 0 초기화 (warm start 는 2회차부터 그대로). snan 빌드 **완주·트랩 0** = 이 경로 미초기화 REAL 읽기 전멸
- **재베이스라인 (PLAN §4-4 의도적 개선 절차 — 사유: 비결정성 제거)**: k1 전부 bitwise 불변, k2 만 변화. its 유일 변화 = ect1_s150 k2 10→9 (개선). **프로덕션 대비 k2 충실도 대폭 개선**: ect1 s1 6.1e-11→**4.6e-16**, s10 1.3e-10→6.3e-13, s30 2.2e-10→2.4e-13, s150 6.1e-8→4.3e-13 (k1 충실도는 전부 불변 — C009 관찰과 정합: 프로덕션 k2 도 사실상 0 근방 시드였던 것). 구골든 3스텝 its 불변. 베이스라인 14파일 재설치 + checksums 갱신 + run_tests s150 ref (9,9)
- **최종 게이트**: 파일 모드 12/12 green. **파일/통신 × direct/mpirun 4조합 res_gate 완전 일치** (레이아웃 의존 소거 입증), 통신 모드 골든 재생 3종 신 베이스라인 bitwise 일치, 통신 모드에서 part###.out·PMG_infor 미생성 물증, mg.in 원복 검증, 프로덕션 빌드 에러 0
- 판정: **PASS** — 증분 ② 완료 + 부수 성과(솔버 결정화·실버그 2건). G2 bitwise 안전망이 이후 증분에서도 성립하는 기반 확보
- 다음: C011-4 — `part_MG###.out`(coarse 레벨 + A_GC Bcast) 통신화. 잔여 유의: 골든 .pre/.post 는 구 바이너리 프로덕션 산출물 (충실도 게이트 green 이므로 유지, 다음 골든 재채취 시 신 바이너리 기준 갱신 권고)
