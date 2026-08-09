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
