# CUPID PMG 솔버 빌드·실행·검증 작업 계획서

> 작성일: 2026-08-09
> 목적: PMG 예조건 BiCGSTAB 논문(`논문 초안 - PMG 예조건화 BiCGSTAB.md`)의 실험 기반 마련
> 대상: `code/Source` (CUPID 소스), `code/2_iSMR_ECT_res1` (iSMR ECT 자연순환 케이스)

---

## 0. 목표와 산출물

| # | 목표 | 산출물 |
|---|------|--------|
| 1 | 기존 Apptainer 이미지(`hpc23.sif`) 재사용 빌드 환경 (신규 빌드 없음) | SIF 연동 검증 결과 (§2-1) |
| 2 | 컨테이너 안팎을 잇는 스크립트 체계 | `scripts/env.sh`, `scripts/in_contain.sh`, `scripts/build.sh`, `scripts/run.sh`, `scripts/prepare_case.sh` |
| 3 | CUPID 전체 빌드 및 iSMR 케이스 스모크 실행 | `Source/cupid.x`, 케이스 실행 로그 |
| 4 | PMG 솔버(압력 수정 방정식) 단독 테스트 하네스 | `pmg_standalone/` (드라이버 + 스텁 모듈 + 테스트) |
| 5 | 골든 회귀 테스트 체계 | `pmg_standalone/golden/` (덤프 데이터 + 기준 출력), `run_tests.sh` |
| 6 | 위 안전망 위에서의 PMG 코드 고도화 — ① 파일 기반 셋업의 np=900 동작 수정 ② file out의 MPI 통신 대체 ③ 코드 간소화·가독성 (§5) | `Source/GMG` 개선 커밋들 + 성능 트래킹 기록 |

### 최종 디렉토리 구조 (목표)

```
code/
├── PLAN.md                      # 이 문서
├── Source/                      # CUPID 원본 소스 (PMG 고도화도 여기서 직접 수행)
├── 2_iSMR_ECT_res1/             # iSMR ECT 실행 케이스
├── scripts/
│   ├── env.sh                   # SIF·METIS_LIB 등 공통 환경변수 (source 용)
│   ├── in_contain.sh            # 임의 명령을 컨테이너 안에서 실행하는 래퍼
│   ├── build.sh                 # 컨테이너 안에서 CUPID 컴파일
│   ├── run.sh                   # 컨테이너 안에서 케이스 실행
│   └── prepare_case.sh          # 케이스 전처리 (7z 해제 등)
└── pmg_standalone/
    ├── stub/                    # CUPID 본체 모듈의 최소 스텁 (MD_matrix 등)
    ├── driver/                  # 테스트 드라이버 (덤프 로더 + 비교기)
    ├── golden/                  # 골든 덤프 + 기준 출력 (대용량은 git 제외)
    ├── makefile
    └── run_tests.sh             # 전체 테스트 러너 (합성 + 골든 회귀)
```

---

## 1. 현황 조사 결과 (as-is)

### 1-1. 호스트 환경

| 항목 | 값 | 시사점 |
|------|-----|--------|
| OS | Ubuntu 24.04, WSL2 | Apptainer는 root 실행이라 `--fakeroot` 불필요 |
| CPU | 12 코어, **AVX2까지만 지원 (AVX512 없음)** | `makefile.in`이 현재 Skylake/AVX512 설정 → **반드시 AVX2로 변경**, 아니면 illegal instruction |
| Apptainer | 1.5.1 설치됨 | 그대로 사용 |
| SIF | `/root/00_apptainer/hpc23.sif` (5.2 GB) **기확보** | oneAPI 2023.2.1 — ifort classic 2021.10.0 + Intel MPI 2021.10, make/git 포함. **신규 컨테이너 빌드 불필요** |
| METIS | 호스트 `/root/00_apptainer/metis-5.0.2/build/Linux-x86_64/libmetis/libmetis.so` 빌드 완료 | 컨테이너 내부엔 METIS 없음 → `/root` 자동 바인드로 컨테이너에서 접근 확인됨. `LIBDIR` 조정 필요 (§3-1) |
| 7z | 호스트에 p7zip-full 설치됨. **컨테이너엔 없음** | `foamGrid.7z` 해제는 호스트에서 수행 |

### 1-2. 빌드 시스템

- 툴체인 요구사항 (`Source/README.txt`, `makefile.in`): **Intel Fortran classic (`mpiifort`)** + MPI + **MKL**(`-mkl`, Pardiso 사용) + **METIS 5.0.2** (`/usr/local/metis-5.0.2/lib -lmetis` 경로 하드코딩)
- 핵심 플래그: `-r8 -fpp -Dmpi_flag -Dmetis_flag -fp-model strict -mcmodel=medium`
  - **`-fp-model strict`는 골든 회귀에 유리** — 동일 플래그·동일 연산 순서면 bitwise 재현 기대 가능
- 빌드 흐름: `Source/`에서 `make` → 서브디렉토리별 병렬 make → 최상위에서 링크 (`makefile:105`) → `Source/cupid.x` 생성
- `makefile.in`의 `p=52`(병렬 잡 수)는 호스트에 맞게 12로 조정 필요
- `ifort` classic은 oneAPI 2025.0에서 제거되었으나, **기확보된 `hpc23.sif`가 oneAPI 2023.2.1(ifort classic 포함)이므로 조건 충족** — `ifx` 전환 불필요
- 단, `makefile.in`이 하드코딩한 `/usr/local/metis-5.0.2`는 컨테이너 안에 존재하지 않음 → `makefile.in.apptainer`에서 호스트 METIS 경로로 교체 (§3-1)

### 1-3. 실행 케이스 `2_iSMR_ECT_res1`

| 파일 | 상태 | 비고 |
|------|------|------|
| `mg.in` | 있음 | `isol_mg = -2` → **MG as preconditioner** (논문의 PMG 예조건 BiCGSTAB 모드) |
| `foamGrid.7z` | 있음 (3 MB) | 내부에 `foamGrid.in` (112 MB) 1개 → 실행 전 해제 필요 |
| `tpfh2o` | 있음 | 증기표 |
| `somaFlow.in` | **없음 (블로커)** | T/H 입력. 없으면 CUPID가 즉시 종료 (`open_files.f90:74`). **사용자 제공 필요** |

실행 방법: 케이스 디렉토리에서 `mpirun -np <N> ./cupid.x`

### 1-4. PMG 솔버 코드 위치와 의존성

- **본체**: `Source/GMG/` — 솔버 `6_solver_pbcg_mg.f90`(`solve_pbcg_mg`), V-cycle `7_SOLVE_GMG.f90`, 셋업 체인 `1_read_input` → `2_Prep_fine_*` → `3_Prep_fine_P` → `4_Prep_MG_GarL` → `5_PREP_GMG*`, 스무더(`poly_smooth`, `Relax_GS*`), coarsening(`coarsening_semi*`) 등 약 40개 파일
- **전용 모듈**: `Source/GMG/module/MD_MG_*.f90` (index, matrix, coord), `MD_MPI.f90`, `MD_geometry.f90`, `MD_matrix.f90`, `MD_parameter.f90`, `MD_connectivity.f90`
- **호출부**: `Source/05_Solver/pressure_solve.f90:185,256`에서 `SOLVE_GMG` 호출 (압력 수정 방정식)
- **CUPID 본체와의 인터페이스** (추출 시 스텁으로 대체할 대상):
  - 행렬: CSR 형식 — `MD_matrix`의 `nnz, ia, ja, ju, au, u, b, alu`
  - 기하: `MD_geometry`의 `nnode, nelem, coord, num_neigh_mg, neigh_mg, imap`
  - 병렬: `MD_MPI`의 `nintf, myrank` 등
  - 본체 모듈: `Zbicg, Zcore, Znode, Zparam` (`00_Module/`) — `ONLY:`로 소수 변수만 참조
- 입력: `mg.in` namelist (레벨 수, coarsening, 스무딩, coarsest 솔버, MPI/OpenMP 옵션)
- **셋업의 파일 기반 데이터 교환**: root 가 `MG_tmp/` 아래 rank별 파일(`part###.out`, `part_MG###.out`, 공유 `PMG_infor`)을 기록(`6_subdomain_infor_mg.f90:131`, unit `50+prc`)하고 각 rank 가 읽어(`2_read_mesh_MPI.f90:89`, unit `10+myrank`) 셋업을 구성. 폴더는 `call system('mkdir MG_tmp')`(`2_Prep_fine_global.f90:83`)로 생성. **고도화 ①·②의 대상** (§0-6, §5)

---

## 2. Phase 1 — 기존 Apptainer 이미지 연동과 스크립트 체계

### 2-1. 기존 SIF 재사용 (신규 컨테이너 빌드 없음)

컨테이너는 새로 만들지 않는다. 기확보된 `/root/00_apptainer/hpc23.sif`를 사용하되, **체계마다 경로가 바뀔 수 있으므로 `env.sh`에 `SIF` 기본값으로만 정의**하고 외부에서 `export SIF=...`로 재정의할 수 있게 한다.

`hpc23.sif` 내용 검증 결과 (2026-08-09 확인):

| 항목 | 상태 |
|------|------|
| ifort classic | 2021.10.0 (oneAPI 2023.2.1) — `%environment` 설정되어 있어 `apptainer exec`만으로 사용 가능 |
| Intel MPI | 2021.10 (`mpiifort`, `mpirun` 즉시 사용 가능) |
| make / git | 있음 |
| MKL | oneAPI 번들 (`-mkl` 링크 가능) |
| METIS | **컨테이너 내부에 없음** → 호스트 빌드본 사용 (아래) |
| 7z | **없음** → `foamGrid.7z` 해제는 호스트에서 (`prepare_case.sh`) |

METIS는 호스트 `/root/00_apptainer/metis-5.0.2/build/Linux-x86_64/libmetis/libmetis.so`(빌드 완료, shared 전용)를 사용한다. Apptainer가 `/root`(홈)를 자동 바인드하므로 컨테이너 안에서 이 경로에 접근 가능함을 확인했다. `makefile.in`의 `/usr/local/metis-5.0.2` 하드코딩과는 다른 경로이므로 `makefile.in.apptainer`에서 `METIS_LIB` 환경변수 기반으로 교체한다 (§3-1).

### 2-2. 스크립트 계약 (인터페이스 정의)

**`scripts/env.sh`** — 모든 스크립트가 source 하는 공통 환경. SIF는 환경변수로 주입받되 기본값 제공:

```bash
# SIF 는 외부에서 export SIF=... 로 재정의 가능 (체계마다 경로가 다를 수 있음)
CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export SIF="${SIF:-/root/00_apptainer/hpc23.sif}"
export METIS_LIB="${METIS_LIB:-/root/00_apptainer/metis-5.0.2/build/Linux-x86_64/libmetis}"
export CUPID_SRC="$CODE_DIR/Source"
export CUPID_CASE="${CUPID_CASE:-$CODE_DIR/2_iSMR_ECT_res1}"
export CUPID_NP="${CUPID_NP:-4}"        # mpirun rank 수
export MAKE_JOBS="${MAKE_JOBS:-12}"     # 병렬 make 잡 수
```

**`scripts/in_contain.sh`** — 임의 명령을 컨테이너 안에서 실행하는 유일한 통로 (빌드·실행·테스트 전부 이 래퍼 경유):

```bash
#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/env.sh"
[[ -f "$SIF" ]] || { echo "SIF not found: $SIF — env.sh 기본값 확인 또는 export SIF=<경로>" >&2; exit 1; }
exec apptainer exec --bind "$CODE_DIR" "$SIF" "$@"
```

바인드 참고: `/root` 자동 바인드로 리포·METIS 모두 접근 가능하지만, 리포가 `/root` 밖에 있는 체계를 대비해 `--bind "$CODE_DIR"`는 유지한다.

**`scripts/build.sh`** — `in_contain.sh` 기반 컴파일 (Phase 2에서 상세)

**`scripts/run.sh`** — `in_contain.sh` 기반 케이스 실행 (Phase 2에서 상세)

### 2-3. 완료 기준 (DoD)

- [ ] `./scripts/in_contain.sh mpiifort --version` → `ifort (IFORT) 2021.10.0` 출력 (사전 검증 완료 — 스크립트 경유로 재확인)
- [ ] `./scripts/in_contain.sh bash -c 'ls "$METIS_LIB"'` 에 `libmetis.so` 존재
- [ ] 잘못된 `SIF` 경로로 실행 시 명확한 에러 메시지 출력

---

## 3. Phase 2 — 컴파일 및 실행

### 3-1. 컴파일 (`scripts/build.sh`)

1. **`makefile.in` 원본 보존**: 컨테이너용 설정을 `Source/makefile.in.apptainer`로 별도 작성하고, `build.sh`가 빌드 시 이를 `makefile.in`으로 복사(또는 백업 후 교체). 변경점:
   - Vector 블록: Skylake/AVX512 → **Broadwell/AVX2** (`VECTOR=core-AVX2`, `ALIGN=AVX2`, `-align array32byte`)
   - `p=52` → `p=12`
   - `LIBDIR = -mkl -L$(METIS_LIB) -lmetis -Wl,-rpath,$(METIS_LIB)` — 컨테이너에 METIS가 없으므로 호스트 빌드본 경로를 `env.sh`의 `METIS_LIB`로 주입 (make는 환경변수를 상속). `libmetis.so`(shared 전용)이므로 **rpath 필수** — 없으면 실행 때마다 `LD_LIBRARY_PATH` 필요. `-mkl`은 ifort 2021.10에서 deprecation 경고만 있고 동작 — 유지, 문제 시 `-qmkl`
2. 빌드 명령: `./scripts/in_contain.sh make -C Source` (내부적으로 서브디렉토리 병렬 make → 최종 링크)
3. 검증: `Source/cupid.x` 생성 확인 + `in_contain.sh ./Source/cupid.x` 단독 실행 시 "lack of somaFlow.in" 메시지(정상적인 입력 체크 도달)까지 확인
4. `clean` 지원: `build.sh clean` → `make clean`
5. **`.gitignore` 추가**: `*.o`, `*.mod`, `*__genmod*`, `cupid.x`, `foamGrid.in`, `*.sif` — 빌드 부산물이 소스 트리에 생기는 구조이므로 필수

### 3-2. 케이스 준비 (`scripts/prepare_case.sh`)

1. `foamGrid.in` 없으면 `7z x foamGrid.7z` 로 해제 (112 MB) — **호스트에서 수행** (컨테이너에 7z 없음)
2. `somaFlow.in` 존재 확인 — **없으면 명확히 안내하고 중단** (§7 미결 사항)
3. `Source/cupid.x` → 케이스 디렉토리로 복사 (README 권장 방식)

### 3-3. 실행 (`scripts/run.sh`)

```bash
# 개념 흐름 (in_contain.sh 경유)
./scripts/prepare_case.sh
./scripts/in_contain.sh bash -c "cd '$CUPID_CASE' && mpirun -np $CUPID_NP ./cupid.x"
```

- rank 수는 `CUPID_NP` 환경변수로 제어 (기본 4, 직렬 검증은 1)
- 로그를 `run_YYYYMMDD_HHMM.log` 형태로 tee 저장 → 이후 골든 기준 채취에 재사용

### 3-4. 완료 기준 (DoD)

- [ ] `build.sh` 1회 실행으로 `cupid.x` 생성 (경고는 허용, 에러 0)
- [ ] `run.sh`로 iSMR 케이스가 최소 수 timestep 진행되고 PMG 솔버(`isol_mg=-2`) 경로 진입 로그 확인
- [ ] 재빌드(clean 후) 시에도 동일하게 성공 — 재현성 확인

---

## 4. Phase 3 — PMG 솔버 추출, 유닛 테스트, 골든 회귀

### 4-0. 전략

- **소스 이원화 금지**: GMG 소스를 복사해 별도 사본을 만들면 원본과 diverge 한다. `Source/GMG/`를 유일한 소스로 유지하고, `pmg_standalone/`은 **스텁 모듈 + 드라이버 + makefile만** 가진 채 `Source/GMG/*.f90`를 직접 컴파일 대상으로 참조한다. 고도화도 `Source/GMG`에서 직접 수행하며, 테스트가 안전망 역할을 한다.
- **테스트 2계층**: ① 합성 문제 유닛 테스트(입력 데이터 불필요 → 즉시 착수 가능), ② 실제 케이스 덤프 기반 골든 회귀(somaFlow.in 확보 후).
- **직렬 우선**: np=1 하네스를 먼저 완성하고, 논문 핵심인 통신 저감(`nlv_glomax`) 검증을 위해 MPI(2, 4 rank) 하네스로 확장.

### 4-1. 스텁 모듈 작성 (`pmg_standalone/stub/`)

GMG가 `USE ... ONLY:`로 참조하는 CUPID 본체 심볼만 최소 구현:

| 스텁 | 제공할 심볼 (조사 결과 기준) |
|------|------------------------------|
| `MD_matrix` | `nnz, ia, ja, ju, au, u, b, alu` (CSR) |
| `MD_geometry` | `nnode, nelem, coord, num_neigh_mg, neigh_mg, imap` |
| `MD_MPI` | `nintf, myrank, myrankt` 등 |
| `MD_parameter` | `maxit, crit, ndom, ndim, teta, ...` |
| `MD_connectivity` | `ia_neigh, ja_neigh, nnz_neigh, ...` |
| `Zbicg, Zcore, Znode, Zparam` | `eps_bicg, myrank, np, nd_max, ns, ndim` 등 `ONLY:` 참조분만 |

- `GMG/module/MD_MG_*.f90`는 GMG 전용이므로 스텁 불필요 — 그대로 컴파일.
- 작성 절차: `grep "^ *USE" Source/GMG/*.f90 | sort -u`로 전체 의존 목록을 기계적으로 추출 → 스텁 뼈대 생성 → 컴파일 에러를 따라 보강 (컴파일러가 명세서 역할).

### 4-2. 덤프 계층 (CUPID 본체 계측)

골든 데이터를 실제 해석에서 채취하기 위한 최소 침습 훅:

1. `Source/05_Solver/pressure_solve.f90`의 `SOLVE_GMG` 호출 직전/직후에 덤프 서브루틴 호출 추가 (신규 파일 `dump_pmg.f90`, 환경변수 `CUPID_PMG_DUMP=<step지정>`일 때만 활성 — 평상시 오버헤드 0)
2. **덤프 내용** (Fortran unformatted stream — `-r8` 정밀도 그대로 보존):
   - 입력: `n, nnz, ia, ja, au, b, u0` + 기하/연결성 배열 + 활성 `mg.in` 파라미터 + (병렬 시) rank별 분할 정보
   - 출력(기준값): 수렴해 `u*`, 반복 수 `its`, 잔차 이력
3. 채취 시점: 초기 과도(행렬 변화 큼) 1개 + 준정상 구간 1개 이상, rank 1/2/4 각각

### 4-3. 테스트 드라이버 (`pmg_standalone/driver/`)

```
driver_pmg.f90 흐름:
  1. 덤프 로드 (또는 합성 문제 생성) → 스텁 모듈 배열에 주입
  2. mg.in 읽기 (1_read_input의 GMG 파트 재사용)
  3. 셋업 체인: Prep_fine → PREP_GMG (MG 계층 구성)
  4. solve_pbcg_mg 호출
  5. 결과 비교 → PASS/FAIL + 지표 출력 (its, 최종잔차, ||u-u*||)
```

- 합성 문제 생성기: 3D 구조 격자 Poisson(등방/고종횡비 이방성 2종) — 격자 무관 수렴성(반복 수 포화), 잔차 단조 감소, 대칭 문제에서의 해 정확도를 검증하는 **유닛 테스트**
- 빌드: `pmg_standalone/makefile`이 `stub/ + driver/ + ../Source/GMG/*.f90`를 컨테이너 안에서 컴파일 (동일 `makefile.in.apptainer` 플래그 사용 — 골든 재현성의 전제)

### 4-4. 골든 회귀 판정 기준

| 변경 유형 | 판정 기준 |
|-----------|-----------|
| 리팩터링(연산 순서 불변: 이름 변경, 모듈화, 죽은 코드 제거) | `u*` **bitwise 일치** + `its` 일치 (`-fp-model strict` 전제) |
| 동등 알고리즘 재배열(루프 교환, 벡터화 등) | `its` 일치(±0) + `‖u−u*‖₂/‖u*‖₂ ≤ 1e-12` |
| 의도적 알고리즘 개선 | 수렴 판정 잔차 충족 + `its` 감소 or 동일 + 물리해 차이 `≤ 1e-8` — **이때 골든 재채취 및 사유 기록** |

- `run_tests.sh`: 합성 유닛 테스트 → 직렬 골든 → MPI 골든 순으로 실행, 요약 표 출력. 모든 실행은 `in_contain.sh` 경유.
- **사이클 단위 운용 판정은 [LOOP.md](LOOP.md) §2 를 따른다** — 1차 지표는 its + 독립 계산 residual(`‖b−Au‖/‖b‖`), 위 표의 bitwise/해 비교는 마일스톤(L3) 검증용.

### 4-5. 완료 기준 (DoD)

- [ ] 합성 문제 3종에서 PMG-BiCGSTAB 수렴 (유닛 테스트 green)
- [ ] iSMR 덤프 재생 시 본체 실행과 `its`/최종잔차/`u*` 일치 (추출 자체의 무결성 검증 — 가장 중요)
- [ ] `run_tests.sh` 1커맨드로 전체 판정 가능

---

## 5. Phase 4 — 코드 고도화 루프

### 5-0. 고도화 3대 목표 (우선순위 순)

| # | 목표 | 배경 |
|---|------|------|
| **G1** | **파일 기반 셋업의 대규모 np 동작 수정** — 현재 코드는 **np=900에서 실행 불가**. 우선 파일 방식을 유지한 채 np=900이 돌도록 고친다 | GMG 셋업이 rank 0에서 `MG_tmp/part###.out`·`part_MG###.out`을 **np개씩 ASCII로 파일 출력**하고 각 rank가 되읽는 구조 (`6_subdomain_infor_mg.f90`, `2_read_mesh_MPI.f90`) |
| **G2** | **file out 제거 — MPI 통신으로 대체** | G1로 동작을 확보한 뒤, 파일 경유 분배 자체를 `MPI_Scatterv` 등 통신으로 교체해 `MG_tmp` 파일 I/O를 없앤다 |
| **G3** | **코드 간소화·가독성 개선** | 현재 코드가 지저분하게 짜여 있어 가독성이 크게 떨어짐. 구조 정리를 독립 목표로 둔다 |

### 5-1. G1 — np=900 실행 가능하게 (파일 방식 유지한 최소 수정)

조사에서 확인된 대규모 np 실패 후보 (수정 전 np=900 재현·원인 확정 필요):

- **유닛 번호 위생**: rank별 유닛을 `10+myrank`(읽기), `50+prc`(쓰기)로 잡음. 소스 검증(2026-08-09) 결과 Fortran 유닛은 **프로세스(rank)-로컬**이라 rank 간 충돌은 성립하지 않고, 쓰기·읽기 모두 반복마다 CLOSE 함(`6_subdomain_infor_mg.f90:188,583`, `2_read_mesh_MPI.f90:205`) — 단독으로는 np=900 실패 원인으로 확정 불가. 단 **동일 rank 내 대역 겹침**은 실재: rank 0은 셋업 중 유닛 51~`50+np`(np=900이면 501·999 포함)를 쓰고, solve 중 `write(501,*)` 잔차 로그(`6_solver_pbcg_mg.f90:280`)·`WRITE(999,*)` 에러 로그(`:271`)와 대역이 겹침(셋업 CLOSE 후라 즉시 크래시는 아니나 잉여 `fort.501` 생성 등 오염) → `NEWUNIT=` 또는 rank 무관 고정 유닛으로 정리
- **파일명 3자리 인코딩**: `CHAR(prc/100+48)` 방식은 np≥1000에서 파일명이 깨짐 (np=900은 가능) → `WRITE(fname,'(I0.4)')` 형태로 일반화
- **rank 0의 직렬 파일 fan-out**: 2×np개 ASCII 파일을 rank 0이 순차 생성 — open-file 한도·파일시스템 부하 등 환경 요인 점검
- 수정 판정: 골든 회귀(§4) green 유지 + np=900 스모크 (클러스터 확보 시; 로컬은 12코어라 oversubscribe 축소 재현으로 대행)

### 5-2. G2 — 파일 경유 분배를 MPI 통신으로 대체

- 대상: `6_subdomain_infor_mg.f90`(쓰기 측) ↔ `2_read_mesh_MPI.f90`(읽기 측)의 파일 왕복 전체
- 방법: 쓰기 측이 파일에 쓰는 배열 시퀀스를 그대로 pack → rank 0에서 `MPI_Scatterv`/`MPI_Send`로 각 rank에 직접 전달 → 읽기 측은 `READ`를 수신 unpack으로 교체. **파일에 쓰던 데이터 목록·순서가 곧 통신 프로토콜 명세**이므로, G1 완료 후 파일 포맷을 기준 삼아 1:1 대응으로 이행
- 이행 안전망: 과도기에 "파일 모드/통신 모드" 스위치를 두고 두 경로의 결과 일치(bitwise)를 확인한 뒤 파일 경로 삭제
- 효과: `MG_tmp` 디렉토리 자체가 불필요해짐 — 대규모 np 실행의 구조적 병목 제거. **셋업 단계 통신 구조 개선으로 논문 서사(통신 저감)와도 연결**

### 5-3. G3 — 간소화·가독성 (동작 불변 리팩터링)

- `IMPLICIT NONE` 전면화, 미사용 변수/죽은 코드 정리 (`mg.in` 주석의 "not use" 파라미터 포함)
- 하드코딩 상수(유닛 번호, 배열 한도 등)의 파라미터화, 모듈 구조 정리
- 골든 회귀의 **bitwise 판정(§4-4 리팩터링 행)** 을 전제로 작은 단위로만 진행

### 5-4. 작업 방식 (공통 루프)

1. **베이스라인 고정**: 골든 green 상태에서 태그 커밋 (`/git_commit` 사용), 성능 기준표 기록 (케이스별 its, 솔버 벽시계 시간)
2. **작은 단위 반복**: 변경 → `build.sh`(또는 standalone make) → `run_tests.sh` → green 확인 → 커밋. red면 즉시 원인 규명 (bitwise 기준이라 회귀 지점이 정확히 드러남)
3. G1→G2→G3 순서를 기본으로 하되, G3 중 유닛 번호 파라미터화처럼 G1과 겹치는 항목은 G1에서 함께 처리
4. 그 외 고도화 후보 (G1~G3 이후): 스무더/coarsest 솔버 성능 개선, 통신 저감 로직(`nlv_glomax`) 개선 — **논문 실험 항목과 직결**
5. 성능 변화는 회귀 판정과 분리해 `pmg_standalone/perf_log.md`에 누적 기록 → 논문 결과 섹션의 원자료

---

## 6. 리스크와 대응

| 리스크 | 영향 | 대응 |
|--------|------|------|
| `somaFlow.in` 부재 | 케이스 실행·골든 채취 불가 | 사용자에게 요청 (§7). 대기 중에도 Phase 1·2(빌드까지)·4-1·4-3(합성 테스트)은 진행 가능 |
| SIF·METIS 절대경로가 체계마다 다름 | 다른 머신에서 스크립트 동작 불가 | `SIF`·`METIS_LIB`를 env.sh 기본값 + 환경변수 재정의 계약으로 관리 (§2-2) |
| ifort↔호스트 CPU 불일치 (AVX512 플래그) | 실행 시 illegal instruction | `makefile.in.apptainer`에서 AVX2 고정 (§3-1) |
| 112 MB `foamGrid.in` 등 대용량 산출물 | git 오염 | `.gitignore` 등록 (SIF는 리포 외부 `/root/00_apptainer`라 무관) |
| 골든이 컴파일 플래그에 민감 | 위양성 회귀 실패 | 플래그를 `makefile.in.apptainer` 한 곳에서만 관리, 골든 메타데이터에 플래그·컴파일러 버전 기록 |
| WSL2에서 Intel MPI 동작 이슈 | 병렬 실행 불안정 | 단일 노드 shm 통신이라 위험 낮음. 문제 시 `I_MPI_FABRICS=shm` 강제 |

---

## 7. 미결 사항 (사용자 확인 필요)

1. **`rv_parameters.in` 확보 (여전히 권장)** — 사용자가 `rv_model=0` 으로 우회해 스모크는 통과했으나(LOG C005-r8), **step 38 에서 과도가 물리적으로 발산** (DP_MAX 지수 성장 후 솔버 한계 도달). 추정 원인이 rv 모델 부재(노심 저항/열원)이므로, 골든 덤프 채취(C009)와 완주 검증에는 rv_model=1 + rv_parameters.in 복원이 필요. (`rv_ht_str=0` 이므로 `ht_str_*.in` 은 불요 추정)
2. **골든 기준 rank 구성** — 1/2/4 rank 제안. 논문 스케일링 실험 계획에 맞춰 조정할지?
3. 논문의 다른 케이스(MSFR 정상상태, iSMR 일일부하추종)도 추후 같은 구조로 추가할지?

> (해소됨) oneAPI 버전 선택 — 기존 `hpc23.sif`(oneAPI 2023.2.1, ifort classic 2021.10.0) 사용으로 확정.
> (해소됨) `somaFlow.in` — 2026-08-09 사용자 제공. 이 소스가 모르는 신형 옵션 4개(`HS_coupling`, `vfporous`, `i_droplet`, `i_fs_temp_intpol`, 전부 0=OFF)는 주석 처리함 (LOG C005-r3, 원본 백업 보관).

---

## 8. 실행 순서 요약 (전체 체크리스트)

```
Phase 1  ☑ scripts/env.sh 작성 (SIF·METIS_LIB 기본값 + 외부 재정의 계약)   [LOG C001]
         ☑ scripts/in_contain.sh 작성                                      [LOG C001]
         ☑ 스크립트 경유 스모크: mpiifort 버전 + $METIS_LIB/libmetis.so 확인 [LOG C001]
Phase 2  ☑ Source/makefile.in.apptainer 작성 (AVX2, p=12)                  [LOG C002]
         ☑ scripts/build.sh 작성 → cupid.x 빌드 성공 (r1: Modules/ mkdir)  [LOG C002]
         ☑ .gitignore 정비                                                 [LOG C003]
         ☑ scripts/prepare_case.sh (7z 해제) + run.sh 작성                 [LOG C003·C005]
         ☑ iSMR 케이스 스모크 실행 — 38스텝·PMG 76솔브 (its 1~3)           [LOG C005-r8]
           (완주 검증은 rv_parameters.in 확보 후 — §7-1)
Phase 3  □ pmg_standalone/stub 모듈 작성 (grep USE 기반)
         □ 합성 문제 생성기 + driver_pmg.f90 → 유닛 테스트 green
         □ 본체에 dump_pmg 훅 추가 (환경변수 게이트)
         □ [실행 가능 시] 골든 덤프 채취 (직렬 → MPI)
         □ run_tests.sh 완성, 추출 무결성 검증 (its/u* 일치)
Phase 4  □ 베이스라인 태그 + 성능 기준표
         □ G1: np=900 실패 원인 확정 (유닛 충돌·파일명·fan-out) → 파일 방식 유지 수정
         □ G2: MG_tmp 파일 왕복 → MPI 통신 대체 (과도기 이중 모드로 bitwise 검증)
         □ G3: 간소화·가독성 리팩터링 (IMPLICIT NONE, 죽은 코드, 상수 파라미터화)
         □ 고도화 반복 (변경→테스트→커밋), perf_log.md 누적
```
