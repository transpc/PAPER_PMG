# PMG 코드 클리닝 인벤토리 — 삭제 중심 (2026-08-14)

목적: PMG(GMG) 관련 루틴의 데드코드를 **파일 삭제 위주**로 정리하기 위한 상세 목록.
이후 단계(폴리싱)에서 파일 내 서브루틴/분기 소거를 수행한다. 본 문서는 그 전 단계의
근거 자료이며, 기존 소스 파일은 아직 변경하지 않았다.

조사 방법: GMG 전체 서브루틴/함수 정의 182개를 추출하고, `Source` 전체 +
`pmg_standalone/driver`·`tests` 에 대해 호출부(case-insensitive `CALL`)를 역추적.
함수형 호출(CALL 없는 사용)이 의심되는 항목은 개별 재검증함 (`f_face`, `f_line`,
`shape_f`, `volm`, `dump_now` 은 함수형으로 살아있음을 확인 → 목록에서 제외).

---

## §0 전제 — 고정된 mg.in 설정 (2_iSMR_ECT_res1, 사실상 fix)

| 파라미터 | 값 | 결과적으로 죽는 경로 |
|---|---|---|
| `isol_mg` | -2 | `SOLVER`(old MG, isol_mg=1) 직접해법 경로 |
| `ihybrid` | 1 | `SOLVER_NEW_MPI` (ihybrid≠1 분기) |
| `smothing` | 'POL' (isth=3) | GAS/JAC/ILU/COG 스무딩 전체 (`Smooth_ILU`, `smooth_CG`, `Relax_GS0P(_BW)`, `ilupcp(_new)`, `lusol0p`) |
| `isemi` | 0 | `coarsening_semi_amg` 계열 |
| `AR_hi` | 0 | 5_PREP_GMG 의 고종횡비 teta 보정 분기 |
| `n_GC` | 1 | 7_SOLVE_GMG 의 n_GC=0 분기 (`PCG_Dig`의 `pcg_dig`, 최심층 분산 GS) |
| `i_dir` | 1 | `SOLVE_EXACT*` 의 CG/Pardiso 분기 (`solve_cg2`, `pardiso_solve`), `pc_ilu1` 셋업 |
| `ipar` | 0 | Pardiso 전체 (`padiso.f`) |
| `id_GS_sym` | 0 | `Relax_GS_SYM` (backward GS) |
| `icommu` | 2 | icommu=1 (root 해석 후 배포) 분기 |
| `igather` | 1 | MPI_REDUCE 취합 분기 |
| `ioplv` | 1 | 수동 `nlevel`/`nlv_glomax` 지정 |

활성 실행 경로 (유일):
`read_input_mg → subdomain_infor_mg(+prep 체인) → assemble_FVM → SOLVE_GMG
→ solve_pbcg_mg(BiCGSTAB) → [예조건자] SOLVER_NEW → V-cycle(poly_smooth POL
+ SOLVE_GC/SOLVE_EXACT i_dir=1 Ainv 직접곱)`

---

## §1 즉시 삭제 가능 파일 — 외부 호출 0 (링크 무해) — ✅ 실행 완료 (2026-08-14)

실행 내역: 3개 파일 `git rm` + 잔여 `.o` 제거(링크가 `./GMG/*.o` 와일드카드라 필수)
+ makefile 2곳 라인 제거. 검증: 컨테이너 CUPID 전체 재빌드 링크 성공,
standalone `run_tests.sh` 12/12 PASS (합성 5 + 골든 iSMR/ECT1 재생 7, 충실도 게이트 포함).

어떤 서브루틴도 파일 밖에서 호출되지 않음. 삭제 시 두 makefile 에서 해당 라인만
제거하면 됨 (이것이 유일한 수반 수정).

| 파일 | 내용물 | 비고 |
|---|---|---|
| `Source/GMG/Relax_GSx.f90` | `Relax_GSS`, `Relax_GS0` | 내부 상호 호출만 존재. 내부에서 부르는 `mt_amux2` 는 SOLVE_GC 도 사용하므로 `mt_amux.f90` 은 유지 |
| `Source/GMG/mt_amuxvr1.f90` | `mt_amuxvr1` + `mt_amuxrr3~99/n` 등 13개 | 전부 내부 연쇄만 존재, 진입점 없음 |
| `Source/GMG/solve_CG_P.f90` | `solve_CG_P`, `mt_amux_P`, `mt_lusol_P` | 병렬 CG 코어스 솔버 구판. 진입점 없음 |

수반 수정 (삭제의 일부):
- `Source/GMG/makefile` 의 FSRCS 에서 3줄 제거
- `pmg_standalone/makefile` 의 GSRC 에서 3개 항목 제거 (두 목록은 항상 동기 유지)

검증: 컨테이너에서 CUPID 재빌드 + standalone `run_tests.sh` 골든 일치.

---

## §2 조건부 데드 파일 — 고정 설정에서 실행 불가, 단 호출 코드가 남아 있음

실행은 절대 안 되지만 호출문이 소스에 존재하므로, **호출부(죽은 분기) 소거 후**
파일 삭제 가능. 폴리싱 단계 첫 대상으로 권장. 괄호는 소거할 호출부.

| 파일 | 죽는 조건 | 소거할 호출부 |
|---|---|---|
| `coarsening_semi_amg.f90` | isemi=1 전용 | `5_PREP_GMG.f90:133`, `5_PREP_GMG_global_coarse.f90` (isemi 분기) |
| `padiso.f` | ipar=1 전용 | `SOLVE_GC.f90:797,875` (SOLVE_EXACT/_MPI 의 ipar 분기). 부수효과: MKL Pardiso 의존 제거 |
| `solve_CG.f90` | `solve_cg2` 는 i_dir∉{1,2} 전용; `solve_cg`/`mt_amux`/`mt_lusol` 은 완전 데드 | `SOLVE_GC.f90:795,873` |
| `mt_precond.f90` | `pc_ilu1` 은 i_dir=0∧ipar=0 셋업 전용; `pc_ilu2` 완전 데드 | `stiffness_GC.f90:281` (+ `stiffness_MG.f90:202` 주석 1줄) |
| `06_solver_pcg_ilu.f90` | 아래 예외 1개 빼고 전부 데드: `ilu0_gpt`/`pcg_ilu_s`/`pbicg_ilu_s` 완전 데드, `ilupcp(_new)` isth=2 전용, `lusol0p` 는 `Smooth_ILU`(isth=2 전용)에서만 호출 | **예외: `amux0P` 는 살아있는 핵심 SpMV** (`6_solver_pbcg_mg` 의 BiCGSTAB 본체가 매 반복 호출). `amux0P` 를 살아있는 파일(예: `mt_amux.f90`)로 이동한 뒤 파일 삭제 |

`PCG_Dig.f90` 은 삭제 불가: `pcg_dig` 자체는 n_GC=0 전용 데드지만, 같은 파일의
`amux0_PCG` 가 `poly_smooth`(POL 스무딩 본체)와 `Relax_GSP` 에서 호출되는 활성
루틴. → 폴리싱 때 `pcg_dig` 서브루틴만 제거.

---

## §3 파일 내 데드 서브루틴 (폴리싱 단계 목록 — 파일은 유지)

### 3-1 완전 데드 (호출부 자체가 없음 — 무조건 삭제 가능)

| 파일 | 서브루틴 |
|---|---|
| `7_SOLVE_GMG.f90` | `matrix_vec` |
| `connectivity_coarse.f90` | `csr_coarse` |
| `linear_algebra.f90` | `matrix_inverse_GS` (무인자 GS 역행렬; `_GS_n`/`linalg_invM` 은 활성) |
| `mt_amux.f90` | `mt_amux1`, `mt_amux2p` |
| `send_receive_mt.f90` | `md_s_r_mt2` + 그 전용 하위 `send_receive_mtc2` (연쇄), `send_receive_mt`(sub) |
| `send_receive.f90` | `send_receive_c` (내부 호출 위치가 죽은 루틴인지 확인 후) |
| `2_Prep_fine_global.f90` / `3_Prep_fine_P.f90` | `csr_fvm` / `csr_fvm_p` — 내부 호출 1곳씩, 호출 지점이 활성인지 확인 후 판단 |

### 3-2 고정 설정상 데드 (분기 가드와 함께 소거)

| 위치 | 대상 | 죽는 조건 |
|---|---|---|
| `7_SOLVE_GMG.f90` | `SOLVER` (구판 V-cycle 반복해법) | isol_mg=1 전용 |
| `7_SOLVE_GMG.f90` | `SOLVER_NEW_MPI` | ihybrid≠1 전용 (`6_solver_pbcg_mg.f90:162,212` 분기) |
| `7_SOLVE_GMG.f90` | n_GC=0 분기 전체 (내부 `PCG_Dig`/`Relax_GSP` 루프) | n_GC=1 고정 |
| `7_SOLVE_GMG.f90:880,1229` | `ilupcp` 셋업 분기 | isth=2 전용 |
| `Relax_GSP.f90` | `Smooth_ILU`, `smooth_CG`, `Relax_GS0P_BW` (`smoothing_fine` 내부 분기 포함) | isth=2/4, id_GS_sym=1 전용 |
| `SOLVE_GC.f90` | `Relax_GS_SYM` 호출 분기 12곳 | id_GS_sym=1 전용 |
| `SOLVE_GC.f90` | `SOLVE_EXACT*` 의 i_dir=2 분기 (분산 Ainv 곱 + ALLREDUCE) | i_dir=1 고정 |
| `PCG_Dig.f90` | `pcg_dig` | n_GC=0 전용 |
| `stiffness_MG.f90:212`, `7_SOLVE_GMG.f90:261,1037,1382` | icommu=1 분기 | icommu=2 고정 |
| 취합 루틴들 | igather=0 (MPI_REDUCE) 분기 | igather=1 고정 |
| `5_PREP_GMG.f90:108-111` | AR_hi=1 teta 보정 | AR_hi=0 고정 |
| `1_read_input.f90` | `isol_start` (namelist 에만 존재, 사용처 0), `icase_MG` (사용처 0), `iter_mg`/`iter_max`("not use" 주석 확인됨) | — |
| `MD_MG_index` 등 모듈 | 위 삭제에 따른 고아 변수 (`icheb(2)`,`icheb(5)` 포함) | — |

### 3-3 보류 (G2 작업과 연동)

- `2_read_mesh_MPI.f90` 의 파일모드(`isetup_comm=0`) 분기: G2 통신화가 완료·검증되면
  파일모드 분기 + `MG_tmp` 산출 코드를 소거. **주의: 현재 mg.in 은 isetup_comm 미지정
  = 파일모드가 기본.** 기본값을 1로 뒤집는 결정이 선행되어야 함.
- `module/dump_pmg.f90`: 골든 회귀 계측 인프라 (env `CUPID_PMG_DUMP` 게이트).
  논문 작업 종료까지 유지.

---

## §4 mg.in 축소 계획과의 연동 (파라미터 기본값 내장)

케이스의 mg.in 은 이미 57줄로 축소 완료(죽은 항목 소거, 유효값 동일 검증됨).
다음 단계로 "핵심 파라미터만 mg.in, 나머지는 코드 기본값" 을 구현하려면
`1_read_input.f90` 에서:

1. namelist READ 이전에 **전 변수 명시적 기본값 대입** (현재는 정적 0 초기화에 암묵 의존).
2. 각 `READ(iu, nml=...)` 에 `iostat` 를 주어 **그룹 부재 허용** (현재는 10개 그룹이
   파일에 순서대로 전부 있어야 함 — 이것이 mg.in 을 더 못 줄이는 유일한 이유).
3. `ncycle_pre`/`crit_pre` 를 지역변수에서 모듈/기본값으로 승격 (현재 파일에서 빠지면
   쓰레기값이 되는 함정).

기본값 내장 시 권장 분류:

- **mg.in 에 남길 튜닝 파라미터**: `teta`, `teta_p`, `alpha`, `itergs`, `icheb(1)`,
  `icheb(3)`, `icheb(4)`, `ip_nmax`, `ip_inter`, `ip_lev`, (`smothing`)
- **코드 기본값으로 고정**: `isol_mg=-2`, `mdf_matrix=1`, `ioplv=1`, `isemi=0`,
  `n_GC=1`, `i_dir=1`, `crit_1=1.d-1`, `icommu=2`, `igather=1`, `ihybrid=1`,
  `nthre=1`, `ncycle_pre=1`, `crit_pre=1.d-1`, `relax=1.d0`, `mxnbne=100`
- **namelist 에서 제거** (§2·§3 삭제 후 자연 소멸): `icase_MG`, `isol_start`,
  `AR_hi`, `id_GS_sym`, `ipar`, `iGS`, `nGS`, `iallocate_c`, `icheb(2)`, `icheb(5)`

동기화 대상: `pmg_standalone/tests/mg.in` (동일 파서 사용).

---

## §5 빌드/런타임 아티팩트 — git 미추적, 디스크 정리 + .gitignore

git 추적 파일 아님(확인됨). 삭제는 `make clean` 성격이며 언제든 재생 가능.

- 소스 트리: `Source/**/*.o`, `*.mod`, `*__genmod.f90`, `Source/Modules/` 산출물
  2,234개, `GMG/module/libcupidCLOSED.a`, `pmg_standalone/build/`,
  `pmg_standalone/pardiso_solve__genmod.*`
- 케이스(`2_iSMR_ECT_res1`): `fort.*`, `log.dat`, `saveout.dat`, `run_*.log`,
  `fort.501.d145c81.bak`, `MG_tmp/`(G2 검증 완료 후), `pmg_dump/`
  — 단 `fort.501` 은 골든 비교 기준으로 참조되므로 재생 절차 확인 후 정리.
- 제안: `.gitignore` 에 `*.o`, `*.mod`, `*__genmod.*`, `build/`, `fort.*`,
  `MG_tmp/`, `pmg_dump/` 추가.

---

## §6 실행 순서 및 검증 절차 제안

1. **커밋 1 (§1)**: 3개 파일 삭제 + makefile 2곳 라인 제거
   → 컨테이너 재빌드 → standalone `run_tests.sh` 골든 일치 → ECT1 러닝 스팟체크.
2. **커밋 2~ (§2)**: 파일별 1커밋 (호출 분기 소거 + 파일 삭제). 매 커밋 골든 검증.
   `06_solver_pcg_ilu.f90` 은 `amux0P` 이동이 선행되므로 마지막 순서 권장.
3. **폴리싱 (§3, §4)**: 파일 내 서브루틴/분기 소거와 read_input 기본값 내장.
   본 문서의 표를 체크리스트로 사용.

주의사항:
- 모든 삭제 후보는 "현 고정 설정" 기준. 설정을 되살릴 가능성이 있는 항목
  (예: GAS 스무딩 비교 실험)이 논문 플랜에 있으면 §2·§3-2 해당 행을 보류할 것.
- genmod 파일은 컴파일 산출물이므로 소스 삭제와 무관하게 정리 가능.
