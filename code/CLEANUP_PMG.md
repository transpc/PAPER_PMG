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

---

# §7 고정 조합 클리닝 플랜 — "현 mg.in 조합 전용" 단일화 (2026-08-19 수립)

전제: **현재 케이스 mg.in 조합만 사용**한다는 결정 (isol_mg=-2, POL, ihybrid=1,
icommu=2, igather=1, n_GC=1, i_dir=1, isemi=0, mdf_matrix=1, ioplv=1).
스무딩 튜닝 파라미터(teta, teta_p, itergs, icheb(1)/(3)/(4), alpha, ip_*)와
신규 강건성 스위치(ieig_pol, il1_gs)는 유지. §2·§3 인벤토리를 실행 플랜으로
구체화한 것이며, **본 절 수립 시점 기준 계획일 뿐 실행 전**이다.

규모: GMG 현재 19,064줄 → 예상 삭제 약 3,600~4,300줄(~20%). §1(794줄)과 별도.

## P0 — 선행 결정 (실행 전 확정 필요)

| # | 결정 | 권고 | 근거 |
|---|---|---|---|
| D1 | `smothing` 스위치(GAS/JAC/ILU/COG) 완전 삭제? | 삭제 | git 이력으로 복원 가능. GAS 는 이번 진단의 A/B 도구였으나 역할 종료. 단 gathered 레벨 GS(SOLVE_COARSE)는 별개로 존치 |
| D2 | `ieig_pol=1` 기본화 + 골든 재베이스라인을 클리닝 전에? | **선행 권고** | 기본화가 먼저면 P3 에서 Lanczos 경로(~200줄) 삭제 가능. 후행이면 스위치·양 경로 유지 |
| D3 | `il1_gs` 유지? | 유지 | np-강건 자산, 비용 ~15줄 |
| D4 | `relax` 처리 | 1.0 하드코딩 | gathered GS 만 사용 — namelist 에서 제거 |
| D5 | MKL 링크 플래그 | 유지 | padiso.f 삭제와 무관하게 CUPID 본체가 -mkl 사용 |
| D6 | isetup_comm 파일 모드 소거 | 이 플랜 범위 밖 | G2 트랙에서 기본값 플립 후 별도 수행 |

## P1 — 조건부 데드 파일 5종 (호출 분기 소거 + 파일 삭제, ~1,500줄)

§2 의 실행. 순서 = 리스크 오름차순, 파일당 1커밋:

1. `coarsening_semi_amg.f90`(374) — isemi 분기 소거 (5_PREP_GMG 2곳)
2. `padiso.f`(127) — ipar 분기 소거 (SOLVE_GC 2곳)
3. `solve_CG.f90`(317) — SOLVE_EXACT* 의 i_dir∉{1,2} 분기 소거
4. `mt_precond.f90`(109) — stiffness_GC:281 분기 소거
5. `06_solver_pcg_ilu.f90`(629) — **선행: amux0P 를 mt_amux.f90 으로 이동**
   (BiCGSTAB 핵심 SpMV), ilupcp 호출부(7_SOLVE_GMG isth==2 셋업 2곳) 소거 후 삭제

## P2 — 솔버 변형 단일화 (~1,700줄, 최대 덩어리)

경로 고정: `SOLVE_GMG → solve_pbcg_mg → SOLVER_NEW → SOLVE_GC_all → SOLVE_COARSE → SOLVE_EXACT(i_dir=1)`

| 파일 | 삭제 대상 | 잔존 |
|---|---|---|
| 7_SOLVE_GMG.f90 (1,527) | `SOLVER`(구판, ~245), `SOLVER_NEW_MPI`(~345), `matrix_vec`, `*_N_MPI`/`residl_MPI`, isol_mg 디스패치 분기, n_GC=0 분기(PCG_Dig 호출), isth==2 셋업, `Dig_mdf_matrix_inv` | `SOLVE_GMG`(디스패치 단순화), `SOLVER_NEW`, `matrix_vec_N`, `residl`, `Dig_mdf_matrix` |
| SOLVE_GC.f90 (1,150) | `SOLVE_GC`(icommu=1), 3번째 변형(ihybrid≠1), `SOLVE_COARSE_MPI`, `SOLVE_EXACT_MPI`, `Relax_GS_MPI`, `Relax_GS_SYM` 호출 분기(id_GS_sym), i_dir=2 분기, igather=0 분기 | `SOLVE_GC_all`, `SOLVE_COARSE`, `SOLVE_EXACT`(i_dir=1 만), `Relax_GS` |
| stiffness_GC.f90 (840) | `stiffness_GC`(icommu=1), igather=0 분기, `STIFF_COARSE2`/`STIFF_*` 미사용 변형, pc_ilu 경로 | `stiffness_GC_all`, `STIFF_COARSE`, `STIFF_EXACT`(i_dir=1) |
| 6_solver_pbcg_mg.f90 (407) | ihybrid 분기(SOLVER_NEW_MPI 호출 2곳), 주석 덩어리 | + breakdown 가드 **복원**(E, 이 사이클에 포함) |
| PCG_Dig.f90 | `pcg_dig` 서브루틴 | `amux0_PCG`(POL 이 사용) |

## P3 — 스무딩 단일화 (~600줄 ± D2 결정)

- `Relax_GSP.f90`(879): `smoothing_fine` 을 POL 전용으로 축약(isth 디스패치 소거),
  `Smooth_ILU`/`smooth_CG`/`Relax_GS0P(_BW)`/`Smooth_GS2_MPI`/`Smooth_GS_BW` 삭제.
  잔존: `smoothing_fine`(축약), `Smooth_GS2`(+l1 diagrc), gathered GS 커널
- `poly_smooth.f90`(423): icheb(3) 메서드-1 분기(~50줄). D2 선행 시
  `lanczos_eig_max`+`compute_eigenvalues`(~200줄)도 삭제 → eig_value = Gershgorin 전용
- 모듈: `isth`, `id_GS_sym`, `icheb(2)/(5)`, `relax`(D4) 소거

## P4 — 완전 데드 서브루틴·고아 변수 (~300줄)

§3-1 목록 실행: `csr_coarse`, `matrix_inverse_GS`, `mt_amux1`/`mt_amux2p`,
`md_s_r_mt2`+`send_receive_mtc2`, `send_receive_mt`(sub), `send_receive_c`(확인 후),
AR_hi/iallocate_c/icase_MG/isol_start/iGS/nGS/ipar/icommu/igather/ihybrid 등
고정된 스위치의 모듈 변수·USE 정리 (읽기만 하고 분기 없어진 것 전부).

## P5 — 입력 계층 최종 정리 (read_input + mg.in)

1. `1_read_input.f90`: 전 파라미터 명시적 기본값 대입 → 각 namelist READ 에
   `iostat` 부여(그룹 부재 허용) → **namelist 그룹 통합**: `&MG_tuning`(teta,
   teta_p, alpha, itergs, icheb, ip_*) + `&MG_options`(ieig_pol, il1_gs,
   isetup_comm, nthre) 2개 수준으로 축약. `ncycle_pre`/`crit_pre` 지역변수 함정 제거
2. mg.in 재작성 (케이스 + tests 동기): ~15줄 수준 목표
3. 사멸 파라미터의 mg.in 잔존 항목 제거는 코드 소거와 같은 커밋으로

## P6 — 범위 밖 (별도 트랙, 이 플랜에서 하지 않음)

- isetup_comm=0 파일 모드 + MG_tmp 산출 코드 (G2 완결 후)
- dump_pmg 계측 (논문 종료 후)
- 2_read_mesh_MPI/6_subdomain_infor_mg 의 셋업 확장성 개선 (O(np·N) 등 — SCALING 이슈)
- .gitignore 정비 (fort.*, build_chk/ 등)

## 검증 프로토콜 (모든 단계 공통)

- **불변식: 기본 설정에서 모든 삭제는 수치 bitwise 중립** — dead code 제거이므로
  게이트가 그대로 안전망. 단계마다: standalone `run_tests.sh` 12/12 (its·md5 불변)
  + CUPID 전체 빌드 링크.
- P2 완료 후 + P5 완료 후: ECT1 케이스 러닝 스팟(np=1/np=4, fort.501 대조).
- 커밋 단위: P1 파일당 1커밋, P2~P5 단계당 1커밋. 각 커밋 메시지에 삭제 줄수 기록.
- 회귀 발견 시: 해당 커밋만 revert 가능하도록 단계 간 의존 최소화(P1→P2→P3 순서 준수).

---

## §7 진행 현황 (2026-08-19 실행)

| 단계 | 상태 | 커밋 | 삭제 규모 |
|---|---|---|---|
| D2 ieig_pol=1 기본화+재베이스라인 | ✅ | 23d4059 | (수치 변경 사이클) |
| P1 조건부 데드 파일 5종 | ✅ | f830434·553d865·e501264 | 1,556줄 + 분기 |
| P2 솔버 변형 단일화 + 가드 복원 | ✅ | 0308815 | -2,004줄 |
| P3 스무딩 단일화 (POL/Gershgorin 전용, Lanczos 삭제) | ✅ | (P3 커밋) | -895줄 |
| P4 완전 데드 서브루틴·AR_hi | ✅ | (P4 커밋) | -824줄 |
| D5 MKL 링크 제거 (LAPACK 보조 4종 로컬) | ✅ | 334a339 | 의존성 제거 |
| P5 입력 계층 재설계 (namelist 2그룹, mg.in 16줄) | ✅ | 7829816 | 고아 13종 |
| **D6 isetup_comm 통신모드 단일화 + MG_tmp 소거** | ⏳ 다음 | — | 분기 56곳 (2_read_mesh 34 + writer 22) |

GMG 총량: 19,064 → **13,766줄 (-28%)**. 전 단계 게이트 12/12 bitwise 동일 + CUPID 링크 green + P5 후 프로덕션 스팟(np=4) its 일치 확인.

D6 주의사항 (착수 시): ELSE(파일모드) 본문이 단문 마커 READ ~ 다층 블록까지 다양 — 정규식 일괄 치환 금지,
사이트별 확인 절제. 완료 후 게이트가 comm 모드를 최초로 상시 검증하게 됨 (rt_ascii shim 이 bitwise 보장,
C011-1 검증). MG_tmp 디렉토리 생성/판독 코드와 iu_prc 유닛 대역도 함께 소거. isetup_comm 변수 자체 제거.

---

# §8 출력·디버그 코드 정리 플랜 (O-시리즈, 2026-08-20 수립 — 계획만)

실측 근거: 프로덕션 np=64 러닝 로그 1,703줄 분석 + GMG 출력 지점 전수조사.

## 인벤토리 (실측)

| 채널 | 실태 | 성격 |
|---|---|---|
| 콘솔 (랭크별 스팸) | `test,pre ilu`·`test,after ilu`(02_IO/read_grid.f90:508 부근), `====>lev_typet`(02_IO/reorder_ilup.f90:94) — **np개씩 반복** | 순수 디버그 잔재. **진원은 GMG 가 아니라 CUPID 초기화부** |
| 콘솔 (GMG, 48곳) | 대부분 메모리부족/정합성 에러 메시지 (+STOP 동반) | 유지 대상 (에러 채널) |
| fort.999 (GMG 31곳) | **전 랭크가 같은 파일에 기록** — np>1 에서 인터리브 오염. 에러(다수)+정보(coarsest info 등) 혼재 | 정책 결정 필요 |
| fort.16 | rank0 이 **매 예조건자 적용마다** icycle 1줄 (러닝당 수천 줄) | 상시 가치 없음 — 삭제 후보 |
| fort.501 | rank0 its 로그 | **골든 게이트 기준 — 절대 불변** |
| fort.101/400/401 | CUPID 측 타이밍 (pressure_solve, rank0 가드) | GMG 범위 밖 — 포함 여부 결정 |
| 주석 디버그 잔해 | GMG 44곳 (write 주석) + system_clock/시간측정 블록 다수 | 일괄 소거 후보 |

## 단계

- **O1 (즉효 — 랭크 스팸 제거)**: 02_IO 의 3종 디버그 프린트 삭제. 이후 np=4 스팟 로그에서
  np회 반복 라인 전수 재확인 → 잔여 동일 처리. 1커밋.
- **O2 (에러 채널 일원화)**: GMG 에러 출력 정책 통일 — 에러는 콘솔(rank 표기)+`stop_mpi` 경유,
  fort.999 는 rank0 전용 진단 로그로 축소(전 랭크 기록 소거) 또는 전폐. ⟵ 결정 (a)
- **O3 (상시 로그 다이어트)**: fort.16 icycle 기록 삭제. fort.501 불변.
  fort.101/400/401 은 perf_log 채취에 쓰이는지 확인 후 결정. ⟵ 결정 (b)
- **O4 (사체 소거)**: 주석 디버그 44곳 + system_clock 잔해 + PAUSE 잔존 확인 일괄 제거 (폴리싱).

## 검증 프로토콜

- 게이트 12/12 (fort.501 bitwise 불변이 핵심 안전망)
- 프로덕션 np=4 스팟: 콘솔 로그 라인 수 before/after 비교 + fort.16 부재 확인
- np>1 에서 fort.999 인터리브 소멸 확인 (O2 후)

## 진행 현황 (2026-08-20 완료)

| 단계 | 상태 | 커밋 | 내용 |
|---|---|---|---|
| O1 랭크 스팸 제거 | ✅ | 10c84ca | 02_IO 디버그 프린트 3종 (test,pre/after ilu·lev_typet) |
| O2 에러 채널 일원화 | ✅ | c6a09f1 | 에러 14곳 fort.999→콘솔(`PMG error` 접두), STOP_MPI 콘솔 사본, PAUSE 3곳 소거 |
| O3 상시 로그 다이어트 | ✅ | 7c8de3c | fort.16 icycle 소거. 타이밍 유닛(101/400/401)은 유지 결정 |
| O4 사체 소거 | ✅ | (O4 커밋) | 디버그 주석·system_clock·미사용 타이밍 선언 75줄 |

최종 상태: GMG 활성 콘솔 출력 = 에러 전용(`PMG error` 접두), fort.999 = rank0 전용 진단 16곳,
fort.501 = 골든 게이트 기준(불변), PAUSE 0, 주석 디버그 사체 0.
검증: 매 단계 게이트 12/12 bitwise 동일 + CUPID 링크 green + np=4 프로덕션 스팟.
※ 게이트 판정은 반드시 **clean 재빌드** 기준 (O3 때 incremental 상태에서 산발 실패 관측 — 코드 무결).

## 결정 기록

- (a) fort.999 정책 → **rank0 전용 진단으로 축소** (에러는 콘솔 일원화). O2 에서 적용
- (b) CUPID 타이밍 유닛(101/400/401) → **유지 결정(O3 시점)**: rank0 가드 확인, 스텝/압력솔버 벽시계 = 논문 성능 계측 채널
- (c) O1 의 02_IO 수정이 "GMG 클리닝" 스코프 밖임을 감안해 별도 커밋으로 분리할지 (권장: 예)
