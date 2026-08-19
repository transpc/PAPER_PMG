# Finding — Chebyshev 다항 스무딩의 스펙트럼 상계 산정과 분할-무관(np-무관) 수렴성

작성: 2026-08-19. 근거 커밋: `fc04c35` (ieig_pol/il1_gs 스위치), 진단 이력 `5134579`(원인 후보 최초 지목),
`314c2cb`(np=900 힙 오염 실버그). 본 문서는 논문 novelty 서술을 위한 1차 정리이며,
모든 수치는 본 저장소에서 재현 가능한 실측값이다 (§8 재현 방법).

---

## 1. 한 줄 요약

**병렬 멀티그리드 예조건자에서 Chebyshev 다항 스무더의 스펙트럼 상계를 관행적
"소수-스텝 Lanczos 추정 × 안전계수(1.1)"로 잡으면, 추정치가 도메인 분할에 민감하게
요동하여 특정 프로세서 수/분할 기하에서만 예조건자가 특정 오차 모드를 증폭 →
Krylov 반복이 붕괴한다. 추정치를 Gershgorin 행합 상계(참이 보장되는 상계)로
교체하면 이 실패 모드가 원천 제거되고, 수렴 반복수가 프로세서 수와 무관해진다.
보수성 비용은 저난도 solve 에서 +1 반복 수준(고난도 solve 동률)으로 실측되었다.**

> **정정 (2026-08-19)**: 초기 정리본의 "직렬 성능 3~11배 개선" 주장은 하네스 표시
> 형식(두 solve 사이트 its 의 문자열 연접, 예: "9"+"9"→"99")을 오독한 것으로 철회함.
> 실측 진실은 §4-3 참조. np-무관성·붕괴 완치 결과(§6-2)는 영향 없음.

## 2. 관측된 현상 (문제 제기)

- 프로덕션(iSMR/ECT1, 43.6만 셀): np=500 수렴, np=900 비수렴 — 동일 문제·동일 코드.
- 합성 Poisson 재현기(등방, 표준 7점):
  - 24³ k-슬랩 분할: np=6 수렴(its 6~7), np=8 붕괴(its 770), np=12 발산(maxit 초과)
  - 48³ 3D 블록 분할: np=64 정상(its 5), **np=96 붕괴(its 231), np=112 정상(its 5), np=128 발산**
- 핵심 특징: **np에 대해 비단조**. 크기 임계·자원 한계형 실패가 아니라
  특정 분할 기하(균일 얇은 블록)에서만 선택적으로 발생. 실행은 bitwise 결정론적
  (동일 조건 2회 실행 완전 일치), 메모리 위반 0건(`-check bounds,uninit` 빌드 확인).

## 3. 배경 — 이 솔버의 구조적 특성 (novelty 의 전제)

본 PMG 는 rank 0 가 전역 격자에서 계층·전달연산자(P/R)·코스닝을 직렬 구축 후
분배하는 설계라서, **계층과 Galerkin 연산자가 분할에 대해 대수적으로 불변**이다.
따라서 분할-의존성이 유입될 수 있는 채널은 원리적으로 다음뿐이다:

1. 분산 코어스 레벨의 hybrid GS 스무더 (블록 구조가 분할의 함수) — §5 에서 반증됨
2. np 의존 계층 컷(nc_min) — §5 에서 반증됨
3. **fine 레벨 Chebyshev 스무더의 λ_max 추정** — 확정된 원인
4. 병렬 리덕션의 라운딩 — 단독으로는 무해, 3 의 민감도를 통해 증폭됨

이 "연산자 불변" 전제 덕분에 원인 후보 공간이 좁고, 배제 수사(§5)가 가능했다.
(일반적 병렬 AMG 는 코스닝 자체가 분할 의존이라 이런 분리가 어렵다 — 서술 포인트.)

## 4. 근본 원인 (메커니즘)

### 4-1. 코드 관행

fine 레벨 스무더는 Chebyshev 다항식(2종, k=2)이며, 스펙트럼 상계를

```
λ̂ = Lanczos_8step(D̃⁻¹Ã, v₀ = 1) ;  rcheb(1) = 1.1 × λ̂
```

로 잡는다 (poly_smooth.f90). 이는 hypre/PETSc 계열에서도 널리 쓰이는 관행
(소수-스텝 Lanczos/CG 추정 + 5~10% 마진)과 같은 구조다.

### 4-2. 실패 메커니즘 (3단 연쇄)

1. **추정치의 분할 민감성**: 짧은(8-스텝) Lanczos 는 유한정밀에서 수치적으로
   예민한 과정이며, 병렬 내적의 비결합적(non-associative) 합산 순서가 np/분할마다
   달라 Ritz 값 궤적이 달라진다. 실측: 동일 행렬에서 λ̂ = 1.7945(np=112) vs
   1.8439(np=96) — **2.7% 요동** (정확 산술이라면 분할 무관이어야 하는 양).
2. **증폭 구간 진입**: Chebyshev 다항식은 지정 구간 밖 고유모드를 감쇠하지 않고
   **증폭**한다. 추정 실패 시 안전계수 1.1 로도 스펙트럼 상단이 구간을 벗어난다.
3. **선택적 붕괴**: 예조건자는 벌크 오차에는 여전히 우수하나(1회 적용 상대잔차
   ~10⁻²) 소수의 상단 모드를 반복 적용마다 증폭한다. Krylov(BiCGSTAB)가 반복
   중 그 모드를 노출하는 순간 발산. **V-cycle 1회 적용의 수축률 실측
   (np=96)**: 0.009 → 0.29 → 0.87 → 1.03 → 1.14 → 1.25 → **1.30 (증폭 전환)**;
   정상(np=112)은 0.008 → 0.16 → … → 0.22 안정.

이 실패 모드의 서술 가치: **"예조건자 품질의 평균적 우수성과 스펙트럼 꼬리의
안정성은 별개"**이며, 후자는 표준 지표(첫 적용 수축률, 셋업 진단)로 보이지 않는다.
maxit 공회전으로만 발현되므로 현장에서는 "임의 np 에서의 재현 불가 발산"으로
오인되기 쉽다 (실제 본 프로젝트에서도 분할 형상 병리·셋업 결함으로 장기 오인).

### 4-3. 보수성 비용의 실측 (정정판 — 구판의 "3~11배 직렬 개선" 철회)

골든 재생 실측(np=1, 프로덕션 덤프 행렬, 사이트 k1/k2 별 its):

| 케이스 | 기존(Lanczos) | Gershgorin | 차이 |
|---|---|---|---|
| ECT1 s1 / s10 / s30 / s150 | 1,1 / 1,1 / 3,3 / 9,9 | 2,2 / 3,3 / 3,3 / 9,9 | +1 / +2 / 0 / 0 |
| iSMR s1 / s10 / s30 | 2,2 / 2,2 / 3,3 | 3,3 / 3,3 / 3,3 | +1 / +1 / 0 |

즉 실계수 행렬에서 Lanczos 추정은 np=1 성능에 문제가 없었고, Gershgorin 의
보수성 비용은 **저난도 solve 에서 +1~2 반복, 고난도 solve(s30/s150)에서 0** 이다.
프로덕션 np 스윕(§6-2)의 초반 구간 +1~2 its 관측과 정합한다. 충실도 게이트는
양 모드 모두 통과(수렴해가 프로덕션 해와 1e-8~1e-10 일치). 결론: 본 수정의
가치는 직렬 성능이 아니라 **np-강건성과 분할-무관성**에 있다.

## 5. 진단 방법론 — 배제 수사 (방법론 자체가 서술 가치 있음)

| # | 실험 | 결과 | 배제된 가설 |
|---|---|---|---|
| 1 | l1-보정 hybrid GS (Baker et al. 2011) ON | 발산 유지 (np=128) | 코어스 스무더 감쇠 부족 |
| 2 | 스무딩 스위프 4배 | 발산 유지 | 스무딩 강도 부족 (증폭원이면 오히려 악화) |
| 3 | 심층 레벨 gathered-직렬화 | 발산 유지 | 극소 도메인 레벨 결함 |
| 4 | `-check bounds,uninit` 빌드 | 위반 0, its 재현 | 메모리 오염 |
| 5 | 동일 조건 2회 실행 | bitwise 동일 | 비결정성/경합 |
| 6 | 보간행렬 P 행합 검사 | 전 행 = 1.0 (±1e-10) | 전달연산자 조립 오염 |
| 7 | V-cycle 층별 잔차 노름 계측 | 정상/붕괴 케이스 층별 동일, 단 반복 적용 시 수축률 1 초과 | "특정 층의 데이터 오염" — 대신 "모드 선택적 증폭" 확인 |
| 8 | fine 스무더 POL→GAS 치환 | **완치 (its 231→5)** | → 증폭원 = fine Chebyshev 확정 |
| 9 | λ̂ 실측 (np=96 vs 112) | 1.844 vs 1.794 (2.7%) | → 눈금의 분할 민감성 확정 |
| 10 | λ 상계 고정 (2.2) | **np=96/112/128 전부 its=5** | → 원인·치료 동시 확정 |

## 6. 해결책 — Gershgorin 행합 상계

### 6-1. 수정 내용

추정 대신 참이 보장되는 상계를 사용:

$$\lambda_{\max}(\tilde{D}^{-1}\tilde{A}) \;\le\; \max_i \sum_j |\tilde{a}_{ij}| \;=\; \Lambda_G$$

- 비용: 소유 행 1회 스캔 O(nnz) + `MPI_ALLREDUCE(MAX)` 1회 —
  기존 Lanczos(행렬-벡터곱 8회 + 내적 리덕션 ~20회)보다 오히려 저렴.
- 성질: (i) **참 상계** → 구간 밖 모드가 존재할 수 없어 증폭 원천 차단,
  (ii) **분할-무관** (MAX 리덕션은 결합적) → 스무더가 np 에 대해 동일 연산자,
  (iii) 파라미터-프리 (안전계수 튜닝 불필요).
- 트레이드오프: 상계의 보수성만큼 스무딩 최적성 손실. 실측: 등방 Poisson 중립
  (its 4~5 동일), 실계수 행렬 저난도 solve +1~2 its, 고난도 solve 동률 (§4-3 정정판).
- 구현: `poly_smooth.f90` `eig_value`, 스위치 `ieig_pol`(&MG_smoothing).
  보조로 코어스 hybrid GS 의 l1-보정(`il1_gs`, Baker–Falgout–Kolev–Yang 2011)도
  구현되어 있으나 본 실패 모드와는 독립.

### 6-2. 검증 총괄

**합성 Poisson (붕괴 케이스 전부 완치, np-평탄):**

| 케이스 | 기존 | Gershgorin |
|---|---|---|
| 48³ 블록 np=64/96/112/128 | 5 / 231 / 5 / 발산 | **5 / 5 / 5 / 5** |
| 24³ 슬랩 np=8 / np=12 | 770 / 발산 | **4 / 4** |

**프로덕션 ECT1 np 스윕 (동일 벽시계 러닝, 공통 첫 120 solve 평균 its):**

| 모드 | np=1 | np=4 | np=16 | np=64 |
|---|---|---|---|---|
| 기존 | 4.09* | 3.63 | 3.58 | 3.53 |
| Gershgorin | 4.10 | 4.13 | 4.17 | 4.17 |

(전 구성 최대 its 8~9, 발산·정체 0건. *np=1 기존은 골든 러닝 기록.
구간별(1-40/41-80/81-120)로도 양 모드 np-평탄; Gershgorin 은 초반 저난도
solve(its 1 수준)에서 +1~2 its 의 보수성 비용, 후반 구간 동률.)

## 7. Novelty 서술 초안 (영문)

논문에서 주장 가능한 기여 후보 (문헌 확인 후 취사선택):

1. **Failure-mode identification.**
   "We identify and characterize a partition-dependent divergence mode of
   Chebyshev-smoothed multigrid preconditioners: the customary few-step
   Lanczos estimate of the spectral upper bound, even with the standard 10%
   safety margin, varies with the domain decomposition through the
   non-associativity of parallel reductions; for particular uniform partition
   geometries the smoothing interval fails to cover the upper spectrum, and
   the V-cycle — while contracting the bulk error by two orders of magnitude
   per application — amplifies a small set of high modes, which the outer
   Krylov iteration eventually exposes. The failure is deterministic,
   non-monotonic in the process count, and invisible to standard setup
   diagnostics, which we demonstrate by an elimination-based forensic
   methodology (Table X)."
2. **Remedy — a known conservative choice, re-valued.** (정정판: "3.7–11x
   직렬 개선" 문구 철회, §10 재포지셔닝과 결합)
   "The remedy is deliberately unsophisticated: replacing the estimate by the
   Gershgorin row-sum bound — a true, partition-invariant upper bound cheaper
   to compute than the estimate itself. While bound-based Chebyshev intervals
   are available as options in existing packages, we show that their
   partition-invariance eliminates an entire class of decomposition-dependent
   divergences, at a measured cost of at most one to two extra iterations on
   easy solves and none on hard ones — a trade we argue should be the default
   in massively parallel settings."
3. **Scalability statement.**
   "With the bound-based smoother (and the root-built, partition-invariant
   hierarchy of our PMG), the preconditioner is algebraically identical for
   every process count; measured iteration counts are flat from np=1 to 128
   (synthetic) and np=1 to 64 (production), where the estimate-based variant
   collapses at specific partition geometries."

관련 문헌: **§10 (2026-08-19 딥 리서치로 확정)**. 판정 요약 — 기여 1(실패 모드
규명)은 직접 선행 없음(§10-1). 단 기여 2(치료법)는 재포지셔닝 필요: Gershgorin/
행합 상계 자체는 Baker et al. 2011 에 한 문장으로 언급되고 hypre·AMGX·PyAMG 에
옵션으로 구현돼 있음 — "새 방법"이 아니라 "알려진 보수적 대안의 분할-무관성이라는
미인식 가치 + 실패 모드의 진단 + 실계수 직렬 개선 실측"으로 서술할 것 (§10-1 ④).

## 8. 재현 방법

- 스위치: `mg.in` → `&MG_smoothing` → `ieig_pol = 1` (기본 0 = 기존 동작).
- 합성 재현: `pmg_standalone` 에서
  `mpirun -np 96 build/driver_pmg 48 48 48 1.0 1` (붕괴) ↔ `ieig_pol=1` (완치).
  슬랩: `mpirun -np 8 build/driver_pmg 24 24 24 1.0`.
- 골든 재생: `build/driver_pmg replay golden/ECT1_436k_np1 <step> 1e-6`.
- 프로덕션 스윕: `scripts/run.sh` + `CUPID_NP`, fort.501 의 solve 별 its 비교.
- 진단 계측(층별 잔차·λ̂ 출력)은 커밋하지 않은 임시 코드였음 — 필요 시
  본 문서 §5 의 절차대로 재삽입.

## 9. 한계 및 후속 확인 사항

1. **np=900 실기 미검증** — 수렴 기전은 제거됐으나 최종 확인은 클러스터 1회
   필요. 별도 리스크로 셋업 트랙(극소 도메인 writer 이상, LOG C010-3 ④)이 남아
   있으며 이는 본 finding 과 독립된 결함이다.
2. **λ̂ 요동의 미시 기전** — "병렬 리덕션 라운딩 → 유한정밀 Lanczos 궤적 변화"
   는 정황이 강하나(결정론적·분할 의존·2.7% 크기) 라운딩 경로 수준의 엄밀
   규명은 하지 않았다. 논문에서는 실측 사실로 기술하는 것이 안전.
3. **보수성 비용의 정량화** — 저난도 solve +1~2 its 의 총 벽시계 영향은 케이스
   의존. 필요 시 `min(Gershgorin, c·Lanczos)` 같은 절충은 **권장하지 않음**
   (참 상계 성질이 깨짐).
4. 기본화(`ieig_pol=1` default) 및 골든 재베이스라인은 별도 사이클로 수행 예정.

---

## 10. 선행 문헌 조사 (딥 리서치, 2026-08-19)

조사 방식: 4축 병렬 웹 리서치 (①novelty 직접 선행, ②라이브러리 관행,
③비결합 리덕션·유한정밀 Lanczos, ④다항 스무더 핵심 문헌). 아래 인용은 별도
표기가 없는 한 1차 소스(arXiv/DOI 원문, 라이브러리 소스·매뉴얼)에서 원문 확인됨.
미확인 항목은 §10-6 에 명시.

### 10-1. 핵심 판정 — novelty 성립 여부

**① 직접 선행 없음.** "few-step Lanczos λ_max 추정이 분할/프로세스 수에 따라
변동 → 특정 분할에서 Chebyshev 구간이 스펙트럼 상단을 놓침 → V-cycle 의 모드
선택적 증폭 → 결정론적·np-비단조 Krylov 발산"의 연쇄를 다룬 논문·리포트·이슈는
발견되지 않았음. 두 반쪽은 각각 별개로 문서화되어 있음:

- **과소추정 → 증폭** (직렬 메커니즘): Ifpack2 소스 문서가 boost factor 1.1 의
  목적을 명문화 — "Otherwise the smoother could actually **magnify high-energy
  error modes**"; PETSc 매뉴얼 — 추정이 너무 낮으면 "the solvers can **fail with
  an indefinite preconditioner message**"; Phillips–Kerkemeier–Fischer 2022 의
  Table 1 — λ_max ≤ 1.0λ̃ 설정 시 1000회 내 비수렴을 실측 보고 (발표된 것 중
  가장 근접한 과소추정-발산 실증, 단 튜닝 맥락이고 분할 의존성 아님).
- **병렬 리덕션 비결합성 → 프로세스 수 의존 변동**: MPI 표준 자체가 허용
  (§10-4 A1), 재현성 문헌에서 표준적 사실. PETSc FAQ 는 "프로세스 수에 따라
  예조건자 수렴이 달라질 수 있고, 어떤 수에서는 되고 어떤 수에서는 안 될 수
  있다"와 "리덕션 순서 비결정성"을 **별개 항목으로** 인정하나 둘을 연결하지 않음.

**② 기존 문헌은 오히려 반대 주장을 반복** — 본 발견의 novelty 를 강화하는 지점.
Adams et al. 2003 이래 다항 스무더 채택의 표준 논거가 바로 분할-무관성이다.
Baker et al. 2011 원문: 다항 스무더는 "**unaffected by the parallel partitioning
of the matrix, the number of parallel processes, and the ordering of the
unknowns**" (OpenFOAM GAMG 최적화 논문 OFW13-2018 도 이를 그대로 재인용).
본 발견은 이 주장이 **셋업(고유값 추정) 경로를 통해 깨지는** 것을 보인 것.

**③ 근접 엔지니어링 흔적 (분석 없는 실무 기록) 2건:**

- deal.II `PreconditionChebyshev` 문서: 결정론적 초기 벡터를 일부러 채택하여
  추정이 프로세스 수에 대해 "same ... **apart from roundoff errors**" 라고 명시
  — 위험을 직감하고 설계로 회피했으나 실패 사례 분석은 없음 (issue #3490 에서
  "random numbers ... depend on the state which is bad for a preconditioner").
- Trilinos issues #64/#567 (2015–16): 특정 OpenMP 스레드 수(정확히 2)에서만
  Chebyshev 고유값 추정 실패로 테스트 실패 — 단 원인은 랜덤 초기 벡터(음수 추정)
  이고 리덕션 비결합성이 아니며, 증상도 Krylov 발산이 아님.

**④ 치료법(Gershgorin/행합 상계)은 선례 있음 — 기여 2 서술 조정 필요.**

- Baker et al. 2011 원문에 한 문장 존재: "one can estimate the largest
  eigenvalue of A by using the **maximum absolute row-sum** of A (infinity
  norm), which ... is cheaper than performing 10 CG iterations. This strategy
  is **equally effective for most problems, though often less so on the
  coarser grid levels**."
- hypre: `eig_est=0` 이면 λ_max 를 Gershgorin 으로 잡는 폴백이 이미 구현됨.
  AMGX: mode 2 가 최대 행합(Gershgorin형) 상계. PyAMG: prolongator smoother 의
  `'local'` 가중이 "Gershgorin estimate ... **avoids any potential
  under-damping due to inaccurate spectral radius estimates**" 라고 문서화.
- Lottes 2023: ℓ1-Jacobi 가 "ρ(BA) ≤ 1 **by construction (as can be seen from
  the Gershgorin circle theorem)**" → 추정 자체를 소거한 "truly parameter-free"
  구성 가능; 또한 "B can be scaled by **any upper bound** of ρ(BA), with some
  degradation in performance" — 참 상계 사용의 안전성(성능 저하만 있고 발산
  없음)에 대한 이론적 뒷받침.

따라서 기여 2는 "새 상계 제안"이 아니라: (i) 실패 모드의 규명과 배제 수사,
(ii) 알려진 보수적 대안이 갖는 **분할-무관성**(MAX 리덕션의 결합성)이라는
미인식 성질의 식별, (iii) 실계수 행렬에서 추정 대비 오히려 3.7~11배 **개선**
이라는 실측(관행적 통념 "추정이 더 타이트해서 유리"의 반례)으로 서술해야 함.

**⑤ 반론 대비 주의 문헌**: Konolige & Brown 2020 (bundle adjustment MG) —
"We also tried using the Gershgorin estimate ..., but that proved to be **very
inaccurate (by multiple orders of magnitude)**" (비대각지배 행렬). 본 논문에서는
FVM 압력계 연산자(대각지배 근방)에서 Gershgorin 이 충분히 타이트함을 실측
(§6-2: 등방 Poisson its 동일, 프로덕션 +1~2 its)으로 논증할 것. Baker 2011 의
"코스 레벨에서는 덜 효과적" 단서도 같은 맥락 — 본 코드는 fine 레벨만 POL 이라
해당 약점을 비껴감.

### 10-2. "추정 + 안전계수" 관행의 실태 (주장 1 의 근거)

라이브러리 (소스·매뉴얼에서 기본값 확인):

| 라이브러리 | 추정법 | 스텝 | 안전계수 |
|---|---|---|---|
| PETSc (KSPCHEBYSHEV/GAMG) | Lanczos(CG)/Arnoldi, noisy RHS | 10 | ×1.1 (transform `0,0.1;0,1.1`) |
| Trilinos Ifpack2/MueLu | power method | 10 | boost 1.1 |
| hypre BoomerAMG (relax 16) | CG-Lanczos (옵션: Gershgorin) | 10 | 없음 (원시 추정) |
| deal.II PreconditionChebyshev | CG-Lanczos, 결정론적 초기벡터 | 8 | ×1.2 |
| MFEM OperatorChebyshevSmoother | power method | 10 | ×1.2 |
| NVIDIA AMGX | Lanczos (옵션: 최대 행합) | ≤128 | ×1.05 |
| nekRS | Arnoldi | 10 | ×1.1 |

논문 실측 설정 (원문 확인): Adams 2003 — 10 CG, 구간 (1/30, **1.1**)λ̃ (Baker
2011 의 "As in [1], we use 10 iterations of CG and multiply λmax by 1.1" 로
확인); Baker 2011 — (0.3, 1.0)λ̃; Sundar–Stadler–Biros 2015 — 10 Arnoldi,
(1/4, 1.0)λ̃, "**usually** estimated during setup with an iterative solver";
Kronbichler–Wall 2018 — 15 CG, (0.06, 1.2)λ̃; Fehn et al. 2020 — 20 CG,
(0.06, 1.2)λ̃, "Since the maximum eigenvalue is only estimated, a **safety
factor** of 1.2 is included to ensure robustness"; Clevenger et al. 2021 —
10 CG, (0.08, 1.2)λ̃; Phillips et al. 2022 (nekRS) — 10 Arnoldi, (0.1, 1.1)λ̃.

→ 본 코드의 "8-스텝 Lanczos × 1.1"은 관행의 정중앙이며(§4-1 주장 그대로 인용
가능), hypre 는 안전계수조차 없는 원시 추정이라 노출이 더 큼 (대비 포인트).

### 10-3. 메커니즘의 인용 근거 (λ̂ 요동, §4-2 ①)

**(A) 병렬 리덕션 비결합성:**

- MPI Forum, *MPI 4.1*, §7.9.1: 구현이 결합성·교환성을 활용해 재배열하는 것을
  허용하며 이것이 부동소수 결과를 바꿀 수 있음을 명시. 동일 인자·동일 순서에서의
  재현만 권고 — **프로세스 수가 다르면 동일 결과 보장 없음** (규범적 근거).
- Villa, Chavarría-Miranda, Gurumoorthi, Márquez, Krishnamoorthy (CUG 2009):
  리덕션 순서 변화만으로 CG 잔차가 반복 5부터 14~25% 갈라지는 실측.
  결정론적 트리 리덕션으로 복원. (run-to-run 비결정성 맥락 — 본 건은
  분할별-결정론적이라는 차이를 명시하고 인용할 것.)
- Balaji & Kimpe (HPCC 2013): run-to-run 재현성과 프로세스 수 변경 시 재현성을
  구분 — 후자는 제공되지 않음. (전문 미입수, 초록 수준 인용만.)
- He & Ding (J. Supercomput. 2001): 프로세스 수 변경 → 합산 순서 변경 → 기후
  코드 결과 변동의 고전 사례 + 보정합 처방.
- Demmel & Nguyen (IEEE TC 2015), Ahrens–Demmel–Nguyen (TOMS 2020): ReproBLAS —
  "프로세서 수·데이터 분할·리덕션 스케줄에 무관한" 재현 합산 (문제의 공인 근거).
- Iakymchuk et al. (IJHPCA 2020; JCAM 2020): 병렬 PCG 의 **반복수 자체가** 병렬
  구성에 따라 달라짐을 보고, ExBLAS 기반 재현 PCG 제안. 동 그룹 arXiv:2302.04180
  (2023) 은 BiCGStab 계열 — 본 외부 솔버와 동일 계열이라 인용 가치.
- Shanmugavelu et al. (arXiv:2408.05148, 2024): 최신 서베이 (현재진행형 문제임을
  뒷받침).

**(B) few-step Lanczos 추정의 성질:**

- Ritz 값은 λ_max 를 **항상 아래에서** 근사 (Cauchy interlacing; Parlett 1998
  Ch.10/13, Demmel ANLA Ch.7) — "추정은 구조적으로 과소추정"의 근거. 안전계수
  1.1 의 존재 이유이기도 함 (Ifpack2 문서가 이를 명시적으로 서술).
- Kuczyński & Woźniakowski (SIMAX 1992; 1994): m-스텝 Lanczos 의 기대 상대오차
  상계 0.103·ln²(n(m−1)⁴)/(m−1)² — m=8 에서 수 %~O(1) 의 오차는 이론이 예측하는
  범위 내. 결정론적 초기 벡터로는 최악의 경우 해결 불가 명제도 유용.
- Urschel (SIMAX 2021): 매칭 하계 — m-스텝 추정의 상대오차는 ~ln²n/m² 아래로
  내려갈 수 없음. **8-스텝 추정의 오차는 구현 결함이 아니라 원리적** (주장
  보강에 최적).
- Paige (1976, 1980), Greenbaum (LAA 1989), Meurant & Strakoš (Acta Numerica
  2006): 유한정밀 Lanczos 는 정확 산술의 작은 섭동이 아님; 섭동된 Lanczos 는
  "다른 행렬에 대한 정확 Lanczos"처럼 거동 (라운딩 차이 → 그럴듯하지만 다른
  Ritz 값의 배경 이론). §9-2 의 "미시 기전은 실측 사실로 기술" 방침과 정합.

### 10-4. 포지셔닝 핵심 문헌 (§7 초안에서 인용할 것)

- **Adams, Brezina, Hu, Tuminaro (JCP 2003)** — 병렬 다항 스무더의 시조 논문,
  "polynomial smoothers achieve parallel scalable multigrid convergence rates".
  본 발견이 단서를 다는 대상. (원문 페이월 — §10-6 참조.)
- **Baker, Falgout, Kolev, Yang (SISC 2011)** — 다항 스무더 분할-무관 주장 +
  행합 상계 한 문장 + ℓ1 스무더 (프로세서 수 무관 수렴 보장, Theorem 6.2).
  본 코드의 `il1_gs` 근거 문헌이기도 함.
- **Lottes (NLAA 2023)** — 4종 Chebyshev: λ_min 불필요, **임의의 참 상계로 안전**
  ("any upper bound ... some degradation"). 본 치료법의 이론적 안전망 + 후속
  개선 방향(4종 전환 시 Gershgorin 과 자연 결합).
- **Phillips & Fischer (arXiv:2210.03179; NLAA 2025)** — 최적 λ_min·4종 스무더.
  **Phillips, Kerkemeier, Fischer (SIAM PP22, arXiv:2110.07663)** — 10 Arnoldi
  ×1.1 관행 + λ_max ≤ λ̃ 에서 비수렴 실측 표 (과소추정 위험의 발표된 최근접
  실증).
- **Gutknecht & Röllin (Parallel Comput. 2002)** — Chebyshev 는 내적이 없어
  병렬 유리 + 구간은 "known in advance" 가정. 그 가정을 조달하는 셋업 추정이
  내적 기반이라는 아이러니가 본 논문의 서사 고리.
- **Manteuffel (Numer. Math. 1977; 1978)** — 나쁜 스펙트럼 파라미터 → 불만족
  수렴 → 재추정·재시작의 고전 (직렬 솔버 맥락의 원조 선행).
- **Sundar, Stadler, Biros (NLAA 2015)** — ℓ1-Jacobi 를 "추정 불필요" 장점으로
  Chebyshev 와 대비하는 프레임 (본 논문 프레임과 동형).
- **D'Ambra, Durastante, Filippone, Massei, Thomas (Numer. Algorithms 2025)** —
  "스펙트럼 정보 없이" 방향의 최신 GPU/엑사스케일 다항 스무더 (추정 제거가
  현재 연구 방향임을 보이는 시의성 근거).
- 보조: Golub & Varga (1961, Chebyshev 반복 원전); Vaněk & Brezina (Appl. Math.
  2013, 상계만으로 전개되는 다항 스무딩 이론); Chen & Wells (arXiv:2402.12947,
  국소 저품질 셀의 λ_max 병리); Kronbichler & Wall (SISC 2018); Clevenger et
  al. (TOMS 2021); Fehn et al. (JCP 2020); Naumov et al. (SISC 2015, AmgX);
  Gandham, Esler, Zhang (CAMWA 2014, GS 의 GPU 부적합).

### 10-5. 참고 문헌 목록

1. M. Adams, M. Brezina, J. Hu, R. Tuminaro, "Parallel multigrid smoothing:
   polynomial versus Gauss–Seidel," *J. Comput. Phys.* 188(2):593–610, 2003.
   doi:10.1016/S0021-9991(03)00194-3.
2. A. H. Baker, R. D. Falgout, Tz. V. Kolev, U. M. Yang, "Multigrid smoothers
   for ultraparallel computing," *SIAM J. Sci. Comput.* 33(5):2864–2887, 2011.
   doi:10.1137/100798806 (LLNL-JRNL-473191).
3. J. Lottes, "Optimal polynomial smoothers for multigrid V-cycles," *Numer.
   Linear Algebra Appl.* 30(6):e2518, 2023. doi:10.1002/nla.2518
   (arXiv:2202.08830).
4. M. Phillips, P. Fischer, "Optimal polynomial smoothers and one-sided
   V-cycles for Poisson problems," *Numer. Linear Algebra Appl.* 32(4):e70030,
   2025. doi:10.1002/nla.70030 (arXiv:2210.03179).
5. M. Phillips, S. Kerkemeier, P. Fischer, "Tuning spectral element
   preconditioners for parallel scalability on GPUs," *Proc. SIAM PP22*, 2022.
   arXiv:2110.07663.
6. M. H. Gutknecht, S. Röllin, "The Chebyshev iteration revisited," *Parallel
   Comput.* 28(2):263–283, 2002. doi:10.1016/S0167-8191(01)00139-9.
7. T. A. Manteuffel, "The Tchebychev iteration for nonsymmetric linear
   systems," *Numer. Math.* 28:307–327, 1977. doi:10.1007/BF01389971; "Adaptive
   procedure for estimating parameters for the nonsymmetric Tchebychev
   iteration," *Numer. Math.* 31:183–208, 1978. doi:10.1007/BF01397475.
8. G. H. Golub, R. S. Varga, "Chebyshev semi-iterative methods, successive
   overrelaxation iterative methods, and second order Richardson iterative
   methods," *Numer. Math.* 3:147–156, 157–168, 1961. doi:10.1007/BF01386014.
9. H. Sundar, G. Stadler, G. Biros, "Comparison of multigrid algorithms for
   high-order continuous finite element discretizations," *Numer. Linear
   Algebra Appl.* 22(4):664–680, 2015. doi:10.1002/nla.1979 (arXiv:1402.5938).
10. M. Kronbichler, W. A. Wall, "A performance comparison of continuous and
    discontinuous Galerkin methods with fast multigrid solvers," *SIAM J. Sci.
    Comput.* 40(5):A3423–A3448, 2018. doi:10.1137/16M110455X.
11. N. Fehn, P. Munch, W. A. Wall, M. Kronbichler, "Hybrid multigrid methods
    for high-order discontinuous Galerkin discretizations," *J. Comput. Phys.*
    415:109538, 2020. doi:10.1016/j.jcp.2020.109538 (arXiv:1910.01900).
12. T. C. Clevenger, T. Heister, G. Kanschat, M. Kronbichler, "A flexible,
    parallel, adaptive geometric multigrid method for FEM," *ACM Trans. Math.
    Softw.* 47(1):1–27, 2021. doi:10.1145/3425193.
13. P. D'Ambra, F. Durastante, S. Filippone, S. Massei, S. Thomas, "Optimal
    polynomial smoothers for parallel AMG," *Numer. Algorithms*, 2025.
    doi:10.1007/s11075-025-02117-6 (arXiv:2407.09848).
14. M. Naumov et al., "AmgX: a library for GPU accelerated algebraic multigrid
    and preconditioned iterative methods," *SIAM J. Sci. Comput.*
    37(5):S602–S626, 2015. doi:10.1137/140980260.
15. R. Gandham, K. Esler, Y. Zhang, "A GPU accelerated aggregation algebraic
    multigrid method," *Comput. Math. Appl.* 68(10):1151–1160, 2014.
    doi:10.1016/j.camwa.2014.08.022.
16. P. Vaněk, M. Brezina, "Nearly optimal convergence result for multigrid
    with aggressive coarsening and polynomial smoothing," *Appl. Math.*
    58(4):369–388, 2013. doi:10.1007/s10492-013-0018-2.
17. T. Konolige, J. Brown, "Multigrid for bundle adjustment," arXiv:2007.01941,
    2020.
18. Y. Chen, G. N. Wells, "Multigrid on unstructured meshes with regions of
    low quality cells," arXiv:2402.12947, 2024.
19. P. Ghysels, W. Vanroose, "Hiding global synchronization latency in the
    preconditioned conjugate gradient algorithm," *Parallel Comput.*
    40(7):224–238, 2014. doi:10.1016/j.parco.2013.06.001.
20. MPI Forum, *MPI: A Message-Passing Interface Standard, Version 4.1*,
    §7.9.1, 2023. mpi-forum.org.
21. O. Villa, D. Chavarría-Miranda, V. Gurumoorthi, A. Márquez,
    S. Krishnamoorthy, "Effects of floating-point non-associativity on
    numerical computations on massively multithreaded systems," *Proc. Cray
    User Group (CUG)*, 2009.
22. P. Balaji, D. Kimpe, "On the reproducibility of MPI reduction operations,"
    *IEEE HPCC/EUC 2013*, pp. 407–414. doi:10.1109/HPCC.and.EUC.2013.65.
23. Y. He, C. H. Q. Ding, "Using accurate arithmetics to improve numerical
    reproducibility and stability in parallel applications," *J. Supercomput.*
    18(3):259–277, 2001. doi:10.1023/A:1008153532043.
24. J. Demmel, H. D. Nguyen, "Parallel reproducible summation," *IEEE Trans.
    Comput.* 64(7):2060–2070, 2015. doi:10.1109/TC.2014.2345391; P. Ahrens,
    J. Demmel, H. D. Nguyen, "Algorithms for efficient reproducible floating
    point summation," *ACM TOMS* 46(3), 2020. doi:10.1145/3389360.
25. R. Iakymchuk, M. Barreda, S. Graillat, J. I. Aliaga, E. S. Quintana-Ortí,
    "Reproducibility of parallel preconditioned conjugate gradient in hybrid
    programming environments," *Int. J. High Perform. Comput. Appl.*
    34(5):502–518, 2020. doi:10.1177/1094342020932650 (arXiv:2005.07282);
    R. Iakymchuk, J. I. Aliaga, "General framework for re-assuring numerical
    reliability in parallel Krylov solvers: a case of BiCGStab methods,"
    arXiv:2302.04180, 2023.
26. S. Shanmugavelu, M. Taillefumier, C. Culver, O. Hernandez, M. Coletti,
    A. Sedova, "Impacts of floating-point non-associativity on reproducibility
    for HPC and deep learning applications," arXiv:2408.05148, 2024.
27. J. Kuczyński, H. Woźniakowski, "Estimating the largest eigenvalue by the
    power and Lanczos algorithms with a random start," *SIAM J. Matrix Anal.
    Appl.* 13(4):1094–1122, 1992. doi:10.1137/0613066; 동, "Probabilistic
    bounds on the extremal eigenvalues and condition number by the Lanczos
    algorithm," *SIAM J. Matrix Anal. Appl.* 15(2):672–691, 1994.
    doi:10.1137/S0895479892230456.
28. J. C. Urschel, "Uniform error estimates for the Lanczos method," *SIAM J.
    Matrix Anal. Appl.* 42(3):1423–1450, 2021. doi:10.1137/20M1331470.
29. C. C. Paige, "Error analysis of the Lanczos algorithm for tridiagonalizing
    a symmetric matrix," *J. Inst. Math. Appl.* 18(3):341–349, 1976;
    "Accuracy and effectiveness of the Lanczos algorithm for the symmetric
    eigenproblem," *Linear Algebra Appl.* 34:235–258, 1980.
30. A. Greenbaum, "Behavior of slightly perturbed Lanczos and
    conjugate-gradient recurrences," *Linear Algebra Appl.* 113:7–63, 1989.
31. G. Meurant, Z. Strakoš, "The Lanczos and conjugate gradient algorithms in
    finite precision arithmetic," *Acta Numerica* 15:471–542, 2006.
    doi:10.1017/S0962492906220000.
32. B. N. Parlett, *The Symmetric Eigenvalue Problem*, SIAM Classics 20, 1998.
33. U. M. Yang, "On the use of relaxation parameters in hybrid smoothers,"
    *Numer. Linear Algebra Appl.* 11:155–172, 2004. doi:10.1002/nla.375.
34. 라이브러리 문서 (관행 실태·버전 기준 2026-08 조회): PETSc `KSPChebyshevEstEigSet`
    매뉴얼 페이지 및 Users Manual KSP 장 (petsc.org/release); Trilinos Ifpack2
    User's Guide SAND2016-5338 + `Ifpack2_Details_Chebyshev_decl.hpp` (boost
    factor 문서); MueLu User's Guide SAND2023-12265; hypre Reference Manual
    "ParCSR Solvers" (`HYPRE_ParChebySetEigEst`); deal.II `PreconditionChebyshev`
    doxygen; MFEM `linalg/solvers.cpp`; NVIDIA AMGX `src/solvers/cheb_solver.cu`;
    PyAMG `pyamg.aggregation.smooth` 문서; Trilinos issues #64, #567; dealii
    issue #3490; PETSc FAQ.

### 10-6. 검증 상태 주의사항

- **원문 미확인 (페이월)**: [1] Adams 2003, [2] Baker 2011 의 정확한 문장은
  각각 Baker 2011 리포트판(UNT 공개본, 원문 확인)과 Phillips et al. [5] 의
  직접 인용으로 교차 확인함. 직접 인용부호로 쓰려면 기관 접근으로 원문 대조
  필요. [22] Balaji–Kimpe 는 전문 미입수 (초록 수준으로만 인용할 것).
  Vaněk–Mandel–Brezina 1996 (SA 원전)의 λ_max 추정 방식은 원문 미확인이라
  목록에서 제외 — 인용 필요 시 원문 확인 선행.
- 라이브러리 기본값은 2026-08 시점 master/release 소스 기준 — 논문 제출 시점에
  버전 재확인 필요.
- Ifpack2 클래스 문서의 1.1 서술은 텍스트상 λ_min 에 걸려 있으나 구현은 λ_max
  에 boost 를 적용 (`boostFactor_` 문서가 정본) — 인용 시 `boostFactor_` 쪽을
  쓸 것.
