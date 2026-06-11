# 논문 초안: 비정렬 격자 PMG 예조건 BiCGSTAB 솔버와 MSFR·iSMR 적용

> 작성일: 2026-06-10 (초안 v0.1)
> 상태: 골격 + 서론/관련연구/방법론 초안. 결과 섹션은 실험 수행 후 채움.
> 관련 노트: [[문헌 조사 - PMG 예조건화 BiCGSTAB]]

---

# 0. 메타 정보

- **가제(영문)**: *A Parallel Semi-Coarsening Geometric Multigrid Preconditioner for BiCGSTAB on Unstructured Grids: Application to MSFR Steady State and iSMR Natural Circulation*
- **가제(국문)**: 비정렬 격자에서의 병렬 반-조밀화 기하적 멀티그리드 예조건 BiCGSTAB 솔버: MSFR 정상상태 및 iSMR 자연순환 적용
- **대상 저널**: *Annals of Nuclear Energy* / *Nuclear Engineering and Design* / *Journal of Computational Physics* (계산 원자력 공학)
- **핵심 기여 (Contributions)**:
  1. 비정렬 격자에 적용 가능한 반-조밀화(semi-coarsening) 기하적 MG를 **BiCGSTAB의 예조건자(preconditioner)** 로 통합. 선행 meshless GMG 연구[@do2024meshless]가 독립 솔버였던 데 비해, 본 연구는 예조건자로 확장하여 강건성과 가속을 동시에 확보.
  2. **통신량 저감** 및 **일정 레벨 이하 serial 계산** 전략으로 병렬 확장성 확보.
  3. **MSFR 정상상태**, **iSMR ECT 자연순환**, **iSMR 노심 일일부하추종 과도** 문제에서 속도 가속·강건성 실증.

---

# 1. 서론

## 1-1. 배경

원자로 열수력 해석에서 비압축성 유동의 압력 보정 방정식(pressure correction / Poisson 형) 풀이는 전체 계산 시간의 지배적 비중을 차지한다. 이 타원형 방정식은 조건수가 크고 격자 세분화에 따라 반복 솔버의 수렴이 급격히 악화되어, 대규모 3차원 해석의 병목이 된다.

KAERI의 CUPID 코드는 비정렬 격자 기반 2유체-3유동장 모델로 원자로 노심·계통을 해석하며, 압력 보정 방정식에 전통적으로 BiConjugate Gradient(BiCG) 계열 솔버를 사용해 왔다[@do2024meshless]. 그러나 BiCG 계열은 격자 크기에 따른 반복 횟수 증가(비-scalable)라는 본질적 한계를 가진다.

멀티그리드(MG)는 이론적으로 격자 크기에 무관한 수렴률을 제공하는 거의 유일한 방법군이다. 선행 연구[@do2024meshless; @do2022highaspect; @ha2024aiaa]에서 저자 그룹은 비정렬 격자에 직접 적용 가능한 **meshless 기하적 멀티그리드(GMG)** 를 개발하여 BiCG 대비 우수한 확장성을 보였다. 본 연구는 이를 **예조건자로 재구성**하고, **병렬 통신 저감** 및 **MSFR·iSMR** 응용으로 확장한다.

## 1-2. 동기: 왜 MG 예조건 BiCGSTAB인가

- **독립 MG 솔버의 한계**: 기하적 MG는 평활화(smoothing)·격자 전이 연산자가 문제에 민감하여, 강한 비등방·불규칙 비정렬 격자에서 단독 사용 시 수렴이 정체될 수 있다.
- **Krylov 가속의 효과**: MG를 BiCGSTAB의 예조건자로 쓰면, MG가 저주파 오차를 효율적으로 제거하고 BiCGSTAB이 MG가 놓친 잔여 스펙트럼 성분을 처리하여 강건성이 크게 향상된다.
- **가변 예조건자 문제**: 대규모 병렬에서 통신을 줄이기 위해 MG 계층을 제한하고 최조밀 격자를 반복법으로 근사하면 예조건자가 매 반복마다 달라지는 **가변(variable) 예조건자**가 된다. 이 경우 표준 BiCGSTAB은 수렴이 깨질 수 있어 **flexible 변형(FBiCGSTAB)** 이 필요하다[@chen2016fbicgstab].

## 1-3. 목표

본 연구는 (1) 비정렬 격자용 semi-coarsening 기하적 MG 예조건자를 BiCGSTAB(필요 시 FBiCGSTAB)에 통합하고, (2) 통신량 저감과 일정 레벨 이하 serial 계산을 구현하며, (3) MSFR 정상상태, iSMR ECT 자연순환, iSMR 노심 일일부하추종 과도 문제에서 가속 성능을 정량적으로 입증하는 것을 목표로 한다.

---

# 2. 관련 연구

## 2-1. 비정렬 격자 기하적 멀티그리드 (저자 그룹 선행 연구)

본 연구가 직접 확장하는 **scalable solver** 계보는 다음과 같다.

- **Meshless GMG (ANE 2024)**[@do2024meshless]: CUPID 코드에 node-coarsening 기반 meshless 기하적 MG를 도입. 정렬·비정렬 격자 모두에서 안정적으로 수렴하며 BiCG 대비 계산 시간을 크게 단축. **본 논문의 출발점이며, 비정렬 격자 적용 전략의 근거.**
- **High aspect ratio GMG (JMST 2022)**[@do2022highaspect]: 종횡비가 큰(신장된) 비등방 격자에서의 meshless GMG. 경계층·박막 유동에 해당하며, semi-coarsening 동기와 직접 연결.
- **Complex geometry cell-coarsening (AIAA J.)**[@ha2024aiaa]: 복잡 형상에 대한 개선된 cell-coarsening 알고리즘. MSFR·iSMR 같은 복잡 형상 적용의 기반.
- **Node-coarsening for linear FE (2021)**[@ha2021nodecoarsening]: 선형 유한요소 이산화에 대한 node-coarsening 알고리즘의 기초.
- **Optimal aggregation level for parallel meshless MG (DDM)**[@ha2025aggregation]: 도메인 분할(DDM) 기반 병렬 meshless MG의 최적 aggregation 레벨 연구 — **본 논문의 "일정 레벨 이하 serial 계산" 및 통신 저감 전략과 직접 연결.**

## 2-2. Semi-coarsening 기하적 멀티그리드

Semi-coarsening은 좌표 방향 하나만 조밀화하여 신장 격자의 비등방 문제에 강건하다[@osti440722]. Mulder[@mulder1989]가 제안한 이래, MSG(Multiple Semi-Coarsened Multigrid)[@dendy2019msg]는 각 좌표 방향으로 조밀화한 추가 coarse grid를 더해 비등방 확산 처리를 개선했다. 비정렬 격자에서 좌표 정렬 방향이 모호하다는 점이 본 연구의 핵심 난제이자 기여 지점이다.

## 2-3. 병렬 MG 통신 저감 및 coarse-level 처리

대규모 병렬 MG의 중심 난제는 coarse-level 통신이다. 주요 전략:

- **Coarse-grid 재분배(incremental agglomeration)**[@reisner2018]: 예측 성능 모델로 프로세스를 점진 축소, 구조 격자 MG를 50만+ 코어로 확장.
- **Node-aware 다단계 통신**[@bienz2020]: inter-node 메시지의 개수·크기를 동시 저감, intra-node 채널 우선 활용.
- **Post-hierarchy sparsification**[@manteuffel2016]: 계층 형성 후 coarse 행렬 엔트리 제거로 통신 저감, 단 과도 시 수렴 저하 trade-off.

이들은 각각 구조 격자/AMG에서 검증되었으며, **비정렬 원자로 CFD 격자에서의 결합 적용은 미해결 영역**이다.

## 2-4. BiCGSTAB과 가변 예조건

BiCGSTAB[@vandervorst1992]은 비대칭계에서 CGS 대비 효율적인 표준 Krylov 솔버다. MG를 가변 예조건자로 쓰는 경우, 계층을 3레벨 수준으로 제한하고 최조밀 격자를 반복법으로 근사한 뒤 **FBiCGSTAB**을 적용하면 대규모에서 Chebyshev 예조건 IBiCGSTAB 대비 2~3배 가속이 보고되었다[@chen2016fbicgstab]. 예조건 연산자 스펙트럼이 소수 점에 클러스터될 때 Krylov 수렴이 가속된다는 이론적 배경도 있다[@lucerolorca2026].

## 2-5. MSFR·SMR 자연순환 수치해석

MSFR 결합 중성자/열수력 해석[@msfr_coupled]에서 연료염은 1차 루프를 약 3~4초에 순환하며, 이는 지연 중성자 선구물질(DNP) 수송과 직결되는 기준 설계 파라미터다. SMR 자연순환은 부력 구동 저레이놀즈·열성층 조건으로 MSFR 강제 난류와 정성적으로 다르며, 별도의 예조건 강건성 검토가 필요하다.

---

# 3. 방법론

## 3-1. 지배 방정식과 이산화

CUPID의 비압축성(또는 저마하) 2유체 모델에서 SIMPLE 계열 압력-속도 연성으로부터 유도되는 압력 보정 방정식:

$$
\nabla \cdot \left( \frac{1}{a_P} \nabla p' \right) = \nabla \cdot \mathbf{u}^*
$$

비정렬 유한체적 이산화로 비대칭(또는 약대칭) 희소 선형계 $A p' = b$ 를 얻는다. $A$는 격자 종횡비·비직교성에 따라 조건수가 악화된다.

## 3-2. Semi-coarsening 기하적 MG 예조건자 (비정렬 격자)

선행 meshless GMG[@do2024meshless; @ha2021nodecoarsening]의 node-coarsening을 기반으로, 비등방 방향을 감지하여 해당 방향으로만 조밀화하는 semi-coarsening 변형을 구성한다.

- **방향 감지**: 비정렬 격자에서 좌표축이 없으므로, 국소 격자 종횡비·연결 강도(행렬 계수 크기)로 "강결합 방향"을 추정한다.
- **격자 전이 연산자**: meshless 보간(restriction/prolongation)으로 정렬·비정렬 혼합 격자에서 작동.
- **평활화(smoother)**: 점/선 완화(point/line relaxation) 또는 다색(multi-color) Gauss-Seidel.

> ⚠️ 미해결: 좌표 정렬 방향이 정의되지 않는 폴리헤드럴 격자에서 semi-coarsening 방향 일관성 확보. [[문헌 조사 - PMG 예조건화 BiCGSTAB]] 4절 열린 질문 1 참조.

## 3-3. 병렬화: 통신 저감과 serial 계산 전환

- **통신 저감**: node-aware 다단계 메시징[@bienz2020] 개념을 적용, halo 교환의 inter-node 메시지 수를 축소. 격자 전이 시 중복 통신 제거.
- **일정 레벨 이하 serial 계산**: aggregation 레벨이 충분히 작아져(프로세스당 미지수 < 임계값) 통신이 계산을 지배하는 지점부터, coarse 문제를 단일(또는 소수) 프로세스로 모아 **serial 직접 해법**으로 처리. 최적 전환 레벨은 선행 DDM aggregation 연구[@ha2025aggregation]의 성능 모델로 결정.
- **대안**: coarse-grid 재분배(incremental agglomeration)[@reisner2018]와의 비교 평가.

## 3-4. BiCGSTAB / FBiCGSTAB 연계

- MG 예조건자가 **고정(fixed)** 인 경우(완전 V-cycle): 표준 BiCGSTAB 적용.
- MG 예조건자가 **가변(variable)** 인 경우(최조밀 격자 반복법 근사): **FBiCGSTAB** 적용[@chen2016fbicgstab]. 내부 MG 정확도와 외부 수렴의 trade-off를 실험적으로 결정.

### 알고리즘 개요 (FBiCGSTAB + PMG 예조건)

```
주어진 A, b, 초기값 x0
r0 = b - A x0,  r̂0 임의 (r̂0·r0 ≠ 0)
for k = 0, 1, 2, ...
    y_k = M_k^{-1} p_k        # M_k: PMG 예조건 (가변 가능)
    v_k = A y_k
    α = (r̂0·r_k) / (r̂0·v_k)
    s_k = r_k - α v_k
    z_k = M_k^{-1} s_k        # PMG 예조건 재적용
    t_k = A z_k
    ω = (t_k·s_k) / (t_k·t_k)
    x_{k+1} = x_k + α y_k + ω z_k
    r_{k+1} = s_k - ω t_k
    수렴 판정: ||r_{k+1}|| / ||b|| < tol
end
```

## 3-5. 검증 및 성능 지표

- **검증**: 제조해(MMS), 표준 벤치마크(lid-driven cavity, 채널 유동)로 정확도 확인.
- **성능 지표**: 반복 횟수의 격자 무관성, 강/약 확장성(strong/weak scaling), BiCG 단독 대비 wall-clock 가속비, 병렬 효율.

---

# 4. 적용 사례 (Case Studies)

세 적용 사례는 모두 **다차원 고속 압력 솔버**의 가치를 실증하기 위해 선정되었다. 핵심 공통 동기는 다음과 같다.

- **긴 물리 시간 척도의 시간전진**: 세 문제 모두 관심 현상이 수 초~수만 초(하루) 규모의 느린 수송·순환·피드백 시간 척도에 지배된다. 정상상태 도달(4-1, 4-2)이든 장시간 운전 과도(4-3)든, 작은 시간 간격으로 수많은 스텝을 시간전진하며 매 스텝마다 압력 보정 방정식을 풀어야 한다. 따라서 **단일 선형 풀이의 가속이 전체 해석 시간에 곱셈적으로 누적**되며, 압력 솔버를 가속하면 총 wall-clock이 직접 단축된다.
- **다차원(3D) 필연성**: 세 문제의 핵심 물리(연료염 내 DNP 이류-확산, 부력 구동 자연순환의 열성층·재순환, 출력 변동 시 노심 유동 재분배)는 본질적으로 3차원 유동장과 결합되어 있어 1D 계통 코드로 환원할 수 없다. 국소 유속·온도 분포가 안전 여유와 직결되므로 다차원 CFD 해석이 필수다.
- **설계·안전 평가를 위한 반복 해석**: 정상상태 해와 운전 과도 해는 설계 변수 스윕, 민감도 분석, 안전 인허가 평가에서 수십~수백 회 반복 계산된다. 고속 솔버는 이 반복 비용을 실용 가능한 수준으로 낮춰 설계 최적화 주기를 단축한다.

## 4-1. MSFR 정상상태 해석

**정상상태 다차원 고속 해석이 중요한 이유**: MSFR은 핵연료가 액체 상태로 1차 루프를 순환하는 독특한 구조로, 연료염이 노심을 약 3~4초에 통과하는 동안 지연 중성자 선구물질(DNP)이 활성 영역 밖으로 수송된다[@msfr_coupled]. 이 DNP의 공간 분포는 노심 형상·유속장에 강하게 의존하므로 점(point) 동특성이나 1D 모델로는 포착할 수 없고, **중성자속·온도·유속이 연성된 3차원 정상상태 분포**를 풀어야 유효 반응도와 출력 분포를 정확히 예측할 수 있다. 정상상태 운전점은 연료염 조성·유량·노심 형상 설계의 기준이 되며, 설계 반복마다 재계산되므로 압력 보정 솔버의 가속이 설계 평가 효율을 좌우한다.

- **대상**: 3000 MWth 기준 MSFR 1차 루프, 정상 난류[@msfr_coupled].
- **물리**: 연료염 순환 3~4초, DNP 수송에 의한 비대칭 이류 항 → 예조건 품질 영향 평가.
- **측정**: PMG-BiCGSTAB vs BiCG 단독의 압력 보정 수렴·총 계산 시간.

## 4-2. iSMR ECT 자연순환

**정상상태 다차원 고속 해석이 중요한 이유**: iSMR의 피동 잔열 제거는 펌프 없이 부력만으로 구동되는 자연순환에 의존하므로, 안전성 입증의 핵심은 자연순환이 충분한 제열 유량을 갖는 **정상 순환 상태로 안정적으로 수렴함**을 보이는 것이다. 자연순환 유량은 루프 전체의 밀도차·열성층·재순환 구조에 의해 결정되는 전역적·다차원적 평형이며, 약한 구동력 탓에 정상상태 도달이 느려 비정상 해석 시 매우 많은 시간 스텝이 필요하다. 또한 부력 구동 저레이놀즈 유동은 수렴이 더디고 불안정해, **강건하면서도 빠른 압력 솔버 없이는 정상상태 해석 자체가 비현실적으로 느려진다**. 따라서 이 문제는 고속·강건 솔버의 실용적 가치를 가장 잘 드러내는 사례다.

- **대상**: 혁신형 SMR(iSMR)의 비상 노심 냉각/제열(ECT) 자연순환 시나리오.
- **물리**: 부력 구동 저레이놀즈·열성층 → 강한 비등방 + 약한 대류. semi-coarsening 강건성 시험의 이상적 케이스.
- **측정**: 자연순환 정상상태 도달까지의 가속비, 강건성(발산 없는 수렴).

## 4-3. iSMR 노심 일일부하추종 유동 분포 과도 해석

**정상상태 단일 해석을 넘어선 "장시간 과도"의 가치**: 앞의 두 사례(4-1, 4-2)가 단일 정상상태 해의 고속화를 다룬다면, 본 사례는 **물리 시간 24시간 규모의 장시간 과도(transient)** 를 시간전진으로 해석하는 문제다. 이는 솔버 가속의 곱셈적 누적 효과가 가장 극적으로 드러나는 구성으로, 단일 압력 풀이의 미세한 가속이 수만~수십만 회의 시간스텝에 걸쳐 누적되어 전체 해석 시간을 결정한다. 따라서 4-1(예조건 품질)·4-2(강건성)에 더해 **장시간 시간전진에서의 누적 가속과 지속적 강건성**이라는 제3의 평가 축을 제공한다.

### 4-3-1. 왜 "하루(24시간)" 규모인가 — 물리적 정당성

1. **부하추종 운전의 본질적 장주기성**: iSMR은 재생에너지 변동성을 보완하는 유연 전원으로 설계되며, 무붕산(soluble-boron-free) 노심과 제어봉 기반 반응도 제어로 **일일부하추종(daily load-following)** 운전을 목표로 한다. 대표적 부하추종 패턴(예: 100 %–50 %–100 % 일일 주기)은 24시간을 한 주기로 하므로, 이 운전 시나리오의 노심 거동을 직접 해석하려면 24시간 물리 시간 전체를 시간전진해야 한다.
2. **느린 피드백 시간상수**: 노심의 핵심 피드백 메커니즘이 수 시간~하루 규모의 시간상수를 가진다.
   - **제논-135 동특성**: I-135(반감기 ≈ 6.6 h) → Xe-135(반감기 ≈ 9.2 h) 붕괴 사슬에 의한 반응도 피드백이 수 시간 규모로 전개되며, 공간 제논 진동의 주기는 약 15~30시간이다. 출력 변동 후 제논 분포의 재평형은 하루 규모에서만 포착된다.
   - **열적 관성**: 노심·구조물·1차계통 냉각재 열용량에 의한 온도 분포 변화가 분~시간 규모로 진행된다.
   - 이러한 피드백은 단일 정상상태 해로 환원되지 않으며, **시간 이력 전체**가 안전·운전성 평가의 대상이다.
3. **출력 변동에 따른 3차원 유동 재분배**: 출력 램프 시 노심 출력 분포가 변하고, 이에 따라 부분 유로(서브채널)별 유량·온도가 재분배된다. 자연순환·혼합대류 성분의 비중이 시간에 따라 변하므로 1D 계통 코드로 환원할 수 없고, 국소 DNBR·온도 여유의 시간 이력을 평가하려면 **3차원 노심 유동장의 과도 해석**이 필요하다.

### 4-3-2. 의미 있는 demonstration 방안

| 측정 항목 | 보이고자 하는 것 | 예상 결과 |
|---|---|---|
| **누적 wall-clock 곡선** | 물리 시간(0→24 h)에 대한 누적 압력 솔버 시간을 PMG-BiCGSTAB vs BiCG 단독으로 비교 | 두 곡선의 격차가 시간에 따라 선형적으로 벌어짐 → 곱셈적 누적의 직관적 제시 ("BiCG로 N일 걸리는 해석을 PMG로 M시간에") |
| **시간스텝별 반복 횟수 안정성** | 출력 램프·제논 진동으로 유동장·조건수가 변하는 구간에서도 반복 횟수가 격자·시간 무관하게 평탄한지 | PMG는 평탄 유지, BiCG 단독은 급램프 구간에서 반복 횟수 급증·발산 위험 → 강건성 실증 |
| **시간스텝 크기 민감도** | 큰 시간스텝에서 대각 우세 약화로 조건수가 악화될 때 예조건 강건성 | PMG가 더 큰 시간스텝을 허용 → 스텝 수 감소 × 스텝당 가속의 이중 가속 |
| **(선택) 계산 자원 관점** | 장시간 HPC 해석의 코어-시간·전력 소비 | 설계 반복(민감도·인허가) 시 누적되는 실질 비용 절감 정량화 |

- **시나리오 정의(예시)**: 100 %–50 %–100 % 일일부하추종(램프율 예: 정격 대비 분당 수 %), 또는 제논 공간 진동을 유발하는 부분 출력 스텝 변화. 시나리오는 노심 유동 분포 변화를 충분히 자극하도록 설계한다.
- **핵심 메시지**: 이 사례는 고속·강건 압력 솔버가 단발성 벤치마크가 아니라 **실제 SMR 운전 해석 워크로드**에서 갖는 가치를 가장 직접적으로 보여준다. 정상상태 가속(4-1)·강건성(4-2)을 실제 운전 시간 척도에서 통합 검증하는 캡스톤 사례에 해당한다.

- **대상**: iSMR 노심(또는 대표 부분 노심)의 일일부하추종 운전 과도, 중성자/열수력 약결합 또는 강결합.
- **물리**: 출력 램프 + 제논 동특성 + 열적 관성 → 시간 변동 유동장. 장시간 시간전진.
- **측정**: 24시간 누적 가속비, 과도 전 구간 반복 횟수 안정성, 허용 시간스텝 확대.

---

# 5. 결과 (TODO — 실험 수행 후 작성)

- [ ] 격자 무관 수렴성 곡선 (반복 횟수 vs 미지수)
- [ ] BiCG 단독 대비 가속비 표 (MSFR / iSMR)
- [ ] Strong/weak scaling (코어 수 vs 효율)
- [ ] 통신 저감 효과 (통신 시간 비중 vs 레벨)
- [ ] serial 전환 레벨 민감도 분석
- [ ] fixed MG vs variable MG(FBiCGSTAB) 비교
- [ ] (4-3) 누적 wall-clock 곡선 (물리 시간 0→24 h vs 누적 솔버 시간)
- [ ] (4-3) 과도 전 구간 시간스텝별 반복 횟수 안정성 (램프·제논 진동 구간)
- [ ] (4-3) 시간스텝 크기 민감도 (PMG 강건성에 의한 큰 스텝 허용)

---

# 6. 결론 (TODO)

비정렬 격자에 적용 가능한 semi-coarsening 기하적 MG를 BiCGSTAB 예조건자로 통합하고, 통신 저감·serial 전환으로 병렬 확장성을 확보하여 MSFR·iSMR 문제에서 가속을 입증한다(예정).

---

# 7. 참고 문헌 (BibTeX)

> 주의: 출판사 페이지(ScienceDirect/Springer/SSRN) 접근이 차단되어, 일부 항목의 **저자 순서·페이지·article number·DOI**는 게재본과 대조 검증 필요(아래 ⚠️ 표시). arXiv·DOI 링크는 하이퍼링크로 연결.

## 7-1. 저자 그룹 선행 연구 (Scalable Solver 계보)

```bibtex
@article{do2024meshless,
  title   = {Highly scalable meshless multigrid solver for 3D thermal-hydraulic analysis of nuclear reactors},
  author  = {Do, Seong Ju and Ha, Sang Truong and Choi, Hyoung Gwon and Yoon, Han Young},
  journal = {Annals of Nuclear Energy},
  volume  = {207},
  pages   = {110713},
  year    = {2024},
  doi     = {10.1016/j.anucene.2024.110713},
  note    = {\url{https://www.sciencedirect.com/science/article/pii/S0306454924003761}}
}
% ⚠️ 저자 순서·article number(110713) 게재본 대조 필요

@article{do2022highaspect,
  title   = {A meshless geometric multigrid method for a grid with a high aspect ratio},
  author  = {Do, Seong Ju and Ha, Sang Truong and Choi, Hyoung Gwon and Yoon, Han Young},
  journal = {Journal of Mechanical Science and Technology},
  year    = {2022},
  doi     = {10.1007/s12206-022-1019-4},
  note    = {\url{https://link.springer.com/article/10.1007/s12206-022-1019-4}}
}
% ⚠️ 저자·권·페이지 대조 필요

@article{ha2024aiaa,
  title   = {Meshless Geometric Multigrid Method for Complex Geometries with Improved Cell Coarsening Algorithm},
  author  = {Ha, Sang Truong and Do, Seong Ju and Choi, Hyoung Gwon and Yoon, Han Young},
  journal = {AIAA Journal},
  year    = {2024},
  doi     = {10.2514/1.J063127},
  note    = {\url{https://arc.aiaa.org/doi/10.2514/1.J063127}}
}
% ⚠️ 저자 순서·권·연도 대조 필요

@article{ha2021nodecoarsening,
  title   = {A meshless geometric multigrid method based on a node-coarsening algorithm for the linear finite element discretization},
  author  = {Ha, Sang Truong and Choi, Hyoung Gwon},
  year    = {2021},
  note    = {\url{https://www.researchgate.net/publication/351771284}}
}
% ⚠️ 저자·저널·DOI 확인 필요 (ResearchGate 등록본)

@misc{ha2025aggregation,
  title        = {Investigation on an Optimal Aggregation Level for a Parallel Meshless Multigrid Method Based on Domain Decomposition Method},
  author       = {Ha, Sang Truong and Yoon, Han Young and Choi, Hyoung Gwon},
  howpublished = {SSRN preprint 5081886},
  year         = {2025},
  note         = {\url{https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5081886}}
}
% ⚠️ 연도·게재 여부 확인 필요
```

## 7-2. Semi-coarsening 멀티그리드

```bibtex
@article{mulder1989,
  title   = {A new multigrid approach to convection problems},
  author  = {Mulder, W. A.},
  journal = {Journal of Computational Physics},
  volume  = {83},
  number  = {2},
  pages   = {303--323},
  year    = {1989},
  doi     = {10.1016/0021-9991(89)90121-6}
}
% ⚠️ semi-coarsening 기초 문헌. 정확한 출처(JCP 1989 vs SIAM J. Numer. Anal.) 확인 필요

@techreport{osti440722,
  title       = {A parallel semicoarsening multigrid algorithm for solving the Reynolds-averaged Navier-Stokes equations},
  institution = {OSTI},
  number      = {440722},
  note        = {\url{https://www.osti.gov/biblio/440722}}
}
% ⚠️ 저자·연도 확인 필요

@article{dendy2019msg,
  title   = {Multiple Semicoarsened Multigrid for anisotropic diffusion problems},
  author  = {Anonymous},
  journal = {arXiv preprint arXiv:1907.12334},
  year    = {2019},
  note    = {\url{https://arxiv.org/abs/1907.12334}}
}
% ⚠️ 저자·정식 게재처 확인 필요
```

## 7-3. 병렬 MG 통신 저감

```bibtex
@article{reisner2018,
  title   = {Scaling Structured Multigrid to 500K+ Cores through Coarse-Grid Redistribution},
  author  = {Reisner, Andrew and Olson, Luke N. and Moulton, J. David},
  journal = {SIAM Journal on Scientific Computing},
  year    = {2018},
  note    = {\url{https://arxiv.org/abs/1803.02481}}
}
% ⚠️ 권·페이지·DOI 확인 필요

@article{bienz2020,
  title   = {Node-Aware Improvements to Allreduce / Reducing Communication in Algebraic Multigrid},
  author  = {Bienz, Amanda and Gropp, William D. and Olson, Luke N.},
  journal = {International Journal of High Performance Computing Applications},
  volume  = {34},
  number  = {5},
  pages   = {547--561},
  year    = {2020},
  note    = {\url{https://arxiv.org/abs/1904.05838}}
}
% ⚠️ 정확한 제목 확인 필요

@article{manteuffel2016,
  title   = {Nonsymmetric Algebraic Multigrid Sparsification / Reducing Communication Costs in AMG},
  author  = {Manteuffel, Thomas A. and others},
  journal = {SIAM Journal on Scientific Computing},
  year    = {2016},
  note    = {\url{https://arxiv.org/abs/1512.04629}}
}
% ⚠️ 정확한 제목·전체 저자·권·페이지 확인 필요
```

## 7-4. BiCGSTAB / 예조건 수렴 이론

```bibtex
@article{vandervorst1992,
  title   = {Bi-CGSTAB: A Fast and Smoothly Converging Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems},
  author  = {van der Vorst, Henk A.},
  journal = {SIAM Journal on Scientific and Statistical Computing},
  volume  = {13},
  number  = {2},
  pages   = {631--644},
  year    = {1992},
  doi     = {10.1137/0913035},
  note    = {\url{https://epubs.siam.org/doi/10.1137/0913035}}
}

@article{chen2016fbicgstab,
  title   = {A Flexible Variant of Bi-CGSTAB with a Multigrid Preconditioner (FBiCGStab)},
  author  = {Chen, Jie and McInnes, Lois Curfman and Zhang, Hong},
  journal = {Journal of Scientific Computing},
  volume  = {68},
  number  = {2},
  pages   = {803--825},
  year    = {2016},
  note    = {\url{https://jiechenjiechen.github.io/pub/fbcgs.pdf}}
}
% ⚠️ 정확한 제목 확인 필요 (PFLOTRAN 적용)

@article{lucerolorca2026,
  title   = {Two-point spectral clustering and Krylov convergence of multigrid-preconditioned solvers},
  author  = {Lucero Lorca, Jose Pablo and McCoid, Conor and Outrata, Michal},
  journal = {arXiv preprint arXiv:2511.12298},
  year    = {2026},
  note    = {\url{https://arxiv.org/abs/2511.12298}; preprint, 동료심사 전}
}
% ⚠️ 제목·저자 확인 필요 (2026년 preprint)
```

## 7-5. MSFR / 원자로 CFD

```bibtex
@article{msfr_coupled,
  title   = {Coupled neutronics and thermal-hydraulics numerical simulations of a Molten Salt Fast Reactor (MSFR)},
  author  = {Anonymous},
  note    = {\url{https://www.academia.edu/70734179/}}
}
% ⚠️ 저자·저널·연도 확인 필요

@article{fnuen2025,
  title   = {Multigrid preconditioning for reactor thermal-hydraulics CFD},
  journal = {Frontiers in Nuclear Engineering},
  year    = {2025},
  note    = {\url{https://www.frontiersin.org/journals/nuclear-engineering/articles/10.3389/fnuen.2025.1597165/full}}
}
% ⚠️ 제목·저자 확인 필요

@inproceedings{nuclearcfd_springer,
  title     = {Parallel multigrid in nuclear CFD},
  booktitle = {Springer LNCS},
  doi       = {10.1007/978-3-030-39647-3_20},
  note      = {\url{https://link.springer.com/chapter/10.1007/978-3-030-39647-3_20}}
}
% ⚠️ 제목·저자·연도 확인 필요
```

---

# 부록: 작성 메모

- 본 초안은 [[문헌 조사 - PMG 예조건화 BiCGSTAB]]의 검증된 발견을 기반으로 함. 기각된 13개 주장(과장·오귀속)은 인용하지 않음.
- BibTeX의 ⚠️ 항목은 출판사 페이지 접근 차단으로 미확정 — 게재본에서 저자 순서·권·페이지·DOI를 직접 대조 후 확정할 것.
- 결과(5절)·결론(6절)은 실제 수치 실험 후 작성.
