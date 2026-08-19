# 05. 2차원 미니어처 데모 — 병렬 PMG 의 각 과정 가시화

> 코드: [demo2d/pmg2d.py](demo2d/pmg2d.py) (알고리즘 본체), [demo2d/make_figs.py](demo2d/make_figs.py) (그림 생성)

## 5.1 목적과 범위

실제 코드(3D 비정렬, Fortran+MPI, ~14k 라인)의 **병렬 구조를 1:1 로 축소 재현**한 파이썬 구현이다. 4개의 "랭크"를 단일 프로세스에서 시뮬레이션하되, 모든 랭크 간 데이터가 **명시적 send/recv 버퍼를 통해서만** 이동하도록 만들어 각 통신 단계를 기록·가시화할 수 있다.

| 실제 코드 | 데모 대응물 | 비고 |
|---|---|---|
| owned(1..nintf) / ghost(nintf+1..nnode) 번호 부여 | `RankLevel.owned / .ghost / .g2l` | 동일한 owned-first 규칙 |
| `nbdom, spt/rpt, sintf/rintf` | `CommPattern.nbdom / .sintf / .rintf` | 수신자 순서로 송신자가 패킹하는 규칙 포함 |
| `send_receive` (pack→ISEND/IRECV→unpack) | `halo_exchange()` + `Traffic` 로거 | 덮어쓰기(비누산) 모델 동일 |
| A/R/P 패턴 (`Neighbor_node_ARP`) | `commA/commR/commP` — A_l·R·P 의 희소패턴에서 유도 | 정의가 코드와 동일 (연산자별 열 의존) |
| Chebyshev 4th-kind (`poly_cheb_smooth`, a=0.3) | `cheb4_smooth()` | 재귀식 그대로 |
| `Smooth_GS2` (D⁻¹ 보정, 고스트 동결) | `gs_smooth()` | 랭크 간 block-Jacobi 특성 재현 |
| `SOLVE_GC_all` (ALLGATHERV + 중복 직접해) | 소유분 수집 → `np.linalg.solve` 를 랭크마다 보유 | gather-to-all 모델 동일 |
| `SOLVER_NEW` V-cycle | `PMG2D.vcycle()` | 단계 순서·통신 위치 동일 |
| `solve_pbcg_mg` BiCGSTAB | `PMG2D.bicgstab()` | 예조건 2회/반복 |

**단순화한 부분** (실제 코드와 다름):
- 격자: 2D 구조 격자 + 기하 2:1 조대화. 실제는 비정렬 + 거리 기반 MIS ([02](02_setup_hierarchy.md) §2.2 그림이 MIS 를 별도 재현).
- R: 행 정규화된 Pᵀ. 실제는 정확히 R = Pᵀ (스케일 없음).
- 랭크는 시뮬레이션이므로 실제 MPI 지연/대역폭 효과는 없음 — 메시지 수·크기만 계수한다.

## 5.2 실행

```bash
cd docs/PMG_algorithm/demo2d
python3 pmg2d.py        # 수렴 확인 (PMG 16회 vs plain 41회)
python3 make_figs.py    # fig07~fig11 재생성
```

## 5.3 데모가 보여주는 것

### (1) 분할과 로컬 번호 부여 — fig07
[04](04_parallelization.md) §4.1 의 interior→interface→ghost(이웃별 블록) 순서를 실제 좌표 위에서 확인. ghost 블록이 연속이라 recv 버퍼가 복사 없이 꽂히는 구조가 눈에 보인다.

### (2) halo 교환 3단계 — fig08
`send_receive` 의 pack(sintf 순서) → 이웃쌍 메시지 → unpack(rintf 덮어쓰기)을 값 추적으로 확인. fine 레벨 5점 스텐실에서는 **대각 이웃 랭크와 통신이 없다**는 점도 드러난다 (코어스 레벨에선 RAP 로 스텐실이 넓어져 대각 통신이 생긴다 — fig11 표에서 레벨 2 의 A 패턴 메시지 수가 상대적으로 많은 이유).

### (3) A/R/P 최소 집합 — fig09
같은 랭크·같은 레벨이라도 연산자별 고스트 집합이 다르다. 합집합 하나로 통신하는 대신 연산별 최소 집합만 교환하는 실제 코드의 설계 의도가 그대로 재현된다.

### (4) GC gather — fig10
소유분 연접(ALLGATHERV) → 순열(`imapgatR`) → 전 랭크 중복 직접해 → 로컬 복원. "코어스 레벨 halo 통신을 집단 통신 1회로 대체"하는 구조.

### (5) 수렴과 통신량 — fig11
- 16×16 Poisson, 4랭크, 3레벨에서 PMG 예조건이 반복 수를 41→16 으로 감소.
- 사이클당 메시지 집계: fine A 패턴이 메시지 수를 지배, 레벨이 내려갈수록 메시지 크기 급감 → 레이턴시 지배 → GC 전환 논리.

## 5.4 코드 구조 요약

```
pmg2d.py
├─ poisson2d / build_hierarchy      # 전역 문제 + P(역거리, ip_nmax=4 유사), R, Galerkin RAP
├─ decompose                        # 소유권(블록) → RankLevel: owned/ghost/g2l
│   └─ build_pattern + 송신측 채움   # 연산자 희소패턴 → commA/commR/commP
├─ halo_exchange (+Traffic)         # pack → "ISEND/IRECV" → unpack, 메시지 로깅
├─ cheb4_smooth / gs_smooth         # 스무더 2종
├─ PMG2D.vcycle                     # SOLVER_NEW 순서 그대로 (통신 위치 포함)
└─ PMG2D.bicgstab                   # 외부 반복, M⁻¹ = V-cycle 1회
```
