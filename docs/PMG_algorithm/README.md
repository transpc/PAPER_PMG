# PMG 예조건자 알고리즘 문서

`code/Source/GMG/` 에 구현된 **PMG(Parallel MultiGrid) 예조건화 BiCGSTAB 압력 솔버**의 알고리즘·병렬화 상세 문서. (작성: 2026-08-19, 클리닝 D5 시점 코드 기준)

## 목차

| 문서 | 내용 | 그림 |
|---|---|---|
| [01_overview.md](01_overview.md) | 전체 구조: SOLVE_GMG 진입, BiCGSTAB 외부 반복, 고정 구성·파라미터 | fig01, fig04 |
| [02_setup_hierarchy.md](02_setup_hierarchy.md) | 계층 구성: MIS 조대화, 보간 P / R=Pᵀ, Galerkin RAP, 연접 저장 구조, GC 구성 | fig05, fig06 |
| [03_vcycle_smoothers.md](03_vcycle_smoothers.md) | V-cycle 단계별 해부, 4차 Chebyshev / GS 스무더, 최조 레벨 해법 | fig02, fig03 |
| [04_parallelization.md](04_parallelization.md) | **병렬화 상세**: 분할·노드 분류·번호 부여, halo 교환, A/R/P 패턴, GC gather, 통신 지도 | fig07–fig11 |
| [05_demo2d.md](05_demo2d.md) | 2D 미니어처 구현(`demo2d/`)으로 각 병렬 과정 가시화 | fig07–fig11 |
| [06_rap_parallel_example.md](06_rap_parallel_example.md) | **병렬 RAP 완전 해부**: 1D 수치 예제로 수식↔코드↔통신 대조, FAQ | fig12, fig13 |

## 30초 요약

1. **알고리즘**: BiCGSTAB 반복당 V-cycle 1회짜리 PMG 예조건자를 2번 적용. 계층은 좌표 거리 기반 MIS C/F 분할(행렬값 미사용) + 역거리 보간 P, **R = Pᵀ 정확한 Galerkin**. fine 스무딩은 4차 Chebyshev(내적 없음 = ALLREDUCE 프리), 코어스는 D⁻¹ 보정 GS.
2. **병렬화**: owner-computes + 덮어쓰기 halo. 셋업(조대화·패턴)은 rank 0 직렬 1회, 수치 RAP 만 병렬. 연산자별(A/R/P) **최소 고스트 집합**을 레벨별로 분리 교환. 전역 노드 ≤ 100 이 되면 **ALLGATHERV 로 전 랭크에 중복 수집 후 통신 없는 직렬 직접해**(gather-to-all)로 코어스 레벨 레이턴시를 제거.

## 그림 재생성

모든 그림은 파이썬으로 생성된다:

```bash
cd docs/PMG_algorithm
python3 scripts/fig01_solver_stack.py     # 개념도 (fig01~fig06)
python3 scripts/fig02_vcycle.py
python3 scripts/fig03_chebyshev.py        # 코드 재귀식 수치 재현
python3 scripts/fig04_bicgstab_flow.py
python3 scripts/fig05_coarsening.py       # coarsening_semi 알고리즘 재현
python3 scripts/fig06_memory_layout.py
python3 scripts/fig12_rap_example.py      # RAP 1D 예제 (수치는 numpy 로 계산·검증)
python3 scripts/fig13_rap_comm.py
python3 demo2d/make_figs.py               # 2D 미니어처 기반 (fig07~fig11)
```

의존성: numpy, scipy, matplotlib. 그림 라벨은 영문(환경에 한글 폰트 부재).

## 폴더 구조

```
docs/PMG_algorithm/
├── README.md            ← 이 파일
├── 01_overview.md … 05_demo2d.md
├── figures/             ← 생성된 PNG (fig01~fig11)
├── scripts/             ← 개념도 생성 스크립트 (style_common.py 공용 스타일)
└── demo2d/              ← 2D 미니어처 PMG 구현 + 데모 그림 생성
    ├── pmg2d.py         (실행: python3 pmg2d.py → 수렴 확인)
    └── make_figs.py
```
