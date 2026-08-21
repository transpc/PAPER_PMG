# PMG 셋업 랭크-로컬화 계획 (S-시리즈) — CUPID `subdomain_info` 방식으로 전환

수립: 2026-08-21. **계획 문서 — 실행 전.** 다른 세션이 단독으로 수행할 수 있도록
배경·좌표·검증까지 자족적으로 기술한다.

목표 한 줄: **PMG 셋업의 `(np × N)` 조밀 배열을 없애서, 2,000만 셀 전처리가
192 GB 노드에서 죽지 않게 한다.** 같은 저장소의 `Closed/subdomain_info.f90` 이
이미 쓰고 있는 패턴을 GMG 셋업에 이식하는 작업이다.

---

## 1. 문제 (실측 근거)

### 1-1. 증상
- 프로덕션 2,000만 셀 해석이 **전처리 단계에서 사망** (192 GB 램 머신).
- 같은 머신에서 CUPID 본체 분할(`Closed/subdomain_info.f90`)은 문제 없음.

### 1-2. 실측 메모리 모델 (2026-08-21 채취)

`pmg_standalone` 합성 케이스로 N·np 스윕 측정 후 회귀:

```
rank 0 피크 ≈ N × (589 B + 18.1 B × np)
```

| 측정점 | rank0 피크 |
|---|---|
| 512k, np=4 | 0.38 GB |
| 512k, np=16 | 0.45 GB |
| 512k, np=64 | 0.94 GB |
| 512k, np=128 | 1.49 GB |
| 2.1M, np=16 | 1.73 GB (모델 예측 1.84 GB, 오차 6%) |

N=20M 대입 (rank 0 **혼자**, 다른 랭크·CUPID 전역배열·OS 제외):

| np | 16 | 64 | 128 | **192** | 256 | 512 |
|---|---|---|---|---|---|---|
| rank0 | 17.6 GB | 35 GB | 58 GB | **81 GB** | 104 GB | 197 GB |

**np 를 늘릴수록 rank 0 부담이 커지는 역설** — 메모리가 모자랄 때 통상적으로
쓰는 "랭크 늘리기" 대응이 PMG 에서는 반대로 작동한다.

### 1-3. 원인 — `(np × N)` 조밀 배열

`Source/GMG/6_subdomain_infor_mg.f90` 은 **rank 0 이 np개 전 랭크 몫을 동시에**
만든다. 그래서 랭크 축과 셀 축을 곱한 2D 배열이 필요하다:

| 위치 | 배열 | 크기 |
|---|---|---|
| `:84` | `iperm(np,nelem)` | np·N |
| `:84` | `jperm(np,nelemt)` | np·nelemt |
| `:85` | `rint(np,nelemt)`, `sint(np,nelemt)` | np·nelemt ×2 |
| `:86` | `lcelem(np,nelemt)` | np·nelemt |
| `:391` | `iperm1(np,nelem1)`, `jperm1(np,nelemt)` | 코스 레벨 반복 |
| `:393` | `rint1(np,nelemt)`, `sint1(np,nelemt)` | 코스 레벨 반복 |
| `:666` | `iperm(np,nelem1)`, `jperm(np,nelemt)` **재할당** | 레벨 승계용 사본 |

**`:665–670` 은 별도 함정**: `DEALLOCATE(iperm,jperm)` → `ALLOCATE(iperm(np,nelem1),
jperm(np,nelemt))` → `iperm = iperm1` 순서라, **복사 시점에 `iperm` 과 `iperm1` 이
동시에 살아 있어 np·N 메모리가 순간 2배**가 된다. `MOVE_ALLOC(iperm1, iperm)` 로
바꾸면 복사도 사라지고 피크도 절반이 된다 — **S2 에 함께 넣을 것 (한 줄 수정, 큰 이득)**.

`nelemt` 는 `Predict_nelemt`(`:760`) 가 정한다. **주의: 피크는 fine 이 아니라
코스 레벨에서 난다.**

```
ilv ≤ 2 : np>50 → nelemt = nelem/np*20      (→ np·nelemt = 20N, N 에 비례)
ilv = 2 : np>50 → nelemt = nelem/np*50
ilv ≥ 3 :         nelemt = nelem            (→ np·nelemt = np·N, 폭발)
```

즉 **레벨 3 이상에서 `nelemt = nelem`** 이라 `iperm1·jperm1·rint1·sint1`
네 개가 전부 np×(레벨 크기)가 된다. 레벨 3 ≈ 0.36N 이므로 np=192, N=20M 에서
이 네 개만으로 ≈ 22 GB.

부수적으로 `Domain_infor_FVM_fine.f90:28` 의 `imark(np,np)` 등 O(np²)
자동배열이 있으나 (np=512 에서 1 MB) 병목은 아니다.

---

## 2. 목표 — CUPID 가 이미 쓰는 패턴

### 2-1. `Closed/subdomain_info.f90` 의 구조 (참조 구현)

- **호출**: `02_IO/read_grid.f90:242` 에서 `IF(np.gt.1)` 아래 — **rank0 가드 없음
  = 모든 랭크가 호출**한다.
- **각 랭크는 자기 서브도메인만 계산**한다. 랭크 축 2D 배열이 **0개**.
- 랭크 간 정보가 필요한 부분은 **전역 O(N) 배열에서 유도**한다:
  - `ia_sub(nelem+1)`, `ja_sub(nelem_sub)` — 전역 분할 CSR (rank0 이 만들어 배포)
  - `flag(nelem)`, `flagt(nelem)`, `ja(nelem)`, `jperms(nelem)`,
    `i_neigh_tmp_nbcon0(nelem+1)` — 랭크별 O(N) 자동배열 (5N 정수 = 20 B/셀)
  - 랭크 축은 전부 **1D**: `iutjp(np)`, `cext(np)`, `irecv_cnt(np)`,
    `jsend_cnt(np)`, `icount(np)`, `ia(np+1)`

### 2-2. 두 방식 비교 (N=20M 기준)

| | CUPID `subdomain_info` | PMG `subdomain_infor_mg` (현재) |
|---|---|---|
| 실행 주체 | 모든 랭크 | rank 0 단독 |
| 랭크당 피크 | ≈ 0.4 GB (전역 O(N) 복제) | 다른 랭크는 작음 |
| rank 0 피크 | ≈ 0.56 GB | **35 GB (np=64) ~ 81 GB (np=192)** |
| 셋업 시간 | 병렬 | **직렬** (`pmg_tser`: 2.1M·np=16 에서 99.5 s) |

### 2-3. 전환 후 기대치

- rank 0 피크: **N × 589 B 수준으로 수렴** (20M → ≈ 12 GB), `18.1 B × np` 항 소멸
- np 확장성: 랭크 늘려도 rank0 메모리 불변
- (Phase 2 까지 가면) 직렬 셋업 시간도 분산

---

## 3. 대상 배열의 의미 (전환 설계의 근거)

`Domain_infor_FVM_fine.f90` 에서 채우는 지점을 역추적한 결과:

| 배열 | 의미 | 채움 위치 | 실제 필요 총량 |
|---|---|---|---|
| `lcelem(proc,k)` | proc 의 k 번째 셀 (전역번호) | `:52`, `:341` | Σ(랭크별 셀수) = **N + 고스트** |
| `jperm(ip,loc)` | 랭크 ip 의 지역번호 loc → 전역번호 | `:176`, `:291`, `:308` | 동일 |
| `iperm(ip,ne)` | 전역 셀 ne → 랭크 ip 안의 지역번호 | `:175`, `:290`, `:307` | 셀당 (소유 1 + 고스트 소수) |
| `rint(prc,·)`, `sint(prc,·)` | prc 의 recv/send 인터페이스 목록 | `:311`, `:294` | **인터페이스 크기** |

**핵심**: 네 배열 모두 "랭크별 가변 길이 리스트"이며 총량은 `np·N` 이 아니라
`N·(1+고스트비)` 다. 조밀 2D 저장이 낭비의 전부다.

---

## 4. 실행 계획 (S1 → S4)

각 단계는 **독립 커밋 + 게이트 통과**를 만족해야 한다. 순서 준수.

### S1 — 계측·기준선 고정 (코드 변경 없음)

1. `§7 검증 프로토콜`의 메모리 측정 스크립트를 재작성(§6-2)하고,
   현재 코드에서 기준선 재채취: (512k, np=4/16/64/128), (2.1M, np=16).
2. 값이 §1-2 표와 ±10% 안에 들어오는지 확인. 벗어나면 머신 상태(다른 작업
   점유)를 의심 — `ps aux --sort=-%cpu | head` 로 확인 후 한산할 때 재측정.
3. 이 기준선이 S2~S4 의 성공 판정 근거가 된다.

### S2 — `Predict_nelemt` 의 코스 레벨 폭발 제거 (선행·저위험)

**가장 싼 큰 이득.** `6_subdomain_infor_mg.f90:760` 의 `Predict_nelemt` 에서
`ilv ≥ 3 → nelemt = nelem` 을 유한 상계로 교체한다.

- 현 상태는 "코스 레벨은 고스트 비율이 커서 안전마진을 무한대로 둔" 것.
- 실제 필요량은 `Domain_infor_FVM_coarse/coarsest` 가 만드는 최대 랭크 크기.
- **권장**: 상계를 `MAX(2000, nelem/np*K)` 형태로 두되 K 를 실측으로 정한다.
  코스 레벨 고스트 비율이 fine 보다 훨씬 크므로 K 는 크게(예: 100) 잡아도
  `np·nelemt = 100N/np·np = 100N` 으로 **N 에 비례**하게 된다.
- **필수 안전장치**: 초과 검출이 **코스 레벨에만** 있다 (`:465` `'nelemt is small'` → STOP). K 를 줄일 때 이 검사가 살아있는지 반드시 확인. 조용한 오버런 금지.
- **동시 수행 권장 (한 줄씩, 저위험)**:
  - `:665–670` 의 복사를 `MOVE_ALLOC` 로 교체 → np·N 순간 2배 해소
  - fine 루프(`:152–260`)에 `nelemt` 초과 검사 추가 (위험표 #2)
- 검증: 게이트 12/12 + np=4/16/64 에서 메모리 재측정 → `18.1 B × np` 계수가
  유의미하게 떨어져야 한다.

> S2 만으로도 20M·np=192 의 rank0 이 81 GB → 대략 30 GB 대로 내려갈 것으로
> 기대된다(레벨3 이상 기여분 소멸). **단독으로 문제를 해결할 수도 있으니
> S2 직후 반드시 재측정하고, 충분하면 S3/S4 는 논문 일정에 따라 보류 가능.**

### S3 — `(np × N)` → CSR 압축 (구조 변경, 알고리즘 불변)

S2 이후에도 남는 `iperm(np,nelem)` 계열을 압축 저장으로 바꾼다.
**알고리즘·통신·수치는 그대로**, 저장 표현만 바꾸는 것이 원칙.

대상과 치환:

| 현재 | 치환 |
|---|---|
| `lcelem(np,nelemt)` | `lc_ptr(np+1)`, `lc_val(총량)` (CSR) |
| `jperm(np,nelemt)` | `jp_ptr(np+1)`, `jp_val(총량)` |
| `rint/sint(np,nelemt)` | `r_ptr/s_ptr(np+1)`, `r_val/s_val(총량)` |
| `iperm(np,nelem)` | 아래 주의 참조 |

**`iperm` 주의**: 이것만 "역방향 조밀 조회"라 CSR 화가 덜 자명하다. 두 선택지:
- (a) 랭크 루프가 `prc` 바깥에 있는 구간에서는 `iperm_1d(nelem)` 하나를 두고
  **현재 prc 것만 채웠다가 지우는** 방식 (rank0 메모리 O(N)).
  → 루프 순서가 prc-major 인지 먼저 확인할 것 (`:152`, `:435` 의 `DO prc=1,np`).
- (b) prc-major 가 아니면 (a) 가 불가하므로, 셀당 (랭크,지역번호) 쌍을
  CSR 로 저장 (`ip_ptr(nelem+1)`, `ip_rank(·)`, `ip_loc(·)`).

**수정 파일**: `6_subdomain_infor_mg.f90`, `Domain_infor_FVM_fine.f90`,
`Domain_infor_FVM_coarse.f90`, `Domain_infor_FVM_coarsest.f90` (시그니처 동반 변경).

**검증**: 게이트 12/12 **bitwise 동일** — 저장 표현만 바꿨으므로 수치가 1비트도
달라지면 안 된다. 이것이 S3 의 안전망이다.

### S4 — 랭크-로컬 계산으로 전환 (CUPID 패턴 본체)

S3 까지 하면 메모리는 해결되지만 **직렬 셋업 시간**은 남는다. 이를 없애려면
CUPID 처럼 각 랭크가 자기 몫만 계산해야 한다.

전환 요지:
1. rank 0 이 만든 **전역 데이터**(레벨별 `celem`(셀→랭크), 연결성, 좌표)를
   전 랭크에 배포한다. **D6 로 이미 MPI 팬아웃 경로가 있으므로** 그 위에 얹는다
   (`MD_MG_index` 의 `stg_*` 스테이징, `2_read_mesh_MPI.f90` 의 SCATTERV/BCAST).
2. `Domain_infor_FVM_*` 의 `do ip=1,np` / `do prc=1,np` 루프를
   **`myrank` 한 값만 처리**하도록 축약한다
   (`Domain_infor_FVM_fine.f90:87,129,168,189,238,259,264,275,282`).
3. 랭크 간 인접 정보(`nbdom`, `nnbdom`)는 CUPID 처럼 **전역 `celem` + 연결성에서
   각 랭크가 유도**한다.
4. 팬아웃 자체가 불필요해진다 — 각 랭크가 자기 데이터를 이미 갖고 있으므로
   `2_read_mesh_MPI.f90` 의 unpack 경로도 대폭 단순화 가능(선택).

**난이도**: S1~S3 보다 크게 높다. 통신 맵 정합성(대칭성)이 깨지면 조용히
오답이 나오므로, 아래 자기검증을 반드시 넣을 것:
- `A 가 B 의 이웃 ⇔ B 가 A 의 이웃` (대칭성)
- `send_count(A→B) == recv_count(B←A)`
- `Σ(랭크별 소유 셀) == 전역 셀 수`, 중복 소유 0

**검증**: 게이트 12/12 bitwise + np=1/4/16 프로덕션 스팟 its 일치 +
`pmg_tser` 가 np 에 반비례해 줄어드는지 확인.

---

## 5. 위험과 함정

| # | 함정 | 대응 |
|---|---|---|
| 1 | **`nelemt` 하한** `IF(nelemt.LE.2000) nelemt=2000` (`:796`) | S2 에서 상계 조정 시 이 하한과 충돌하지 않는지 확인 |
| 2 | **조용한 오버런 — 확인됨**: `nelemt` 초과 검사가 **코스 레벨(`:465`)에만 있고 fine 레벨(`:152–260`)에는 없다** | `:465` 검사를 절대 제거하지 말 것 + **fine 루프에도 동일 검사를 추가**할 것 (S2 착수 시 첫 작업 권장). 검사 없이 `nelemt` 상계를 줄이면 조용한 메모리 오염이 된다 |
| 3 | **`ioplv` 재진입** — `GOTO 500`(`:76`)으로 셋업이 **2회 실행**될 수 있음 | 새로 만드는 배열도 재할당 가드 필요 (`stg_mg_init` 가 이미 이 패턴, `MD_MG_index`) |
| 4 | **INTEGER(4) 오버플로** — np·N > 2³¹ (np≥108, N=20M) | S3 로 np·N 배열이 사라지면 자동 해소. 그 전까지는 큰 np 금지 |
| 5 | **극소 도메인 병리 (미해결 기존 결함)** | 랭크당 ~15셀에서 `MPI_Isend` count 음수 → SIGSEGV (24³ np=900 재현기). S4 가 이 코드를 건드리므로 **동시에 고쳐질 수도, 가려질 수도** 있음. 별도 추적 |
| 6 | **게이트 판정은 clean 재빌드 기준** | incremental 상태에서 산발 실패 관측 이력 있음 (`CLEANUP_PMG.md §8`) |
| 7 | **머신 부하** | 다른 작업(예: `IBM_unstructured/build/ns`)이 CPU 를 점유하면 게이트가 흔들린다. 실패 시 `ps aux --sort=-%cpu` 확인 후 재실행 |

---

## 6. 검증 프로토콜

### 6-1. 기능 (매 단계 필수)

```bash
cd code/pmg_standalone
../scripts/in_contain.sh make clean && ../scripts/in_contain.sh make driver
./run_tests.sh          # 12/12 PASS + bitwise baseline 일치
/home/sjdo/PAPER_PMG/code/scripts/build.sh   # CUPID 링크
```

- S3 까지는 **bitwise 동일**이 요구 조건 (저장 표현만 변경).
- S4 도 원칙적으로 bitwise 동일해야 한다 (같은 분할 → 같은 연산자).
  달라지면 분할 결과가 바뀐 것이므로 **버그로 간주하고 추적**할 것.
- 프로덕션 스팟: `CUPID_NP=4 code/scripts/run.sh` 후 `fort.501` its 대조.

### 6-2. 메모리 측정 스크립트

```bash
#!/usr/bin/env bash
# 전처리 피크: rank0 최대 RSS + 전체 합.  사용: mem.sh <nx> <np>
HERE=/home/sjdo/PAPER_PMG/code/pmg_standalone
IC=/home/sjdo/PAPER_PMG/code/scripts/in_contain.sh
n=$1; np=$2
cd "$HERE"; rm -f tests/fort.501
"$IC" bash -c "ulimit -s unlimited && cd '$HERE/tests' && \
  I_MPI_FABRICS=shm mpirun -np $np ../build/driver_pmg $n $n $n 1.0 1" >/dev/null 2>&1 &
runpid=$!; maxsum=0; maxone=0
while kill -0 $runpid 2>/dev/null; do
  read sum one < <(ps -o rss= -C driver_pmg 2>/dev/null | awk '{s+=$1; if($1>m)m=$1} END{print s+0, m+0}')
  (( sum > maxsum )) && maxsum=$sum; (( one > maxone )) && maxone=$one
  sleep 0.2
done
wait $runpid
printf "N=%-9s np=%-4s rank0peak=%7.2fGB total=%7.2fGB its=%s\n" "$((n*n*n))" "$np" \
  "$(echo "$maxone/1048576"|bc -l)" "$(echo "$maxsum/1048576"|bc -l)" \
  "$(head -1 tests/fort.501 2>/dev/null | awk '{print $1}')"
```

**주의**: `driver_pmg` 합성 드라이버는 **모든 랭크가 전역 메시를 복제**한다
(`driver/driver_pmg.f90` 의 `ALLOCATE(celem(ncell))` 등이 rank0 가드 없음).
따라서 `total` 은 CUPID 실제보다 과대. **`rank0peak` 만 지표로 쓸 것.**

### 6-3. 성공 판정

| 단계 | 판정 기준 |
|---|---|
| S2 | `18.1 B × np` 계수가 크게 감소 (np=128 에서 rank0 피크 30% 이상 감소) |
| S3 | np 계수가 **거의 0** — np=4 와 np=128 의 rank0 피크가 같은 자릿수 |
| S4 | 위 + `pmg_tser` 가 np 에 반비례 |
| 최종 | N=20M 환산 rank0 ≤ 15 GB, 프로덕션 20M 전처리 완주 |

---

## 7. 참고 좌표 요약

| 대상 | 파일:라인 |
|---|---|
| PMG 셋업 본체 (rank0 단독) | `Source/GMG/6_subdomain_infor_mg.f90` |
| ├ `(np×N)` 할당 (fine) | `:82–86` |
| ├ fine 팬아웃 루프 | `:152–260` |
| ├ 코스 레벨 루프 | `:303–682` |
| ├ `(np×N)` 할당 (coarse) | `:391–393` |
| ├ 코스 팬아웃 루프 | `:435–637` (버퍼 초과 검사 `:465`) |
| ├ 레벨 승계 재할당 (피크 2배 지점) | `:665–670` |
| └ `Predict_nelemt` | `:760` |
| 분할 커널 (np 루프 내장) | `Domain_infor_FVM_{fine,coarse,coarsest}.f90` |
| 팬아웃 수신·언팩 | `Source/GMG/2_read_mesh_MPI.f90` |
| 스테이징 버퍼 정의 | `Source/GMG/module/MD_MG_index.f90` (`stg_*`) |
| **참조 구현 (CUPID)** | `Source/Closed/subdomain_info.f90` |
| └ 호출부 (전 랭크) | `Source/02_IO/read_grid.f90:242` |
| PMG 셋업 호출부 (rank0) | `Source/02_IO/read_grid.f90:248–253` |

관련 문서: `code/CLEANUP_PMG.md`(클리닝 이력·검증 관행), `code/LOG.md`,
`Finding - Chebyshev 스펙트럼 상계와 np-무관 수렴성.md`(수렴성 트랙 — 본 계획과 독립).

---

## 8. 착수 시 첫 3개 명령

```bash
# 1) 현재 기준선 재확인 (머신 한산할 때)
bash mem.sh 80 4 ; bash mem.sh 80 128       # §6-2 스크립트 저장 후

# 2) Predict_nelemt 현재 로직 확인
sed -n '760,790p' code/Source/GMG/6_subdomain_infor_mg.f90

# 3) 게이트 green 확인 (작업 시작 전 기준선)
cd code/pmg_standalone && ./run_tests.sh
```
