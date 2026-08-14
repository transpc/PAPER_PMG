#!/usr/bin/env bash
# PMG 테스트 러너 — 합성 유닛 테스트 (LOOP C007) + 골든 재생 회귀 (LOOP C009)
# 판정 (LOOP.md §2 L2):
#   합성: 드라이버의 독립 residual 판정(PASS/FAIL) + its == ref (다르면 WARN)
#   골든: 드라이버 게이트(수렴+충실도) + its == 하네스 ref + 베이스라인 bitwise(md5)
#         ※ its/bitwise 기준은 "하네스 베이스라인"이다 — 프로덕션 its 와는 다를 수
#           있음 (예조건자 시드 상태 차이, LOG C009). 프로덕션 대비는 충실도 게이트가 담당.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IC="$HERE/../scripts/in_contain.sh"

"$IC" make -C "$HERE" driver >/dev/null 2>&1 || { echo "BUILD FAIL"; exit 1; }

# name nx ny nz aspect ref_its   (ref_its 출처: LOOP.md §2 기준표, C007 채취)
CASES=(
  "iso24    24 24 24   1.0  4"
  "aniso24  24 24 24 100.0  4"
  "iso48    48 48 48   1.0  5"
)

fail=0
printf "%-10s %-9s %-5s %-12s %s\n" "case" "ncell" "its" "res_gate" "verdict"
for c in "${CASES[@]}"; do
  set -- $c; name=$1; nx=$2; ny=$3; nz=$4; asp=$5; ref=$6
  out=$("$IC" bash -c "ulimit -s unlimited && cd '$HERE/tests' && rm -rf MG_tmp fort.* && ../build/driver_pmg $nx $ny $nz $asp" 2>&1)
  rc=$?
  its=$(awk '{print $1; exit}' "$HERE/tests/fort.501" 2>/dev/null)
  res=$(echo "$out" | grep -o "res_gate= *[^ ]*" | head -1 | awk -F= '{print $2}' | tr -d ' ')
  verdict=FAIL
  echo "$out" | grep -q "VERDICT PASS" && [ "$rc" -eq 0 ] && verdict=PASS
  if [ "$verdict" = PASS ] && [ "${its:-x}" != "$ref" ]; then
    verdict="WARN(its=${its:-?},ref=$ref)"
  fi
  [ "$verdict" = FAIL ] && fail=1
  printf "%-10s %-9s %-5s %-12s %s\n" "$name" "$((nx*ny*nz))" "${its:-?}" "${res:-?}" "$verdict"
done

# ---- MPI 합성 (C010-2) — np 스케일 회귀. np≥6 수렴 붕괴는 미해결 (LOG C010-2) ----
# name np ref_its   (ref: C010-2 채취, iso24 슬랩 분할)
MCASES=(
  "iso24_np2 2 4"
  "iso24_np4 4 4"
)
for c in "${MCASES[@]}"; do
  set -- $c; name=$1; n=$2; ref=$3
  out=$("$IC" bash -c "ulimit -s unlimited && cd '$HERE/tests' && rm -rf MG_tmp fort.* && I_MPI_FABRICS=shm mpirun -np $n ../build/driver_pmg 24 24 24 1.0" 2>&1)
  rc=$?
  its=$(awk '{print $1; exit}' "$HERE/tests/fort.501" 2>/dev/null)
  res=$(echo "$out" | grep -o "res_gate= *[^ ]*" | head -1 | awk -F= '{print $2}' | tr -d ' ')
  verdict=FAIL
  echo "$out" | grep -q "VERDICT PASS" && [ "$rc" -eq 0 ] && verdict=PASS
  if [ "$verdict" = PASS ] && [ "${its:-x}" != "$ref" ]; then
    verdict="WARN(its=${its:-?},ref=$ref)"
  fi
  [ "$verdict" = FAIL ] && fail=1
  printf "%-10s %-9s %-5s %-12s %s\n" "$name" "13824" "${its:-?}" "${res:-?}" "$verdict"
done

# ---- 골든 재생 회귀 (C009, ECT1 추가 C017) — 바이너리가 있을 때만 (fresh clone 은 SKIP) ----
# 세트별 fid_gate: 프로덕션·하네스 its 차이(예조건자 시드 문맥, LOG C009)의 발현 크기가
# 케이스마다 달라서 — 구입력 관측 ≤1e-10 → 1e-9 유지, ECT1 관측 ≤2.4e-7 → 1e-6 (LOG C017)
run_golden() {  # $1=디렉토리명 $2=접두어 $3=fid_gate, 이후 인자 = "step r1 r2" 목록
  local dir="$1" pfx="$2" fgate="$3"; shift 3
  local GOLD="$HERE/golden/$dir"
  if [ ! -f "$GOLD/setup_r0.bin" ]; then
    printf "%-10s %-9s %-5s %-12s %s\n" "$pfx" "-" "-" "-" "SKIP(채취 필요 — golden/$dir/meta.md)"
    return
  fi
  local c s r1 r2 name out rc its verdict k fid
  for c in "$@"; do
    set -- $c; s=$1; r1=$2; r2=$3; name="${pfx}_s${s}"
    if [ ! -f "$GOLD/s${s}_k1_r0_c1.pre" ]; then
      printf "%-10s %-9s %-5s %-12s %s\n" "$name" "-" "-" "-" "SKIP(덤프 없음)"; continue
    fi
    out=$("$IC" bash -c "ulimit -s unlimited && cd '$HERE/tests' && rm -rf MG_tmp fort.501 && ../build/driver_pmg replay ../golden/$dir $s $fgate" 2>&1)
    rc=$?
    its=$(awk '{printf "%s ", $1}' "$HERE/tests/fort.501" 2>/dev/null | tr -d ' \n')
    verdict=FAIL
    echo "$out" | grep -q "VERDICT PASS" && [ "$rc" -eq 0 ] && verdict=PASS
    if [ "$verdict" = PASS ] && [ "$its" != "${r1}${r2}" ]; then
      verdict="WARN(its=$its,ref=${r1}${r2})"
    fi
    # 베이스라인 bitwise (있을 때만)
    if [ "$verdict" = PASS ] && [ -f "$GOLD/baseline_s${s}_k1.u" ]; then
      for k in 1 2; do
        cmp -s "$HERE/tests/replay_s${s}_k${k}.u" "$GOLD/baseline_s${s}_k${k}.u" \
          || verdict="FAIL(bitwise k$k)"
      done
    fi
    case "$verdict" in FAIL*) fail=1;; esac
    fid=$(echo "$out" | grep -o "fidelity= *[^ ]*" | tail -1 | awk -F= '{print $2}' | tr -d ' ')
    printf "%-10s %-9s %-5s %-12s %s\n" "$name" "436136" "${its:-?}" "fid=${fid:-?}" "$verdict"
  done
}

# 구입력(d145c81) 골든 — 하네스 ref its: meta.md (C009/C014)
run_golden iSMR436k_np1 gold 1e-9  "1 2 2"  "10 2 2"  "30 3 3"
# ECT1 신입력(29ee238) 골든 — 하네스 ref its: meta.md (C017)
# s150 = dt 포화(0.05s)·정착 상태의 어려운 행렬 (its 9~10) — 초기 과도(s1~30)보다 감지력 높음
run_golden ECT1_436k_np1 ect1 1e-6  "1 1 1"  "10 1 1"  "30 3 3"  "150 9 9"

exit $fail
