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

# ---- 골든 재생 회귀 (C009) — 바이너리가 있을 때만 (fresh clone 은 SKIP) ----
GOLD="$HERE/golden/iSMR436k_np1"
# name step ref_its_k1 ref_its_k2   (하네스 기준, meta.md 참조)
GCASES=(
  "gold_s1   1  2 2"
  "gold_s10 10  2 2"
  "gold_s30 30  3 3"
)
if [ -f "$GOLD/setup_r0.bin" ]; then
  for c in "${GCASES[@]}"; do
    set -- $c; name=$1; s=$2; r1=$3; r2=$4
    if [ ! -f "$GOLD/s${s}_k1_r0_c1.pre" ]; then
      printf "%-10s %-9s %-5s %-12s %s\n" "$name" "-" "-" "-" "SKIP(덤프 없음)"; continue
    fi
    out=$("$IC" bash -c "ulimit -s unlimited && cd '$HERE/tests' && rm -rf MG_tmp fort.501 && ../build/driver_pmg replay ../golden/iSMR436k_np1 $s" 2>&1)
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
else
  printf "%-10s %-9s %-5s %-12s %s\n" "golden" "-" "-" "-" "SKIP(채취 필요 — golden/iSMR436k_np1/meta.md)"
fi
exit $fail
