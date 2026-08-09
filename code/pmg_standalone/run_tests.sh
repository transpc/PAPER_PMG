#!/usr/bin/env bash
# PMG 합성 유닛 테스트 러너 (LOOP C007) — 골든 회귀는 C009 에서 이 스크립트를 확장
# 판정 (LOOP.md §2 L2): 드라이버의 독립 residual 판정(PASS/FAIL) + its == ref_its (다르면 WARN)
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
printf "%-10s %-9s %-5s %-12s %s\n" "case" "ncell" "its" "rel_res" "verdict"
for c in "${CASES[@]}"; do
  set -- $c; name=$1; nx=$2; ny=$3; nz=$4; asp=$5; ref=$6
  out=$("$IC" bash -c "ulimit -s unlimited && cd '$HERE/tests' && rm -rf MG_tmp fort.* && ../build/driver_pmg $nx $ny $nz $asp" 2>&1)
  rc=$?
  its=$(awk '{print $1; exit}' "$HERE/tests/fort.501" 2>/dev/null)
  res=$(echo "$out" | grep -o "rel_res= *[^ ]*" | head -1 | awk -F= '{print $2}' | tr -d ' ')
  verdict=FAIL
  echo "$out" | grep -q "VERDICT PASS" && [ "$rc" -eq 0 ] && verdict=PASS
  if [ "$verdict" = PASS ] && [ "${its:-x}" != "$ref" ]; then
    verdict="WARN(its=${its:-?},ref=$ref)"
  fi
  [ "$verdict" = FAIL ] && fail=1
  printf "%-10s %-9s %-5s %-12s %s\n" "$name" "$((nx*ny*nz))" "${its:-?}" "${res:-?}" "$verdict"
done
exit $fail
