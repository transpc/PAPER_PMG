#!/usr/bin/env bash
# np 스케일 전용 측정 — PMG vs hypre (2026-08-21)
#
# 고정 문제(기본 48³)를 np 를 바꿔가며 풀어 강한 확장성(strong scaling)을 본다.
#   its   : 분할에 따른 수렴 저하 여부 (논문의 np-무관 주장 대조군)
#   time  : 실제 병렬 효율
#
# 주의 — 코어 가용량: 이 머신은 20코어지만 다른 작업이 상주할 수 있다.
#   np + (타 작업 랭크 수) > nproc 이면 오버서브스크립션으로 시간이 급증한다.
#   그 구간의 수치는 솔버 특성이 아니므로 NP_MAX 로 상한을 두고 측정할 것.
#   실행 전 `uptime` 부하와 `nproc` 을 확인하고 상한을 정한다.
#
# 통계: REPEAT 회 중앙값 (0 이하 표본은 클럭 이상치로 배제).
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IC="$HERE/../scripts/in_contain.sh"
REPEAT="${REPEAT:-3}"
NPS="${NPS:-1 2 4 8}"
GRID="${GRID:-48 48 48 1.0}"
IPART="${IPART:-0}"          # 0=k-슬랩, 1=3D 블록
OUT="${OUT:-$HERE/tests/np_scan.tsv}"
RAW="$HERE/tests/.np_raw"

"$IC" make -C "$HERE" driver_hypre >/dev/null 2>&1 || { echo "BUILD FAIL"; exit 1; }
: > "$OUT"
printf "np\tipart\tpmg_its\tpmg_t\thyp_its\thyp_tset\thyp_tsol\tnsample\n" >> "$OUT"

echo "== np 스캔: grid=($GRID) ipart=$IPART repeat=$REPEAT =="
printf "%-4s %-8s %-9s %-8s %-9s %-9s %s\n" "np" "pmg_its" "pmg_t" "hyp_its" "hyp_tset" "hyp_tsol" "n"
for np in $NPS; do
  : > "$RAW"
  for r in $(seq 1 "$REPEAT"); do
    cmd="../build/driver_hypre $GRID $IPART"
    [ "$np" -gt 1 ] && cmd="I_MPI_FABRICS=shm mpirun -np $np $cmd"
    out=$("$IC" bash -c "ulimit -s unlimited && cd '$HERE/tests' && rm -rf MG_tmp fort.* && $cmd" 2>&1)
    echo "$out" | grep -q "VERDICT PASS" || { echo "  !! np=$np repeat=$r 실패"; continue; }
    echo "$out" | grep "^ RESULT" | while IFS= read -r l; do
      printf "%s\t%s\t%s\t%s\t%s\n" \
        "$(sed -n 's/.*pmg_its= *\([^ ]*\).*/\1/p' <<<"$l")" \
        "$(sed -n 's/.* pmg_t= *\([^ ]*\).*/\1/p'  <<<"$l")" \
        "$(sed -n 's/.*hyp_its= *\([^ ]*\).*/\1/p' <<<"$l")" \
        "$(sed -n 's/.*hyp_tset= *\([^ ]*\).*/\1/p' <<<"$l")" \
        "$(sed -n 's/.*hyp_tsol= *\([^ ]*\).*/\1/p' <<<"$l")" >> "$RAW"
    done
  done
  [ -s "$RAW" ] || { printf "%-4s %s\n" "$np" "표본 없음"; continue; }
  awk -F'\t' -v np="$np" -v ip="$IPART" -v out="$OUT" '
    function med(a, c,   i, j, t, m) {
      if (c == 0) return -1
      for (i = 1; i <= c; i++) t[i] = a[i]
      for (i = 2; i <= c; i++) { m = t[i]; j = i-1
        while (j > 0 && t[j] > m) { t[j+1] = t[j]; j-- }; t[j+1] = m }
      return (c % 2) ? t[(c+1)/2] : (t[c/2] + t[c/2+1])/2.0 }
    { pit = $1; hit = $3
      if ($2+0 > 0) pt[++cp] = $2+0
      if ($4+0 > 0) hs[++ch] = $4+0
      if ($5+0 > 0) ho[++co] = $5+0 }
    END { printf "%s\t%s\t%s\t%.4f\t%s\t%.4f\t%.4f\t%d\n",
                 np, ip, pit, med(pt,cp), hit, med(hs,ch), med(ho,co), cp >> out
          printf "%-4s %-8s %-9.4f %-8s %-9.4f %-9.4f %d\n",
                 np, pit, med(pt,cp), hit, med(hs,ch), med(ho,co), cp }
  ' "$RAW"
done
echo; echo "결과 TSV: $OUT"
