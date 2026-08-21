#!/usr/bin/env bash
# PMG vs hypre 광범위 스윕 (2026-08-20) — run_hypre_compare.sh 의 확장판
#
# 축:
#   A) 격자 크기 (np=1)         : 약수렴 상수의 문제크기 의존성
#   B) 이방성 aspect (np=1)     : semi-coarsening(PMG) vs 대수적 강연결(AMG)
#   C) np 스케일 (동일 문제)    : 분할 강건성 — 논문의 np-무관 수렴 주장 대조군
#   D) 실제 CUPID 행렬 (재생)   : 합성 Poisson 이 놓치는 현실 조건수
#   E) hypre 설정 변형          : default(ℓ1-GS 13/14) vs sgs(hybrid SGS 6)
#                                 vs agg(aggressive coarsening — 셋업 절감)
#
# 타이밍은 REPEAT 회 반복 후 최솟값 채택 (WSL2 노이즈가 커서 1회 측정은 무의미).
# 결과는 TSV 로 $OUT 에 누적 — 표 요약은 stdout.
#
# 사용: ./sweep_hypre.sh [축...]   (인자 없으면 A B C D E 전부)
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IC="$HERE/../scripts/in_contain.sh"
OUT="${OUT:-$HERE/tests/sweep_results.tsv}"
REPEAT="${REPEAT:-3}"

"$IC" make -C "$HERE" driver_hypre >/dev/null 2>&1 || { echo "BUILD FAIL"; exit 1; }
: > "$OUT"
printf "axis\tcase\tcfg\ttag\tncell\tnp\tpmg_its\tpmg_t\thyp_its\thyp_tset\thyp_tsol\tfid\tverdict\n" >> "$OUT"

# run <축> <이름> <cfg> <ncell> <np> <드라이버 인자...>
run() {
  local axis="$1" name="$2" cfg="$3" nc="$4" np="$5"; shift 5
  local envs="HYPRE_CFG=$cfg"
  case "$cfg" in
    sgs)     envs="$envs HYPRE_RELAX=6" ;;
    agg)     envs="$envs HYPRE_AGG=1" ;;
    default) ;;
  esac
  local cmd="../build/driver_hypre $*"
  [ "$np" -gt 1 ] && cmd="I_MPI_FABRICS=shm mpirun -np $np $cmd"
  local r out rc best_line
  declare -A bp bh1 bh2   # 최솟값 누적 (tag 별)
  local tags=""
  for r in $(seq 1 "$REPEAT"); do
    out=$("$IC" bash -c "ulimit -s unlimited && cd '$HERE/tests' && rm -rf MG_tmp fort.* && $envs $cmd" 2>&1)
    rc=$?
    if ! echo "$out" | grep -q "VERDICT"; then
      [ "$r" -eq "$REPEAT" ] && {
        printf "%s\t%s\t%s\t-\t%s\t%s\t?\t?\t?\t?\t?\t?\tFAIL(no VERDICT)\n" \
          "$axis" "$name" "$cfg" "$nc" "$np" >> "$OUT"
        return 1
      }
      sleep 1; continue
    fi
    local verdict=FAIL
    echo "$out" | grep -q "VERDICT PASS" && [ "$rc" -eq 0 ] && verdict=PASS
    # DIVERGE 판정: 드라이버가 FAIL 이거나 hypre maxiter 도달
    while IFS= read -r line; do
      local tag pits pt hits hts hso fid
      tag=$(echo  "$line" | awk '{print $2}')
      pits=$(echo "$line" | sed -n 's/.*pmg_its= *\([^ ]*\).*/\1/p')
      pt=$(echo   "$line" | sed -n 's/.*pmg_t= *\([^ ]*\).*/\1/p')
      hits=$(echo "$line" | sed -n 's/.*hyp_its= *\([^ ]*\).*/\1/p')
      hts=$(echo  "$line" | sed -n 's/.*hyp_tset= *\([^ ]*\).*/\1/p')
      hso=$(echo  "$line" | sed -n 's/.*hyp_tsol= *\([^ ]*\).*/\1/p')
      fid=$(echo  "$line" | sed -n 's/.* fid= *\([^ ]*\).*/\1/p')
      echo "$tags" | grep -qw "$tag" || tags="$tags $tag"
      # 최솟값 갱신 (bash 산술은 정수뿐 → awk 로 비교)
      local key="$tag"
      if [ -z "${bp[$key]:-}" ]; then
        bp[$key]="$pits|$pt"; bh1[$key]="$hits|$hts"; bh2[$key]="$hso|$fid|$verdict"
      else
        local op="${bp[$key]#*|}" oh1="${bh1[$key]#*|}" oh2="${bh2[$key]%%|*}"
        awk -v a="$pt"  -v b="$op"  'BEGIN{exit !(a<b)}' && bp[$key]="$pits|$pt"
        awk -v a="$hts" -v b="$oh1" 'BEGIN{exit !(a<b)}' && bh1[$key]="$hits|$hts"
        awk -v a="$hso" -v b="$oh2" 'BEGIN{exit !(a<b)}' && bh2[$key]="$hso|$fid|$verdict"
      fi
    done < <(echo "$out" | grep "^ RESULT")
  done
  local t
  for t in $tags; do
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$axis" "$name" "$cfg" "$t" "$nc" "$np" \
      "${bp[$t]%%|*}" "${bp[$t]#*|}" "${bh1[$t]%%|*}" "${bh1[$t]#*|}" \
      "${bh2[$t]%%|*}" "$(echo "${bh2[$t]}" | cut -d'|' -f2)" \
      "$(echo "${bh2[$t]}" | cut -d'|' -f3)" >> "$OUT"
    printf "  %-8s %-10s %-8s %-5s pmg(its=%s t=%s)  hyp(its=%s set=%s sol=%s)  %s\n" \
      "$axis" "$name" "$cfg" "$t" "${bp[$t]%%|*}" "${bp[$t]#*|}" \
      "${bh1[$t]%%|*}" "${bh1[$t]#*|}" "${bh2[$t]%%|*}" \
      "$(echo "${bh2[$t]}" | cut -d'|' -f3)"
  done
}

AXES="${*:-A B C D E}"
has() { echo " $AXES " | grep -q " $1 "; }

# ---- A) 격자 크기 (np=1, 등방) ----
if has A; then
  echo "== A) 격자 크기 스케일 (np=1, 등방 Poisson) =="
  run A g24  default 13824   1  24 24 24 1.0
  run A g32  default 32768   1  32 32 32 1.0
  run A g48  default 110592  1  48 48 48 1.0
  run A g64  default 262144  1  64 64 64 1.0
  run A g80  default 512000  1  80 80 80 1.0
fi

# ---- B) 이방성 (np=1, 32³) ----
if has B; then
  echo "== B) 이방성 aspect (np=1, 32^3) =="
  run B a1     default 32768 1  32 32 32 1.0
  run B a10    default 32768 1  32 32 32 10.0
  run B a100   default 32768 1  32 32 32 100.0
  run B a1000  default 32768 1  32 32 32 1000.0
fi

# ---- C) np 스케일 (48³ 고정 문제, 슬랩 분할) ----
# LOG C010-2: PMG 는 np>=6 에서 수렴 붕괴 미해결 — 그 구간이 그대로 드러난다.
if has C; then
  echo "== C) np 스케일 (48^3 고정, 슬랩 분할) =="
  for n in 1 2 4 8 16; do
    run C np$n default 110592 "$n"  48 48 48 1.0
  done
  echo "== C') np 스케일 (48^3, 3D 블록 분할 ipart=1) =="
  for n in 2 4 8; do
    run Cb np$n default 110592 "$n"  48 48 48 1.0 1
  done
fi

# ---- D) 실제 CUPID 행렬 (골든 재생, np=1) ----
if has D; then
  echo "== D) 실제 행렬 (골든 재생, np=1) =="
  gold() {
    local dir="$1" pfx="$2"; shift 2
    local G="$HERE/golden/$dir" s
    [ -f "$G/setup_r0.bin" ] || { echo "  SKIP $pfx (골든 없음)"; return; }
    for s in "$@"; do
      [ -f "$G/s${s}_k1_r0_c1.pre" ] || continue
      run D "${pfx}_s$s" default 436136 1  replay "../golden/$dir" "$s"
    done
  }
  gold iSMR436k_np1  gold 1 10 30
  gold ECT1_436k_np1 ect1 1 10 30 150
fi

# ---- E) hypre 설정 변형 (대표 케이스에서) ----
if has E; then
  echo "== E) hypre 설정 변형 =="
  for c in default sgs agg; do
    run E g48    "$c" 110592 1  48 48 48 1.0
    run E a100   "$c" 32768  1  32 32 32 100.0
    run E np8    "$c" 110592 8  48 48 48 1.0
    [ -f "$HERE/golden/ECT1_436k_np1/setup_r0.bin" ] && \
      run E ect1_s150 "$c" 436136 1  replay ../golden/ECT1_436k_np1 150
  done
fi

echo
echo "결과 TSV: $OUT"
