#!/usr/bin/env bash
# 시간 비교 전용 — 최적화 플래그를 맞춘 빌드로 재측정 (2026-08-20)
#
# 왜 필요한가: 기본 빌드는 양쪽 최적화 수준이 다르다.
#   PMG   : makefile.in FOPTFLAGS = -O -fp-model strict ... -no-simd -no-vec
#           (골든 bitwise 재현을 위한 빌드 — 벡터화 끔, FP 재결합 금지)
#   hypre : configure CFLAGS=-O2 (자동 벡터화 켜짐, 기본 fp-model)
# 이 상태의 벽시계 비교는 PMG 에 불리하다. 반복수(its)는 산술 차이와 무관하므로
# 기본 빌드 결과를 그대로 쓰면 되지만, 시간은 이 스크립트로 다시 재야 한다.
#
# 두 빌드를 모두 재고 나란히 보고한다:
#   strict = build/       (기본, -fp-model strict -no-vec)
#   fast   = build_fast/  (-O3 -xcore-AVX2, hypre 도 install-fast 로 동일 수준)
#
# 통계: REPEAT 회 반복의 **중앙값**. 최솟값이 아닌 이유 — WSL2 에서 클럭 이상치가
# 관측되어(2026-08-20) 최솟값은 이상치를 그대로 채택한다. 0 이하 표본은 버린다.
#
# 주의: build_fast 는 bitwise 재현 빌드가 아니다. 골든 게이트(run_tests.sh)에는
#       절대 쓰지 말 것 — 그쪽은 build/ 의 strict 빌드를 유지한다.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IC="$HERE/../scripts/in_contain.sh"
REPEAT="${REPEAT:-5}"
OUT="${OUT:-$HERE/tests/fair_timing.tsv}"
RAW="$HERE/tests/.fair_raw"

echo "== hypre fast 변형 빌드 =="
HYPRE_VARIANT=fast "$HERE/../scripts/build_hypre.sh" || exit 1
echo "== PMG strict 빌드 =="
"$IC" make -C "$HERE" driver_hypre >/dev/null 2>&1 || { echo "BUILD FAIL(strict)"; exit 1; }
echo "== PMG fast 빌드 (build_fast/) =="
"$IC" make -C "$HERE" driver_hypre_fast >/dev/null 2>&1 || { echo "BUILD FAIL(fast)"; exit 1; }

: > "$OUT"
printf "build\tcase\ttag\tncell\tnp\tpmg_its\tpmg_tprep\tpmg_t\thyp_its\thyp_tasm\thyp_tamg\thyp_tsol\tnsample\n" >> "$OUT"

# measure <빌드라벨> <바이너리디렉토리> <이름> <ncell> <np> <드라이버 인자...>
measure() {
  local blab="$1" bdir="$2" name="$3" nc="$4" np="$5"; shift 5
  local cmd="../$bdir/driver_hypre $*"
  [ "$np" -gt 1 ] && cmd="I_MPI_FABRICS=shm mpirun -np $np $cmd"
  : > "$RAW"
  local r out l
  for r in $(seq 1 "$REPEAT"); do
    out=$("$IC" bash -c "ulimit -s unlimited && cd '$HERE/tests' && rm -rf MG_tmp fort.* && HYPRE_CFG=$blab $cmd" 2>&1)
    echo "$out" | grep -q "VERDICT PASS" || { echo "  !! $blab/$name 실행 실패 (repeat $r)"; continue; }
    while IFS= read -r l; do
      # tag its_p tprep t_pmg its_h tasm tset tsol
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$(sed -n 's/.*RESULT \([^ ]*\).*/\1/p' <<<"$l")" \
        "$(sed -n 's/.*pmg_its= *\([^ ]*\).*/\1/p'   <<<"$l")" \
        "$(sed -n 's/.*pmg_tprep= *\([^ ]*\).*/\1/p' <<<"$l")" \
        "$(sed -n 's/.* pmg_t= *\([^ ]*\).*/\1/p'    <<<"$l")" \
        "$(sed -n 's/.*hyp_its= *\([^ ]*\).*/\1/p'   <<<"$l")" \
        "$(sed -n 's/.*hyp_tasm= *\([^ ]*\).*/\1/p'  <<<"$l")" \
        "$(sed -n 's/.*hyp_tset= *\([^ ]*\).*/\1/p'  <<<"$l")" \
        "$(sed -n 's/.*hyp_tsol= *\([^ ]*\).*/\1/p'  <<<"$l")" >> "$RAW"
    done < <(echo "$out" | grep "^ RESULT")
  done
  [ -s "$RAW" ] || { echo "  !! $blab/$name 표본 없음"; return 1; }
  # tag 별 중앙값 (0 이하 표본 제외). hyp_tamg = hyp_tset − hyp_tasm (순수 AMG 셋업)
  awk -F'\t' -v b="$blab" -v n="$name" -v nc="$nc" -v np="$np" -v out="$OUT" '
    function med(arr, cnt,   i, tmp, m) {
      if (cnt == 0) return -1
      for (i = 1; i <= cnt; i++) tmp[i] = arr[i]
      # 삽입정렬 (표본 수가 작다)
      for (i = 2; i <= cnt; i++) { m = tmp[i]; j = i-1
        while (j > 0 && tmp[j] > m) { tmp[j+1] = tmp[j]; j-- }; tmp[j+1] = m }
      return (cnt % 2) ? tmp[(cnt+1)/2] : (tmp[cnt/2] + tmp[cnt/2+1])/2.0
    }
    { t = $1
      if (!(t in pit)) { pit[t] = $2; hit[t] = $5 }
      # 시간 표본은 양수만 채택 (클럭 이상치 배제). tamg 는 0 이어도 유효(k2 재사용)
      if ($3+0 >  0) { pp[t][++cp[t]] = $3+0 }
      if ($4+0 >  0) { pt[t][++ct[t]] = $4+0 }
      if ($6+0 >= 0 && $7+0 >= 0) { ha[t][++ca[t]] = $6+0; hg[t][++cg[t]] = ($7+0)-($6+0) }
      if ($8+0 >  0) { ho[t][++co[t]] = $8+0 }
    }
    END { for (t in pit)
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%.4f\t%.4f\t%s\t%.4f\t%.4f\t%.4f\t%d\n",
                   b, n, t, nc, np, pit[t], med(pp[t],cp[t]), med(pt[t],ct[t]),
                   hit[t], med(ha[t],ca[t]), med(hg[t],cg[t]), med(ho[t],co[t]),
                   ct[t] >> out
    }
  ' "$RAW"
  awk -F'\t' -v b="$blab" -v n="$name" '$1==b && $2==n {
      printf "  %-6s %-10s %-4s  pmg(its=%s prep=%s solve=%s)  hyp(its=%s asm=%s amg=%s sol=%s)  n=%s\n",
             $1,$2,$3,$6,$7,$8,$9,$10,$11,$12,$13 }' "$OUT"
}

CASES_RUN="${CASES:-g48 g80 np4 ect1_s150}"
for pair in "strict build" "fast build_fast"; do
  set -- $pair; blab=$1; bdir=$2
  echo "== $blab 빌드 =="
  for c in $CASES_RUN; do
    case "$c" in
      g48)  measure "$blab" "$bdir" g48  110592 1  48 48 48 1.0 ;;
      g80)  measure "$blab" "$bdir" g80  512000 1  80 80 80 1.0 ;;
      np4)  measure "$blab" "$bdir" np4  110592 4  48 48 48 1.0 ;;
      ect1_s150)
        [ -f "$HERE/golden/ECT1_436k_np1/setup_r0.bin" ] && \
          measure "$blab" "$bdir" ect1_s150 436136 1  replay ../golden/ECT1_436k_np1 150 ;;
    esac
  done
done

echo; echo "결과 TSV: $OUT"
