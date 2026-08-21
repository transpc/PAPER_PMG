#!/usr/bin/env bash
# 프로덕션 관점 시간 재측정 (2026-08-21)
#
# 무엇을 비교하는가 — CUPID 는 타임스텝마다 압력을 다시 푼다. 따라서 반복
# 계산에서 실제로 매번 지불하는 비용끼리 견주는 것이 맞다:
#     PMG   : RAP (SOLVE_GMG(icase=1) 내부 Galerkin 재계산)
#     hypre : BoomerAMG 셋업 전체 (강연결+조대화+보간+RAP)
# PMG 의 계층/전달연산자 구축(subdomain_infor_MG 등)은 1회성이므로 별도 계상
# 한다. 단 그것은 rank 0 직렬이라 np 가 커지면 무시할 수 없다 (run_np_scan 참조).
#
# 계측 전제 (앞선 오측정의 교훈):
#   - 시계: driver 가 clock_gettime(CLOCK_MONOTONIC) 을 직접 바인딩한다.
#           MPI_WTIME / ifort SYSTEM_CLOCK 은 이 환경에서 CLOCK_REALTIME 이라
#           WSL2 시각 동기화 때 역행한다 (음수 경과시간 실측).
#   - 격리: tests_ab/ 에서 실행한다. tests/ 는 run_tests.sh 가 쓰므로 겹치면
#           양쪽 결과가 모두 오염된다 (2026-08-21 실제 충돌).
#   - 부하: 다른 작업이 코어를 점유하면 2~3배까지 부풀려진다. 실행 전
#           `uptime` 으로 확인할 것.
#   - 통계: REPEAT 회 중앙값, 0 이하 표본 배제.
#
# hypre 레벨 수도 함께 기록한다 — Num levels=1 이면 조대화 실패(=AMG 아님)라
# 비교의 의미가 완전히 달라진다. 실제로 쉬운 골든 케이스들이 그랬다.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IC="$HERE/../scripts/in_contain.sh"
WD="$HERE/tests_ab"
REPEAT="${REPEAT:-5}"
OUT="${OUT:-$WD/prod_timing.tsv}"

mkdir -p "$WD"; [ -f "$WD/mg.in" ] || cp "$HERE/tests/mg.in" "$WD/"
"$IC" make -C "$HERE" driver_hypre >/dev/null 2>&1 || { echo "BUILD FAIL"; exit 1; }
: > "$OUT"
printf "case\tlev\tpmg_its\tpmg_rap\tpmg_solve\tpmg_ser1x\thyp_its\thyp_amg\thyp_solve\tn\n" >> "$OUT"

echo "부하 확인: $(uptime | sed 's/.*load average/load/')"
printf "%-12s %-4s %-8s %-9s %-10s %-10s %-8s %-9s %-10s\n" \
  "case" "lev" "pmg_its" "PMGrap" "PMGsolve" "PMG1x(ser)" "hyp_its" "hypAMG" "hypSolve"

measure() {   # $1=이름  나머지=드라이버 인자
  local name="$1"; shift
  local raw="$WD/.prod_raw"; : > "$raw"
  local lev="?" out r
  # 레벨 수는 1회만 (반복해도 같다)
  out=$("$IC" bash -c "ulimit -s unlimited && cd '$WD' && rm -rf MG_tmp fort.* && HYPRE_PRINT=3 ../build/driver_hypre $*" 2>&1)
  lev=$(sed -n 's/.*Num levels = *\([0-9]*\).*/\1/p' <<<"$out" | head -1)
  for r in $(seq 1 "$REPEAT"); do
    out=$("$IC" bash -c "ulimit -s unlimited && cd '$WD' && rm -rf MG_tmp fort.* && ../build/driver_hypre $*" 2>&1)
    echo "$out" | grep -q "VERDICT PASS" || continue
    # k1(=재생) 또는 syn(=합성) 행만: RAP 이 분리 계측된 행
    echo "$out" | grep -E "^ RESULT (k1|syn)" | while IFS= read -r l; do
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$(sed -n 's/.*pmg_its= *\([^ ]*\).*/\1/p'   <<<"$l")" \
        "$(sed -n 's/.*pmg_trap= *\([^ ]*\).*/\1/p'  <<<"$l")" \
        "$(sed -n 's/.* pmg_t= *\([^ ]*\).*/\1/p'    <<<"$l")" \
        "$(sed -n 's/.*pmg_tser= *\([^ ]*\).*/\1/p'  <<<"$l")" \
        "$(sed -n 's/.*hyp_its= *\([^ ]*\).*/\1/p'   <<<"$l")" \
        "$(sed -n 's/.*hyp_tset= *\([^ ]*\).*/\1/p'  <<<"$l")" \
        "$(sed -n 's/.*hyp_tsol= *\([^ ]*\).*/\1/p'  <<<"$l")" >> "$raw"
    done
  done
  [ -s "$raw" ] || { printf "%-12s %s\n" "$name" "표본 없음"; return; }
  # hyp_amg = hyp_tset − hyp_tasm 는 tasm 이 tset 에 포함돼 있어 따로 빼야 하나,
  # 여기서는 tset 을 그대로 쓴다 (IJ 조립은 하네스 고유 비용이라 아래 주석 참조)
  awk -F'\t' -v nm="$name" -v lv="${lev:-?}" -v out="$OUT" '
    function med(a,c,  i,j,t,m){if(c==0)return -1; for(i=1;i<=c;i++)t[i]=a[i];
      for(i=2;i<=c;i++){m=t[i];j=i-1;while(j>0&&t[j]>m){t[j+1]=t[j];j--}t[j+1]=m}
      return (c%2)?t[(c+1)/2]:(t[c/2]+t[c/2+1])/2.0}
    { pit=$1; hit=$5
      if($2+0>0) rp[++cr]=$2+0
      if($3+0>0) sv[++cs]=$3+0
      if($4+0>0) se[++ce]=$4+0
      if($6+0>0) hs[++ch]=$6+0
      if($7+0>0) ho[++co]=$7+0 }
    END{ printf "%s\t%s\t%s\t%.4f\t%.4f\t%.4f\t%s\t%.4f\t%.4f\t%d\n",
                nm,lv,pit,med(rp,cr),med(sv,cs),med(se,ce),hit,med(hs,ch),med(ho,co),cs >> out
         printf "%-12s %-4s %-8s %-9.4f %-10.4f %-10.4f %-8s %-9.4f %-10.4f\n",
                nm,lv,pit,med(rp,cr),med(sv,cs),med(se,ce),hit,med(hs,ch),med(ho,co) }
  ' "$raw"
}

measure syn48  48 48 48 1.0
for s in 1 10 30;      do measure "gold_s$s" replay ../golden/iSMR436k_np1  $s; done
for s in 1 10 30 150;  do measure "ect1_s$s" replay ../golden/ECT1_436k_np1 $s; done

echo; echo "결과 TSV: $OUT"
echo "주의: hypAMG 는 IJ 조립(하네스 고유 비용)을 포함한다 — RESULT 행의"
echo "      hyp_tasm 을 빼면 순수 BoomerAMG 셋업이다."
