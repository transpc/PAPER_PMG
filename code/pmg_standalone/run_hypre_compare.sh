#!/usr/bin/env bash
# PMG vs hypre(BiCGSTAB+BoomerAMG) 비교 러너 (2026-08-20)
# 선행: ../scripts/build_hypre.sh (hypre v2.33.0 — code/external/ 에 설치)
#
# 케이스 구성은 run_tests.sh 와 동일 (합성 iso/aniso + np 스케일 + 골든 재생).
# 드라이버(driver_hypre)가 한 프로세스에서 두 솔버를 같은 행렬로 순차 실행하고
# 'RESULT <tag> pmg_its= .. pmg_t= .. hyp_its= .. hyp_tset= .. hyp_tsol= ..
#  fid= .. reldiff= ..' 라인을 출력 — 여기서 표로 집계한다.
#   pmg_t   : SOLVE_GMG 벽시계 (icase=1 은 RAP 등 per-solve 셋업 포함)
#   hyp_tset: hypre IJ 조립 + BoomerAMG setup / hyp_tsol: BiCGSTAB 반복
#   fid     : mean(Δ) 제거 후 max|u_hypre−u_pmg| / (1+max|u_pmg,k1|)
#   판정    : 드라이버 VERDICT (합성 = 두 솔버 독립잔차 게이트,
#             재생 = hypre 자기보고 relres ≤ eps) + fid > 1e-4 는 드라이버가 WARN
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IC="$HERE/../scripts/in_contain.sh"

"$IC" make -C "$HERE" driver_hypre >/dev/null 2>&1 || { echo "BUILD FAIL"; exit 1; }

fail=0
printf "%-12s %-4s %-9s %-8s %-8s %-8s %-9s %-9s %-10s %s\n" \
  "case" "tag" "ncell" "pmg_its" "pmg_t" "hyp_its" "hyp_tset" "hyp_tsol" "fid" "verdict"

# 실행 + 집계. WSL2/apptainer 바인드에서 MG_tmp 쓰기→읽기 사이 간헐 지연이
# 관측되어 (2026-08-20, driver 는 동일한데 재실행 시 정상) 1회 재시도한다.
run_case() {  # $1=이름 $2=ncell 표기, 이후 = 드라이버 실행 커맨드 (tests/ 기준)
  local name="$1" nc="$2"; shift 2
  local out rc try
  for try in 1 2; do
    out=$("$IC" bash -c "ulimit -s unlimited && cd '$HERE/tests' && rm -rf MG_tmp fort.* && $*" 2>&1)
    rc=$?
    echo "$out" | grep -q "VERDICT" && break
    sleep 1   # 일시 실패 (MG_tmp 지연 등) — 재시도
  done
  local verdict=FAIL
  echo "$out" | grep -q "VERDICT PASS" && [ "$rc" -eq 0 ] && verdict=PASS
  echo "$out" | grep -q "WARN\[" && verdict="$verdict+WARN(fid)"
  [ "${verdict%%+*}" = FAIL ] && fail=1
  local n=0
  while IFS= read -r line; do
    n=$((n+1))
    # 필드 파싱 ('라벨= 값' 쌍)
    local tag pmg_its pmg_t hyp_its hyp_tset hyp_tsol fid
    tag=$(echo "$line" | awk '{print $2}')
    pmg_its=$(echo "$line" | sed -n 's/.*pmg_its= *\([^ ]*\).*/\1/p')
    pmg_t=$(echo "$line"   | sed -n 's/.*pmg_t= *\([^ ]*\).*/\1/p')
    hyp_its=$(echo "$line" | sed -n 's/.*hyp_its= *\([^ ]*\).*/\1/p')
    hyp_tset=$(echo "$line"| sed -n 's/.*hyp_tset= *\([^ ]*\).*/\1/p')
    hyp_tsol=$(echo "$line"| sed -n 's/.*hyp_tsol= *\([^ ]*\).*/\1/p')
    fid=$(echo "$line"     | sed -n 's/.* fid= *\([^ ]*\).*/\1/p')
    printf "%-12s %-4s %-9s %-8s %-8s %-8s %-9s %-9s %-10s %s\n" \
      "$name" "$tag" "$nc" "${pmg_its:-?}" "${pmg_t:-?}" "${hyp_its:-?}" \
      "${hyp_tset:-?}" "${hyp_tsol:-?}" "${fid:-?}" "$verdict"
  done < <(echo "$out" | grep "^ RESULT")
  if [ "$n" -eq 0 ]; then
    printf "%-12s %-4s %-9s %-8s %-8s %-8s %-9s %-9s %-10s %s\n" \
      "$name" "-" "$nc" "?" "?" "?" "?" "?" "?" "FAIL(no RESULT)"
    fail=1
  fi
}

# ---- 합성 (run_tests.sh 의 CASES/MCASES 와 동일 구성) ----
run_case iso24     13824  "../build/driver_hypre 24 24 24 1.0"
run_case aniso24   13824  "../build/driver_hypre 24 24 24 100.0"
run_case iso48     110592 "../build/driver_hypre 48 48 48 1.0"
run_case iso24_np2 13824  "I_MPI_FABRICS=shm mpirun -np 2 ../build/driver_hypre 24 24 24 1.0"
run_case iso24_np4 13824  "I_MPI_FABRICS=shm mpirun -np 4 ../build/driver_hypre 24 24 24 1.0"

# ---- 골든 재생 (바이너리가 있을 때만 — fresh clone 은 SKIP) ----
run_golden() {  # $1=디렉토리명 $2=접두어, 이후 = step 목록
  local dir="$1" pfx="$2"; shift 2
  local GOLD="$HERE/golden/$dir" s
  if [ ! -f "$GOLD/setup_r0.bin" ]; then
    printf "%-12s %-4s %-9s %s\n" "$pfx" "-" "-" "SKIP(채취 필요 — golden/$dir/meta.md)"
    return
  fi
  for s in "$@"; do
    if [ ! -f "$GOLD/s${s}_k1_r0_c1.pre" ]; then
      printf "%-12s %-4s %-9s %s\n" "${pfx}_s${s}" "-" "-" "SKIP(덤프 없음)"; continue
    fi
    run_case "${pfx}_s${s}" 436136 "../build/driver_hypre replay ../golden/$dir $s"
  done
}

run_golden iSMR436k_np1 gold 1 10 30
run_golden ECT1_436k_np1 ect1 1 10 30 150

exit $fail
