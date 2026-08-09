#!/usr/bin/env bash
# 케이스 실행 — 컨테이너 안에서 mpirun (LOOP C005, 계약: PLAN §3-3)
# 사용: run.sh            (CUPID_NP rank 로 실행, 기본 4)
#       CUPID_NP=1 run.sh (직렬 검증)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/env.sh"

"$HERE/prepare_case.sh"

ts="$(date +%Y%m%d_%H%M)"
log="$CUPID_CASE/run_${ts}_np${CUPID_NP}.log"
echo "run: np=$CUPID_NP  case=$CUPID_CASE  log=$log"
# ulimit: -auto 컴파일(로컬 전부 스택) + 대형 자동 배열 → 기본 8MB 스택에서 SIGSEGV (LOG C005-r4/r5 로 확정)
# I_MPI_FABRICS: WSL2 단일 노드 공유메모리 통신 강제 (PLAN §6)
"$HERE/in_contain.sh" bash -c "ulimit -s unlimited && cd '$CUPID_CASE' && I_MPI_FABRICS=shm mpirun -np $CUPID_NP ./cupid.x" 2>&1 | tee "$log"
