#!/usr/bin/env bash
# PMG 전처리 피크 메모리 측정 (PLAN_SETUP_RANKLOCAL §6-2)
# 사용: bash mem.sh <nx> <np>   — rank0peak(=단일 프로세스 최대 RSS)만 지표로 사용
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
