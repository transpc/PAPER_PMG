#!/usr/bin/env bash
# 임의 명령을 컨테이너 안에서 실행하는 유일한 통로 (LOOP C001, 계약: PLAN §2-2)
# 사용: ./in_contain.sh <명령> [인자...]   예) ./in_contain.sh mpiifort --version
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
[[ -f "$SIF" ]] || { echo "SIF not found: $SIF — env.sh 기본값 확인 또는 export SIF=<경로>" >&2; exit 1; }
# /root 는 자동 바인드되지만, 리포가 /root 밖에 있는 체계 대비로 CODE_DIR 명시 바인드 유지
exec apptainer exec --bind "$CODE_DIR" "$SIF" "$@"
