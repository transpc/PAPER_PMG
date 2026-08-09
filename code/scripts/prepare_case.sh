#!/usr/bin/env bash
# 케이스 전처리 (LOOP C003, 계약: PLAN §3-2)
# 1) foamGrid.in 해제(호스트 7z — 컨테이너에는 7z 없음)  2) somaFlow.in 확인  3) cupid.x 복사
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
cd "$CUPID_CASE"

# 1) foamGrid.in — 없으면 7z 해제 (호스트에서 수행)
if [[ ! -f foamGrid.in ]]; then
  command -v 7z >/dev/null || { echo "7z not found on host — apt install p7zip-full" >&2; exit 1; }
  7z x -y foamGrid.7z
fi

# 2) somaFlow.in — 없으면 실행 불가 (PLAN §7-1 블로커)
[[ -f somaFlow.in ]] || {
  echo "BLOCKED: somaFlow.in 없음 — 원본 케이스에서 가져와 $CUPID_CASE 에 두세요 (PLAN §7-1)" >&2
  exit 2
}

# 3) 실행 파일 복사 (README 권장 방식)
[[ -x "$CUPID_SRC/cupid.x" ]] || { echo "cupid.x 없음 — scripts/build.sh 먼저 실행" >&2; exit 1; }
cp -f "$CUPID_SRC/cupid.x" .
echo "case ready: $CUPID_CASE"
