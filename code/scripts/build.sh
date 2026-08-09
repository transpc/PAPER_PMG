#!/usr/bin/env bash
# CUPID 전체 빌드 — 컨테이너 안에서 make (LOOP C002, 계약: PLAN §3-1)
# 사용: build.sh [clean]
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/env.sh"

# makefile.in 교체: 원본 1회 백업 후 apptainer 변형 적용 (PLAN §3-1)
[[ -f "$CUPID_SRC/makefile.in.orig" ]] || cp "$CUPID_SRC/makefile.in" "$CUPID_SRC/makefile.in.orig"
cp "$CUPID_SRC/makefile.in.apptainer" "$CUPID_SRC/makefile.in"

# .mod 출력 디렉토리 — 전 서브 makefile 이 -module (상대경로)/Modules 사용, git 에는 없음 (C002-r1)
mkdir -p "$CUPID_SRC/Modules"

if [[ "${1:-}" == "clean" ]]; then
  exec "$HERE/in_contain.sh" make -C "$CUPID_SRC" clean
fi
exec "$HERE/in_contain.sh" make -C "$CUPID_SRC"
