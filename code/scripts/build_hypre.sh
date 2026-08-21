#!/usr/bin/env bash
# hypre 빌드 — PMG vs hypre 비교 하네스(pmg_standalone/driver_hypre)용
#
# 소스와 산출물은 CUPID 본체와 **한 곳을 공유**한다: code/Source/HYPRE/
#   - 버전 고정 tarball 이 리포에 들어 있어(hypre-<VER>.tar.gz) 네트워크가
#     필요 없다. 컨테이너/HPC 네이티브/오프라인 모두 같은 경로.
#   - 기본(base) 변형은 Source/HYPRE/makefile 이 CUPID `make` 중에 자동으로
#     세운다. 이 스크립트는 그것과 동일한 산출물을 만들거나(없을 때),
#     시간 비교용 fast 변형을 추가로 세운다.
#   ※ 두 벌을 따로 두면 버전·플래그가 조용히 갈라져 스탠드얼론 측정치와
#     프로덕션 측정치를 같은 표에 올릴 수 없다 — 그래서 일원화한다.
#
# 산출물: code/Source/HYPRE/hypre-<VER>/{install,install-fast}/{lib,include}
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/env.sh"

VER="${HYPRE_VER:-2.33.0}"
HDIR="$CODE_DIR/Source/HYPRE"
SRC="$HDIR/hypre-$VER"
TARBALL="$HDIR/hypre-$VER.tar.gz"
# 빌드 변형: 기본(install, -O2) / fast(install-fast, -O3+AVX2)
# — PMG 는 makefile.in 이 -fp-model strict -no-simd -no-vec (bitwise 재현 빌드) 라
#   시간 비교를 하려면 양쪽 최적화 수준을 맞춘 변형이 필요하다 (2026-08-20).
VARIANT="${HYPRE_VARIANT:-base}"
if [ "$VARIANT" = fast ]; then
  PREFIX="$SRC/install-fast"
  BFLAGS="${HYPRE_CFLAGS:--O3 -xcore-AVX2}"
  BDIR="$SRC/src-fast"
else
  PREFIX="$SRC/install"
  BFLAGS="${HYPRE_CFLAGS:--O2}"
  BDIR="$SRC/src"
fi

if [ ! -d "$SRC/src" ]; then
  [ -f "$TARBALL" ] || { echo "tarball 없음: $TARBALL"; exit 1; }
  echo "== hypre $VER 전개 (in-tree tarball) =="
  tar xzf "$TARBALL" -C "$HDIR"
fi

if [ -f "$PREFIX/lib/libHYPRE.a" ]; then
  echo "== 이미 빌드됨: $PREFIX/lib/libHYPRE.a =="
  exit 0
fi

# fast 변형은 소스 트리를 따로 두어 base 빌드의 오브젝트와 섞이지 않게 한다
if [ "$BDIR" != "$SRC/src" ] && [ ! -d "$BDIR" ]; then
  cp -a "$SRC/src" "$BDIR"
  "$HERE/in_contain.sh" bash -c "cd '$BDIR' && make distclean >/dev/null 2>&1" || true
fi

echo "== hypre $VER [$VARIANT] configure + make (mpiicc, CFLAGS='$BFLAGS') =="
"$HERE/in_contain.sh" bash -c "
  cd '$BDIR' &&
  ./configure CC=mpiicc CFLAGS='$BFLAGS' F77=mpiifort FFLAGS='$BFLAGS' \
              --prefix='$PREFIX' --without-openmp --disable-shared \
  > configure.log 2>&1 &&
  make -j${MAKE_JOBS} > make.log 2>&1 &&
  make install > install.log 2>&1
" || { echo "BUILD FAIL — $BDIR/{configure,make,install}.log 확인"; exit 1; }

echo "== OK: $PREFIX/lib/libHYPRE.a =="
