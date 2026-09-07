#!/usr/bin/env bash
set -euo pipefail

VERSION="10.0.0"
COMMIT="9faf13dd8f1128dd75ed233a9627027b0422c3f7"
ROOT="${CHRONO_ROOT:-/home/ubuntu/third_party}"
SRC="$ROOT/chrono-$VERSION"
BUILD="$ROOT/chrono-build-$VERSION"
PREFIX="$ROOT/chrono-install-$VERSION"

if [[ ! -d "$SRC/.git" ]]; then
  git clone --depth 1 --branch "$VERSION" https://github.com/projectchrono/chrono.git "$SRC"
fi
[[ "$(git -C "$SRC" rev-parse HEAD)" == "$COMMIT" ]] || { echo "CHRONO_SOURCE_COMMIT_MISMATCH" >&2; exit 2; }

cmake -S "$SRC" -B "$BUILD" -G 'Unix Makefiles' \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DCH_ENABLE_MODULE_PYTHON=ON \
  -DCH_ENABLE_MODULE_CASCADE=OFF \
  -DCH_ENABLE_MODULE_IRRLICHT=OFF \
  -DCH_ENABLE_MODULE_VSG=OFF \
  -DCH_ENABLE_MODULE_VEHICLE=OFF \
  -DCH_ENABLE_MODULE_FEA=OFF \
  -DCH_ENABLE_MODULE_FSI=OFF \
  -DCH_ENABLE_MODULE_SENSOR=OFF \
  -DCH_ENABLE_MODULE_POSTPROCESS=OFF \
  -DCH_ENABLE_MODULE_PARDISO_MKL=OFF
cmake --build "$BUILD" --parallel "${CHRONO_JOBS:-4}"

echo "PYCHRONO_BUILD=PASS"
echo "PYCHRONO_SOURCE_COMMIT=$COMMIT"
echo "PYCHRONO_PYTHONPATH=$BUILD/bin"
echo "PYCHRONO_LD_LIBRARY_PATH=$BUILD/lib"
