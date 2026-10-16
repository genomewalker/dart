#!/usr/bin/env bash
# libtaph consumer contract for DART: configure against <libtaph_dir>, build, run ctest.
# Called by libtaph/scripts/bump-libtaph.sh; touches nothing outside <build_dir>.
# Extra configure flags (compilers, CMAKE_PREFIX_PATH) come from DART_CMAKE_ARGS.
set -euo pipefail
[ $# -eq 2 ] || { echo "usage: $0 <build_dir> <libtaph_dir>" >&2; exit 2; }
B=$1
L=$2
SRC=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$B"
# shellcheck disable=SC2086
cmake -S "$SRC" -B "$B" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON \
    -DLIBTAPH_LOCAL_PATH="$L" ${DART_CMAKE_ARGS:-} > "$B/configure.log" 2>&1 \
    || { tail -20 "$B/configure.log" >&2; exit 1; }
cmake --build "$B" -j "${JOBS:-16}" > "$B/build.log" 2>&1 || { grep -m5 ' error ' "$B/build.log" >&2; exit 1; }
(cd "$B" && ctest -j "${JOBS:-8}" --output-on-failure > ctest.log 2>&1) || { tail -30 "$B/ctest.log" >&2; exit 1; }
echo "dart ctest: $(grep -oE '[0-9]+% tests passed, [0-9]+ tests failed out of [0-9]+' "$B/ctest.log")"
"$B/dart" --version | grep '^libtaph '
