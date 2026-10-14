#!/usr/bin/env bash
# Build gravity_core.dll (or .so) and the physics test harness.
set -e
cd "$(dirname "$0")"
mkdir -p bin
CXX=${CXX:-g++}

if [[ "$OSTYPE" == msys* || "$OSTYPE" == cygwin* || "$OSTYPE" == win32 ]]; then
  "$CXX" -O2 -std=c++17 -shared -static -static-libgcc -static-libstdc++ \
    -o bin/gravity_core.dll core/api.cpp core/expr.cpp core/sim.cpp
else
  "$CXX" -O2 -std=c++17 -shared -fPIC -o bin/gravity_core.so \
    core/api.cpp core/expr.cpp core/sim.cpp
fi
"$CXX" -O2 -std=c++17 -static-libgcc -static-libstdc++ \
  -o bin/gs_test core/test_main.cpp core/expr.cpp core/sim.cpp
"$CXX" -O2 -std=c++17 -static-libgcc -static-libstdc++ \
  -o bin/gs_bench core/bench_main.cpp core/expr.cpp core/sim.cpp
echo "Build OK. Run tests with: bin/gs_test"
