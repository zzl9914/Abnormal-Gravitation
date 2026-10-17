#!/usr/bin/env bash
# Build gravity_core.dll (or .so), the physics test harness and the benchmark.
# wasm 构建（前端直连内核）用 web/tools/build-kernel.js（需要 clang + wasm32）。
set -e
cd "$(dirname "$0")"
mkdir -p bin
CXX=${CXX:-g++}
SRC="src/core/expr.cpp src/core/sim.cpp src/core/gsmath.cpp"

if [[ "$OSTYPE" == msys* || "$OSTYPE" == cygwin* || "$OSTYPE" == win32 ]]; then
  # Simulation 含固定容量数组成员（约 1MB），测试里是栈对象，Windows 上必须加大栈
  "$CXX" -O2 -std=c++17 -shared -static -static-libgcc -static-libstdc++ \
    -o bin/gravity_core.dll src/core/api.cpp $SRC
  "$CXX" -O2 -std=c++17 -static-libgcc -static-libstdc++ -Wl,--stack,67108864 \
    -o bin/gs_test.exe src/core/test_main.cpp $SRC
  "$CXX" -O2 -std=c++17 -static-libgcc -static-libstdc++ -Wl,--stack,67108864 \
    -o bin/gs_bench.exe src/core/bench_main.cpp $SRC
else
  "$CXX" -O2 -std=c++17 -shared -fPIC -o bin/gravity_core.so src/core/api.cpp $SRC
  "$CXX" -O2 -std=c++17 -o bin/gs_test src/core/test_main.cpp $SRC
  "$CXX" -O2 -std=c++17 -o bin/gs_bench src/core/bench_main.cpp $SRC
fi
echo "Build OK. Run tests with: bin/gs_test"
