@echo off
rem Build gravity_core.dll, the physics test harness and the benchmark
rem (requires MSYS2 g++).
cd /d "%~dp0"
set CXX=g++
%CXX% -O2 -std=c++17 -shared -static -static-libgcc -static-libstdc++ -o bin\gravity_core.dll core\api.cpp core\expr.cpp core\sim.cpp
if errorlevel 1 exit /b 1
%CXX% -O2 -std=c++17 -static-libgcc -static-libstdc++ -o bin\gs_test.exe core\test_main.cpp core\expr.cpp core\sim.cpp
if errorlevel 1 exit /b 1
%CXX% -O2 -std=c++17 -static-libgcc -static-libstdc++ -o bin\gs_bench.exe core\bench_main.cpp core\expr.cpp core\sim.cpp
if errorlevel 1 exit /b 1
echo Build OK: bin\gravity_core.dll, bin\gs_test.exe, bin\gs_bench.exe
