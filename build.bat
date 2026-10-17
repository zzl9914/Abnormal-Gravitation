@echo off
chcp 65001 >nul
cd /d "%~dp0"
set CXX=g++
set SRC=src\core\expr.cpp src\core\sim.cpp src\core\gsmath.cpp
rem Simulation 含固定容量数组成员（约 1MB），测试里是栈对象，必须加大栈
set STACK=-Wl,--stack,67108864
%CXX% -O2 -std=c++17 -shared -static -static-libgcc -static-libstdc++ -o bin\gravity_core.dll src\core\api.cpp %SRC%
if errorlevel 1 exit /b 1
%CXX% -O2 -std=c++17 -static-libgcc -static-libstdc++ %STACK% -o bin\gs_test.exe src\core\test_main.cpp %SRC%
if errorlevel 1 exit /b 1
%CXX% -O2 -std=c++17 -static-libgcc -static-libstdc++ %STACK% -o bin\gs_bench.exe src\core\bench_main.cpp %SRC%
if errorlevel 1 exit /b 1
echo Build OK: bin\gravity_core.dll, bin\gs_test.exe, bin\gs_bench.exe
