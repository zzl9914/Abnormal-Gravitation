@echo off
chcp 65001 >nul
cd /d "%~dp0"
rem 编译 wasm 内核并生成 src/web/gravity_wasm.js（base64 内联）
node src\web\tools\build-kernel.js
if errorlevel 1 exit /b 1
echo Build OK: src\web\gravity_wasm.js
