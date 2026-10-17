@echo off
chcp 65001 >nul
cd /d "%~dp0"
rem 启动静态服务器（前端直连 wasm 内核，无需 Python）
node web\tools\serve.js
if errorlevel 1 pause
