@echo off
rem Launch the gravity simulation server (managed Python venv preferred).
cd /d "%~dp0server"
set PY="C:\Users\flyin\.workbuddy\binaries\python\envs\default\Scripts\python.exe"
if not exist %PY% set PY=python
%PY% server.py --port 8123
