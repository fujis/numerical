@echo off
rem ---------------------------------------------------------------------
rem  Launcher for the numerical-methods samples (Python version)
rem  Note: keep this file ASCII-only. Japanese messages are printed by
rem        numerical_menu.py so that cmd.exe never has to parse them.
rem ---------------------------------------------------------------------
chcp 65001 >nul
setlocal
set "PYTHONIOENCODING=utf-8"
set "PYTHONUTF8=1"
cd /d "%~dp0"

rem Python command (change here if "python" is not on your PATH)
set "PY=python"

%PY% "%~dp0glviewer_mpl.py"
echo.
pause

endlocal
