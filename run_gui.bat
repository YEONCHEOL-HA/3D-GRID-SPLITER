@echo off
REM 빌드 없이 GUI 를 바로 실행합니다. Python 이 없으면 자동으로 설치합니다.
chcp 65001 >nul
cd /d "%~dp0"

if not exist "%~dp0_python_setup.bat" goto :NOSETUP
call "%~dp0_python_setup.bat"
if errorlevel 1 goto :FAIL
if not defined PY goto :FAIL

echo [*] 필요한 패키지 확인 중...
%PY% -m pip install --quiet --upgrade pip
%PY% -m pip install --quiet numpy matplotlib
%PY% gui.py
if errorlevel 1 pause
exit /b 0

:NOSETUP
echo [!] _python_setup.bat 이 같은 폴더에 없습니다.
pause
exit /b 1

:FAIL
pause
exit /b 1
