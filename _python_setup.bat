@echo off
REM ===================================================================
REM  _python_setup.bat
REM  Python 을 찾고, 없으면 자동으로 설치한다.
REM  호출한 배치에 PY 변수를 남긴다.  사용 예)  %PY% -m pip install ...
REM  자동 설치를 원하지 않으면 호출 전에  set PYSETUP_NOINSTALL=1
REM ===================================================================
set "PYVER=3.13.15"

call :FINDPY
if defined PY exit /b 0
if defined PYSETUP_NOINSTALL goto :GIVEUP

echo.
echo [*] Python 이 설치돼 있지 않습니다. 지금 자동으로 설치합니다.
echo     사용자 계정에만 설치하므로 관리자 권한은 필요 없습니다.
echo.

REM ---------- 1차 시도: winget ----------
where winget >nul 2>nul
if errorlevel 1 goto :DIRECT
echo [*] winget 으로 Python %PYVER% 설치 중...
winget install -e --id Python.Python.3.13 --scope user --silent --accept-source-agreements --accept-package-agreements
call :FINDPY
if defined PY goto :INSTALLED

REM ---------- 2차 시도: python.org 설치 파일 직접 내려받기 ----------
:DIRECT
set "SUF=-amd64"
set "ARCHX=%PROCESSOR_ARCHITECTURE%"
if defined PROCESSOR_ARCHITEW6432 set "ARCHX=%PROCESSOR_ARCHITEW6432%"
if /i "%ARCHX%"=="ARM64" set "SUF=-arm64"
if /i "%ARCHX%"=="x86" set "SUF="

set "URL=https://www.python.org/ftp/python/%PYVER%/python-%PYVER%%SUF%.exe"
set "DL=%TEMP%\python-%PYVER%%SUF%.exe"

echo [*] 설치 파일을 내려받는 중 (약 30MB)
echo     %URL%
if exist "%DL%" del /q "%DL%"
where curl.exe >nul 2>nul
if errorlevel 1 goto :USEPS
curl.exe -L --fail --silent --show-error -o "%DL%" "%URL%"
if exist "%DL%" goto :CHECKDL

:USEPS
echo [*] PowerShell 로 내려받는 중...
powershell -NoProfile -ExecutionPolicy Bypass -Command "$ProgressPreference='SilentlyContinue';[Net.ServicePointManager]::SecurityProtocol=[Net.SecurityProtocolType]::Tls12;Invoke-WebRequest -Uri '%URL%' -OutFile '%DL%' -UseBasicParsing"

:CHECKDL
if not exist "%DL%" goto :DLFAIL
set "TOOSMALL="
for %%A in ("%DL%") do if %%~zA LSS 5000000 set "TOOSMALL=1"
if defined TOOSMALL goto :DLBROKEN

echo [*] 설치 중... 2~3분 걸릴 수 있습니다. 창이 멈춘 것처럼 보여도 그대로 두세요.
start /wait "" "%DL%" /quiet InstallAllUsers=0 PrependPath=1 Include_pip=1 Include_launcher=1 Include_test=0
del /q "%DL%" 2>nul
call :FINDPY

:INSTALLED
if not defined PY goto :GIVEUP
echo [*] Python 준비 완료.
echo.
exit /b 0

:DLFAIL
echo [!] 설치 파일을 내려받지 못했습니다. 인터넷 연결이나 회사 방화벽을 확인하세요.
goto :GIVEUP

:DLBROKEN
echo [!] 내려받은 파일이 손상된 것 같습니다. 잠시 후 다시 실행해 보세요.
del /q "%DL%" 2>nul
goto :GIVEUP


REM ---------------- Python 찾기 ----------------
REM  Microsoft Store 의 가짜 python.exe(앱 실행 별칭)는 걸러낸다.
:FINDPY
set "PY="
call py -3 -c "import sys" >nul 2>nul
if not errorlevel 1 set "PY=py -3"
if defined PY goto :EOF

set "CAND="
for /f "delims=" %%p in ('python -c "import sys;print(sys.executable)" 2^>nul') do set "CAND=%%p"
if not defined CAND goto :FINDDIR
REM  경로에 WindowsApps 가 들어 있으면 Store 별칭이므로 쓰지 않는다
if "%CAND%"=="%CAND:WindowsApps=%" set "PY=python"
if defined PY goto :EOF

:FINDDIR
for %%D in ("%LOCALAPPDATA%\Programs\Python" "%ProgramFiles%" "C:\Program Files" "C:") do if not defined PY call :CHECKDIR "%%~D"
goto :EOF

:CHECKDIR
for /f "delims=" %%P in ('dir /b /ad "%~1\Python3*" 2^>nul') do if exist "%~1\%%P\python.exe" set PY="%~1\%%P\python.exe"
goto :EOF


:GIVEUP
echo.
echo  [!] Python 을 자동으로 준비하지 못했습니다. 직접 설치해 주세요.
echo.
echo     https://www.python.org/downloads/windows/   에서 3.13 설치
echo     설치 첫 화면의 "Add python.exe to PATH" 를 반드시 체크하세요.
echo.
exit /b 1
