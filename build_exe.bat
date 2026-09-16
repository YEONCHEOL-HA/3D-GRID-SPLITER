@echo off
REM ===================================================================
REM  3D GRID SPLITTER - Windows EXE 빌드 (Python 자동 설치 포함)
REM  gui.py / grid_splitter.py 와 같은 폴더에 두고 더블클릭하세요.
REM ===================================================================
chcp 65001 >nul
cd /d "%~dp0"

REM SLIM=1 : scipy 를 빼고 빌드 (용량 약 40MB 감소, 기능 동일)
set SLIM=1
REM ONEFILE=1 : exe 한 개로 배포 / 0 : 폴더 배포(실행이 빠름)
set ONEFILE=1
set APPNAME=3D-GRID-SPLITTER

echo ===================================================
echo    3D GRID SPLITTER  -  EXE 빌드
echo ===================================================
echo.

if not exist "%~dp0_python_setup.bat" goto :NOSETUP
call "%~dp0_python_setup.bat"
if errorlevel 1 goto :FAIL
if not defined PY goto :FAIL

echo [*] 사용할 Python: %PY%
%PY% -c "import sys;print('    version',sys.version.split()[0]);print('   ',sys.executable)"
echo.

echo [1/3] 필요한 패키지 설치 중...
%PY% -m pip install --upgrade pip
%PY% -m pip install numpy matplotlib pyinstaller
if errorlevel 1 (
  echo [!] 패키지 설치 실패. 인터넷 연결 또는 회사 프록시 설정을 확인하세요.
  pause
  exit /b 1
)
if "%SLIM%"=="0" %PY% -m pip install scipy
echo.

echo [2/3] EXE 빌드 중... 몇 분 걸릴 수 있습니다.
set EXTRA=
if "%SLIM%"=="1" set EXTRA=--exclude-module scipy
set MODE=--onedir
if "%ONEFILE%"=="1" set MODE=--onefile

%PY% -m PyInstaller --noconfirm --clean %MODE% --windowed ^
  --name "%APPNAME%" ^
  --collect-data matplotlib ^
  --hidden-import matplotlib.backends.backend_agg ^
  --hidden-import mpl_toolkits.mplot3d ^
  --exclude-module PyQt5 --exclude-module PyQt6 ^
  --exclude-module PySide2 --exclude-module PySide6 ^
  --exclude-module IPython --exclude-module notebook ^
  --exclude-module pandas --exclude-module pytest ^
  %EXTRA% gui.py

if errorlevel 1 (
  echo [!] 빌드 실패. 위 메시지를 확인하세요.
  pause
  exit /b 1
)

echo.
echo [3/3] 완료!
if "%ONEFILE%"=="1" (
  echo    dist\%APPNAME%.exe  를 실행하세요.
) else (
  echo    dist\%APPNAME%\%APPNAME%.exe  를 실행하세요. 폴더 전체를 함께 복사해야 합니다.
)
echo.
pause
exit /b 0

:NOSETUP
echo [!] _python_setup.bat 이 같은 폴더에 없습니다.
echo     압축을 푼 파일 전체를 한 폴더에 두고 다시 실행하세요.
pause
exit /b 1

:FAIL
pause
exit /b 1
