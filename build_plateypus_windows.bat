@echo off
setlocal enabledelayedexpansion

REM ============================================================
REM Build Windows executable for plateypus using the active conda env
REM ============================================================

REM Change this to your actual conda environment name
set CONDA_ENV_NAME=plateypus-py310

REM Repo root is assumed to be the directory containing this .bat file
set REPO_ROOT=%~dp0

echo.
echo ============================================================
echo Activating conda environment: %CONDA_ENV_NAME%
echo ============================================================
call conda activate %CONDA_ENV_NAME%

if errorlevel 1 (
    echo.
    echo ERROR: Could not activate conda environment: %CONDA_ENV_NAME%
    echo Make sure this is being run from Anaconda Prompt or Miniconda Prompt.
    exit /b 1
)

echo.
echo ============================================================
echo Moving to repo root
echo ============================================================
cd /d "%REPO_ROOT%"

echo.
echo ============================================================
echo Python being used:
echo ============================================================
where python
python --version

echo.
echo ============================================================
echo Installing local plateypus in editable mode
echo ============================================================
python -m pip install -e .

if errorlevel 1 (
    echo.
    echo ERROR: Failed to install plateypus in editable mode.
    exit /b 1
)

echo.
echo ============================================================
echo Installing/updating PyInstaller
echo ============================================================
python -m pip install --upgrade pyinstaller

if errorlevel 1 (
    echo.
    echo ERROR: Failed to install PyInstaller.
    exit /b 1
)

echo.
echo ============================================================
echo Checking PyInstaller version
echo ============================================================
python -m PyInstaller --version

echo.
echo ============================================================
echo Cleaning old build artifacts
echo ============================================================
if exist build rmdir /s /q build
if exist dist rmdir /s /q dist
if exist plateypus.spec del /q plateypus.spec

echo.
echo ============================================================
echo Ensuring launcher exists
echo ============================================================
if not exist packaging mkdir packaging

(
echo from plateypus.__main__ import main
echo.
echo.
echo if __name__ == "__main__":
echo     main()
) > packaging\plateypus_launcher.py

echo.
echo ============================================================
echo Building plateypus executable folder
echo ============================================================
python -m PyInstaller ^
  --noconfirm ^
  --clean ^
  --onedir ^
  --windowed ^
  --name plateypus ^
  --collect-data plateypus ^
  --collect-submodules plateypus ^
  --add-data "plateypus\plateypusLogo.png;plateypus" ^
  --add-data "plateypus\templates;plateypus\templates" ^
  packaging\plateypus_launcher.py

if errorlevel 1 (
    echo.
    echo ERROR: PyInstaller build failed.
    exit /b 1
)

echo.
echo ============================================================
echo Build complete.
echo Executable is here:
echo %REPO_ROOT%dist\plateypus\plateypus.exe
echo.
echo Test it by running:
echo dist\plateypus\plateypus.exe
echo ============================================================

endlocal
