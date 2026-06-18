@echo off
setlocal enabledelayedexpansion

REM ============================================================
REM plateypus release script
REM ============================================================

REM CHANGE THESE EACH RELEASE
set NEW_VERSION=0.6.12
set CONDA_ENV_NAME=plateypus-py310
set COMMIT_MESSAGE=Release %NEW_VERSION%

REM Usually you do not need to change this if this .bat file is in repo root
set REPO_ROOT=%~dp0

REM Inno Setup compiler path. Adjust if yours is installed somewhere else.
set INNO_COMPILER=C:\Program Files (x86)\Inno Setup 6\ISCC.exe

echo.
echo ============================================================
echo plateypus release %NEW_VERSION%
echo Repo root: %REPO_ROOT%
echo Conda env: %CONDA_ENV_NAME%
echo ============================================================
echo.

cd /d "%REPO_ROOT%"

if errorlevel 1 (
    echo ERROR: Could not cd to repo root.
    exit /b 1
)

REM ============================================================
REM 1. Activate conda environment
REM ============================================================

echo.
echo ============================================================
echo Activating conda environment
echo ============================================================
call conda activate %CONDA_ENV_NAME%

if errorlevel 1 (
    echo ERROR: Could not activate conda environment: %CONDA_ENV_NAME%
    echo Try running this from Anaconda Prompt.
    exit /b 1
)

echo.
echo Python path:
where python
python --version

REM ============================================================
REM 2. Update version in setup.py and Inno script
REM ============================================================

echo.
echo ============================================================
echo Updating version numbers
echo ============================================================

python -c "from pathlib import Path; import re; p=Path('setup.py'); s=p.read_text(); s=re.sub(r'version\s*=\s*\"[^\"]+\"', 'version=\"%NEW_VERSION%\"', s); p.write_text(s)"

if errorlevel 1 (
    echo ERROR: Failed to update setup.py version.
    exit /b 1
)

python -c "from pathlib import Path; import re; p=Path('installer/plateypus_inno_setup.iss'); s=p.read_text(); s=re.sub(r'#define MyAppVersion \"[^\"]+\"', '#define MyAppVersion \"%NEW_VERSION%\"', s); p.write_text(s)"

if errorlevel 1 (
    echo ERROR: Failed to update Inno Setup version.
    exit /b 1
)

REM ============================================================
REM 3. Install local editable package and test imports
REM ============================================================

echo.
echo ============================================================
echo Installing local editable plateypus
echo ============================================================

python -m pip install -e .

if errorlevel 1 (
    echo ERROR: Editable install failed.
    exit /b 1
)

echo.
echo ============================================================
echo Testing plateypus import
echo ============================================================

python -c "import plateypus; print('plateypus imported from:', plateypus.__file__)"

if errorlevel 1 (
    echo ERROR: plateypus import failed.
    exit /b 1
)

python -c "from plateypus.__main__ import main; print('main function found:', main)"

if errorlevel 1 (
    echo ERROR: Could not import plateypus.__main__.main.
    exit /b 1
)

REM ============================================================
REM 4. Commit and push to GitHub
REM ============================================================

echo.
echo ============================================================
echo Git status
echo ============================================================

git status

echo.
echo ============================================================
echo Adding, committing, and pushing
echo ============================================================

git add .

git commit -m "%COMMIT_MESSAGE%"

if errorlevel 1 (
    echo.
    echo WARNING: git commit failed.
    echo This may simply mean there were no changes to commit.
    echo Continuing to push/tag/build may still be OK, but check the output above.
    echo.
)

git push origin main

if errorlevel 1 (
    echo ERROR: git push failed.
    exit /b 1
)

REM ============================================================
REM 5. Create and push Git tag
REM ============================================================

echo.
echo ============================================================
echo Creating Git tag v%NEW_VERSION%
echo ============================================================

git tag v%NEW_VERSION%

if errorlevel 1 (
    echo.
    echo WARNING: Could not create tag v%NEW_VERSION%.
    echo It may already exist.
    echo.
)

git push origin v%NEW_VERSION%

if errorlevel 1 (
    echo.
    echo WARNING: Could not push tag v%NEW_VERSION%.
    echo It may already exist remotely.
    echo.
)

REM ============================================================
REM 6. Build and upload PyPI package
REM ============================================================

echo.
echo ============================================================
echo Installing/updating PyPI build tools
echo ============================================================

python -m pip install --upgrade build twine

if errorlevel 1 (
    echo ERROR: Could not install/update build and twine.
    exit /b 1
)

echo.
echo ============================================================
echo Cleaning old PyPI build outputs
echo ============================================================

if exist pypi_dist rmdir /s /q pypi_dist
if exist plateypus.egg-info rmdir /s /q plateypus.egg-info
if exist build rmdir /s /q build

echo.
echo ============================================================
echo Building PyPI package
echo ============================================================

python -m build --outdir pypi_dist

if errorlevel 1 (
    echo ERROR: PyPI build failed.
    exit /b 1
)

echo.
echo ============================================================
echo Checking PyPI package
echo ============================================================

python -m twine check pypi_dist\*

if errorlevel 1 (
    echo ERROR: twine check failed.
    exit /b 1
)

echo.
echo ============================================================
echo Uploading to PyPI
echo ============================================================

python -m twine upload pypi_dist\*

if errorlevel 1 (
    echo ERROR: PyPI upload failed.
    echo Possible causes:
    echo - Version %NEW_VERSION% already exists on PyPI
    echo - Missing/incorrect PyPI token
    echo - Network/authentication issue
    exit /b 1
)

REM ============================================================
REM 7. Build PyInstaller onedir executable
REM ============================================================

echo.
echo ============================================================
echo Installing/updating PyInstaller
echo ============================================================

python -m pip install --upgrade pyinstaller

if errorlevel 1 (
    echo ERROR: Could not install/update PyInstaller.
    exit /b 1
)

echo.
echo ============================================================
echo Cleaning old PyInstaller outputs
echo ============================================================

if exist pyinstaller_build rmdir /s /q pyinstaller_build
if exist pyinstaller_dist rmdir /s /q pyinstaller_dist

echo.
echo ============================================================
echo Building PyInstaller onedir executable
echo ============================================================

python -m PyInstaller ^
  --noconfirm ^
  --clean ^
  --onedir ^
  --windowed ^
  --name plateypus ^
  --icon "plateypus\plateypus.ico" ^
  --distpath pyinstaller_dist ^
  --workpath pyinstaller_build ^
  --collect-data plateypus ^
  --collect-submodules plateypus ^
  --add-data "plateypus\plateypusLogo.png;plateypus" ^
  --add-data "plateypus\plateypus.ico;plateypus" ^
  --add-data "plateypus\templates;plateypus\templates" ^
  packaging\plateypus_launcher.py

if errorlevel 1 (
    echo ERROR: PyInstaller build failed.
    exit /b 1
)

if not exist pyinstaller_dist\plateypus\plateypus.exe (
    echo ERROR: Expected executable was not created:
    echo pyinstaller_dist\plateypus\plateypus.exe
    exit /b 1
)

REM ============================================================
REM 8. Compile Inno Setup installer
REM ============================================================

echo.
echo ============================================================
echo Compiling Inno Setup installer
echo ============================================================

if not exist "%INNO_COMPILER%" (
    echo ERROR: Inno Setup compiler not found at:
    echo %INNO_COMPILER%
    echo.
    echo Edit INNO_COMPILER in this .bat file.
    exit /b 1
)

if exist installer_output rmdir /s /q installer_output

"%INNO_COMPILER%" "%REPO_ROOT%installer\plateypus_inno_setup.iss"

if errorlevel 1 (
    echo ERROR: Inno Setup compile failed.
    exit /b 1
)

REM ============================================================
REM 9. Done
REM ============================================================

echo.
echo ============================================================
echo Release complete.
echo.
echo PyPI package uploaded for version:
echo %NEW_VERSION%
echo.
echo Windows executable folder:
echo %REPO_ROOT%pyinstaller_dist\plateypus
echo.
echo Installer should be here:
echo %REPO_ROOT%installer_output\plateypus-installer-%NEW_VERSION%.exe
echo.
echo Next recommended test:
echo Install the generated installer in Windows Sandbox or another Windows user account.
echo ============================================================
echo.

endlocal
