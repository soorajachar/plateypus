conda activate plateypus-py310
cd C:\Users\User\Documents\clonedRepos\plateypus

REM ============================================================
REM 1. After editing code, test local editable package
REM ============================================================
python -m pip install -e .
plateypus

REM ============================================================
REM 2. Commit and push source to GitHub
REM ============================================================
git status
git add .
git commit -m "Release 0.6.12"
git push origin main

git tag v0.6.12
git push origin v0.6.12

REM ============================================================
REM 3. Build and upload PyPI package
REM ============================================================
python -m pip install --upgrade build twine

if exist pypi_dist rmdir /s /q pypi_dist
if exist plateypus.egg-info rmdir /s /q plateypus.egg-info

python -m build --outdir pypi_dist
python -m twine check pypi_dist\*
python -m twine upload pypi_dist\*

REM ============================================================
REM 4. Build PyInstaller onedir executable
REM ============================================================
python -m pip install --upgrade pyinstaller

if exist pyinstaller_build rmdir /s /q pyinstaller_build
if exist pyinstaller_dist rmdir /s /q pyinstaller_dist

python -m PyInstaller ^
  --noconfirm ^
  --clean ^
  --onedir ^
  --windowed ^
  --name plateypus ^
  --icon "plateypus\plateypus.ico" ^
  --distpath pyinstaller_dist ^
  --workpath pyinstaller_build ^
  --specpath packaging ^
  --collect-data plateypus ^
  --collect-submodules plateypus ^
  --add-data "plateypus\plateypusLogo.png;plateypus" ^
  --add-data "plateypus\plateypus.ico;plateypus" ^
  --add-data "plateypus\templates;plateypus\templates" ^
  packaging\plateypus_launcher.py

REM ============================================================
REM 5. Test executable
REM ============================================================
pyinstaller_dist\plateypus\plateypus.exe

REM ============================================================
REM 6. Compile installer
REM ============================================================
"C:\Program Files (x86)\Inno Setup 6\ISCC.exe" "C:\path\to\plateypus\installer\plateypus_inno_setup.iss"
