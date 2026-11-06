@echo off
setlocal
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
  py -3 -m venv .venv
  if errorlevel 1 goto error
)
.venv\Scripts\python.exe -m pip install --disable-pip-version-check -r requirements.txt
if errorlevel 1 goto error
.venv\Scripts\python.exe -m streamlit run app.py
if errorlevel 1 goto error
exit /b 0
:error
echo See the error above and the installation instructions in README.md.
pause
exit /b 1
