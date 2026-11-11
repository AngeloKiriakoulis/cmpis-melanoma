@echo off
REM Stage 14: redraw the CMPIS weight and drug-effect figures from saved tables (under 1 min).
cd /d "%~dp0"
set "RS="
where Rscript >nul 2>nul && set "RS=Rscript"
if not defined RS for /d %%D in ("C:\Program Files\R\R-*") do if exist "%%D\bin\Rscript.exe" set "RS=%%D\bin\Rscript.exe"
if not defined RS ( echo Could not find Rscript.exe. & pause & exit /b 1 )
"%RS%" -e "source('config.R'); source('R/14_Manuscript_Figures.R')"
pause
