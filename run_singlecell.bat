@echo off
REM Stage 13: single-cell validation (Fig 7) - downloads Tirosh, Jerby-Arnon and Sade-Feldman from GEO into data\sc (about 1 GB; files already there are reused), then runs the analysis. 20-60 min. Double-click.
cd /d "%~dp0"
set "RS="
where Rscript >nul 2>nul && set "RS=Rscript"
if not defined RS for /d %%D in ("C:\Program Files\R\R-*") do if exist "%%D\bin\Rscript.exe" set "RS=%%D\bin\Rscript.exe"
if not defined RS ( echo Could not find Rscript.exe. & pause & exit /b 1 )
echo ================= 13_SingleCell_Validation  (download + analysis; this takes a while) =================
"%RS%" -e "source('config.R'); source('R/13_SingleCell_Validation.R')"
echo.
echo ==================== FINISHED - review output above ====================
pause
