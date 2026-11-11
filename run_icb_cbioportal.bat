@echo off
REM ====================================================================
REM  Stage 15: real immunotherapy validation (Liu/Van Allen, Riaz, Hugo,
REM  Gide, Van Allen 2015) from cBioPortal. Needs internet on first run;
REM  data cached in data\icb\cbio_*.rds. Log: results\stage15_log.txt
REM ====================================================================
cd /d "%~dp0"
set "RS="
where Rscript >nul 2>nul && set "RS=Rscript"
if not defined RS for /d %%D in ("C:\Program Files\R\R-*") do if exist "%%D\bin\Rscript.exe" set "RS=%%D\bin\Rscript.exe"
if not defined RS ( echo Could not find Rscript.exe. & pause & exit /b 1 )
echo Using Rscript: "%RS%"
echo Running Stage 15 (downloads ~5 cohorts, a few minutes)...
"%RS%" -e "source('config.R'); source('R/15_ICB_cBioPortal.R')" > results\stage15_log.txt 2>&1
type results\stage15_log.txt
echo.
echo ==================== FINISHED - tell Claude it is done ====================
pause
