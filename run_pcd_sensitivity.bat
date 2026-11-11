@echo off
REM ====================================================================
REM  Run the PCD gene-set sensitivity stage.
REM  Runs the stage script DIRECTLY (source config + script) so it does
REM  not require the {targets}/callr machinery to be installed yet.
REM  Reuses existing caches; writes to results\.
REM ====================================================================
cd /d "%~dp0"

set "RS="
where Rscript >nul 2>nul && set "RS=Rscript"
if not defined RS for /d %%D in ("C:\Program Files\R\R-*") do if exist "%%D\bin\Rscript.exe" set "RS=%%D\bin\Rscript.exe"
if not defined RS for /d %%D in ("C:\Program Files\Microsoft\R Open\R-*") do if exist "%%D\bin\Rscript.exe" set "RS=%%D\bin\Rscript.exe"
if not defined RS (
  echo Could not find Rscript.exe. Edit this file and set RS to its full path.
  pause & exit /b 1
)
echo Using Rscript: "%RS%"
echo.
"%RS%" -e "source('config.R'); source('R/06_PCD_Sensitivity.R')"
echo.
echo ==================== FINISHED - review results\ above ====================
pause
