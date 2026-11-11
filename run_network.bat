@echo off
REM Run the mechanistic network figure stage, directly.
cd /d "%~dp0"
set "RS="
where Rscript >nul 2>nul && set "RS=Rscript"
if not defined RS for /d %%D in ("C:\Program Files\R\R-*") do if exist "%%D\bin\Rscript.exe" set "RS=%%D\bin\Rscript.exe"
if not defined RS ( echo Could not find Rscript.exe. & pause & exit /b 1 )
echo Using Rscript: "%RS%"
echo.
"%RS%" -e "source('config.R'); source('R/09_Mechanistic_Network.R')"
echo.
echo ==================== FINISHED - review results\ above ====================
pause
