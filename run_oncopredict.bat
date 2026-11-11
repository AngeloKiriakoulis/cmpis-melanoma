@echo off
REM Stage 12: re-predict a 12-drug panel with the oncoPredict/pRRophetic ridge pipeline and compare with Stage 11. ~15-30 min. Double-click.
cd /d "%~dp0"
set "RS="
where Rscript >nul 2>nul && set "RS=Rscript"
if not defined RS for /d %%D in ("C:\Program Files\R\R-*") do if exist "%%D\bin\Rscript.exe" set "RS=%%D\bin\Rscript.exe"
if not defined RS ( echo Could not find Rscript.exe. & pause & exit /b 1 )
echo ================= 12_Drug_oncoPredict_Check  (this takes a while) =================
"%RS%" -e "source('config.R'); source('R/12_Drug_oncoPredict_Check.R')"
echo.
echo ==================== FINISHED - review output above ====================
pause
