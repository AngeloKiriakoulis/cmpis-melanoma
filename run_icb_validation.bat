@echo off
REM ====================================================================
REM  Run the anti-PD-1 validation stage, directly.
REM  Downloads GSE78220 + GSE91061 via GEOquery on FIRST run (needs
REM  internet), then analyses. Writes to results\. Re-runs reuse the
REM  GEO cache. To add a cohort manually, drop <cohort>_expr.txt and
REM  <cohort>_meta.txt in data\icb\ (see R\08_ICB_Validation.R header).
REM ====================================================================
cd /d "%~dp0"
set "RS="
where Rscript >nul 2>nul && set "RS=Rscript"
if not defined RS for /d %%D in ("C:\Program Files\R\R-*") do if exist "%%D\bin\Rscript.exe" set "RS=%%D\bin\Rscript.exe"
if not defined RS ( echo Could not find Rscript.exe. & pause & exit /b 1 )
echo Using Rscript: "%RS%"
echo.
"%RS%" -e "source('config.R'); source('R/08_ICB_Validation.R')"
echo.
echo ==================== FINISHED - review results\ above ====================
pause
