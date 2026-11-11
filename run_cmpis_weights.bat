@echo off
REM Rebuild CMPIS weights (Stage 10), then re-run the stages that use them (7, 8, 9).
cd /d "%~dp0"
set "RS="
where Rscript >nul 2>nul && set "RS=Rscript"
if not defined RS for /d %%D in ("C:\Program Files\R\R-*") do if exist "%%D\bin\Rscript.exe" set "RS=%%D\bin\Rscript.exe"
if not defined RS ( echo Could not find Rscript.exe. & pause & exit /b 1 )
for %%S in (10_CMPIS_Reweight 07_Classifier_Benchmark 08_ICB_Validation 09_Mechanistic_Network) do (
  echo.
  echo ================= %%S =================
  "%RS%" -e "source('config.R'); source('R/%%S.R')"
)
echo.
echo ==================== FINISHED - review output above ====================
pause
