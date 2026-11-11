@echo off
REM ====================================================================
REM  Run the analysis stages 6-15 in dependency order (after Stages 1-5;
REM  see README). Logs: logs\stage_NN.txt. Takes about 1 hour.
REM ====================================================================
cd /d "%~dp0"
set "RS="
where Rscript >nul 2>nul && set "RS=Rscript"
if not defined RS for /d %%D in ("C:\Program Files\R\R-*") do if exist "%%D\bin\Rscript.exe" set "RS=%%D\bin\Rscript.exe"
if not defined RS ( echo Could not find Rscript.exe. & pause & exit /b 1 )
if not exist logs mkdir logs
for %%S in (10_CMPIS_Reweight 06_PCD_Sensitivity 07_Classifier_Benchmark 08_ICB_Validation 11_CMPIS_Clinical 09_Mechanistic_Network 12_Drug_oncoPredict_Check 13_SingleCell_Validation 14_Manuscript_Figures 15_ICB_cBioPortal) do (
  echo [%time%] Stage %%S ...
  "%RS%" -e "source('config.R'); source('R/%%S.R')" > "logs\stage_%%S.txt" 2>&1
  if errorlevel 1 ( echo    FAILED - see logs\stage_%%S.txt ) else ( echo    ok )
)
echo.
echo ==================== FINISHED - tell Claude it is done ====================
pause
