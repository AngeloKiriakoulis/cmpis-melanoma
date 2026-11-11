@echo off
REM Stage 9 (redraw network figure) + Stage 11 (CMPIS clinical panels). Double-click.
cd /d "%~dp0"
set "RS="
where Rscript >nul 2>nul && set "RS=Rscript"
if not defined RS for /d %%D in ("C:\Program Files\R\R-*") do if exist "%%D\bin\Rscript.exe" set "RS=%%D\bin\Rscript.exe"
if not defined RS ( echo Could not find Rscript.exe. & pause & exit /b 1 )
for %%S in (09_Mechanistic_Network 11_CMPIS_Clinical) do (
  echo.
  echo ================= %%S =================
  "%RS%" -e "source('config.R'); source('R/%%S.R')"
)
echo.
echo ==================== FINISHED - review output above ====================
pause
