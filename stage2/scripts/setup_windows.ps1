# Stage 2 — Windows local setup and optional execution
# Run from PowerShell in the cloned/extracted repository:
#   powershell -ExecutionPolicy Bypass -File .\stage2\scripts\setup_windows.ps1
#   powershell -ExecutionPolicy Bypass -File .\stage2\scripts\setup_windows.ps1 -DataDir "D:\path\to\GEO files" -RunPreflight
#   powershell -ExecutionPolicy Bypass -File .\stage2\scripts\setup_windows.ps1 -RunAnalysis
# No automatic network downloads or R package installations.
param(
    [string]$DataDir = "",
    [string]$RscriptPath = "",
    [switch]$RunPreflight,
    [switch]$RunAnalysis,
    [switch]$SkipPathways
)
$ErrorActionPreference = "Stop"
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$dataTarget = Join-Path $repoRoot "stage2\data"
$resultTarget = Join-Path $repoRoot "stage2\results"
New-Item -ItemType Directory -Force -Path $dataTarget, $resultTarget | Out-Null
$names = @("GSE256068_processed_data.csv.gz", "GSE256068_series_matrix.txt.gz")
Write-Host "Stage 2 repository: $repoRoot"
Write-Host "Data folder: $dataTarget"
if ($DataDir -ne "") {
    if (-not (Test-Path -LiteralPath $DataDir -PathType Container)) {
        throw "DataDir does not exist: $DataDir"
    }
    foreach ($name in $names) {
        $source = Join-Path $DataDir $name
        $target = Join-Path $dataTarget $name
        if (-not (Test-Path -LiteralPath $source -PathType Leaf)) {
            throw "Required GEO file missing: $source"
        }
        if ((Resolve-Path -LiteralPath $source).Path -ne $target) {
            Copy-Item -LiteralPath $source -Destination $target -Force
            Write-Host "Copied: $name"
        }
    }
}
$missing = @($names | Where-Object { -not (Test-Path -LiteralPath (Join-Path $dataTarget $_) -PathType Leaf) })
if ($missing.Count -gt 0) {
    Write-Host "Missing GEO inputs: $($missing -join ', ')" -ForegroundColor Yellow
    Write-Host "Place the original .gz files in stage2\data or rerun with -DataDir 'path'."
    if ($RunPreflight -or $RunAnalysis) { throw "Cannot run without the two original GEO files." }
    exit 0
}
foreach ($name in $names) {
    $file = Get-Item -LiteralPath (Join-Path $dataTarget $name)
    Write-Host ("Found: {0} ({1:N1} MB)" -f $name, ($file.Length / 1MB))
}
if (-not ($RunPreflight -or $RunAnalysis)) {
    Write-Host "Local folders are ready. Add -RunPreflight or -RunAnalysis to continue."
    exit 0
}
if ($RscriptPath -ne "") {
    if (-not (Test-Path -LiteralPath $RscriptPath -PathType Leaf)) { throw "RscriptPath not found: $RscriptPath" }
    $rscript = $RscriptPath
} else {
    $cmd = Get-Command Rscript.exe -ErrorAction SilentlyContinue
    if ($cmd) { $rscript = $cmd.Source }
    else {
        $candidates = @(Get-ChildItem -Path "C:\Program Files\R" -Filter "Rscript.exe" -Recurse -ErrorAction SilentlyContinue |
            Sort-Object FullName -Descending)
        if ($candidates.Count -eq 0) { throw "Rscript.exe not found. Install R, or pass -RscriptPath 'C:\...\Rscript.exe'." }
        $rscript = $candidates[0].FullName
    }
}
Write-Host "Using Rscript: $rscript"
Push-Location $repoRoot
try {
    & $rscript "stage2/scripts/00_preflight.R"
    if ($LASTEXITCODE -ne 0) { throw "Preflight failed (exit $LASTEXITCODE)" }
    Write-Host "Review stage2\results\preflight_report.txt before interpreting any results."
    $preflightText = Get-Content -LiteralPath (Join-Path $resultTarget "preflight_report.txt") -Raw
    if ($preflightText -notmatch "TLE-HS selected:\s*64" -or
        $preflightText -notmatch "Hippocampal controls selected:\s*13" -or
        $preflightText -notmatch "Selected expression IDs matched:\s*77\s*/\s*77") {
        throw "Preflight did not confirm 64 TLE-HS, 13 controls, and 77/77 matched sample IDs. Inspect preflight_report.txt; analysis aborted."
    }
    if ($RunAnalysis) {
        $arguments = @("stage2/scripts/01_differential_expression.R")
        if ($SkipPathways) { $arguments += "--skip-pathways" }
        & $rscript @arguments
        if ($LASTEXITCODE -ne 0) { throw "Differential expression failed (exit $LASTEXITCODE)" }
        & $rscript "stage2/scripts/02_review_results.R"
        if ($LASTEXITCODE -ne 0) { throw "Statistical review failed (exit $LASTEXITCODE)" }
        Write-Host "Stage 2 output: stage2\results\review_summary.txt" -ForegroundColor Green
    }
} finally { Pop-Location }
