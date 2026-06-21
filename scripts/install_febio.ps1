# install_febio.ps1 - Automatic FEBio Solver installer
# Downloads febio4-Windows-X64 from GitHub Actions CI artifacts (46MB)
#
# Usage:
#   powershell -ExecutionPolicy Bypass -File install_febio.ps1

$ErrorActionPreference = "Stop"

$REPO = "febiosoftware/FEBio"
$RUN_ID = "27706397336"
$ARTIFACT = "febio4-Windows-X64"
$INSTALL_DIR = "$env:USERPROFILE\.febio\bin"
$TARGET = "$INSTALL_DIR\febio4.exe"

Write-Host "=== FEBio Solver Auto Install ===" -ForegroundColor Cyan
Write-Host "Target: $TARGET" -ForegroundColor Cyan
Write-Host ""

# 1. Check/install gh CLI
$gh = Get-Command gh -ErrorAction SilentlyContinue
if (-not $gh) {
    Write-Host "[1/4] Installing gh (GitHub CLI)..." -ForegroundColor Yellow
    winget install GitHub.cli --silent --accept-package-agreements 2>&1 | Out-Null
    $env:Path = [Environment]::GetEnvironmentVariable("Path", "User") + ";$env:Path"
    $gh = Get-Command gh -ErrorAction SilentlyContinue
    if (-not $gh) {
        Write-Host "[1/4] FAILED: gh install failed. Install manually from https://cli.github.com/" -ForegroundColor Red
        exit 1
    }
    Write-Host "[1/4] gh CLI installed!" -ForegroundColor Green
} else {
    Write-Host "[1/4] gh CLI already installed: $($gh.Source)" -ForegroundColor Green
}

# 2. Check authentication
Write-Host ""
$authCheck = gh auth status 2>&1
if ($LASTEXITCODE -ne 0) {
    Write-Host "[2/4] GitHub authentication required." -ForegroundColor Yellow
    Write-Host "  A browser window will open for GitHub login." -ForegroundColor Yellow
    Write-Host "  Press Enter to continue..." -ForegroundColor Yellow
    Read-Host
    gh auth login --web
    if ($LASTEXITCODE -ne 0) {
        Write-Host "[2/4] Auth failed. Run 'gh auth login' manually." -ForegroundColor Red
        exit 1
    }
    Write-Host "[2/4] GitHub authentication OK!" -ForegroundColor Green
} else {
    Write-Host "[2/4] GitHub already authenticated" -ForegroundColor Green
}

# 3. Download artifact
Write-Host ""
$DOWNLOAD_DIR = "$env:TEMP\febio_download"
if (Test-Path $DOWNLOAD_DIR) { Remove-Item -Recurse -Force $DOWNLOAD_DIR }
New-Item -ItemType Directory -Path $DOWNLOAD_DIR -Force | Out-Null

Write-Host "[3/4] Downloading febio4 artifact (46MB)..." -ForegroundColor Yellow
gh run download $RUN_ID --repo $REPO --name $ARTIFACT --dir $DOWNLOAD_DIR
if ($LASTEXITCODE -ne 0) {
    Write-Host "[3/4] Download failed. Trying latest runs..." -ForegroundColor Yellow
    gh run download --repo $REPO --name $ARTIFACT --dir $DOWNLOAD_DIR -L 1
    if ($LASTEXITCODE -ne 0) {
        Write-Host "[3/4] Download failed. Check run ID or network." -ForegroundColor Red
        exit 1
    }
}
Write-Host "[3/4] Download complete!" -ForegroundColor Green

# 4. Install (extract febio4.exe)
Write-Host ""
Write-Host "[4/4] Installing..." -ForegroundColor Yellow
if (-not (Test-Path $INSTALL_DIR)) {
    New-Item -ItemType Directory -Path $INSTALL_DIR -Force | Out-Null
}

# Unzip or copy files
if (Test-Path "$DOWNLOAD_DIR\$ARTIFACT.zip") {
    Expand-Archive -Path "$DOWNLOAD_DIR\$ARTIFACT.zip" -DestinationPath "$env:TEMP\febio_extract" -Force
    $febioExe = Get-ChildItem -Path "$env:TEMP\febio_extract" -Recurse -Filter "febio4.exe" | Select-Object -First 1
    if ($febioExe) {
        Copy-Item $febioExe.FullName $TARGET -Force
    }
    Remove-Item -Recurse -Force "$env:TEMP\febio_extract" -ErrorAction SilentlyContinue
} else {
    Copy-Item "$DOWNLOAD_DIR\*" $INSTALL_DIR -Recurse -Force
}

# 5. Verify
if (Test-Path $TARGET) {
    $ver = & $TARGET --version 2>&1 | Select-Object -First 1
    Write-Host ""
    Write-Host "=== INSTALL COMPLETE ===" -ForegroundColor Green
    Write-Host "  Path: $TARGET" -ForegroundColor Green
    Write-Host "  Version: $ver" -ForegroundColor Green
    Write-Host ""
    Write-Host "Run solver:" -ForegroundColor Cyan
    Write-Host "  python run_folding.py --solve" -ForegroundColor White
} else {
    Write-Host "ERROR: febio4.exe not found after extraction." -ForegroundColor Red
    exit 1
}
