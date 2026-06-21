#!/usr/bin/env python3
"""
FEBio Foldable Display - Runner & Post-processor

1. generate .feb file
2. run FEBio solver
3. convert results (xplt -> hdf5)
4. visualize

Usage:
  python scripts/run_folding.py                    # .feb generation only
  python scripts/run_folding.py --solve            # generate + solve
  python scripts/run_folding.py --solve --view     # generate + solve + visualize
"""

import os
import sys
import subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FEB_FILE = os.path.join(ROOT, "display_folding.feb")
XPLT_FILE = os.path.join(ROOT, "display_folding.xplt")
HDF5_FILE = os.path.join(ROOT, "display_folding_results.h5")


def find_febio():
    """Locate febio4.exe."""
    candidates = [
        os.path.join(os.path.expanduser("~"), ".febio", "bin", "febio4.exe"),
        os.path.join(os.environ.get("ProgramFiles", "C:/Program Files"), "FEBio Studio", "bin", "febio4.exe"),
        os.path.join(os.environ.get("ProgramFiles", "C:/Program Files"), "FEBio", "bin", "febio4.exe"),
        os.path.join(os.environ.get("ProgramFiles(x86)", "C:/Program Files (x86)"), "FEBio", "bin", "febio4.exe"),
        "C:/FEBio/febio4.exe",
    ]
    for p in os.environ.get("PATH", "").split(os.pathsep):
        candidates.append(os.path.join(p, "febio4.exe"))
    for cand in candidates:
        if os.path.isfile(cand):
            return os.path.normpath(cand)
    return None


def run_solver(feb_path, febio_path=None):
    if febio_path is None:
        febio_path = find_febio()
    if febio_path is None:
        print("[FAIL] febio4.exe not found.")
        print("  Run: powershell -ExecutionPolicy Bypass -File scripts/install_febio.ps1")
        return False
    if not os.path.isfile(feb_path):
        print(f"[FAIL] .feb not found: {feb_path}")
        return False

    print(f"[Solver] {febio_path}")
    print(f"[Input]  {feb_path}")
    print("[Run] Solving... (may take minutes)")
    try:
        result = subprocess.run(
            [febio_path, feb_path],
            cwd=os.path.dirname(feb_path),
            capture_output=True, text=True, timeout=3600,
        )
        print(f"[Done] Exit code: {result.returncode}")
        return result.returncode == 0
    except subprocess.TimeoutExpired:
        print("[FAIL] Timeout (1h)")
        return False
    except Exception as e:
        print(f"[FAIL] {e}")
        return False


def convert_xplt(xplt_path, hdf5_path):
    try:
        from pyfebio import xplt
        print(f"[XPLT] Converting: {xplt_path}")
        xplt.to_hdf5(xplt_path, hdf5_path)
        print(f"[XPLT] Saved: {hdf5_path}")
        return True
    except ImportError:
        print("[XPLT] pip install pyfebio")
        return False
    except Exception as e:
        print(f"[XPLT] Failed: {e}")
        return False


def generate_feb():
    ret = os.system(f'"{sys.executable}" -m src --solver febio')
    return ret == 0


def main():
    args = set(sys.argv[1:])
    do_solve = "--solve" in args
    do_view = "--view" in args

    if not generate_feb():
        print("[FAIL] .feb generation failed")
        sys.exit(1)

    if do_solve:
        if not run_solver(FEB_FILE):
            print("[SKIP] Solve failed")
            return

    if do_solve and os.path.exists(XPLT_FILE):
        convert_xplt(XPLT_FILE, HDF5_FILE)
        if do_view:
            print(f"[View] Open: {HDF5_FILE}")
            print("  Run: python scripts/view_results.py")

    print("[Done]")


if __name__ == "__main__":
    main()
