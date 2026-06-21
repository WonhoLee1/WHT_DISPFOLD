#!/usr/bin/env python3
"""
Batch generate all solver input files.

Usage:
  python scripts/run_all.py
  python scripts/run_all.py --gmsh   # force GMSH mesh
"""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PYTHON = sys.executable
SRC = os.path.join(ROOT, "src")


def run(solver: str, gmsh: bool = False):
    cmd = f'"{PYTHON}" -m src --solver {solver}'
    if gmsh:
        cmd += ' --gmsh'
    print(f"\n{'=' * 60}")
    print(f"  Generating: {solver}")
    print(f"{'=' * 60}")
    ret = os.system(cmd)
    if ret != 0:
        print(f"  [FAIL] {solver}")
    else:
        print(f"  [OK] {solver}")
    return ret


def main():
    gmsh = "--gmsh" in sys.argv[1:]
    failures = 0
    for solver in ['opti2d', 'opti3d', 'ccx', 'febio']:
        if run(solver, gmsh) != 0:
            failures += 1
    print(f"\n{'=' * 60}")
    if failures:
        print(f"  {4 - failures}/4 OK, {failures} failed")
    else:
        print("  All 4 solver files generated successfully!")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
