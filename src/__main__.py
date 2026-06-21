#!/usr/bin/env python3
"""
Foldable Display Folding — Multi-Solver Pre-Processor (CLI entry).

Generates FEM input files for 4 solvers:
  - OptiStruct 2D (plane strain, CQPSTN+PPLANE+MATHE)
  - OptiStruct 3D (solid, CHEXA+PSOLID+MATHE)
  - CalculiX   3D (solid, C3D8R+*HYPERELASTIC)
  - FEBio      3D (hex8+neo-Hookean+viscoelastic)

Usage:
  python -m src --solver ccx
  python -m src --solver opti2d
  python -m src --solver opti3d
  python -m src --solver febio
  python -m src --solver all          # all four in one go
"""

import argparse
import os
import sys
import numpy as np

from .materials import DEFAULT_PARAMS
from .mesh import build_mesh, extrude_to_3d


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Foldable display folding - multi-solver pre-processor")
    parser.add_argument('--solver', '-s', default='ccx',
                        choices=['opti2d', 'opti3d', 'ccx', 'febio', 'all'],
                        help="Target solver (default: ccx)")
    parser.add_argument('--outdir', '-o', default=None,
                        help="Output directory (default: current dir)")
    parser.add_argument('--gmsh', action='store_true',
                        help="Force GMSH mesh generation")
    parser.add_argument('--no-gmsh', dest='gmsh', action='store_false',
                        help="Force pure Python mesh")
    parser.set_defaults(gmsh=None)
    return parser.parse_args()


def main():
    p = DEFAULT_PARAMS
    args = parse_args()
    nlyr = p['n_layers']
    t = p['layer_thick']
    y_total = nlyr * t
    layers_str = '-'.join(['A' if i % 2 == 0 else 'B' for i in range(nlyr)])

    # Output directory
    out_dir = args.outdir
    if out_dir is None:
        out_dir = os.getcwd()
    os.makedirs(out_dir, exist_ok=True)

    # Determine mesh backend
    use_gmsh = args.gmsh
    if use_gmsh is None:
        use_gmsh = False
        try:
            import gmsh
            gmsh.initialize()
            gmsh.finalize()
            use_gmsh = True
        except (ImportError, OSError):
            pass
        except Exception:
            pass

    solver_list = ['opti2d', 'opti3d', 'ccx', 'febio'] if args.solver == 'all' else [args.solver]

    for solver in solver_list:
        solver_name = {'opti2d': 'OptiStruct 2D (plane strain)',
                       'opti3d': 'OptiStruct 3D (solid)',
                       'ccx':    'CalculiX 3D (solid)',
                       'febio':  'FEBio 3D (solid)'}[solver]

        print("=" * 65)
        print(f"  Foldable Display Folding - {solver_name}")
        print("=" * 65)
        print(f"  Layers: {layers_str} ({nlyr} layers x {t * 1000:.0f}um = {y_total:.3f}mm)")
        print(f"  Hinge:  X={p['hinge_L']}mm to X={p['hinge_R']}mm")
        print(f"  Fold:   {p['fold_angle']}deg/side, {p['fold_time']}s, {p['nsteps']} steps")
        print()

        if solver == 'febio':
            # FEBio uses its own full 3D mesh (left wing + hinge + right wing)
            from .writers.febio import build_mesh_3d, write as write_febio
            febio_p = {**p}
            nodes_3d, elems_3d, nsets_3d = build_mesh_3d(febio_p)
            path = os.path.join(out_dir, 'display_folding.feb')
            write_febio(febio_p, nodes_3d, elems_3d, nsets_3d, path)

        else:
            # Common path: 2D quad mesh → optionally extrude to 3D
            print(f"  Mesh backend: {'GMSH' if use_gmsh else 'Python'}")

            # Inject ny_per_layer for mesh builder
            mesh_p = {**p}
            mesh_p.setdefault('ny_per_layer', 3)

            nodes_2d, elems_2d, nsets_2d = build_mesh(mesh_p, use_gmsh=use_gmsh)
            print(f"  2D mesh: {len(nodes_2d)} nodes, {len(elems_2d)} quads")

            if solver == 'opti2d':
                from .writers.optistruct_2d import write as write_opti2d
                path = os.path.join(out_dir, 'display_folding_2d.fem')
                write_opti2d(p, nodes_2d, elems_2d, nsets_2d, path)

            elif solver == 'opti3d':
                from .writers.optistruct_3d import write as write_opti3d
                nodes_3d, elems_3d, nsets_3d = extrude_to_3d(
                    nodes_2d, elems_2d, nsets_2d, p['depth'])
                print(f"  3D extrusion: {len(nodes_3d)} nodes, {len(elems_3d)} hexes")
                path = os.path.join(out_dir, 'display_folding_3d.fem')
                write_opti3d(p, nodes_3d, elems_3d, nsets_3d, path)

            elif solver == 'ccx':
                from .writers.calculix import write as write_ccx
                nodes_3d, elems_3d, nsets_3d = extrude_to_3d(
                    nodes_2d, elems_2d, nsets_2d, p['depth'])
                print(f"  3D extrusion: {len(nodes_3d)} nodes, {len(elems_3d)} hexes")
                path = os.path.join(out_dir, 'display_folding.inp')
                write_ccx(p, nodes_3d, elems_3d, nsets_3d, path)

        print()

    # Print next steps
    if args.solver in ('opti2d', 'opti3d'):
        print("  [Next steps on OptiStruct PC]")
        print(f"    1. Copy files to the target PC")
        print(f"    2. Run: optistruct display_folding_*.fem")
        print(f"    3. View results in HyperView")
    elif args.solver == 'ccx':
        print("  [Next steps]")
        print(f"    1. Run: ccx display_folding")
        print(f"    2. View results in cgx or ParaView")
    elif args.solver == 'febio':
        print("  [Next steps]")
        print(f"    1. Run: febio4.exe display_folding.feb")
        print(f"    2. View results in FEBioStudio or ParaView")
    elif args.solver == 'all':
        print("  [All files generated]")
        print("    Copy to target PC(s) and run with respective solvers.")
    print("=" * 65)


if __name__ == '__main__':
    main()
