#!/usr/bin/env python3
"""
Foldable Display Folding — Multi-Solver Pre-Processor (CLI entry).

Generates FEM input files for 4 solvers with RBE hinge control:
  - OptiStruct 2D (plane strain, CQPSTN+PPLANE+RBE2)
  - OptiStruct 3D (solid, CHEXA+PSOLID+RBE2)
  - CalculiX   3D (solid, C3D8R+KINEMATIC COUPLING)
  - FEBio      3D (hex8+neo-Hookean+rigid body)

Usage:
  python -m src --solver ccx
  python -m src --solver opti2d
  python -m src --solver opti3d
  python -m src --solver febio
  python -m src --solver all
"""

import argparse
import os
import sys

from .materials import DEFAULT_PARAMS, RBE_MASTER_IDS
from .mesh import build_mesh, extrude_to_3d, build_hinge_cylinders


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

    out_dir = args.outdir
    if out_dir is None:
        out_dir = os.getcwd()
    os.makedirs(out_dir, exist_ok=True)

    # Mesh backend detection
    use_gmsh = args.gmsh
    if use_gmsh is None:
        use_gmsh = False
        try:
            import gmsh
            gmsh.initialize()
            gmsh.finalize()
            use_gmsh = True
        except Exception:
            pass

    solvers = ['opti2d', 'opti3d', 'ccx', 'febio'] if args.solver == 'all' else [args.solver]

    for solver in solvers:
        solver_name = {
            'opti2d': 'OptiStruct 2D (plane strain, RBE2)',
            'opti3d': 'OptiStruct 3D (solid, RBE2)',
            'ccx':    'CalculiX 3D (direct nodal displacement)',
            'febio':  'FEBio 3D (rigid body)',
        }[solver]

        print("=" * 65)
        print(f"  Foldable Display Folding - {solver_name}")
        print("=" * 65)
        print(f"  Layers: {layers_str} ({nlyr} layers x {t*1000:.0f}um"
              f" = {y_total:.3f}mm)")
        W = p['total_width']
        rbe = p['rbe_region']
        print(f"  Width:  {W:.0f}mm ({-W/2:.0f}..{W/2:.0f}),"
              f" RBE region: |X| > {rbe:.0f}mm")
        print(f"  Hinges: X={p['hinge_L']:.0f}mm, X={p['hinge_R']:.0f}mm")
        print(f"  Fold:   {p['fold_angle']}deg/side, {p['fold_time']}s,"
              f" {p['nsteps']} steps")
        print()

        if solver == 'febio':
            from .writers.febio import build_mesh_3d, write as write_febio
            nodes_3d, elems_3d, nsets_3d = build_mesh_3d(p)
            path = os.path.join(out_dir, 'display_folding.feb')
            write_febio(p, nodes_3d, elems_3d, nsets_3d, path)

        else:
            print(f"  Mesh: {'GMSH' if use_gmsh else 'Python'} structured quad")
            nodes_2d, elems_2d, nsets_2d = build_mesh(p, use_gmsh=use_gmsh)
            print(f"  2D mesh: {len(nodes_2d)} nodes, {len(elems_2d)} quads")
            print(f"  RBE slaves: L={len(nsets_2d['BOTTOM_L'])}"
                  f" R={len(nsets_2d['BOTTOM_R'])}")

            if solver == 'opti2d':
                from .writers.optistruct_2d import write as write_o2d
                path = os.path.join(out_dir, 'display_folding_2d.fem')
                write_o2d(p, nodes_2d, elems_2d, nsets_2d, path)

            elif solver == 'opti3d':
                from .writers.optistruct_3d import write as write_o3d
                nodes_3d, elems_3d, nsets_3d = extrude_to_3d(
                    nodes_2d, elems_2d, nsets_2d, p['depth'])
                print(f"  Extruded: {len(nodes_3d)} nodes, {len(elems_3d)} hexes")
                path = os.path.join(out_dir, 'display_folding_3d.fem')
                write_o3d(p, nodes_3d, elems_3d, nsets_3d, path)

            elif solver == 'ccx':
                from .writers.calculix import write as write_ccx
                nodes_3d, elems_3d, nsets_3d = extrude_to_3d(
                    nodes_2d, elems_2d, nsets_2d, p['depth'])
                print(f"  Extruded: {len(nodes_3d)} nodes, {len(elems_3d)} hexes")
                cyl_depth = p.get('hinge_cylinder_depth', p['depth'])
                nodes_3d, elems_3d, nsets_3d = build_hinge_cylinders(
                    nodes_3d, elems_3d, nsets_3d, p, cyl_depth)
                path = os.path.join(out_dir, 'display_folding.inp')
                write_ccx(p, nodes_3d, elems_3d, nsets_3d, path)

        print()

    # Next steps
    s = args.solver
    if s == 'opti2d':
        print("  Run: optistruct display_folding_2d.fem  (OptiStruct PC)")
    elif s == 'opti3d':
        print("  Run: optistruct display_folding_3d.fem  (OptiStruct PC)")
    elif s == 'ccx':
        print("  Run: ccx display_folding")
    elif s == 'febio':
        print("  Run: febio4.exe display_folding.feb")
    elif s == 'all':
        print("  All 4 files generated. Copy to respective solver PC(s).")
    print("=" * 65)


if __name__ == '__main__':
    main()
