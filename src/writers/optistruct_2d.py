"""
OptiStruct 2D .fem writer -- plane strain with RBE2 hinge control.

Generates an OptiStruct input file using:
  - CQPSTN (quad plane strain) + PPLANE + MATHE (Neo-Hookean)
  - RBE2 elements: bottom RBE region nodes -> hinge master nodes
  - SPCD rotation on hinge master nodes (left -90deg, right +90deg)
"""

import math
import numpy as np
from ..materials import nh_c10, nh_d1, RBE_MASTER_IDS


def write(p: dict, nodes: np.ndarray, elems: list, node_sets: dict,
          filename: str) -> None:
    """Write OptiStruct 2D plane strain .fem with RBE hinge control."""
    nlyr = p['n_layers']
    x_hinge_L = p['hinge_L']
    x_hinge_R = p['hinge_R']
    mid_L = RBE_MASTER_IDS['L']
    mid_R = RBE_MASTER_IDS['R']

    lines = []
    L = lines.append

    # --- Header ---
    L("SOL 600")
    L("CEND")
    L("$")
    L("$ OptiStruct 2D Plane Strain -- RBE hinge control")
    L("$")
    L("SUBCASE 1")
    L("  LABEL = Display_Folding")
    L("  TITLE = Foldable display, RBE-controlled hinge folding")
    L("  ANALYSIS = NLSTAT")
    L("  SPC = 1")
    L("  LOAD = 1")
    L("  NLPARM = 10")
    L("  NLOUT = 1")
    L("$")
    L("BEGIN BULK")
    L("$")
    L("PARAM, LGDISP, 1")
    L("$")

    # --- Nodes (display mesh) ---
    L("$ DISPLAY NODES")
    for i in range(len(nodes)):
        L(f"GRID, {i+1},, {nodes[i,0]:.8f}, {nodes[i,1]:.8f}, 0.0")

    # --- Hinge master nodes ---
    L("$")
    L("$ HINGE MASTER REFERENCE NODES")
    L("$")
    L(f"GRID, {mid_L},, {x_hinge_L}, 0.0, 0.0")
    L(f"GRID, {mid_R},, {x_hinge_R}, 0.0, 0.0")
    L("$")

    # --- Materials ---
    L("$")
    L("$ MATERIALS (Neo-Hookean hyperelastic)")
    L("$")
    a_c10 = nh_c10(p['A_E'], p['A_nu'])
    a_d1 = nh_d1(p['A_E'], p['A_nu'])
    L(f"MATHE, 1, NEOHOOKE, {a_c10:.4f}, {a_d1:.6f}")
    L("$")
    b_c10 = nh_c10(p['B_E'], p['B_nu'])
    b_d1 = nh_d1(p['B_E'], p['B_nu'])
    L(f"MATHE, 2, NEOHOOKE, {b_c10:.4f}, {b_d1:.6f}")
    L("$")

    # --- Properties ---
    L("$ PROPERTIES (PPLANE -- plane strain)")
    for li in range(nlyr):
        mid = 1 if li % 2 == 0 else 2
        L(f"PPLANE, {li+1}, {mid}, {p['depth']:.4f}")
    L("$")

    # --- Elements ---
    L("$ ELEMENTS (CQPSTN -- quad plane strain)")
    for eid, (layer_idx, conn) in enumerate(elems, 1):
        pid = layer_idx + 1
        L(f"CQPSTN, {eid}, {pid}, {conn[0]}, {conn[1]}, {conn[2]}, {conn[3]}")

    # --- RBE2 elements: bottom RBE region -> hinge masters ---
    L("$")
    L("$ RBE2 -- bottom RBE region nodes to hinge masters")
    L("$   DOFs 1(UX) 2(UY) 6(RZ)")
    L("$")
    n_slaves_L = len(node_sets['BOTTOM_L'])
    n_slaves_R = len(node_sets['BOTTOM_R'])

    # RBE2 for left hinge (max 10 slaves per card, use continuation)
    L(f"$ LEFT hinge: {n_slaves_L} slave nodes")
    _write_rbe2(lines, 900001, mid_L, [1, 2, 6], node_sets['BOTTOM_L'])
    L("$")
    # RBE2 for right hinge
    L(f"$ RIGHT hinge: {n_slaves_R} slave nodes")
    _write_rbe2(lines, 900002, mid_R, [1, 2, 6], node_sets['BOTTOM_R'])
    L("$")

    # --- SPC: fix RBE master nodes ---
    # The master nodes have their UX/UY/RZ prescribed by SPCD,
    # so we fix the remaining DOFs (UZ only in 2D).
    L("$ SPC -- fix RBE master out-of-plane")
    L(f"SPC, 1, {mid_L}, 3, 0.0")
    L(f"SPC, 1, {mid_R}, 3, 0.0")
    L("$")

    # --- SPCD rotation on hinge masters ---
    angle_rad_L = -math.radians(p['fold_angle'])  # left: -90deg (RZ)
    angle_rad_R = +math.radians(p['fold_angle'])  # right: +90deg (RZ)
    L("$ SPCD -- hinge rotation (DOF 6 = RZ)")
    L(f"SPCD, 2, {mid_L}, 6, {angle_rad_L:.8e}")
    L(f"SPCD, 2, {mid_R}, 6, {angle_rad_R:.8e}")
    L("$")

    # --- TABLED1 ramp ---
    L("$")
    L("TABLED1, 1, , FOLDING_RAMP")
    for i in range(11):
        t = p['fold_time'] * i / 10
        scale = 1.0 * i / 10
        L(f"+,      {t:.6f}, {scale:.6f}")
    L("+,      ENDT")
    L("$")
    L("TLOAD1, 1, , , 2, 1")
    L("$")

    # --- NLPARM / NLOUT ---
    dt = p['fold_time'] / p['nsteps']
    L(f"NLPARM, 10, {p['nsteps']}, {dt:.6f}, UPW")
    L("$")
    L("NLOUT, 1, 1")
    L("$")

    L("ENDDATA")

    with open(filename, 'w', encoding='ascii') as f:
        f.write('\n'.join(lines))

    print(f"[OK] OptiStruct 2D: {filename}")
    print(f"     Nodes: {len(nodes)+2}, Elements: {len(elems)}"
          f", RBE slaves: L={n_slaves_L} R={n_slaves_R}")
    print(f"     Hinge rotation: L={p['fold_angle']:+.0f}deg,"
          f" R={p['fold_angle']:+.0f}deg")


def _write_rbe2(lines, eid: int, master: int, dofs: list[int],
                slaves: list[int]):
    """Write RBE2 card(s) with continuation lines (10 slaves/line)."""
    dof_str = ''.join(str(d) for d in dofs)
    first = True
    for i in range(0, len(slaves), 10):
        chunk = slaves[i:i + 10]
        if first:
            ids_str = ', '.join(str(n) for n in chunk)
            lines.append(f"RBE2, {eid}, {master}, {dof_str}, {ids_str}")
            first = False
        else:
            ids_str = ', '.join(str(n) for n in chunk)
            lines.append(f"+, {ids_str}")
