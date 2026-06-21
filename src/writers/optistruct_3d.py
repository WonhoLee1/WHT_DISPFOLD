"""
OptiStruct 3D .fem writer -- solid CHEXA with RBE2 hinge control.

Generates an OptiStruct input file using:
  - CHEXA (8-node hex) + PSOLID + MATHE (Neo-Hookean)
  - RBE2 elements: bottom RBE region nodes -> hinge master nodes
  - SPCD rotation on hinge master nodes (left -90deg, right +90deg)
"""

import math
import numpy as np
from ..materials import nh_c10, nh_d1, RBE_MASTER_IDS


def write(p: dict, nodes: np.ndarray, elems: list, node_sets: dict,
          filename: str) -> None:
    """Write OptiStruct 3D solid .fem with RBE hinge control."""
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
    L("$ OptiStruct 3D Solid -- RBE hinge control")
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

    # --- Nodes ---
    L("$ DISPLAY NODES (x, y, z)")
    for i in range(len(nodes)):
        L(f"GRID, {i+1},, {nodes[i,0]:.8f}, {nodes[i,1]:.8f}, {nodes[i,2]:.8f}")

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
    L("$ PROPERTIES (PSOLID -- solid)")
    for li in range(nlyr):
        mid = 1 if li % 2 == 0 else 2
        L(f"PSOLID, {li+1}, {mid}")
    L("$")

    # --- Elements ---
    L("$ ELEMENTS (CHEXA -- 8-node hex)")
    for eid, (layer_idx, conn) in enumerate(elems, 1):
        pid = layer_idx + 1
        L(f"CHEXA, {eid}, {pid}, {conn[0]}, {conn[1]}, {conn[2]}, {conn[3]}, "
          f"{conn[4]}, {conn[5]}, {conn[6]}, {conn[7]}")

    # --- Symmetry BC: Y=0 face → UY=0 (DOF 2) ---
    # Coord system: X=width, Y=depth(0..1mm,sym), Z=thickness(0..0.15mm)
    L("$")
    L("$ SPC -- symmetry plane Y=0: UY=0 (DOF 2)")
    L("$")
    eps = 1e-12
    for i in range(len(nodes)):
        if abs(nodes[i, 1]) < eps:      # Y=0 face
            L(f"SPC, 1, {i+1}, 2, 0.0")
    L("$")

    # --- RBE2 elements: Z=0 panel bottom nodes → hinge masters ---
    # Solid elements have only translational DOFs → list DOFs 1,2,3 only.
    # Master RY (DOF 5) rotation drives slave translations via rigid-body kinematics.
    L("$")
    L("$ RBE2 -- Z=0 bottom-surface nodes to hinge masters (DOFs 1,2,3)")
    L("$")
    nL = len(node_sets['BOTTOM_L'])
    nR = len(node_sets['BOTTOM_R'])
    L(f"$ LEFT hinge: {nL} slave nodes")
    _write_rbe2(lines, 900001, mid_L, [1, 2, 3], node_sets['BOTTOM_L'])
    L("$")
    L(f"$ RIGHT hinge: {nR} slave nodes")
    _write_rbe2(lines, 900002, mid_R, [1, 2, 3], node_sets['BOTTOM_R'])
    L("$")

    # --- SPC: fix all translations on hinge masters (DOF 1-3) ---
    L("$")
    L("$ SPC -- fix master node translations (position anchored)")
    L(f"SPC, 1, {mid_L}, 123, 0.0")
    L(f"SPC, 1, {mid_R}, 123, 0.0")
    L("$")

    # --- SPCD rotation about Y-axis (DOF 5 = RY) ---
    # Left: -90 deg (CW from +Y), Right: +90 deg (CCW from +Y) → U-shape
    angle_L = -math.radians(p['fold_angle'])
    angle_R = +math.radians(p['fold_angle'])
    L("$ SPCD -- hinge rotation (DOF 5 = RY)")
    L(f"SPCD, 2, {mid_L}, 5, {angle_L:.8e}")
    L(f"SPCD, 2, {mid_R}, 5, {angle_R:.8e}")
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

    print(f"[OK] OptiStruct 3D: {filename}")
    print(f"     Nodes: {len(nodes)+2}, Elements: {len(elems)},"
          f" RBE slaves: L={nL} R={nR}")
    print(f"     Hinge rotation: L={p['fold_angle']:+.0f}deg,"
          f" R={p['fold_angle']:+.0f}deg")


def _write_rbe2(lines, eid: int, master: int, dofs: list[int],
                slaves: list[int]):
    """Write RBE2 card(s) with continuation lines (10 slaves/line)."""
    dof_str = ''.join(str(d) for d in dofs)
    for i in range(0, len(slaves), 10):
        chunk = slaves[i:i + 10]
        ids_str = ', '.join(str(n) for n in chunk)
        if i == 0:
            lines.append(f"RBE2, {eid}, {master}, {dof_str}, {ids_str}")
        else:
            lines.append(f"+, {ids_str}")
