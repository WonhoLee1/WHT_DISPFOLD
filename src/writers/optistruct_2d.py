"""
OptiStruct 2D .fem writer — plane strain (CQPSTN + PPLANE + MATHE).

Generates an OptiStruct input file for the foldable display problem using
2D plane strain quad elements with Neo-Hookean hyperelastic material.
"""

import os
import numpy as np
from ..materials import nh_c10, nh_d1
from ..kinematics import compute_rotation_displacement


def write(p: dict, nodes: np.ndarray, elems: list, node_sets: dict,
          filename: str) -> int:
    """
    Write OptiStruct 2D plane strain .fem file.

    Returns number of SPCD entries written.
    """
    nlyr = p['n_layers']
    x0, x1 = p['hinge_L'], p['hinge_R']
    lines = []
    L = lines.append

    # Header
    L("SOL 600")
    L("CEND")
    L("$")
    L("$ OptiStruct 2D Plane Strain - Foldable Display Folding")
    L("$")
    L("SUBCASE 1")
    L("  LABEL = Display_Folding")
    L("  TITLE = Multi-layer foldable display, hinge region folding")
    L("  ANALYSIS = NLSTAT")
    L("  SPC = 1")
    L("  LOAD = 1")
    L("  NLPARM = 10")
    L("  NLOUT = 1")
    L("$")
    L("BEGIN BULK")
    L("$")
    L("$ PARAMETERS")
    L("$")
    L("PARAM, LGDISP, 1")
    L("$")

    # Nodes
    L("$ NODES")
    for i in range(len(nodes)):
        L(f"GRID, {i + 1},, {nodes[i, 0]:.8f}, {nodes[i, 1]:.8f}, 0.0")

    # Materials (Neo-Hookean)
    L("$")
    L("$ MATERIALS (Neo-Hookean hyperelastic)")
    L("$")
    a_c10 = nh_c10(p['A_E'], p['A_nu'])
    a_d1 = nh_d1(p['A_E'], p['A_nu'])
    L(f"$ A-layer: C10={a_c10:.4f}, D1={a_d1:.6f}  (E={p['A_E']}MPa, nu={p['A_nu']})")
    L(f"MATHE, 1, NEOHOOKE, {a_c10:.4f}, {a_d1:.6f}")
    L("$")
    b_c10 = nh_c10(p['B_E'], p['B_nu'])
    b_d1 = nh_d1(p['B_E'], p['B_nu'])
    L(f"$ B-layer (PSA): C10={b_c10:.4f}, D1={b_d1:.6f}  (E={p['B_E']}MPa, nu={p['B_nu']})")
    L(f"MATHE, 2, NEOHOOKE, {b_c10:.4f}, {b_d1:.6f}")
    L("$")

    # Properties (PPLANE - plane strain)
    L("$ PROPERTIES (PPLANE - plane strain)")
    for li in range(nlyr):
        mid = 1 if li % 2 == 0 else 2
        L(f"PPLANE, {li + 1}, {mid}, {p['depth']:.4f}")
    L("$")

    # Elements (CQPSTN)
    L("$ ELEMENTS (CQPSTN - quad plane strain)")
    for eid, (layer_idx, conn) in enumerate(elems, 1):
        pid = layer_idx + 1
        L(f"CQPSTN, {eid}, {pid}, {conn[0]}, {conn[1]}, {conn[2]}, {conn[3]}")

    # SPC
    L("$")
    L("$ SPC SET 1: Z-translation only (hinge pivots)")
    L("$   UX/UY are prescribed by rotation SPCD")
    L("$")
    for nid in node_sets['BOTTOM']:
        x = nodes[nid - 1, 0]
        if abs(x - x0) < 1e-10 or abs(x - x1) < 1e-10:
            L(f"SPC, 1, {nid}, 3, 0.0")
    L("$")

    # SPCD rotation
    L("$")
    L(f"$ SPCD SET 2: Rotation about hinges (band = {p.get('hinge_band_elements', 0)} elements)")
    L("$")
    hinge_y = nlyr * p['layer_thick'] / 2.0
    band_width = p.get('hinge_band_elements', 0) * (x1 - x0) / p['nx_hinge']
    spcd_left = compute_rotation_displacement(
        p, nodes, node_sets['EDGE_L_BAND'], x0, hinge_y, +1, band_width)
    spcd_right = compute_rotation_displacement(
        p, nodes, node_sets['EDGE_R_BAND'], x1, hinge_y, -1, band_width)
    spcd_all = spcd_left + spcd_right
    for nid, dof, val in spcd_all:
        L(f"SPCD, 2, {nid}, {dof}, {val:.8e}")
    L("$")

    # TABLED1 ramp
    L("$")
    L("$ TABLED1: Linear ramp 0->1")
    L("$")
    L("TABLED1, 1, , FOLDING_RAMP")
    for i in range(11):
        t = p['fold_time'] * i / 10
        scale = 1.0 * i / 10
        L(f"+,      {t:.6f}, {scale:.6f}")
    L("+,      ENDT")
    L("$")
    L("$ TLOAD1: SPCD set 2 scaled by TABLED1 1")
    L("$")
    L("TLOAD1, 1, , , 2, 1")
    L("$")

    # NLPARM / NLOUT
    dt = p['fold_time'] / p['nsteps']
    L("$")
    L("$ NLPARM: Nonlinear static")
    L("$")
    L(f"NLPARM, 10, {p['nsteps']}, {dt:.6f}, UPW")
    L("$")
    L("$ NLOUT: Output at all steps")
    L("$")
    L("NLOUT, 1, 1")
    L("$")

    L("ENDDATA")

    with open(filename, 'w', encoding='ascii') as f:
        f.write('\n'.join(lines))

    print(f"[OK] OptiStruct 2D: {filename}")
    print(f"     Nodes: {len(nodes)}, Elements: {len(elems)}, SPCD: {len(spcd_all)}")
    return len(spcd_all)
