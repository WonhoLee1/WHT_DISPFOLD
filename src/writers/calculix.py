"""
CalculiX .inp writer — solid (C3D8R + *HYPERELASTIC).

Generates a CalculiX input file using 8-node reduced-integration hex elements
with Neo-Hookean hyperelastic material and AMPLITUDE ramp.
"""

import os
import numpy as np
from ..materials import nh_c10, nh_d1
from ..kinematics import compute_rotation_displacement


def write(p: dict, nodes: np.ndarray, elems: list, node_sets: dict,
          filename: str) -> int:
    """Write CalculiX 3D .inp file. Returns SPCD count."""
    nlyr = p['n_layers']
    x0, x1 = p['hinge_L'], p['hinge_R']
    lines = []
    L = lines.append

    L("*HEADING")
    L("Foldable Display Folding - CalculiX 3D")
    L("5-layer A-B-A-B-A, hinge 35-65mm, Neo-Hookean hyperelastic")
    L("")

    # Nodes
    L("*NODE")
    for i in range(len(nodes)):
        L(f"{i+1}, {nodes[i,0]:.8e}, {nodes[i,1]:.8e}, {nodes[i,2]:.8e}")

    # Element sets per layer
    for li in range(nlyr):
        elset_name = f"EL_LAYER{li+1}"
        mat_name = "MAT_A" if li % 2 == 0 else "MAT_B"
        L(f"*ELEMENT, TYPE=C3D8R, ELSET={elset_name}")
        for eid, (layer_idx, conn) in enumerate(elems, 1):
            if layer_idx == li:
                L(f"{eid}, {conn[0]}, {conn[1]}, {conn[2]}, {conn[3]}, "
                  f"{conn[4]}, {conn[5]}, {conn[6]}, {conn[7]}")
        L(f"*SOLID SECTION, ELSET={elset_name}, MATERIAL={mat_name}")
    L("")

    # Materials
    L("*MATERIAL, NAME=MAT_A")
    a_c10 = nh_c10(p['A_E'], p['A_nu'])
    a_d1 = nh_d1(p['A_E'], p['A_nu'])
    L("*HYPERELASTIC, NEO HOOKE")
    L(f"{a_c10:.6e}, {a_d1:.6e}")
    L("")

    L("*MATERIAL, NAME=MAT_B")
    b_c10 = nh_c10(p['B_E'], p['B_nu'])
    b_d1 = nh_d1(p['B_E'], p['B_nu'])
    L("*HYPERELASTIC, NEO HOOKE")
    L(f"{b_c10:.6e}, {b_d1:.6e}")
    L("")

    # Node sets
    def _write_nset(name, ids):
        L(f"*NSET, NSET={name}")
        for i in range(0, len(ids), 16):
            chunk = ids[i:i+16]
            L(", ".join(str(n) for n in chunk))

    _write_nset("PIVOT_L",
                [n for n in node_sets['BOTTOM']
                 if abs(nodes[n-1, 0] - x0) < 1e-10])
    _write_nset("PIVOT_R",
                [n for n in node_sets['BOTTOM']
                 if abs(nodes[n-1, 0] - x1) < 1e-10])
    _write_nset("EDGE_L", node_sets['EDGE_L'])
    _write_nset("EDGE_R", node_sets['EDGE_R'])
    L("")

    # SPC: fix Z at hinge pivots
    L("*BOUNDARY")
    L("PIVOT_L, 3, 3, 0.0")
    L("PIVOT_R, 3, 3, 0.0")
    L("")

    # Amplitude ramp
    L("*AMPLITUDE, NAME=RAMP, TIME=TOTAL TIME")
    L(f"0.0, 0.0, {p['fold_time']:.6f}, 1.0")
    L("")

    # Prescribed displacement (rotation)
    hinge_y = nlyr * p['layer_thick'] / 2.0
    bw = p.get('hinge_band_elements', 0) * (x1 - x0) / p['nx_hinge']
    spcd_left = compute_rotation_displacement(
        p, nodes, node_sets['EDGE_L_BAND'], x0, hinge_y, +1, bw)
    spcd_right = compute_rotation_displacement(
        p, nodes, node_sets['EDGE_R_BAND'], x1, hinge_y, -1, bw)
    spcd_all = spcd_left + spcd_right

    L("*BOUNDARY, AMPLITUDE=RAMP")
    for nid, dof, val in spcd_all:
        L(f"{nid}, {dof}, {dof}, {val:.8e}")
    L("")

    # Step
    L("*STEP, NLGEOM")
    dt = p['fold_time'] / p['nsteps']
    L("*STATIC")
    L(f"{dt:.6e}, {p['fold_time']:.6e}")
    L("")
    L("*NODE FILE, OUTPUT=2")
    L("U")
    L("*EL FILE, OUTPUT=2")
    L("S, E")
    L("")
    L("*END STEP")

    with open(filename, 'w', encoding='ascii') as f:
        f.write('\n'.join(lines))

    print(f"[OK] CalculiX 3D: {filename}")
    print(f"     Nodes: {len(nodes)}, Elements: {len(elems)}, SPCD: {len(spcd_all)}")
    return len(spcd_all)
