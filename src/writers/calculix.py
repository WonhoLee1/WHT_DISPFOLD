"""
CalculiX .inp writer -- solid C3D8R, U-fold about Y-axis.

Coordinate system: X=width(-80..+80), Y=depth(0..1mm, sym Y=0), Z=thickness(0..0.15mm)
Pivot axes run along Y at (X=hinge_L, Z=0) and (X=hinge_R, Z=0).

Folding motion (90 deg):
  Left panel  (X < hinge_L): rotates -90 deg about Y at X=hinge_L  (CW from +Y)
  Right panel (X > hinge_R): rotates +90 deg about Y at X=hinge_R  (CCW from +Y)
  Result: U-shape with arms pointing up (+Z), orange (top layer) inside the curve.

Arc-following AMPLITUDE tables (smooth step s(t)=t^2(3-2t), theta(t)=theta_final*s(t)):
    AMP_UX(t) = 1 - cos(theta(t))   (peaks at 1 when theta=pi/2)
    AMP_UZ(t) = sin(theta(t))        (peaks at 1 when theta=pi/2)
    AMP_ROT(t) = s(t)               (linear smooth ramp 0->1 for cylinder DOF5)
"""

import math
import numpy as np
from ..materials import nh_c10, nh_d1, RBE_MASTER_IDS


def _rotation_displacement(x: float, z: float, x_h: float, theta: float) -> tuple:
    """Rotation in XZ plane about Y-axis pivot at (x_h, Z=0). Returns (UX, UZ)."""
    c = math.cos(theta)
    s = math.sin(theta)
    rx = x - x_h
    ux = rx * (c - 1.0) - z * s
    uz = rx * s + z * (c - 1.0)
    return ux, uz


def write(p: dict, nodes: np.ndarray, elems: list, node_sets: dict,
          filename: str) -> None:
    nlyr = p['n_layers']
    x_hinge_L = p['hinge_L']
    x_hinge_R = p['hinge_R']
    rbe = p['rbe_region']

    angle_deg = p['fold_angle']
    # Left: -90 deg (CW from +Y), Right: +90 deg (CCW from +Y) → U-shape arms up
    theta_L = -math.radians(angle_deg)
    theta_R = +math.radians(angle_deg)

    # Arc-following AMPLITUDE tables (smooth step s(t)=t^2(3-2t))
    n_amp = 21
    ux_amp_vals, uz_amp_vals, rot_amp_vals = [], [], []
    for i in range(n_amp):
        t = i / (n_amp - 1)
        s = t * t * (3.0 - 2.0 * t)
        th = 0.5 * math.pi * s
        ux_amp_vals.append((t, 1.0 - math.cos(th)))
        uz_amp_vals.append((t, math.sin(th)))
        rot_amp_vals.append((t, s))

    def _format_amplitude(name, vals, lines_target):
        lines_target.append(f"*AMPLITUDE, NAME={name}, TIME=TOTAL TIME")
        for chunk_start in range(0, len(vals), 8):
            chunk = vals[chunk_start:chunk_start + 8]
            line_parts = []
            for tv in chunk:
                line_parts.append(f"{tv[0]:.6e}")
                line_parts.append(f"{tv[1]:.12e}")
            lines_target.append(", ".join(line_parts))
        lines_target.append("")

    lines = []
    L = lines.append

    L("*HEADING")
    L("Foldable Display Folding - CalculiX 3D with Hinge Cylinders")
    L(f"5-layer A-B-A-B-A, hinge cylinders at X={x_hinge_L:.0f}, X={x_hinge_R:.0f}")
    L("")

    # --- Nodes ---
    L("*NODE, NSET=ALL")
    for i in range(len(nodes)):
        L(f"{i+1}, {nodes[i,0]:.8e}, {nodes[i,1]:.8e}, {nodes[i,2]:.8e}")
    for mid, hx in [(RBE_MASTER_IDS['L'], p['hinge_L']),
                    (RBE_MASTER_IDS['R'], p['hinge_R'])]:
        L(f"{mid}, {hx:.8e}, 0.00000000e+00, 0.00000000e+00")
    L("")

    # --- Display element sets per layer ---
    for li in range(nlyr):
        elset_name = f"EL_LAYER{li+1}"
        mat_name = "MAT_A" if li % 2 == 0 else "MAT_B"
        L(f"*ELEMENT, TYPE=C3D8R, ELSET={elset_name}")
        count = 0
        for eid, (layer_idx, conn) in enumerate(elems, 1):
            if layer_idx == li:
                L(f"{eid}, {conn[0]}, {conn[1]}, {conn[2]}, {conn[3]}, "
                  f"{conn[4]}, {conn[5]}, {conn[6]}, {conn[7]}")
                count += 1
        if count > 0:
            L(f"*SOLID SECTION, ELSET={elset_name}, MATERIAL={mat_name}")
    L("")

    # --- Hinge cylinder element sets ---
    for side_key, elset_name, layer_val in [
            ('L', 'CYL_L', -1), ('R', 'CYL_R', -2)]:
        cyl_elems = [(eid, conn) for eid, (layer_idx, conn)
                     in enumerate(elems, 1) if layer_idx == layer_val]
        if not cyl_elems:
            continue
        L(f"*ELEMENT, TYPE=C3D8R, ELSET={elset_name}")
        for eid, conn in cyl_elems:
            L(f"{eid}, {conn[0]}, {conn[1]}, {conn[2]}, {conn[3]}, "
              f"{conn[4]}, {conn[5]}, {conn[6]}, {conn[7]}")
        L(f"*SOLID SECTION, ELSET={elset_name}, MATERIAL=MAT_CYL")
    L("")

    # --- Materials ---
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

    # Cylinder material (steel-like linear elastic; *RIGID BODY overrides)
    L("*MATERIAL, NAME=MAT_CYL")
    L("*ELASTIC")
    L("2.100000e+05, 3.000000e-01")
    L("")

    # --- Hinge cylinder rigid bodies ---
    for side_key, elset_name, master_id in [
            ('L', 'CYL_L', RBE_MASTER_IDS['L']),
            ('R', 'CYL_R', RBE_MASTER_IDS['R'])]:
        cyl_elems = [eid for eid, (layer_idx, _)
                     in enumerate(elems, 1) if layer_idx == (-1 if side_key == 'L' else -2)]
        if not cyl_elems:
            continue
        L(f"*RIGID BODY, REF NODE={master_id}, ELSET={elset_name}")
    L("")

    # --- Fixed boundary conditions (pre-step) ---
    eps = 1e-12

    # Symmetry BC: Y=0 face (depth symmetry plane) → UY=0 (DOF 2)
    # Exclude cylinder element nodes (layer_idx < 0) to avoid MPC/SPC conflict
    # with *RIGID BODY at nodes where cylinder ring passes through Z=0.
    cyl_node_ids = set()
    for layer_idx, conn in elems:
        if layer_idx < 0:
            cyl_node_ids.update(conn)
    L("*BOUNDARY")
    for i in range(len(nodes)):
        nid = i + 1
        if nid in cyl_node_ids:
            continue
        if abs(nodes[i, 1]) < eps:       # Y=0 symmetry plane
            L(f"{nid}, 2, 2, 0.0")
    L("")

    # Center-line constraint: X=0, Z=0 nodes → UX=0 (symmetry in X)
    # Prevents rigid-body drift of the free hinge zone during large fold.
    L("*BOUNDARY")
    for i in range(len(nodes)):
        nid = i + 1
        if nid in cyl_node_ids:
            continue
        x, y, z = nodes[i, 0], nodes[i, 1], nodes[i, 2]
        if abs(x) < eps and abs(z) < eps:
            L(f"{nid}, 1, 1, 0.0")
    L("")

    # Cylinder REF nodes: fix all translations (UX=UY=UZ=0), DOF 1-3
    L("*BOUNDARY")
    for master_id in [RBE_MASTER_IDS['L'], RBE_MASTER_IDS['R']]:
        L(f"{master_id}, 1, 3, 0.0")
    L("")

    # --- Amplitude definitions ---
    _format_amplitude("AMP_UX",  ux_amp_vals,  lines)
    _format_amplitude("AMP_UZ",  uz_amp_vals,  lines)
    _format_amplitude("AMP_ROT", rot_amp_vals, lines)

    # --- Step ---
    dt = p['fold_time'] / p['nsteps']
    L("*STEP, NLGEOM, INC=2000")
    L("*STATIC")
    L(f"{dt:.6e}, {p['fold_time']:.6e}, 1.000000e-06, {dt:.6e}")
    L("*CONTROLS, PARAMETERS=TIME INCREMENTATION")
    L("4, 8, 9, 16, 10, 4, 12, 20, 6, 3, 3")
    L("")

    # --- Prescribed UX on Z=0 panel nodes (DOF 1) ---
    # Skip cylinder nodes (already constrained via *RIGID BODY MPC).
    L("*BOUNDARY, AMPLITUDE=AMP_UX")
    nL = nR = 0
    for i in range(len(nodes)):
        nid = i + 1
        if nid in cyl_node_ids:
            continue
        x, y, z = nodes[i, 0], nodes[i, 1], nodes[i, 2]
        if abs(z) > eps:
            continue
        if x < -rbe + eps:
            ux, _ = _rotation_displacement(x, z, x_hinge_L, theta_L)
            L(f"{nid}, 1, 1, {ux:.8e}")
            nL += 1
        elif x > rbe - eps:
            ux, _ = _rotation_displacement(x, z, x_hinge_R, theta_R)
            L(f"{nid}, 1, 1, {ux:.8e}")
            nR += 1
    L("")

    # --- Prescribed UZ on Z=0 panel nodes (DOF 3) ---
    L("*BOUNDARY, AMPLITUDE=AMP_UZ")
    for i in range(len(nodes)):
        nid = i + 1
        if nid in cyl_node_ids:
            continue
        x, y, z = nodes[i, 0], nodes[i, 1], nodes[i, 2]
        if abs(z) > eps:
            continue
        if x < -rbe + eps:
            _, uz = _rotation_displacement(x, z, x_hinge_L, theta_L)
            L(f"{nid}, 3, 3, {uz:.8e}")
        elif x > rbe - eps:
            _, uz = _rotation_displacement(x, z, x_hinge_R, theta_R)
            L(f"{nid}, 3, 3, {uz:.8e}")
    L("")

    # --- Cylinder rigid body Y-rotation (DOF 5), same angle as panel BCs ---
    L("*BOUNDARY, AMPLITUDE=AMP_ROT")
    L(f"{RBE_MASTER_IDS['L']}, 5, 5, {theta_L:.8e}")
    L(f"{RBE_MASTER_IDS['R']}, 5, 5, {theta_R:.8e}")
    L("")

    L("*NODE FILE")
    L("U")
    L("*EL FILE")
    L("S, E")
    L("")
    L("*END STEP")

    with open(filename, 'w', encoding='ascii') as f:
        f.write('\n'.join(lines))

    n_display = sum(1 for _, li in [(e, li) for e, (li, _) in enumerate(elems)] if li >= 0)
    n_cyl = sum(1 for _, li in [(e, li) for e, (li, _) in enumerate(elems)] if li < 0)
    print(f"[OK] CalculiX 3D: {filename}")
    print(f"     Nodes: {len(nodes)}, Elements: {len(elems)} "
          f"(display={n_display}, hinge_cyl={n_cyl})")
    n_slaves_L = len(node_sets.get('HINGE_BOTTOM_L', []))
    n_slaves_R = len(node_sets.get('HINGE_BOTTOM_R', []))
    print(f"     Hinge RBE2 slaves (KINEMATIC COUPLING): L={n_slaves_L} R={n_slaves_R}")
    print(f"     Cylinder ref nodes: L={RBE_MASTER_IDS['L']} R={RBE_MASTER_IDS['R']} (UX,UY fixed)")
    print(f"     Prescribed displacement on |X|>{rbe:.0f}: L={nL} R={nR}")
    print(f"     Hinge rotation: L={-angle_deg:+.0f}deg, R={+angle_deg:+.0f}deg")
