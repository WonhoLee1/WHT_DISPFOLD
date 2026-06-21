"""
FEBio .feb writer — solid (hex8 + neo-Hookean/viscoelastic).

Generates a FEBio input file with:
  - Full 3D mesh (left wing + hinge + right wing)
  - Neo-Hookean structural layers (A)
  - Viscoelastic PSA layers (B: Neo-Hookean + Prony series)
  - Rigid body controlled wing rotation
"""

import os
import math
import numpy as np
from ..materials import DEFAULT_PARAMS, FEBIO_VISCO_PARAMS


# ---------------------------------------------------------------------------
# FEBio mesh generator (full 3D, left-wing + hinge + right-wing)
# ---------------------------------------------------------------------------

def build_mesh_3d(p: dict) -> tuple:
    """
    Build full 3D hex mesh for FEBio (left wing + hinge + right wing).

    Returns
    -------
    nodes : (N, 3) ndarray
    elems : dict of {part_name: list_of_connectivity}
    node_sets : dict of {name: [nid, ...]}
    """
    L_tot = 100.0  # full display width [mm]
    L_h = p['hinge_L']
    R_h = p['hinge_R']
    W = 1.0         # Z-depth [mm]

    layers = [(f"L{i+1}", 'A' if i % 2 == 0 else 'B') for i in range(p['n_layers'])]
    nlyr = len(layers)
    t_lyr = p['layer_thick']
    ny_per = FEBIO_VISCO_PARAMS['ny_per_layer']
    nx_l = FEBIO_VISCO_PARAMS['nx_left']
    nx_h = FEBIO_VISCO_PARAMS['nx_hinge']
    nx_r = FEBIO_VISCO_PARAMS['nx_right']
    nz = FEBIO_VISCO_PARAMS['nz']

    # Y-coordinates
    y_edges = [0.0]
    for li in range(nlyr):
        for _ in range(ny_per):
            y_edges.append(y_edges[-1] + t_lyr / ny_per)
    Y_nodes = np.array(y_edges)
    ny_tot = len(Y_nodes) - 1

    # X-coordinates (3 zones)
    x_l = np.linspace(0, L_h, nx_l + 1)
    x_h = np.linspace(L_h, R_h, nx_h + 1)
    x_r = np.linspace(R_h, L_tot, nx_r + 1)
    x_edges = np.concatenate([x_l[:-1], x_h[:-1], x_r])
    nx_tot = len(x_edges) - 1

    z_edges = np.array([0.0, W])
    nnx, nny, nnz = nx_tot + 1, ny_tot + 1, nz + 1

    def nid(ix, iy, iz):
        return iz * (nnx * nny) + iy * nnx + ix + 1

    # Nodes
    nodes = []
    for iz in range(nnz):
        for iy in range(nny):
            for ix in range(nnx):
                nodes.append([float(x_edges[ix]), float(Y_nodes[iy]), float(z_edges[iz])])
    nodes = np.array(nodes)

    # Elements
    ix_hinge_start = nx_l
    ix_hinge_end = nx_l + nx_h

    elems = {name: [] for name, _ in layers}
    elems['LEFT_WING'] = []
    elems['RIGHT_WING'] = []

    for iz in range(nz):
        for iy in range(ny_tot):
            layer_idx = min(iy // ny_per, nlyr - 1)
            lname, _ = layers[layer_idx]
            for ix in range(nx_tot):
                is_left = ix < ix_hinge_start
                is_right = ix >= ix_hinge_end
                conn = [
                    nid(ix, iy, iz), nid(ix + 1, iy, iz),
                    nid(ix + 1, iy + 1, iz), nid(ix, iy + 1, iz),
                    nid(ix, iy, iz + 1), nid(ix + 1, iy, iz + 1),
                    nid(ix + 1, iy + 1, iz + 1), nid(ix, iy + 1, iz + 1),
                ]
                if is_left:
                    elems['LEFT_WING'].append(conn)
                elif is_right:
                    elems['RIGHT_WING'].append(conn)
                else:
                    elems[lname].append(conn)

    # Node sets
    NBOT, NTOP, NHL, NHR = [], [], [], []
    for iz in range(nnz):
        for iy in range(nny):
            for ix in range(nnx):
                nd = nid(ix, iy, iz)
                xy = x_edges[min(ix, nx_tot)]
                if iy == 0:
                    NBOT.append(nd)
                if iy == nny - 1:
                    NTOP.append(nd)
                if abs(xy - L_h) < 1e-10:
                    NHL.append(nd)
                if abs(xy - R_h) < 1e-10:
                    NHR.append(nd)

    node_sets = {
        'ALL':     sorted(range(1, len(nodes) + 1)),
        'BOTTOM':  sorted(set(NBOT)),
        'TOP':     sorted(set(NTOP)),
        'HINGE_L': sorted(set(NHL)),
        'HINGE_R': sorted(set(NHR)),
    }

    return nodes, elems, node_sets


# ---------------------------------------------------------------------------
# .feb file writer
# ---------------------------------------------------------------------------

def write(p: dict, nodes: np.ndarray, elems: dict, node_sets: dict,
          filename: str) -> None:
    """Write FEBio .feb file with viscoelastic PSA and rigid wing rotation."""
    MID_A, MID_B = 1, 2
    MID_LW, MID_RW = 3, 4

    L_h = p['hinge_L']
    R_h = p['hinge_R']
    angle = math.radians(p['fold_angle'])
    t_fold = p['fold_time']
    nsteps = p['nsteps']

    def nh_c1(E, nu): return E / (2.0 * (1.0 + nu))
    def nh_k(E, nu): return E / (3.0 * (1.0 - 2.0 * nu))

    AC1 = nh_c1(p['A_E'], p['A_nu'])
    AK = nh_k(p['A_E'], p['A_nu'])
    BC1 = nh_c1(p['B_E'], p['B_nu'])
    BK = nh_k(p['B_E'], p['B_nu'])

    lines = []
    L = lines.append

    L('<?xml version="1.0" encoding="ISO-8859-1"?>')
    L('<febio_spec version="4.0">')
    L('  <Module type="solid"/>')

    # Control
    L('  <Control>')
    L(f'    <time_steps>{nsteps}</time_steps>')
    L(f'    <step_size>{t_fold / nsteps:.6e}</step_size>')
    L('    <solver type="solid">')
    L('      <max_refs>25</max_refs>')
    L('      <diverge_reform>1</diverge_reform>')
    L('      <reform_each_time_step>1</reform_each_time_step>')
    L('      <qn_method type="1">')
    L('        <max_ups>10</max_ups>')
    L('      </qn_method>')
    L('    </solver>')
    L('    <time_stepper>')
    L('      <dtmin>1e-8</dtmin>')
    L('      <dtmax>0.01</dtmax>')
    L('      <max_retries>15</max_retries>')
    L('      <opt_iter>5</opt_iter>')
    L('    </time_stepper>')
    L('  </Control>')

    # Materials
    L('  <Material>')
    L(f'    <material id="1" name="A_layer" type="neo-Hookean">')
    L(f'      <c1>{AC1:.6f}</c1>')
    L(f'      <k>{AK:.6f}</k>')
    L('    </material>')
    # B-layer with viscoelasticity (if Prony params exist)
    L(f'    <material id="2" name="B_PSA" type="viscoelastic">')
    L('      <elastic type="neo-Hookean">')
    L(f'        <c1>{BC1:.6f}</c1>')
    L(f'        <k>{BK:.6f}</k>')
    L('      </elastic>')
    L('      <relaxation>')
    L(f'        <g1>{FEBIO_VISCO_PARAMS["B_g1"]}</g1> <t1>{FEBIO_VISCO_PARAMS["B_t1"]}</t1>')
    L(f'        <g2>{FEBIO_VISCO_PARAMS["B_g2"]}</g2> <t2>{FEBIO_VISCO_PARAMS["B_t2"]}</t2>')
    L('      </relaxation>')
    L('    </material>')
    # Rigid body wings
    L(f'    <material id="3" name="LeftWing" type="rigid body">')
    L('      <density>1e-9</density>')
    L(f'      <center_of_mass>{L_h}, 0, 0</center_of_mass>')
    L('    </material>')
    L(f'    <material id="4" name="RightWing" type="rigid body">')
    L('      <density>1e-9</density>')
    L(f'      <center_of_mass>{R_h}, 0, 0</center_of_mass>')
    L('    </material>')
    L('  </Material>')

    # Geometry
    L('  <Geometry>')
    L('    <Nodes name="all">')
    for i, (x, y, z) in enumerate(nodes):
        L(f'      <node id="{i + 1}">{x:.8f}, {y:.8f}, {z:.8f}</node>')
    L('    </Nodes>')

    # Deformable layers (A / B)
    mat_of = {'A': MID_A, 'B': MID_B}
    for name, mtype in [("L1", 'A'), ("L2", 'B'), ("L3", 'A'), ("L4", 'B'), ("L5", 'A')]:
        el = elems.get(name, [])
        if not el:
            continue
        mid = mat_of[mtype]
        L(f'    <Elements type="hex8" mat="{mid}" elset="{name}">')
        for eid, conn in enumerate(el, 1):
            L(f'      <elem id="{eid}"> {",".join(str(n) for n in conn)} </elem>')
        L('    </Elements>')

    # Rigid wing parts
    for part, mid in [('LEFT_WING', 3), ('RIGHT_WING', 4)]:
        el = elems[part]
        L(f'    <Elements type="hex8" mat="{mid}" elset="{part}">')
        for eid, conn in enumerate(el, 1):
            L(f'      <elem id="{eid}"> {",".join(str(n) for n in conn)} </elem>')
        L('    </Elements>')

    for name in ['ALL', 'BOTTOM', 'TOP', 'HINGE_L', 'HINGE_R']:
        ids = node_sets[name]
        L(f'    <NodeSet name="{name}">')
        for nd in ids:
            L(f'      <node id="{nd}"/>')
        L('    </NodeSet>')
    L('  </Geometry>')

    # Boundary: plane strain (u_z=0 on all nodes)
    L('  <Boundary>')
    L('    <bc type="zero_displacement">')
    L('      <node_set>ALL</node_set>')
    L('      <z_dof>1</z_dof>')
    L('    </bc>')
    L('  </Boundary>')

    # Rigid body control
    L('  <Rigid>')
    for rb_id in [3, 4]:
        L(f'    <rigid_bc type="rigid_fixed">')
        L(f'      <rb>{rb_id}</rb>')
        L('      <x_dof>1</x_dof> <y_dof>1</y_dof> <z_dof>1</z_dof>')
        L('      <Ru_dof>1</Ru_dof> <Rv_dof>1</Rv_dof>')
        L('    </rigid_bc>')
    L('    <rigid_bc type="rigid_rotation">')
    L('      <rb>3</rb>')
    L('      <dof>Rz</dof>')
    L(f'      <value lc="1">{angle:.6f}</value>')
    L('    </rigid_bc>')
    L('    <rigid_bc type="rigid_rotation">')
    L('      <rb>4</rb>')
    L('      <dof>Rz</dof>')
    L(f'      <value lc="2">{angle:.6f}</value>')
    L('    </rigid_bc>')
    L('  </Rigid>')

    # Load curves
    L('  <LoadData>')
    L('    <loadcurve id="1" type="linear">')
    L('      <point>0, 0</point>')
    L(f'      <point>{t_fold}, 1</point>')
    L('    </loadcurve>')
    L('    <loadcurve id="2" type="linear">')
    L('      <point>0, 0</point>')
    L(f'      <point>{t_fold}, -1</point>')
    L('    </loadcurve>')
    L('  </LoadData>')

    # Output
    L('  <Output>')
    L('    <plotfile type="xplt">')
    L('      <file>display_folding.xplt</file>')
    L('      <data_selection><data>ALL_DATA</data></data_selection>')
    L('      <compression>0</compression>')
    L('    </plotfile>')
    L('    <logfile>')
    L('      <file>display_folding_log.txt</file>')
    L('      <node_data data="u;R"><delimiter>,</delimiter></node_data>')
    L('      <element_data data="s1;s2;s3;s4;s5;s6;e1;e2;e3;e4;e5;e6">'
      '<delimiter>,</delimiter></element_data>')
    L('    </logfile>')
    L('  </Output>')
    L('</febio_spec>')

    with open(filename, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))

    total = sum(len(v) for v in elems.values())
    print(f"[OK] FEBio: {filename}")
    print(f"     Nodes: {len(nodes)}, Elements: {total}")
    print(f"     Left wing (rigid): {len(elems['LEFT_WING'])} elems")
    total_def = sum(len(elems[n]) for n in ['L1', 'L2', 'L3', 'L4', 'L5'])
    print(f"     Hinge (deform):    {total_def} elems")
    print(f"     Right wing (rigid): {len(elems['RIGHT_WING'])} elems")
