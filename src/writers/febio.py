"""
FEBio .feb writer -- hex8 solid with rigid-body hinge control.

Generates a FEBio input file for the full display (-60..+60mm):
  - Full 3D mesh with structured hex8 elements
  - Neo-Hookean structural layers (A) + viscoelastic PSA (B)
  - Rigid body controlled hinge rotation (left -90deg, right +90deg)
  - Bottom RBE-region nodes are attached to rigid bodies
"""

import math
import numpy as np
from ..materials import DEFAULT_PARAMS, FEBIO_VISCO_PARAMS


def build_mesh_3d(p: dict) -> tuple:
    """
    Build full 3D hex mesh for FEBio (X=0 symmetric, full 120mm width).

    Returns
    -------
    nodes : (N, 3) ndarray
    elems : dict of {part_name: list_of_connectivity}
    node_sets : dict of {name: [nid, ...]}
    """
    x0 = -p['total_width'] / 2.0
    x1 = +p['total_width'] / 2.0
    W = p['depth']
    rbe = p['rbe_region']

    nlyr = p['n_layers']
    layers = [(f"L{i+1}", 'A' if i % 2 == 0 else 'B') for i in range(nlyr)]
    t_lyr = p['layer_thick']
    ny_per = p['ny_per_layer']
    nz = FEBIO_VISCO_PARAMS['nz']

    # Mesh divisions (finer in hinge, coarser in wings/RBE)
    # RBE left : x0 .. -rbe
    # Left wing: -rbe .. hinge_L
    # Hinge    : hinge_L .. hinge_R
    # Right wing: hinge_R .. rbe
    # RBE right : rbe .. x1
    hinge_L = p['hinge_L']
    hinge_R = p['hinge_R']

    # Use finer resolution: ~20/mm in hinge, ~10/mm in wings/RBE
    def _n_elems(length, density):
        return max(2, int(round(length * density)))

    d_hinge = 20.0   # elements/mm in hinge
    d_other = 10.0   # elements/mm elsewhere

    nx_Lrbe  = _n_elems(-rbe - x0, d_other)       # left RBE region
    nx_Lwing = _n_elems(hinge_L - (-rbe), d_other) # left wing
    nx_hinge = _n_elems(hinge_R - hinge_L, d_hinge) # hinge
    nx_Rwing = _n_elems(rbe - hinge_R, d_other)    # right wing
    nx_Rrbe  = _n_elems(x1 - rbe, d_other)         # right RBE region

    # X coordinates
    x_edges = np.concatenate([
        np.linspace(x0, -rbe, nx_Lrbe + 1)[:-1],
        np.linspace(-rbe, hinge_L, nx_Lwing + 1)[:-1],
        np.linspace(hinge_L, hinge_R, nx_hinge + 1)[:-1],
        np.linspace(hinge_R, rbe, nx_Rwing + 1)[:-1],
        np.linspace(rbe, x1, nx_Rrbe + 1),
    ])
    nx_tot = len(x_edges) - 1

    # Y coordinates
    y_edges = [0.0]
    for li in range(nlyr):
        for _ in range(ny_per):
            y_edges.append(y_edges[-1] + t_lyr / ny_per)
    Y_nodes = np.array(y_edges)
    ny_tot = len(Y_nodes) - 1

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

    # Elements by region
    def _find_ix_range(x_start, x_end):
        """Return (i_start, i_end) covering elements in [x_start, x_end]."""
        i_start = int(np.searchsorted(x_edges, x_start, side='right'))
        i_end = int(np.searchsorted(x_edges, x_end, side='left'))
        return max(0, i_start), min(nx_tot, i_end)

    i_Lwing_s, i_Lwing_e = _find_ix_range(-rbe, hinge_L)
    i_hinge_s, i_hinge_e = _find_ix_range(hinge_L, hinge_R)
    i_Rwing_s, i_Rwing_e = _find_ix_range(hinge_R, rbe)
    i_Rrbe_s, _ = _find_ix_range(rbe, x1)
    _, i_Lrbe_e = _find_ix_range(x0, -rbe)

    def _elements_in_ix_range(ix_start, ix_end, layer_filter=None):
        """Return list of connectivity for elements in given ix range."""
        result = []
        for iz in range(nz):
            for iy in range(ny_tot):
                li = min(iy // ny_per, nlyr - 1)
                if layer_filter is not None and li != layer_filter:
                    continue
                for ix in range(ix_start, ix_end):
                    conn = [
                        nid(ix, iy, iz), nid(ix + 1, iy, iz),
                        nid(ix + 1, iy + 1, iz), nid(ix, iy + 1, iz),
                        nid(ix, iy, iz + 1), nid(ix + 1, iy, iz + 1),
                        nid(ix + 1, iy + 1, iz + 1), nid(ix, iy + 1, iz + 1),
                    ]
                    result.append(conn)
        return result

    # Deformable: left wing + hinge + right wing
    elems = {}
    for li in range(nlyr):
        lname = f"L{li+1}"
        wing_elems = _elements_in_ix_range(i_Lwing_s, i_Lwing_e, layer_filter=li)
        hinge_elems = _elements_in_ix_range(i_hinge_s, i_hinge_e, layer_filter=li)
        wing_elems += _elements_in_ix_range(i_Rwing_s, i_Rwing_e, layer_filter=li)
        elems[lname] = wing_elems + hinge_elems

    # RBE regions (all layers, rigid body)
    elems['LEFT_RBE'] = _elements_in_ix_range(0, i_Lrbe_e) if i_Lrbe_e > 0 else []
    elems['RIGHT_RBE'] = _elements_in_ix_range(i_Rrbe_s, nx_tot) if i_Rrbe_s < nx_tot else []

    # Node sets
    eps = 1e-12
    NBOT, NTOP = [], []
    for iz in range(nnz):
        for iy in range(nny):
            for ix in range(nnx):
                nd = nid(ix, iy, iz)
                y = Y_nodes[min(iy, ny_tot)]
                if abs(y) < eps:
                    NBOT.append(nd)
                if abs(y - nlyr * t_lyr) < eps:
                    NTOP.append(nd)

    node_sets = {
        'ALL':     sorted(range(1, len(nodes) + 1)),
        'BOTTOM':  sorted(set(NBOT)),
        'TOP':     sorted(set(NTOP)),
    }

    return nodes, elems, node_sets


def write(p: dict, nodes: np.ndarray, elems: dict, node_sets: dict,
          filename: str) -> None:
    """Write FEBio .feb file with rigid-body hinge control."""
    MID_A, MID_B = 1, 2
    MID_LW, MID_RW = 3, 4  # left/right rigid wing bodies

    hinge_L = p['hinge_L']
    hinge_R = p['hinge_R']
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
    L('      <qn_method type="1"><max_ups>10</max_ups></qn_method>')
    L('    </solver>')
    L('    <time_stepper>')
    L('      <dtmin>1e-8</dtmin><dtmax>0.01</dtmax>')
    L('      <max_retries>15</max_retries><opt_iter>5</opt_iter>')
    L('    </time_stepper>')
    L('  </Control>')

    # Materials
    L('  <Material>')
    # A: neo-Hookean
    L(f'    <material id="1" name="A_layer" type="neo-Hookean">')
    L(f'      <c1>{AC1:.6f}</c1><k>{AK:.6f}</k>')
    L('    </material>')
    # B: viscoelastic
    L(f'    <material id="2" name="B_PSA" type="viscoelastic">')
    L('      <elastic type="neo-Hookean">')
    L(f'        <c1>{BC1:.6f}</c1><k>{BK:.6f}</k>')
    L('      </elastic>')
    L('      <relaxation>')
    L(f'        <g1>{FEBIO_VISCO_PARAMS["B_g1"]}</g1>'
      f' <t1>{FEBIO_VISCO_PARAMS["B_t1"]}</t1>')
    L(f'        <g2>{FEBIO_VISCO_PARAMS["B_g2"]}</g2>'
      f' <t2>{FEBIO_VISCO_PARAMS["B_t2"]}</t2>')
    L('      </relaxation>')
    L('    </material>')
    # Rigid body wings
    L(f'    <material id="3" name="LeftAssy" type="rigid body">')
    L('      <density>1e-9</density>')
    L(f'      <center_of_mass>{hinge_L}, 0, 0</center_of_mass>')
    L('    </material>')
    L(f'    <material id="4" name="RightAssy" type="rigid body">')
    L('      <density>1e-9</density>')
    L(f'      <center_of_mass>{hinge_R}, 0, 0</center_of_mass>')
    L('    </material>')
    L('  </Material>')

    # Geometry
    L('  <Geometry>')
    L('    <Nodes name="all">')
    for i, (x, y, z) in enumerate(nodes):
        L(f'      <node id="{i+1}">{x:.8f}, {y:.8f}, {z:.8f}</node>')
    L('    </Nodes>')

    # Deformable layers (A / B)
    mat_of = {'A': MID_A, 'B': MID_B}
    for li in range(p['n_layers']):
        lname = f"L{li+1}"
        mtype = 'A' if li % 2 == 0 else 'B'
        el = elems.get(lname, [])
        if not el:
            continue
        mid = mat_of[mtype]
        L(f'    <Elements type="hex8" mat="{mid}" elset="{lname}">')
        for eid, conn in enumerate(el, 1):
            L(f'      <elem id="{eid}"> {",".join(str(n) for n in conn)} </elem>')
        L('    </Elements>')

    # Rigid RBE regions
    for part, mid in [('LEFT_RBE', 3), ('RIGHT_RBE', 4)]:
        el = elems.get(part, [])
        if not el:
            continue
        L(f'    <Elements type="hex8" mat="{mid}" elset="{part}">')
        for eid, conn in enumerate(el, 1):
            L(f'      <elem id="{eid}"> {",".join(str(n) for n in conn)} </elem>')
        L('    </Elements>')

    # Node sets
    for name in ['ALL', 'BOTTOM', 'TOP']:
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
    # Left: Rz = +angle (fold inward)
    L('    <rigid_bc type="rigid_rotation">')
    L('      <rb>3</rb><dof>Rz</dof>')
    L(f'      <value lc="1">{angle:.6f}</value>')
    L('    </rigid_bc>')
    # Right: Rz = -angle (fold inward)
    L('    <rigid_bc type="rigid_rotation">')
    L('      <rb>4</rb><dof>Rz</dof>')
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
    print(f"     Hinge rotation: L=+{p['fold_angle']:.0f}deg,"
          f" R=-{p['fold_angle']:.0f}deg")
