"""
Mesh generation for foldable display (X=0 symmetric).

Geometry:
  Total width: 120mm (-60 .. +60), Y total = n_layers * layer_thick
  Hinges:  X = -15mm, +15mm  (reference points, Y=0)
  RBE:     bottom-surface nodes with |X| > 35mm

Two backends:
  [1] GMSH (structured quad, transfinite)
  [2] Pure NumPy structured quad

Output:
  nodes:    ndarray (N, 2)  — XY coordinates
  elems:    list of (layer_idx, [n1, n2, n3, n4])
  node_sets: {name: [nid, ...]}
    - BOTTOM       : all Y=0 nodes
    - BOTTOM_L     : Y=0 nodes with X < -rbe_region  (left RBE slaves)
    - BOTTOM_R     : Y=0 nodes with X > +rbe_region  (right RBE slaves)
    - TOP          : Y = total_thickness
    - ALL          : all node IDs
"""

from typing import Optional
import math
import numpy as np

try:
    import gmsh
    _HAS_GMSH = True
except ImportError:
    _HAS_GMSH = False


def build_mesh(p: dict, use_gmsh: Optional[bool] = None) -> tuple:
    """
    Build structured 2D quad mesh for the full display cross-section.

    Parameters
    ----------
    p : dict
        Keys: total_width, hinge_L, hinge_R, rbe_region, layer_thick,
              n_layers, ny_per_layer, nx_total.
    use_gmsh : bool or None
        True=force GMSH, False=force Python, None=auto-detect.

    Returns
    -------
    nodes, elems, node_sets
    """
    if use_gmsh is None:
        use_gmsh = _HAS_GMSH
    if use_gmsh and _HAS_GMSH:
        return _build_mesh_gmsh(p)
    return _build_mesh_python(p)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _extract_node_sets(nodes, p):
    """Build node sets from coordinate inspection."""
    x0 = -p['total_width'] / 2.0
    x1 = +p['total_width'] / 2.0
    rbe = p['rbe_region']
    y_tot = p['n_layers'] * p['layer_thick']
    eps = 1e-12

    bottom, bottom_l, bottom_r, top = [], [], [], []
    for nid, (x, y) in enumerate(nodes, 1):
        if abs(y) < eps:
            bottom.append(nid)
            if x < -rbe + eps:
                bottom_l.append(nid)
            elif x > rbe - eps:
                bottom_r.append(nid)
        if abs(y - y_tot) < eps:
            top.append(nid)

    return {
        'ALL':      list(range(1, len(nodes) + 1)),
        'BOTTOM':   sorted(set(bottom)),
        'BOTTOM_L': sorted(set(bottom_l)),
        'BOTTOM_R': sorted(set(bottom_r)),
        'TOP':      sorted(set(top)),
    }


# ---------------------------------------------------------------------------
# GMSH backend
# ---------------------------------------------------------------------------

def _build_mesh_gmsh(p: dict) -> tuple:
    import gmsh

    x0 = -p['total_width'] / 2.0
    x1 = +p['total_width'] / 2.0
    t = p['layer_thick']
    nlyr = p['n_layers']
    ny_per = p['ny_per_layer']
    nx = p['nx_total']
    y_total = nlyr * t

    gmsh.initialize()
    gmsh.model.add("display_folding")

    rect = gmsh.model.occ.addRectangle(x0, 0, 0, x1 - x0, y_total)
    gmsh.model.occ.synchronize()

    surfaces = gmsh.model.getEntities(2)
    surf_tag = surfaces[0][1]
    curves = [c[1] for c in gmsh.model.getBoundary([(2, surf_tag)], oriented=False)]

    y_mid = y_total / 2.0
    x_mid = 0.0

    curve_bottom = curve_right = curve_top = curve_left = None
    for c in curves:
        xmin, ymin, _, xmax, ymax, _ = gmsh.model.getBoundingBox(1, c)
        dx, dy = xmax - xmin, ymax - ymin
        if abs(dx) > abs(dy):
            if ymax < y_mid:
                curve_bottom = c
            else:
                curve_top = c
        else:
            if xmax < x_mid:
                curve_left = c
            else:
                curve_right = c

    ny = nlyr * ny_per
    gmsh.model.mesh.setTransfiniteSurface(surf_tag)
    gmsh.model.mesh.setRecombine(2, surf_tag)

    if curve_bottom is not None:
        gmsh.model.mesh.setTransfiniteCurve(curve_bottom, nx + 1)
    if curve_top is not None:
        gmsh.model.mesh.setTransfiniteCurve(curve_top, nx + 1)
    if curve_left is not None:
        gmsh.model.mesh.setTransfiniteCurve(curve_left, ny + 1)
    if curve_right is not None:
        gmsh.model.mesh.setTransfiniteCurve(curve_right, ny + 1)

    gmsh.model.mesh.generate(2)

    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    node_map = {tag: i for i, tag in enumerate(node_tags)}
    nodes = np.column_stack([node_coords[0::3], node_coords[1::3]])

    elem_types, elem_tags, elem_conn = gmsh.model.mesh.getElements(dim=2)
    elems = []
    for etype, etags, econn in zip(elem_types, elem_tags, elem_conn):
        if etype not in (3, 16):
            continue
        npe = 4 if etype == 3 else 9
        for i in range(len(etags)):
            conn = list(econn[i * npe:(i + 1) * npe])
            cy = np.mean([nodes[node_map[n] - 1][1] for n in conn[:4]])
            layer_idx = int(cy // t)
            layer_idx = min(layer_idx, nlyr - 1)
            elems.append((layer_idx, conn[:4]))

    # Re-index nodes to 1..N
    nid_map = {old: new for new, old in enumerate(sorted(node_tags), 1)}
    old_to_new = {old: i + 1 for i, old in enumerate(sorted(node_tags))}

    # Remap node coordinates
    nodes_sorted = np.zeros((len(node_tags), 2))
    for old_tag, (x, y) in zip(sorted(node_tags), [nodes[node_map[t]] for t in sorted(node_tags)]):
        idx = old_to_new[old_tag] - 1
        nodes_sorted[idx] = [x, y]
    nodes = nodes_sorted

    # Remap element connectivity
    elems = [(li, [old_to_new[n] for n in conn]) for li, conn in elems]

    gmsh.finalize()

    # Build node sets from re-indexed nodes
    node_sets = _extract_node_sets(nodes, p)
    return nodes, elems, node_sets


# ---------------------------------------------------------------------------
# Pure Python backend
# ---------------------------------------------------------------------------

def _build_mesh_python(p: dict) -> tuple:
    x0 = -p['total_width'] / 2.0
    x1 = +p['total_width'] / 2.0
    t = p['layer_thick']
    nlyr = p['n_layers']
    ny_per = p['ny_per_layer']
    nx = p['nx_total']

    ny = nlyr * ny_per
    x_edges = np.linspace(x0, x1, nx + 1)

    y_edges = [0.0]
    for li in range(nlyr):
        for _ in range(ny_per):
            y_edges.append(y_edges[-1] + t / ny_per)
    y_edges = np.array(y_edges)

    nnx, nny = nx + 1, ny + 1

    def nid(ix, iy):
        return iy * nnx + ix + 1

    nodes = np.zeros((nnx * nny, 2))
    for iy in range(nny):
        for ix in range(nnx):
            nodes[iy * nnx + ix] = [x_edges[ix], y_edges[iy]]

    elems = []
    for iy in range(ny):
        layer_idx = min(iy // ny_per, nlyr - 1)
        for ix in range(nx):
            conn = [nid(ix, iy), nid(ix + 1, iy),
                    nid(ix + 1, iy + 1), nid(ix, iy + 1)]
            elems.append((layer_idx, conn))

    node_sets = _extract_node_sets(nodes, p)
    return nodes, elems, node_sets


# ---------------------------------------------------------------------------
# 3D Extrusion
# ---------------------------------------------------------------------------

def extrude_to_3d(nodes_2d: np.ndarray,
                  elems_2d: list,
                  node_sets_2d: dict,
                  depth: float) -> tuple:
    """
    Extrude 2D quad mesh to 3D hex mesh by replicating in Z.

    Returns
    -------
    nodes_3d : (2N, 3) XYZ (bottom then top)
    elems_3d : list of (layer_idx, [h1..h8])
    node_sets_3d : same keys, IDs extend to both Z layers
    """
    n2 = len(nodes_2d)
    # Coordinate system: X=width, Y=depth(extrusion), Z=thickness(layer-stacking)
    # 2D mesh: col0=X(width), col1=Z(thickness)
    xyz_front = np.column_stack([nodes_2d[:, 0], np.zeros(n2),         nodes_2d[:, 1]])
    xyz_back  = np.column_stack([nodes_2d[:, 0], np.full(n2, depth),   nodes_2d[:, 1]])
    nodes_3d = np.vstack([xyz_front, xyz_back])

    elems_3d = []
    for layer_idx, qconn in elems_2d:
        n1, n2id, n3, n4 = qconn
        # Swap Y=depth face to nodes 1-4 so local-ζ points in -Y direction:
        # det(J)=Δx*Δy*Δz > 0  (was negative when front face was nodes 1-4)
        h = [n1 + n2, n2id + n2, n3 + n2, n4 + n2, n1, n2id, n3, n4]
        elems_3d.append((layer_idx, h))

    node_sets_3d = {}
    for name, ids in node_sets_2d.items():
        combined = sorted(set(ids) | set(nid + n2 for nid in ids))
        node_sets_3d[name] = combined

    return nodes_3d, elems_3d, node_sets_3d


# ---------------------------------------------------------------------------
# Hinge cylinder mesh generation
# ---------------------------------------------------------------------------

def build_hinge_cylinders(nodes_3d: np.ndarray,
                          elems_3d: list,
                          node_sets_3d: dict,
                          p: dict,
                          depth: float) -> tuple:
    """
    Build hollow cylinder meshes at hinge positions and append to the
    3D display mesh.

    Each hinge is a hollow cylinder (OD x ID) centered at (hinge_x, 0, 0)
    extruded in Z.  The cylinder mesh is appended after the display nodes.

    Parameters
    ----------
    nodes_3d : (N, 3) ndarray — existing extruded display nodes
    elems_3d : list of (layer_idx, [8-node connectivity]) — existing hexes
    node_sets_3d : dict
    p : dict with keys: hinge_L, hinge_R, hinge_cylinder_od,
        hinge_cylinder_id, hinge_cylinder_ntheta
    depth : float — Z-depth of the cylinder (= display depth)

    Returns
    -------
    nodes_3d, elems_3d, node_sets_3d — modified with cylinder geometry
    """
    n_display_nodes = len(nodes_3d)
    n_display_elems = len(elems_3d)
    od = p['hinge_cylinder_od']   # outer diameter [mm]
    id_ = p['hinge_cylinder_id']  # inner diameter [mm]
    ro = od / 2.0
    ri = id_ / 2.0
    ntheta = p['hinge_cylinder_ntheta']
    eps = 1e-12

    hinge_positions = [
        ('L', p['hinge_L']),
        ('R', p['hinge_R']),
    ]

    # In new coords: X=width, Y=depth, Z=thickness(Z=0 = bottom/hinge surface)
    # Bottom-surface nodes: Z=0 (nodes_3d[:,2]); front face: Y=0 (nodes_3d[:,1])
    bot_ids_L = [nid for nid in range(1, n_display_nodes + 1)
                 if abs(nodes_3d[nid - 1, 2]) < eps
                 and abs(nodes_3d[nid - 1, 1]) < eps
                 and abs(nodes_3d[nid - 1, 0] - p['hinge_L']) < ro - eps]
    bot_ids_R = [nid for nid in range(1, n_display_nodes + 1)
                 if abs(nodes_3d[nid - 1, 2]) < eps
                 and abs(nodes_3d[nid - 1, 1]) < eps
                 and abs(nodes_3d[nid - 1, 0] - p['hinge_R']) < ro - eps]

    node_sets_3d['HINGE_BOTTOM_L'] = bot_ids_L
    node_sets_3d['HINGE_BOTTOM_R'] = bot_ids_R

    node_sets_3d['HINGE_CYL_L'] = []
    node_sets_3d['HINGE_CYL_R'] = []
    elset_L = []
    elset_R = []

    next_nid = n_display_nodes + 1
    next_eid = n_display_elems + 1

    for side_key, hinge_x in hinge_positions:
        # Node layout per cylinder (local indexed as below):
        #   idx 0..ntheta-1          : inner ring, Z=0
        #   idx ntheta..2*ntheta-1   : outer ring, Z=0
        #   idx 2*ntheta..3*ntheta-1 : inner ring, Z=depth
        #   idx 3*ntheta..4*ntheta-1 : outer ring, Z=depth
        # Cylinder axis along Y; rings in XZ plane, center at (hinge_x, Z=-ro)
        # → cylinder top (Z=0) is tangent to display bottom surface.
        cyl_z_center = -ro
        cyl_nids = []
        for iy in range(2):
            y = 0.0 if iy == 0 else depth
            for radius in (ri, ro):
                for i_th in range(ntheta):
                    th = 2.0 * math.pi * i_th / ntheta
                    x = hinge_x + radius * math.cos(th)
                    z = cyl_z_center + radius * math.sin(th)
                    nodes_3d = np.vstack([nodes_3d, [x, y, z]])
                    cyl_nids.append(next_nid)
                    next_nid += 1

        # Hex connectivity wraps around the ring (i_th -> i_next = i_th+1 mod ntheta)
        for i_th in range(ntheta):
            i_next = (i_th + 1) % ntheta
            bi = cyl_nids[i_th]
            bj = cyl_nids[i_next]
            bo = cyl_nids[ntheta + i_next]
            bo_i = cyl_nids[ntheta + i_th]
            ti = cyl_nids[2 * ntheta + i_th]
            tj = cyl_nids[2 * ntheta + i_next]
            to = cyl_nids[3 * ntheta + i_next]
            to_i = cyl_nids[3 * ntheta + i_th]

            # Swap Y=depth face to nodes 1-4 (same sign fix as display hexes)
            conn = [ti, to_i, to, tj, bi, bo_i, bo, bj]
            if side_key == 'L':
                elset_L.append(next_eid)
                elems_3d.append((-1, conn))
            else:
                elset_R.append(next_eid)
                elems_3d.append((-2, conn))
            next_eid += 1

    node_sets_3d['HINGE_CYL_L'] = sorted(elset_L)
    node_sets_3d['HINGE_CYL_R'] = sorted(elset_R)

    print(f"     Hinge cylinders: {ntheta} elems each (@ {p['hinge_L']:.0f}, {p['hinge_R']:.0f})")
    print(f"     Hinge RBE slaves: L={len(bot_ids_L)} R={len(bot_ids_R)} "
          f"(display bottom y=0 nodes inside cylinder OD)")

    return nodes_3d, elems_3d, node_sets_3d
