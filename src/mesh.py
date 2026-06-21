"""
Mesh generation for foldable display hinge region.

Two backends:
  [1] GMSH (structured quad, transfinite) — requires ``gmsh`` package
  [2] Pure NumPy structured quad — no external dependency

Output format (both backends):
  nodes:    ndarray (N, 2)  — XY coordinates
  elems:    list of (layer_idx, [n1, n2, n3, n4])
  node_sets: {name: [nid, ...]}
    - EDGE_L      : nodes at X = hinge_L
    - EDGE_R      : nodes at X = hinge_R
    - EDGE_L_BAND : nodes near hinge_L (band_width)
    - EDGE_R_BAND : nodes near hinge_R (band_width)
    - BOTTOM      : nodes at Y = 0
    - TOP         : nodes at Y = total_thickness
"""

from typing import Sequence
import numpy as np

try:
    import gmsh
    _HAS_GMSH = True
except ImportError:
    _HAS_GMSH = False


def build_mesh(p: dict, use_gmsh: bool | None = None) -> tuple:
    """
    Build structured 2D quad mesh for the hinge region.

    Parameters
    ----------
    p : dict
        Parameter dictionary with keys: hinge_L, hinge_R, layer_thick,
        n_layers, ny_per_layer (default 3), nx_hinge, hinge_band_elements.
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
    else:
        return _build_mesh_python(p)


def _build_mesh_gmsh(p: dict) -> tuple:
    """GMSH structured quad mesh for hinge region."""
    import gmsh

    x0 = p['hinge_L']
    x1 = p['hinge_R']
    t = p['layer_thick']
    nlyr = p['n_layers']
    ny_per = p.get('ny_per_layer', 3)
    nx = p['nx_hinge']
    y_total = nlyr * t

    gmsh.initialize()
    gmsh.model.add("display_folding")

    # Full rectangle
    rect_tag = gmsh.model.occ.addRectangle(x0, 0, 0, x1 - x0, y_total)
    gmsh.model.occ.synchronize()

    surfaces = gmsh.model.getEntities(2)
    surf_tag = surfaces[0][1]

    # Identify curves by bounding box
    curves = gmsh.model.getBoundary([(2, surf_tag)], oriented=False)
    curves = [c[1] for c in curves]

    y_mid = y_total * 0.5
    x_mid = (x0 + x1) * 0.5

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

    # Extract node coordinates
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    node_map = {tag: i for i, tag in enumerate(node_tags)}
    nodes = np.column_stack([node_coords[0::3], node_coords[1::3]])

    # Extract elements (quad only)
    elem_types, elem_tags, elem_conn = gmsh.model.mesh.getElements(dim=2)
    elems = []
    for etype, etags, econn in zip(elem_types, elem_tags, elem_conn):
        if etype not in (3, 16):
            continue
        npe = 4 if etype == 3 else 9
        for i in range(len(etags)):
            conn = list(econn[i * npe:(i + 1) * npe])
            cy = np.mean([nodes[node_map[n] - 1][1] for n in conn[:4]])
            layer_idx = int(cy // p['layer_thick'])
            layer_idx = min(layer_idx, nlyr - 1)
            elems.append((layer_idx, conn[:4]))

    # Node sets
    eps = 1e-10
    edge_L, edge_R, bottom, top = [], [], [], []
    for tag, coord in zip(node_tags, nodes):
        x, y = coord[0], coord[1]
        if abs(x - x0) < eps:
            edge_L.append(tag)
        if abs(x - x1) < eps:
            edge_R.append(tag)
        if abs(y) < eps:
            bottom.append(tag)
        if abs(y - y_total) < eps:
            top.append(tag)

    band_x = p.get('hinge_band_elements', 5) * (x1 - x0) / nx

    def _band(tags, pts, x_min, x_max):
        ids = []
        for tag, pt in zip(tags, pts):
            if x_min - eps <= pt[0] <= x_max + eps:
                ids.append(tag)
        return sorted(set(ids))

    node_sets = {
        'EDGE_L':      sorted(set(edge_L)),
        'EDGE_R':      sorted(set(edge_R)),
        'EDGE_L_BAND': _band(node_tags, nodes, x0, x0 + band_x),
        'EDGE_R_BAND': _band(node_tags, nodes, x1 - band_x, x1),
        'BOTTOM':      sorted(set(bottom)),
        'TOP':         sorted(set(top)),
    }

    gmsh.finalize()
    return nodes, elems, node_sets


def _build_mesh_python(p: dict) -> tuple:
    """Pure NumPy structured quad mesh, no GMSH needed."""
    x0 = p['hinge_L']
    x1 = p['hinge_R']
    t = p['layer_thick']
    nlyr = p['n_layers']
    ny_per = p.get('ny_per_layer', 3)
    nx = p['nx_hinge']

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

    # Nodes
    nodes = np.zeros((nnx * nny, 2))
    for iy in range(nny):
        for ix in range(nnx):
            nodes[iy * nnx + ix] = [x_edges[ix], y_edges[iy]]

    # Elements
    elems = []
    for iy in range(ny):
        layer_idx = min(iy // ny_per, nlyr - 1)
        for ix in range(nx):
            conn = [nid(ix, iy), nid(ix + 1, iy),
                    nid(ix + 1, iy + 1), nid(ix, iy + 1)]
            elems.append((layer_idx, conn))

    # Node sets
    eps = 1e-10
    edge_L = [nid(0, iy) for iy in range(nny)]
    edge_R = [nid(nx, iy) for iy in range(nny)]
    bottom = [nid(ix, 0) for ix in range(nnx)]
    top = [nid(ix, ny) for ix in range(nnx)]

    n_band = max(1, p.get('hinge_band_elements', 5))
    n_band_x = min(nx, int(round(n_band)))
    L_band = [nid(ix, iy) for ix in range(n_band_x + 1) for iy in range(nny)]
    R_band = [nid(nx - ix, iy) for ix in range(n_band_x + 1) for iy in range(nny)]

    node_sets = {
        'EDGE_L':      edge_L,
        'EDGE_R':      edge_R,
        'EDGE_L_BAND': L_band,
        'EDGE_R_BAND': R_band,
        'BOTTOM':      bottom,
        'TOP':         top,
    }
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
    elems_3d : list of (layer_idx, [h1..h8]) — C3D8/CHEXA ordering
    node_sets_3d : same keys, IDs extend to both Z layers
    """
    n2 = len(nodes_2d)
    xyz_bot = np.column_stack([nodes_2d[:, 0], nodes_2d[:, 1], np.zeros(n2)])
    xyz_top = np.column_stack([nodes_2d[:, 0], nodes_2d[:, 1], np.full(n2, depth)])
    nodes_3d = np.vstack([xyz_bot, xyz_top])

    elems_3d = []
    for layer_idx, qconn in elems_2d:
        n1, n2id, n3, n4 = qconn
        h = [n1, n2id, n3, n4, n1 + n2, n2id + n2, n3 + n2, n4 + n2]
        elems_3d.append((layer_idx, h))

    node_sets_3d = {}
    for name, ids in node_sets_2d.items():
        combined = sorted(set(ids) | set(nid + n2 for nid in ids))
        node_sets_3d[name] = combined

    return nodes_3d, elems_3d, node_sets_3d


def nid_map_2d(p: dict) -> tuple:
    """Return (nnx, nny, nid_func) for the 2D Python mesh indexing scheme."""
    nx = p['nx_hinge']
    ny = p['n_layers'] * p.get('ny_per_layer', 3)
    nnx, nny = nx + 1, ny + 1
    return nnx, nny, lambda ix, iy: iy * nnx + ix + 1
