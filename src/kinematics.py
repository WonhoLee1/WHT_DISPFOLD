"""
Folding rotation displacement calculator.

Computes the prescribed displacement for nodes rotating about a hinge axis,
with optional linear taper over a bandwidth to prevent edge-localized shear.
"""

import math
import numpy as np


def compute_rotation_displacement(
    p: dict,
    nodes: np.ndarray,
    node_ids: list[int],
    hinge_x: float,
    hinge_y: float,
    angle_sign: int,
    band_width: float = 0.0,
) -> list[tuple[int, int, float]]:
    """
    Compute displacement for nodes rotating about (hinge_x, hinge_y).

    If ``band_width > 0``, the rotation angle varies linearly from theta_max
    at the hinge edge to zero at ``hinge_x +/- band_width``.

    Parameters
    ----------
    p : dict
        Must contain ``fold_angle`` [deg].
    nodes : (N, 2) ndarray
        Node XY coordinates.
    node_ids : list[int]
        Node IDs (1-based) to prescribe.
    hinge_x, hinge_y : float
        Hinge axis coordinate.
    angle_sign : int
        +1 for CCW, -1 for CW.
    band_width : float
        Linear taper zone [mm]; 0 = constant theta.

    Returns
    -------
    list of (nid, dof, value)
        dof=1 (UX), dof=2 (UY). Only non-zero entries included.
    """
    theta_max = math.radians(p['fold_angle']) * angle_sign
    entries = []

    for nid in node_ids:
        x0 = nodes[nid - 1, 0]
        y0 = nodes[nid - 1, 1]
        dx = x0 - hinge_x
        dy = y0 - hinge_y

        if band_width > 1e-15:
            xi = abs(dx) / band_width
            alpha = 1.0 - min(1.0, xi)
            theta = theta_max * alpha
        else:
            theta = theta_max

        cos_t = math.cos(theta)
        sin_t = math.sin(theta)

        ux = dx * (cos_t - 1.0) - dy * sin_t
        uy = dx * sin_t + dy * (cos_t - 1.0)

        if abs(ux) > 1e-15:
            entries.append((nid, 1, ux))
        if abs(uy) > 1e-15:
            entries.append((nid, 2, uy))

    return entries
