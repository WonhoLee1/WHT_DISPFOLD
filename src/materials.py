"""
Neo-Hookean material constants for multi-layer foldable display.

Layer stack (bottom to top, A-B-A-B-A):
  A: neo-Hookean (E=2000MPa, nu=0.34) - stiff structural layer
  B: neo-Hookean (E=1-50MPa, nu=0.45) - soft PSA layer

Both layers use the same Neo-Hookean form:
  C10 = E / (4 * (1 + nu))
  D1  = 6 * (1 - 2*nu) / E
"""

from typing import Tuple


def nh_c10(E: float, nu: float) -> float:
    """Neo-Hookean C10 constant from Young's modulus and Poisson ratio."""
    return E / (4.0 * (1.0 + nu))


def nh_d1(E: float, nu: float) -> float:
    """Neo-Hookean D1 (incompressibility) constant."""
    return 6.0 * (1.0 - 2.0 * nu) / E


def layer_stack(n_layers: int) -> list[Tuple[str, str]]:
    """
    Return list of (name, type) pairs in stacking order (bottom to top).

    Example for n_layers=5:
      [('L1', 'A'), ('L2', 'B'), ('L3', 'A'), ('L4', 'B'), ('L5', 'A')]
    """
    return [(f"L{i+1}", 'A' if i % 2 == 0 else 'B')
            for i in range(n_layers)]


def layer_mid(layer_idx: int, mid_a: int = 1, mid_b: int = 2) -> int:
    """Return material ID for a given layer index (0-based)."""
    return mid_a if layer_idx % 2 == 0 else mid_b


DEFAULT_PARAMS = {
    # Geometry [mm]
    'hinge_L':       35.0,      # left hinge X-coordinate
    'hinge_R':       65.0,      # right hinge X-coordinate
    'depth':          1.0,      # Z-depth for 3D [mm]

    # Layer stack
    'layer_thick':    0.030,    # each layer thickness [mm] (30 um)
    'n_layers':       5,        # A-B-A-B-A

    # Mesh X-division (hinge region only)
    'nx_hinge':      600,       # elements in hinge region

    # Bandwidth: prescribe displacement for a band of nodes near each hinge edge
    'hinge_band_elements': 5,

    # Material A: stiff structural layer (Neo-Hookean)
    'A_E':           2000.0,    # Young's modulus [MPa]
    'A_nu':           0.34,

    # Material B: soft PSA layer (Neo-Hookean)
    'B_E':            50.0,     # Young's modulus [MPa]
    'B_nu':            0.45,

    # Folding
    'fold_angle':     90.0,     # [deg] per wing
    'fold_time':       1.0,     # [s] folding duration
    'nsteps':        200,       # solver increments
}

# FEBio-specific viscoelastic parameters (for febio solver only)
FEBIO_VISCO_PARAMS = {
    'B_g1':  0.5,     # Prony g1
    'B_t1':  0.1,     # Prony tau1 [s]
    'B_g2':  0.3,     # Prony g2
    'B_t2':  2.0,     # Prony tau2 [s]
    'nx_left':   28,  # left wing elements (0~35mm)
    'nx_hinge':  60,  # hinge elements (35~65mm)
    'nx_right':  28,  # right wing elements (65~100mm)
    'ny_per_layer': 2,
    'nz':           1,
}
