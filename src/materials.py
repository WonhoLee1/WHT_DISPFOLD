"""
Neo-Hookean material constants for multi-layer foldable display.

Layer stack (bottom to top, A-B-A-B-A):
  A: neo-Hookean (E=2000MPa, nu=0.34) - stiff structural layer
  B: neo-Hookean (E=1-50MPa, nu=0.45) - soft PSA layer

Both layers use the same Neo-Hookean form:
  C10 = E / (4 * (1 + nu))
  D1  = 6 * (1 - 2*nu) / E

Geometry (X=0 symmetric):
  Total width: 120mm (-60 .. +60)
  Hinges:      X = -15mm, +15mm (reference points at Y=0)
  RBE region:  X < -35mm (left), X > +35mm (right) — bottom surface nodes
"""

from typing import Tuple, List


def nh_c10(E: float, nu: float) -> float:
    """Neo-Hookean C10 constant from Young's modulus and Poisson ratio."""
    return E / (4.0 * (1.0 + nu))


def nh_d1(E: float, nu: float) -> float:
    """Neo-Hookean D1 (incompressibility) constant."""
    return 6.0 * (1.0 - 2.0 * nu) / E


def layer_stack(n_layers: int) -> List[Tuple[str, str]]:
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
    # Geometry [mm] (X=0 symmetric)
    # Coordinate axes: X=width(-80..+80), Y=depth(0..1mm, sym at Y=0), Z=thickness(0..0.15mm)
    'total_width':    160.0,     # full display width [mm] (-80 to +80)
    'hinge_L':       -10.0,      # left hinge X position [mm] (cylinder axis at Z=0)
    'hinge_R':        10.0,      # right hinge X position [mm]
    'rbe_region':     10.0,      # |X| > 10mm bottom-surface nodes get rotation BCs
    'depth':           1.0,      # Y-depth [mm] (symmetry at Y=0)

    # Layer stack (Z direction: Z=0 bottom/hinge surface, Z=0.15mm top/orange)
    'layer_thick':     0.030,    # each layer thickness [mm] (30 um)
    'n_layers':        5,        # A-B-A-B-A
    'ny_per_layer':    3,        # elements per layer through-thickness

    # Mesh division
    'nx_total':      1600,       # total X-direction elements (160mm @ 10/mm)

    # Hinge cylinder geometry (hollow, axis along Y, center at X=hinge, Z=0)
    'hinge_cylinder_od':     10.0,   # outer diameter [mm]
    'hinge_cylinder_id':      4.0,   # inner diameter [mm]
    'hinge_cylinder_ntheta':  16,    # circumferential elements
    'hinge_cylinder_depth':    1.0,  # Y-depth [mm] (same as display depth)

    # Material A: stiff structural layer (Neo-Hookean)
    'A_E':           2000.0,     # Young's modulus [MPa]
    'A_nu':            0.34,

    # Material B: soft PSA layer (Neo-Hookean)
    'B_E':             50.0,     # Young's modulus [MPa]
    'B_nu':             0.45,

    # Folding
    'fold_angle':      90.0,     # [deg] per wing
    'fold_time':        1.0,     # [s] folding duration
    'nsteps':         200,       # solver increments
}

# Master node IDs for RBE (reserved high IDs to avoid mesh collision)
RBE_MASTER_IDS = {
    'L': 999997,      # left hinge master node ID
    'R': 999998,      # right hinge master node ID
}

# FEBio-specific viscoelastic parameters
FEBIO_VISCO_PARAMS = {
    'B_g1':  0.5,     # Prony g1
    'B_t1':  0.1,     # Prony tau1 [s]
    'B_g2':  0.3,     # Prony g2
    'B_t2':  2.0,     # Prony tau2 [s]
    'nz':     1,      # Z-division
}
