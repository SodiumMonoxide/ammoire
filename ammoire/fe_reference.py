"""Domain operations of Sun et al. (doi:10.1038/s41467-026-69635-2).

For bilayer MnPS3 the z axis is the layer normal and the xz plane is
the initial mirror plane. Bulk BiFeO3 uses the pseudocubic Cartesian axes.
"""

import numpy as np

from .fe_symmetry import SpinLaueGroup


def _case(group, material, order, path, initial, final, label, operation):
    return {
        "material": material,
        "order": order,
        "path": path,
        "initial_state": initial,
        "final_state": final,
        "operation": label,
        "operation_cartesian": operation.tolist(),
        "momentum_action": group.momentum_action(operation).tolist(),
        "response_type": group.classify(operation),
    }


def literature_cases():
    """Return the six MnPS3 and two BiFeO3 switching-path classifications.

    The real-space operations in Eqs. (11)–(14), Supplementary Fig. S7,
    and Fig. 5 act on Cartesian columns, with the rightmost factor first.
    The three bilayer paths start from P1+ and are evaluated separately
    for each interlayer magnetic order. The bulk paths differ in whether
    the octahedral rotations switch together with the [111] polarization.
    """
    mirror_y = np.diag([1, -1, 1])
    mirror_z = np.diag([1, 1, -1])
    s = np.sqrt(3) / 2
    c3z = np.array([[-0.5, -s, 0], [s, -0.5, 0], [0, 0, 1]])
    layer = SpinLaueGroup(2, [(mirror_y, -1)])
    cases = []
    for order, label, base in (("L1", "Mz", mirror_z), ("L2", "Mz My", mirror_z @ mirror_y)):
        paths = ((label, base), (label + " C3z", base @ c3z), (label + " C3z^2", base @ c3z @ c3z))
        for path, (operation_label, operation) in enumerate(paths, 1):
            cases.append(_case(layer, "MnPS3", order, path, "P1+", f"P{path}-", operation_label, operation))

    # x → y → z about [111]. The twofold about [-110] interchanges
    # the magnetic sublattices and sends the polarization to -[111].
    c3_111 = np.eye(3)[[2, 0, 1]]
    c2_m110 = np.array([[0, -1, 0], [-1, 0, 0], [0, 0, -1]])
    bulk = SpinLaueGroup(3, [(c3_111, 1), (c2_m110, -1)])
    for path, (label, operation) in enumerate((("I", -np.eye(3)), ("C2[-110]", c2_m110)), 1):
        cases.append(_case(bulk, "BiFeO3", "G-type", path, "+P[111]", "-P[111]", label, operation))
    return cases


def transform_conductivity(tensors, operation):
    """Map spin-resolved charge conductivities between fixed-order domains.

    Supply (2,d,d) tensors ordered up/down, with j_s = sigma_s @ E for
    Cartesian columns. For an active orthogonal domain map O (d×d),
    the result is O @ sigma_s @ O.T for each conserved spin s. Use the
    in-plane action for 2D transport. Units and spin labels are preserved.
    The input tensors may be complex and need not be symmetric.
    """
    values = np.asarray(tensors)
    if values.shape not in ((2, 2, 2), (2, 3, 3)):
        raise ValueError("tensors must have shape (2,d,d) for d=2 or 3")
    if not np.issubdtype(values.dtype, np.number) or not np.isfinite(values).all():
        raise ValueError("tensors must contain finite numerical values")
    if np.iscomplexobj(operation):
        raise ValueError("the Cartesian operation must be real")
    matrix = np.asarray(operation, dtype=float)
    if matrix.shape != values.shape[1:] or not np.isfinite(matrix).all():
        raise ValueError("the Cartesian operation must be finite and match the tensor dimension")
    if not np.allclose(matrix.T @ matrix, np.eye(len(matrix)), atol=1e-8, rtol=0):
        raise ValueError("the Cartesian operation must be orthogonal")
    return matrix @ values @ matrix.T
