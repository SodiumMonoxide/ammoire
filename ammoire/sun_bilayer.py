"""Representative no-SOC honeycomb bilayer for Sun et al., Supplement Eq. (10).

Each sheet has A=(0,0), B=(0,-a/√3) with Bravais vectors a(1,0) and
a(1/2,√3/2). Basis order is bottom A,B then top A,B, each with physical
spin. The top sheet is shifted by a/3 at angle stacking × 60 degrees.
The relative shift and shortest interlayer bond stencil are explicit model
choices informed by the schematic in Fig. S1. Exchange magnitudes and
hopping ratios use its illustrative values. No polarization is calculated.

The moment signs (+,−,−,+) multiply the common Néel direction in the
reference domain. `layer_sign` tracks this pattern when the real-space
operation permutes sheets or sublattices in the common bottom-to-top basis.
Up to a change of origin, Mz maps the reference (stacking=0, layer_sign=+1)
to (3,−1), Mz My to (3,+1), and Mz C3z to (5,−1). My reverses y, C3z is
an active 120-degree rotation. All act trivially on transported spins.
"""

import numpy as np


SOURCE = "https://doi.org/10.1038/s41467-026-69635-2"
# Cell translations to B sites nearest A, in Cartesian units of a.
_AB_CELLS = np.array([[0.0, 0.0], [0.5, np.sqrt(3) / 2], [-0.5, np.sqrt(3) / 2]])
_AB_BONDS = _AB_CELLS + [0.0, -1 / np.sqrt(3)]


def evaluate(q, p, spin):
    """Return the eight-band cell-gauge Hamiltonian at Cartesian q=ak."""
    stacking = p['stacking']
    if stacking not in range(6) or p['layer_sign'] not in (-1, 1):
        raise ValueError('stacking must be an integer from 0 to 5 and layer_sign must be ±1')
    angle = stacking * np.pi / 3
    shift = np.array([np.cos(angle), np.sin(angle)]) / 3

    # Equal-sublattice bonds join the same cell. The two cross-sublattice
    # neighbours are chosen in real space before any k-dependent calculation.
    ab_cell = _AB_CELLS[np.argmin(np.linalg.norm(shift + _AB_BONDS, axis=1))]
    ba_cell = -_AB_CELLS[np.argmin(np.linalg.norm(shift - _AB_BONDS, axis=1))]
    intra = p['t'] * np.exp(1j * (q @ _AB_CELLS.T)).sum(axis=-1)
    hop = np.zeros(q.shape[:-1] + (4, 4), dtype=complex)
    hop[..., 0, 1] = hop[..., 2, 3] = intra
    hop[..., 0, 2] = hop[..., 1, 3] = p['t_inter']
    hop[..., 0, 3] = p['t_inter'] * np.exp(1j * (q @ ab_cell))
    hop[..., 1, 2] = p['t_inter'] * np.exp(1j * (q @ ba_cell))
    hop += hop.conj().swapaxes(-2, -1)
    moments = p['layer_sign'] * np.diag([1, -1, -1, 1])
    return (np.kron(hop, np.eye(2)) - p['j'] * np.kron(moments, spin)
            - p['mu'] * np.eye(8))
