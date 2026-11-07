"""Cairo pentagonal electronic model, Li et al., npj QM 10, 83 (2025).

Source: https://arxiv.org/html/2412.16857v1, Eqs. (S5)-(S7).
Adapted from the 2026-11-07 local daily implementation, without its legacy
imports. Six orbitals (A,B,alpha,beta,gamma,eta), each followed by physical
spin up/down. k is dimensionless ka in this evaluator; the registry handles
the lattice length. Energies use the paper's |t| convention, not eV.

The removable d-dependent embedding is set to zero. By default return a
reciprocal-periodic cell gauge, U H_source U^dagger. Transform any imported
orbital-dependent observables or tunneling to this same gauge.
"""

import numpy as np


SOURCE = "https://doi.org/10.1038/s41535-025-00793-0"
# Gauge offsets for the d=0 source matrix; not relaxed atomic positions.
GAUGE_OFFSETS = np.array([[0, 0], [.5, .5], [.5, 0], [-.5, 0], [0, .5], [0, -.5]])
GAUGE_OFFSETS.setflags(write=False)


def cell_gauge(k):
    """Diagonal of U(k), duplicated over physical spin."""
    return np.repeat(np.exp(1j * (np.asarray(k) @ GAUGE_OFFSETS.T)), 2, axis=-1)


def evaluate(k, p, spin, *, periodic=True):
    """Registry evaluator. j is the signed JS exchange energy in the source.

    periodic=False exposes the d=0 matrix for formula/gauge cross-checks.
    Non-collinear uniform Neel directions are spin rotations of the no-SOC
    source model; the published quantization direction is (0,0,1).
    """
    kx, ky = np.moveaxis(k, -1, 0)
    x, y = np.exp(.5j * kx), np.exp(.5j * ky)
    plus, minus = 1 + p['delta'], 1 - p['delta']
    ab = -2 * p['tm'] * (plus * np.cos((kx + ky) / 2) + minus * np.cos((kx - ky) / 2))
    ag = -p['tn2'] * x.conj() * (plus * y.conj() + minus * y)
    ae = -p['tn2'] * y.conj() * (plus * x.conj() + minus * x)
    t = p['t']
    bonds = ((0, 1, ab), (0, 2, -t * minus * x), (0, 3, -t * minus * x.conj()),
             (0, 4, -t * plus * y), (0, 5, -t * plus * y.conj()),
             (1, 2, -t * minus * y), (1, 3, -t * minus * y.conj()),
             (1, 4, -t * plus * x.conj()), (1, 5, -t * plus * x),
             (2, 3, -p['tn1'] * plus), (2, 4, ag), (2, 5, ae),
             (3, 4, ae.conj()), (3, 5, ag.conj()), (4, 5, -p['tn1'] * minus))
    orbital = np.zeros(kx.shape + (6, 6), complex)
    diag = np.arange(6)
    orbital[..., diag, diag] = np.array([p['eps_m']] * 2 + [p['eps_nm']] * 4) - p['mu']
    for i, j, value in bonds:
        orbital[..., i, j] = value
        orbital[..., j, i] = np.conjugate(value)
    h = np.kron(orbital, np.eye(2)) + p['j'] * np.kron(np.diag([-1, 1, 0, 0, 0, 0]), spin)
    if periodic:
        u = cell_gauge(k)
        h = u[..., :, None] * h * u[..., None, :].conj()
    return h
