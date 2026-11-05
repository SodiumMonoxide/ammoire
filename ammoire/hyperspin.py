"""Projected HAM block of Ma et al., PRL 137, 116704 (2026), Eqs. (4), (6).

The six periodic coefficients below form an illustrative square-lattice
parameterization. The four states are AFM subsystem/label pairs (A+, A-,
B+, B-). The label Pauli matrices do not represent physical spin.
"""

from dataclasses import dataclass

import numpy as np

I2 = np.eye(2, dtype=complex)
SY = np.array([[0, -1j], [1j, 0]])
SZ = np.diag([1, -1])


def projected_block(h0, h1, h2, h3, h4, h5):
    """Assemble Eq. (6) in its original order, with broadcast coefficients."""
    h0, h1, h2, h3, h4, h5 = np.broadcast_arrays(h0, h1, h2, h3, h4, h5)
    if any(np.any(x.imag) for x in (h0, h1, h4, h5)):
        raise ValueError("h0, h1, h4, h5 must be real")
    out = np.empty(h0.shape + (4, 4), complex)
    out[..., :2, :2] = (h0 + h1)[..., None, None] * I2 + h4[..., None, None] * SY
    out[..., 2:, 2:] = (h0 - h1)[..., None, None] * I2 + h5[..., None, None] * SY
    hopping = h2[..., None, None] * I2 + h3[..., None, None] * SY
    out[..., :2, 2:] = hopping
    out[..., 2:, :2] = hopping.conj().swapaxes(-1, -2)
    return out


@dataclass(frozen=True)
class SquareHAM:
    """Odd-wave demo with H(-k)=H(k)* and opposite Q sectors at ±k."""

    a: float = 1.0
    t: float = -0.35
    imbalance: float = 0.12
    hybrid: float = 0.17
    odd_mix: float = 0.20
    odd_a: float = 0.45
    odd_b: float = 0.30
    offset: float = 0.0

    name = "square-ham"
    dimension = 2
    size = 4
    geometry = "square"

    def __post_init__(self):
        if not np.isfinite(self.a) or self.a <= 0:
            raise ValueError("a must be positive and finite")

    @property
    def lattice(self):
        return self.a * np.eye(2)

    @property
    def reciprocal(self):
        return (2 * np.pi / self.a) * np.eye(2)

    @property
    def charge(self):
        """Dimensionless Q = 2 J_y / ℏ."""
        return np.kron(I2, SY)

    def coefficients(self, k):
        k = np.asarray(k, dtype=float)
        if k.ndim == 0 or k.shape[-1] != 2:
            raise ValueError("k must have shape (..., 2)")
        x, y = np.moveaxis(k * self.a, -1, 0)
        return (
            self.offset + 2 * self.t * (np.cos(x) + np.cos(y)),
            self.imbalance,
            self.hybrid,
            self.odd_mix * np.sin(x),
            self.odd_a * np.sin(x),
            self.odd_b * np.sin(x),
        )

    def hamiltonian(self, k):
        return projected_block(*self.coefficients(k))

    def sector_bands(self, k):
        """Independent analytic E[..., Q index (+, -), branch (low, high)]."""
        h0, h1, h2, h3, h4, h5 = np.broadcast_arrays(*self.coefficients(k))
        q = np.array([1, -1])
        center = h0[..., None] + (h4 + h5)[..., None] * q / 2
        mass = h1[..., None] + (h4 - h5)[..., None] * q / 2
        radius = np.sqrt(mass**2 + np.abs(h2[..., None] + h3[..., None] * q)**2)
        return center[..., None] + radius[..., None] * np.array([-1, 1])


def charge_channels(tunneling, q_top, q_bottom):
    """Split T into same- and opposite-charge channels for Hermitian Q²=I."""
    tunneling = np.asarray(tunneling, complex)
    q_top, q_bottom = (np.asarray(q, complex) for q in (q_top, q_bottom))
    if tunneling.ndim != 2:
        raise ValueError("T must be a matrix")
    for q, size in zip((q_top, q_bottom), tunneling.shape):
        identity = np.eye(size)
        if (q.shape != identity.shape or not np.allclose(q, q.conj().T)
                or not np.allclose(q @ q, identity)):
            raise ValueError("Q must be a Hermitian involution in its layer basis")
    conjugated = q_top @ tunneling @ q_bottom
    return (tunneling + conjugated) / 2, (tunneling - conjugated) / 2
