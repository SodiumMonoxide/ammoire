"""Periodic real-space TB for eight validated commensurate triangular cells.

This is the geometry demonstration used in the 2026-11-06 calculations,
ported to a standalone NumPy module. It is motivated by twisted A-type
antiferromagnets (https://arxiv.org/abs/2404.17146 and
https://arxiv.org/abs/2407.15097), rather than a fitted NiCl2 Hamiltonian or
an equation-by-equation reproduction of either paper.

The default (m,n)=(2,1) cell contains seven sites per layer, with one orbital
and physical up/down spin per site: 28 states in total. Other validated
cells contain 7, 13, or 19 sites per layer (28, 52, or 76 states).
Basis order is lower
layer sites first, upper layer sites second, with up/down interleaved at
each site. Energies are eV; real vectors are in units of the monolayer
lattice constant a and Cartesian momenta in inverse a. Real and reciprocal
vectors are rows. The cell Bloch gauge is exactly reciprocal-periodic.

The ``twisted=False`` zero-twist control retains the default seven-site
supercell and lower layer's orientation. The angle-scan endpoints (2,2)
and (1,0) are also folded into seven-site cells, as in the subsequent
2026-11-06 angle study; intermediate cells use the factor-three reduction
when m-n is divisible by three. Only those eight previously checked
integer pairs are supported, rather than arbitrary incommensurate angles.
Magnetism is a fixed collinear opposite exchange field on the two layers;
no magnetic-order optimization, SOC, or atomic relaxation is performed.
"""

from itertools import product

import numpy as np


COMMENSURATE_ANGLE_DEG = 21.786789298261812
COMMENSURATE_PAIRS = ((2, 2), (3, 2), (2, 1), (5, 2), (3, 1), (4, 1), (7, 1), (1, 0))
INTERLAYER_CUTOFF = 1.75
_SZ = np.diag([1.0, -1.0])


def _rotation(angle):
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s], [s, c]])


def _geometry(twisted, m, n):
    unit = np.array([[1.0, 0.0], [0.5, np.sqrt(3) / 2]])
    if not twisted:
        lower = unit @ _rotation(-np.deg2rad(COMMENSURATE_ANGLE_DEG / 2)).T
        return 0.0, np.array([[1, 2], [-2, 3]]) @ lower, (lower, lower.copy())
    if m == n:
        return 0.0, np.array([[1, 2], [-2, 3]]) @ unit, (unit.copy(), unit.copy())
    area = m * m + m * n + n * n
    theta = np.arccos((m * m + 4 * m * n + n * n) / (2 * area))
    layers = (unit @ _rotation(-theta / 2).T, unit @ _rotation(theta / 2).T)
    if n == 0:
        return theta, np.array([[1, 2], [-2, 3]]) @ layers[0], layers
    lower_map = np.array([[n, m], [-m, m + n]])
    upper_map = np.array([[m, n], [-n, m + n]])
    if (m - n) % 3 == 0:
        reduction = np.array([[1, 1], [-1, 2]])
        if np.any((reduction @ lower_map) % 3) or np.any((reduction @ upper_map) % 3):
            raise RuntimeError("Invalid integer commensurate-cell reduction")
        lower_map, upper_map = reduction @ lower_map // 3, reduction @ upper_map // 3
    cell = lower_map @ layers[0]
    if not np.allclose(cell, upper_map @ layers[1], atol=2e-13, rtol=0):
        raise RuntimeError("The two layer lattices must share a commensurate cell")
    return theta, cell, layers


def _cell_sites(cell, lattice):
    corners = np.array([[0, 0], [1, 0], [0, 1], [1, 1]]) @ cell @ np.linalg.inv(lattice)
    lo = np.floor(corners.min(axis=0)).astype(int)
    hi = np.ceil(corners.max(axis=0)).astype(int)
    indices = np.array(list(product(*(range(a, b + 1) for a, b in zip(lo, hi)))))
    points = indices @ lattice
    fractional = points @ np.linalg.inv(cell)
    fractional[np.abs(fractional) < 1e-12] = 0
    inside = ((fractional >= 0) & (fractional < 1 - 1e-12)).all(axis=1)
    sites = points[inside]
    if len(sites) != round(np.linalg.det(cell) / np.linalg.det(lattice)):
        raise RuntimeError("Commensurate site count must equal the cell area")
    return sites


def _momentum(k):
    points = np.asarray(k, dtype=float)
    if points.ndim == 0 or points.shape[-1] != 2 or not np.all(np.isfinite(points)):
        raise ValueError("k must be finite Cartesian momenta with shape (..., 2)")
    return points


class CommensurateTriangularBilayer:
    """Validated commensurate triangular bilayer or the original AA control.

    Parameters are illustrative: nearest-neighbor intralayer hopping
    -0.18 eV, interlayer hopping -0.07 exp[-d_parallel^2/(2 width^2)] eV
    truncated at d_parallel < 1.75 a, and exchange J*layer*sigma_z with
    J=0.60 eV. ``width`` is in a. The integer pair (m,n) selects one of
    ``COMMENSURATE_PAIRS``: approximately 0, 13.17, 21.79, 27.80, 32.20,
    38.21, 46.83, or 60 degrees. ``twisted=False`` preserves the original
    rotated seven-site AA control and requires default m=2, n=1.
    """

    dimension = 2
    geometry = "triangular-supercell"
    energy_unit = "eV"
    length_unit = "a"

    def __init__(self, twisted=True, *, m=2, n=1, hopping=-0.18, interlayer=-0.07,
                 exchange=0.60, width=0.42):
        if not isinstance(twisted, (bool, np.bool_)):
            raise ValueError("twisted must be a boolean selecting a commensurate cell or the AA control")
        if (not isinstance(m, (int, np.integer)) or isinstance(m, (bool, np.bool_))
                or not isinstance(n, (int, np.integer)) or isinstance(n, (bool, np.bool_))
                or (m, n) not in COMMENSURATE_PAIRS):
            raise ValueError(f"(m,n) must be one of the validated pairs: {COMMENSURATE_PAIRS}")
        if not twisted and (m, n) != (2, 1):
            raise ValueError("twisted=False is the original (2,1)-area AA control; use (2,2) for the scan endpoint")
        parameters = np.asarray([hopping, interlayer, exchange, width], dtype=float)
        if parameters.shape != (4,) or not np.all(np.isfinite(parameters)) or width <= 0:
            raise ValueError("hopping, interlayer, exchange must be finite; width must be positive")
        self.twisted = bool(twisted)
        self.m, self.n = int(m), int(n)
        self.hopping, self.interlayer, self.exchange, self.width = map(float, parameters)
        theta, self.lattice, self.primitive_layers = _geometry(self.twisted, self.m, self.n)
        self.twist_deg = float(np.rad2deg(theta))
        self.reciprocal = 2 * np.pi * np.linalg.inv(self.lattice).T
        self.sites = np.vstack([_cell_sites(self.lattice, layer)
                                for layer in self.primitive_layers])
        self.nsites = len(self.sites)
        self.sites_per_layer = self.nsites // 2
        self.size = 2 * self.nsites
        self.layer = np.repeat([-1.0, 1.0], self.sites_per_layer)
        self.exchange_operator = np.kron(np.diag(self.layer), _SZ)
        translations = np.array(list(product(range(-3, 4), repeat=2))) @ self.lattice
        displacements = self.sites[None, :, :] - self.sites[:, None, :]
        distances = np.linalg.norm(displacements + translations[:, None, None, :], axis=-1)
        same_layer = self.layer[:, None] == self.layer[None, :]
        intra = np.isclose(distances, 1, atol=1e-12, rtol=0) & same_layer
        if not np.array_equal(intra.sum(axis=(0, 2)), np.full(self.nsites, 6)):
            raise RuntimeError("Each triangular site must have six nearest neighbors")
        inter = (distances < INTERLAYER_CUTOFF) & ~same_layer
        hopping_matrices = (self.hopping * intra + self.interlayer
                            * np.exp(-distances**2 / (2 * self.width**2)) * inter)
        active = np.any(hopping_matrices, axis=(1, 2))
        self.translations = translations[active]
        self.hoppings = hopping_matrices[active].reshape(-1, self.nsites**2)
        # Geometry and stored matrix coefficients describe one immutable basis.
        for array in (self.lattice, self.reciprocal, self.sites, self.layer,
                      self.exchange_operator, self.translations, self.hoppings,
                      *self.primitive_layers):
            array.setflags(write=False)

    def hamiltonian(self, k):
        """Return H(k) in eV, shape (..., size, size); k is Cartesian in 1/a."""
        momentum = _momentum(k)
        phases = np.exp(1j * (momentum @ self.translations.T))
        spatial = (phases @ self.hoppings).reshape(*momentum.shape[:-1], self.nsites, self.nsites)
        return np.kron(spatial, np.eye(2)) + self.exchange * self.exchange_operator

    def spin_bands(self, k):
        """Ascending eigenvalues per physical spin, shape (..., 2, nsites).

        Channel 0 is up, channel 1 down. The final index is the local
        energy rank within a spin channel, not a tracked orbital label.
        """
        h = self.hamiltonian(k)
        return np.stack([np.linalg.eigvalsh(h[..., ::2, ::2]),
                         np.linalg.eigvalsh(h[..., 1::2, 1::2])], axis=-2)

    def splitting(self, k):
        """E_up - E_down for equal within-spin energy ranks, in eV."""
        bands = self.spin_bands(k)
        return bands[..., 0, :] - bands[..., 1, :]

    def primitive_dispersion(self, k, layer):
        """Unfolded spin-independent nearest-neighbor dispersion, in eV.

        ``layer`` is 0 (lower) or 1 (upper). The difference of these two
        dispersions is a same-primitive-momentum reference; at a general
        supercell k it need not equal the lowest sorted spin-band split.
        They coincide in the decoupled, small-k isolated-edge regime.
        """
        if (not isinstance(layer, (int, np.integer)) or isinstance(layer, (bool, np.bool_))
                or layer not in (0, 1)):
            raise ValueError("layer must be 0 (lower) or 1 (upper)")
        lattice = self.primitive_layers[layer]
        bonds = np.vstack([lattice, lattice[1] - lattice[0]])
        return 2 * self.hopping * np.cos(_momentum(k) @ bonds.T).sum(axis=-1)
