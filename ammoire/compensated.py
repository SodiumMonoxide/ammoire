"""Commensurate TB control for a compensated buckled triangular parent.

Each monolayer contains two triangular sublayers with opposite exchange
fields. The model uses representative parameters from the MnSe-inspired
control in the local twist project. The material prediction that motivated
the control is https://doi.org/10.1103/PhysRevMaterials.8.L051401.

There is one orbital per site. The basis lists four sublayers from bottom
to top, with physical up/down spin interleaved at each site. The row lattice
is in units of the monolayer lattice constant a, and Cartesian k is in 1/a.
Energies and hopping amplitudes are in eV. The model uses fixed collinear
exchange and ideal positions with no SOC or structural relaxation.
"""

from itertools import product

import numpy as np

from .commensurate import COMMENSURATE_PAIRS, INTERLAYER_CUTOFF, _cell_sites, _geometry, _momentum


def _shifted_sublayer(cell, primitive, unshifted, sign):
    shift = sign * primitive.sum(axis=0) / 3
    fractional = (unshifted + shift) @ np.linalg.inv(cell)
    return (fractional - np.floor(fractional + 1e-12)) @ cell


class CompensatedTriangularBilayer:
    """Buckled antiferromagnetic parents with independently tunable interface terms.

    ``interlayer`` scales the Gaussian hopping between the facing
    sublayers. ``facing_potential`` scales the onsite shift from summing
    the same Gaussian kernel. Its mean within each facing sublayer and
    its zero-mean registry part have separate ``uniform_scale`` and
    ``registry_scale`` multipliers. Set ``interlayer=facing_potential=0``
    for the decoupled-parent control, or ``exchange=0`` for the spin
    degeneracy check with both interface terms present.

    The isolated monolayers are spectrally spin degenerate. An interface
    term can lift this degeneracy after twisting. The effect sizes depend
    on the representative parameters, so this model provides an ablation
    of the two mechanisms rather than a material fit.
    """

    dimension = 2
    geometry = "compensated-triangular-supercell"
    energy_unit = "eV"
    length_unit = "a"

    def __init__(self, *, m=2, n=1, planar_hopping=-0.16, buckled_hopping=-0.13,
                 interlayer=-0.065, facing_potential=-0.04, exchange=0.62,
                 uniform_scale=1.0, registry_scale=1.0, width=0.42):
        if (not isinstance(m, (int, np.integer)) or isinstance(m, (bool, np.bool_))
                or not isinstance(n, (int, np.integer)) or isinstance(n, (bool, np.bool_))
                or (m, n) not in COMMENSURATE_PAIRS):
            raise ValueError(f"(m,n) must be one of the validated pairs: {COMMENSURATE_PAIRS}")
        values = np.asarray([planar_hopping, buckled_hopping, interlayer, facing_potential,
                             exchange, uniform_scale, registry_scale, width], dtype=float)
        if values.shape != (8,) or not np.all(np.isfinite(values)) or values[-1] <= 0:
            raise ValueError("all model parameters must be finite and width must be positive")
        self.m, self.n = int(m), int(n)
        (self.planar_hopping, self.buckled_hopping, self.interlayer, self.facing_potential,
         self.exchange, self.uniform_scale, self.registry_scale, self.width) = map(float, values)
        theta, self.lattice, self.primitive_layers = _geometry(True, self.m, self.n)
        self.twist_deg = float(np.rad2deg(theta))
        self.reciprocal = 2 * np.pi * np.linalg.inv(self.lattice).T
        unshifted = [_cell_sites(self.lattice, primitive) for primitive in self.primitive_layers]
        self.sites_per_sublayer = len(unshifted[0])
        self.sites = np.vstack([
            _shifted_sublayer(self.lattice, primitive, sites, sign)
            for primitive, sites in zip(self.primitive_layers, unshifted) for sign in (-1, 1)
        ])
        self.nsites = len(self.sites)
        self.size = 2 * self.nsites
        self.sublayer = np.repeat(np.arange(4), self.sites_per_sublayer)
        self._exchange_matrix = np.diag(self.exchange * np.where(self.sublayer % 2 == 0, -1, 1))

        translations = np.array(list(product(range(-2, 3), repeat=2))) @ self.lattice
        relative = self.sites[None, :, :] - self.sites[:, None, :]
        distances = np.linalg.norm(relative + translations[:, None, None, :], axis=-1)
        same_plane = self.sublayer[:, None] == self.sublayer[None, :]
        monolayer = self.sublayer // 2
        same_mono = monolayer[:, None] == monolayer[None, :]
        facing = (np.abs(self.sublayer[:, None] - self.sublayer[None, :]) == 1) & ~same_mono
        planar = np.isclose(distances, 1.0, atol=1e-12, rtol=0) & same_plane
        buckled = np.isclose(distances, 1 / np.sqrt(3), atol=1e-12, rtol=0) & same_mono & ~same_plane
        if np.any(planar.sum(axis=(0, 2)) != 6) or np.any(buckled.sum(axis=(0, 2)) != 3):
            raise RuntimeError("each site must have six planar and three buckled nearest neighbors")
        interface = facing & (distances < INTERLAYER_CUTOFF)
        kernel = np.exp(-distances**2 / (2 * self.width**2)) * interface
        active = np.any(planar | buckled | interface, axis=(1, 2))
        self._translations = translations[active]
        hopping = (self.planar_hopping * planar[active] + self.buckled_hopping * buckled[active]
                   + self.interlayer * kernel[active])
        self._hoppings = hopping.reshape(-1, self.nsites**2)
        potential = self.facing_potential * kernel.sum(axis=(0, 2)).reshape(4, self.sites_per_sublayer)
        mean = potential.mean(axis=1, keepdims=True)
        potential = self.uniform_scale * mean + self.registry_scale * (potential - mean)
        self._onsite = np.diag(potential.ravel())
        for array in (self.lattice, self.reciprocal, self.sites, self.sublayer, *self.primitive_layers):
            array.setflags(write=False)

    def _spatial(self, k):
        points = _momentum(k)
        phase = np.exp(1j * (points @ self._translations.T))
        hopping = (phase @ self._hoppings).reshape(*points.shape[:-1], self.nsites, self.nsites)
        return hopping + self._onsite

    def hamiltonian(self, k):
        """Return periodic-gauge H(k) in eV with shape (..., size, size)."""
        spatial = self._spatial(k)
        result = np.zeros(spatial.shape[:-2] + (self.size, self.size), dtype=complex)
        result[..., ::2, ::2] = spatial + self._exchange_matrix
        result[..., 1::2, 1::2] = spatial - self._exchange_matrix
        return result

    def spin_bands(self, k):
        """Ascending within-spin energies, with up and down on the penultimate axis."""
        spatial = self._spatial(k)
        return np.stack((np.linalg.eigvalsh(spatial + self._exchange_matrix),
                         np.linalg.eigvalsh(spatial - self._exchange_matrix)), axis=-2)

    def splitting(self, k):
        """Return E_up minus E_down in eV at equal within-spin energy rank."""
        energies = self.spin_bands(k)
        return energies[..., 0, :] - energies[..., 1, :]
