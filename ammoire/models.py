"""Bloch tight-binding models with explicit provenance.

Momenta are Cartesian, in inverse units of the real-space lattice. Matrix
indices run over site/orbital first and spin (up, down) second. Hamiltonians
use the reciprocal-periodic cell gauge. The built-in Pauli spin operator has
eigenvalues ±1; multiply by ℏ/2 for angular momentum.
"""

from dataclasses import dataclass, replace
from types import MappingProxyType
from typing import Callable, Mapping

import numpy as np


IDENTITY = np.eye(2, dtype=complex)
PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]])
PAULI.setflags(write=False)
TAU_Z = np.kron(PAULI[2], IDENTITY)
SQRT3 = np.sqrt(3.0)
HEX_NEIGHBORS = np.array([[1, 0], [-0.5, SQRT3 / 2], [-0.5, -SQRT3 / 2]])

EZAWA = "https://doi.org/10.1103/PhysRevB.111.125420"
ANTONENKO = "https://doi.org/10.1103/PhysRevLett.134.096703"
BREKKE = "https://doi.org/10.1103/PhysRevB.108.224421"
ROIG = "https://doi.org/10.1103/PhysRevB.110.144412"
GONDOLF = "https://doi.org/10.1103/PhysRevB.111.174436"


@dataclass(frozen=True)
class ModelSpec:
    dimension: int
    geometry: str
    size: int
    form: str
    source: str
    implementation: str
    defaults: Mapping[str, float]
    evaluate: Callable


@dataclass(frozen=True)
class BlochModel:
    name: str
    spec: ModelSpec
    params: Mapping[str, float]
    neel: tuple[float, float, float]

    @property
    def dimension(self):
        return self.spec.dimension

    @property
    def geometry(self):
        return self.spec.geometry

    @property
    def size(self):
        return self.spec.size

    @property
    def source(self):
        return self.spec.source

    @property
    def _lengths(self):
        lengths = [self.params['a'], self.params.get('b', self.params['a'])]
        if self.dimension == 3:
            lengths.append(self.params.get('c', self.params['a']))
        return np.asarray(lengths)

    @property
    def lattice(self):
        """Primitive direct vectors as rows, with a1 along x."""
        basis = np.eye(self.dimension)
        if self.geometry == 'hexagonal':
            basis[1, :2] = [0.5, SQRT3 / 2]
        return basis * self._lengths

    @property
    def reciprocal(self):
        """Primitive reciprocal vectors as rows, including 2π."""
        return 2 * np.pi * np.linalg.inv(self.lattice).T

    @property
    def spin_operators(self):
        """Three physical Pauli-spin matrices in the model basis."""
        orbital_eye = np.eye(self.size // 2)
        return np.array([np.kron(orbital_eye, s) for s in PAULI])

    def with_n(self, neel):
        """Return a copy with a normalized Néel direction."""
        return replace(self, neel=_unit_vector(neel))

    def hamiltonian(self, k):
        """Return H(k) with shape (..., size, size) from k of shape (..., d)."""
        momentum = _momentum(k, self.dimension)
        spin = np.einsum('i,ijk->jk', self.neel, PAULI)
        return self.spec.evaluate(momentum * self._lengths, self.params, spin)

    def slice_kz(self, kz):
        """Hold Cartesian kz fixed in a 3D model to obtain a 2D section."""
        kz = float(kz)
        if self.dimension != 3 or not np.isfinite(kz):
            raise ValueError('a finite kz and a 3D model are required')
        return MomentumSlice(self, kz)


@dataclass(frozen=True)
class MomentumSlice:
    """Fixed-kz bulk section, with the parent orbitals retained."""

    parent: BlochModel
    kz: float

    dimension = 2

    @property
    def name(self):
        return f'{self.parent.name}@kz={self.kz:g}'

    @property
    def geometry(self):
        return {'cubic': 'square', 'tetragonal': 'square', 'orthorhombic': 'rectangular'}.get(
            self.parent.geometry, self.parent.geometry
        )

    @property
    def size(self):
        return self.parent.size

    @property
    def neel(self):
        return self.parent.neel

    @property
    def source(self):
        return self.parent.source

    @property
    def lattice(self):
        return self.parent.lattice[:2, :2]

    @property
    def reciprocal(self):
        return self.parent.reciprocal[:2, :2]

    @property
    def spin_operators(self):
        return self.parent.spin_operators

    def with_n(self, neel):
        return replace(self, parent=self.parent.with_n(neel))

    def hamiltonian(self, k):
        momentum = _momentum(k, self.dimension)
        kz = np.full(momentum.shape[:-1] + (1,), self.kz)
        return self.parent.hamiltonian(np.concatenate([momentum, kz], axis=-1))


def _unit_vector(direction):
    n = np.asarray(direction)
    if n.shape != (3,) or np.iscomplexobj(n):
        raise ValueError('neel must be a finite nonzero 3-vector')
    n = n.astype(float)
    if not np.all(np.isfinite(n)) or not np.any(n):
        raise ValueError('neel must be a finite nonzero 3-vector')
    n = n / np.max(np.abs(n))
    return tuple(n / np.linalg.norm(n))


def _momentum(k, dimension):
    momentum = np.asarray(k)
    if np.iscomplexobj(momentum) or momentum.ndim == 0 or momentum.shape[-1] != dimension:
        raise ValueError(f'expected real momentum with shape (..., {dimension})')
    momentum = np.asarray(momentum, dtype=float)
    if not np.all(np.isfinite(momentum)):
        raise ValueError('momentum must be finite')
    return momentum


def _square_effective(harmonic):
    def evaluate(q, p, spin):
        x, y = np.moveaxis(q, -1, 0)
        energy = p['t'] * (2 - np.cos(x) - np.cos(y)) - p['mu']
        splitting = p['j'] * harmonic(x, y)
        return energy[..., None, None] * IDENTITY + splitting[..., None, None] * spin

    return evaluate


def _hex_i(q, p, spin):
    # Two sine products on triangular-lattice bond stars. Their product is
    # inversion-even and has a sixth-order Γ expansion with six nodal lines.
    nearest = q @ HEX_NEIGHBORS.T
    perpendicular = SQRT3 * q @ HEX_NEIGHBORS[:, ::-1].T
    energy = (2 / 3) * p['t'] * (3 - np.cos(nearest).sum(axis=-1)) - p['mu']
    # This prefactor matches the continuum coefficient j of Eq. (111).
    splitting = (-16 / (3 * SQRT3)) * p['j'] * np.sin(nearest).prod(axis=-1)
    splitting *= np.sin(perpendicular).prod(axis=-1)
    return energy[..., None, None] * IDENTITY + splitting[..., None, None] * spin


def _hex_g(q, p, spin):
    nearest = q[..., :2] @ HEX_NEIGHBORS.T
    z = q[..., 2]
    energy = (2 / 3) * p['t'] * (3 - np.cos(nearest).sum(axis=-1))
    energy += p['tz'] * (1 - np.cos(z)) - p['mu']
    splitting = -4 * p['j'] * np.sin(z) * np.sin(nearest).prod(axis=-1)
    return energy[..., None, None] * IDENTITY + splitting[..., None, None] * spin


def _cubic_i(q, p, spin):
    cosine = np.cos(q)
    cx, cy, cz = np.moveaxis(cosine, -1, 0)
    energy = p['t'] * (3 - cosine.sum(axis=-1)) - p['mu']
    splitting = 8 * p['j'] * (cy - cx) * (cz - cy) * (cx - cz)
    return energy[..., None, None] * IDENTITY + splitting[..., None, None] * spin


def _two_site(energy, tx, tz, phase, exchange, spin, soc=None):
    """Two sites with opposite exchange, in the cell gauge.

    For embedded-basis tx τx + τy λ·σ, the A–B block in the cell
    gauge is exp(-ik·(rB-rA)) (tx I - i λ·σ).
    """
    result = energy[..., None, None] * np.eye(4) + tz[..., None, None] * TAU_Z
    result = result + exchange * np.kron(PAULI[2], spin)
    intersite = tx[..., None, None] * IDENTITY
    if soc is not None:
        intersite -= 1j * np.einsum('...i,ijk->...jk', soc, PAULI)
    intersite *= np.exp(-1j * phase)[..., None, None]
    result[..., :2, 2:] = intersite
    result[..., 2:, :2] = np.swapaxes(intersite.conj(), -1, -2)
    return result


def _lieb(q, p, spin):
    x, y = np.moveaxis(q, -1, 0)
    energy = -2 * p['t2'] * (np.cos(x) + np.cos(y)) - p['mu']
    tx = -4 * p['t1'] * np.cos(x / 2) * np.cos(y / 2)
    tz = -2 * p['td'] * (np.cos(x) - np.cos(y))
    soc = np.zeros(q.shape[:-1] + (3,))
    soc[..., 2] = p['soc'] * np.sin(x / 2) * np.sin(y / 2)
    return _two_site(energy, tx, tz, (x + y) / 2, p['j'], spin, soc)


def _brekke_lieb(q, p, spin):
    # Orbital order (A, nonmagnetic, B), at (1/2,0), (0,0), (0,1/2).
    # q measures the primitive cell, twice the nearest-neighbor distance.
    x, y = np.moveaxis(q / 2, -1, 0)
    cx, cy = np.cos(x), np.cos(y)
    h = np.zeros(q.shape[:-1] + (6, 6), dtype=complex)
    h[..., :2, :2] = -p['mu'] * IDENTITY - p['j'] * spin
    h[..., 2:4, 2:4] = (p['eps_nm'] - p['mu']) * IDENTITY
    h[..., 4:, 4:] = -p['mu'] * IDENTITY + p['j'] * spin
    h[..., :2, 2:4] = (2 * p['t'] * cx * np.exp(1j * x))[..., None, None] * IDENTITY
    h[..., 2:4, 4:] = (2 * p['t'] * cy * np.exp(-1j * y))[..., None, None] * IDENTITY
    h[..., :2, 4:] = (4 * p['t2'] * cx * cy * np.exp(1j * (x - y)))[..., None, None] * IDENTITY
    h[..., 2:4, :2] = h[..., :2, 2:4].conj()
    h[..., 4:, 2:4] = h[..., 2:4, 4:].conj()
    h[..., 4:, :2] = h[..., :2, 4:].conj()
    return h


def _gondolf_dxy(q, p, spin):
    x, y = np.moveaxis(q, -1, 0)
    cx, cy = np.cos(x), np.cos(y)
    energy = -2 * p['t2'] * (cx + cy) - 4 * p['t3'] * cx * cy - p['mu']
    tx = -4 * p['t1'] * np.cos(x / 2) * np.cos(y / 2)
    tz = 4 * p['t4'] * np.sin(x) * np.sin(y)
    return _two_site(energy, tx, tz, (x + y) / 2, p['j'], spin)


def _roig_223(q, p, spin):
    cx, cy, cz = np.moveaxis(np.cos(q), -1, 0)
    hx, hy, hz = np.moveaxis(np.cos(q / 2), -1, 0)
    sx, sy, sz = np.moveaxis(np.sin(q / 2), -1, 0)
    energy = p['t1'] * (cx + cy + cz) - p['mu']
    tx = p['t2'] * hx * hy * hz
    tz = p['t3'] * (cx - cy) * (cx - cz) * (cy - cz)
    soc = p['soc'] * np.stack([hx * sy * sz, sx * hy * sz, sx * sy * hz], axis=-1)
    return _two_site(energy, tx, tz, q.sum(axis=-1) / 2, p['j'], spin, soc)


def _roig_194(q, p, spin):
    x, y, z = np.moveaxis(q, -1, 0)
    half_y = SQRT3 * y / 2
    fx = np.sin(x) + np.sin(x / 2) * np.cos(half_y)
    fy = SQRT3 * np.cos(x / 2) * np.sin(half_y)
    energy = p['t1'] * (np.cos(x) + 2 * np.cos(x / 2) * np.cos(half_y))
    energy += p['t2'] * np.cos(z) - p['mu']
    tx = p['t3'] * np.cos(z / 2)
    tz = p['t4'] * np.sin(z) * fy * (fy**2 - 3 * fx**2)
    in_plane_soc = p['soc'] * np.cos(z / 2)
    soc = np.stack([
        in_plane_soc * (fx**2 - fy**2),
        -2 * in_plane_soc * fx * fy,
        p['soc_z'] * np.sin(z / 2) * fx * (fx**2 - 3 * fy**2),
    ], axis=-1)
    return _two_site(energy, tx, tz, z / 2, p['j'], spin, soc)


def _roig_rutile(q, p, spin):
    cx, cy, cz = np.moveaxis(np.cos(q), -1, 0)
    hx, hy, hz = np.moveaxis(np.cos(q / 2), -1, 0)
    sx, sy, sz = np.moveaxis(np.sin(q / 2), -1, 0)
    energy = p['t1'] * (cx + cy) + p['t2'] * cz + p['t3'] * cx * cy - p['mu']
    energy += (p['t4'] * (cx + cy) + p['t5'] * cx * cy) * cz
    tx = p['t8'] * hx * hy * hz
    tz = (p['t6'] + p['t7'] * cz) * 4 * sx * hx * sy * hy
    soc = np.stack([
        p['soc'] * sz * sx * hy,
        -p['soc'] * sz * sy * hx,
        p['soc_z'] * hx * hy * hz * (cx - cy),
    ], axis=-1)
    return _two_site(energy, tx, tz, q.sum(axis=-1) / 2, p['j'], spin, soc)


def _roig_fesb2(q, p, spin):
    cx, cy, cz = np.moveaxis(np.cos(q), -1, 0)
    hx, hy, hz = np.moveaxis(np.cos(q / 2), -1, 0)
    sx, sy, sz = np.moveaxis(np.sin(q / 2), -1, 0)
    energy = p['t1x'] * cx + p['t1y'] * cy + p['t2'] * cz + p['t3'] * cx * cy - p['mu']
    energy += (p['t4x'] * cx + p['t4y'] * cy + p['t5'] * cx * cy) * cz
    tx = p['t8'] * hx * hy * hz
    tz = (p['t6'] + p['t7'] * cz) * 4 * sx * hx * sy * hy
    soc = np.stack([
        p['soc_x'] * sz * sx * hy,
        p['soc_y'] * sz * sy * hx,
        p['soc_z'] * hx * hy * hz,
    ], axis=-1)
    return _two_site(energy, tx, tz, q.sum(axis=-1) / 2, p['j'], spin, soc)


def _roig_kappa_cl(q, p, spin):
    x, y = np.moveaxis(q, -1, 0)
    cx, cy = np.cos(x), np.cos(y)
    energy = p['t1x'] * cx + p['t1y'] * cy + p['t3'] * cx * cy - p['mu']
    tx = p['t8'] * np.cos(x / 2) * np.cos(y / 2)
    tz = p['t6'] * np.sin(x) * np.sin(y)
    return _two_site(energy, tx, tz, (x + y) / 2, p['j'], spin)


def _layer_group(offset, form):
    """Representative spin-independent dispersion with Table II harmonics."""
    def evaluate(q, p, spin):
        x, y = np.moveaxis(q, -1, 0)
        cx, cy = np.cos(x), np.cos(y)
        phase = q @ np.asarray(offset)
        energy = -2 * p['t'] * (cx + cy) - p['mu']
        tx = p['tx'] * np.cos(x * offset[0]) * np.cos(y * offset[1])
        dx2y2, dxy = cx - cy, np.sin(x) * np.sin(y)
        harmonics = {'dxy': dxy, 'dx2-y2': dx2y2, 'g': dxy * dx2y2, 'd-mixed': dx2y2}
        tz = p['td'] * harmonics[form]
        if form == 'd-mixed':
            tz += p['td_xy'] * dxy
        return _two_site(energy, tx, tz, phase, p['j'], spin)

    return evaluate


def _spec(dimension, geometry, size, form, source, implementation, evaluate, **defaults):
    parameters = MappingProxyType({'a': 1.0, 'mu': 0.0, 'j': 0.2, **defaults})
    return ModelSpec(dimension, geometry, size, form, source, implementation, parameters, evaluate)


_MODELS = {
    'square-dxy': _spec(2, 'square', 2, 'dxy', EZAWA, 'Eq. (68), configurable scales',
                       _square_effective(lambda x, y: np.sin(x) * np.sin(y)), t=1.0),
    'square-dx2-y2': _spec(2, 'square', 2, 'dx2-y2', ANTONENKO,
                          'Representative single-band harmonic',
                          _square_effective(lambda x, y: np.cos(x) - np.cos(y)), t=1.0),
    'square-g': _spec(2, 'square', 2, 'g', EZAWA, 'Eq. (101), configurable scales',
                     _square_effective(lambda x, y: np.sin(x) * np.sin(y) * (np.cos(x) - np.cos(y))),
                     t=1.0),
    'hex-i': _spec(2, 'hexagonal', 2, 'i', EZAWA,
                   'Planar adaptation of Eqs. (116)–(117), normalized to Eq. (111) near Γ',
                   _hex_i, t=1.0),
    'hex-g-3d': _spec(3, 'hexagonal', 2, 'g', EZAWA, 'Eqs. (108)–(109), configurable scales',
                      _hex_g, t=1.0, tz=1.0, c=1.0),
    'cubic-i': _spec(3, 'cubic', 2, 'i', EZAWA, 'Eq. (124), configurable scales',
                     _cubic_i, t=1.0),
    'lieb-dx2-y2': _spec(2, 'square', 4, 'dx2-y2', ANTONENKO,
                        'Eqs. (S8), (S38)–(S41) with J′ and net moment zero',
                        _lieb, t1=0.5, t2=0.25, td=1.0, j=2.0, soc=0.0),
    'brekke-lieb': _spec(2, 'square', 6, 'dx2-y2', BREKKE,
                        'Eq. (21) in periodic cell gauge, Fig. 2 defaults with a=2a_paper',
                        _brekke_lieb, a=2.0, t=1.0, t2=0.0, j=0.4, eps_nm=0.0),
    'gondolf-dxy': _spec(2, 'square', 4, 'dxy', GONDOLF, 'Appendix A and Table I defaults in eV',
                         _gondolf_dxy, t1=0.425, t2=0.05, t3=-0.025, t4=-0.075, mu=0.3),
    'roig-223': _spec(3, 'cubic', 4, 'i', ROIG, 'Eqs. (35)–(36), representative scales',
                      _roig_223, t1=-0.3, t2=1.0, t3=0.5, soc=0.0),
    'roig-194': _spec(3, 'hexagonal', 4, 'g', ROIG,
                      'Eqs. (32)–(34), CrSb/MnTe-inspired with representative scales',
                      _roig_194, t1=-0.3, t2=-0.2, t3=1.0, t4=0.15, soc=0.0, soc_z=0.0, c=1.0),
    'roig-ruo2': _spec(3, 'tetragonal', 4, 'dxy', ROIG,
                       'One-orbital Eqs. (3), (5)–(6), Fig. 2 and Table 5 energy defaults in eV',
                       _roig_rutile, t1=-0.05, t2=0.7, t3=0.5, t4=-0.15, t5=-0.4,
                       t6=-0.6, t7=0.3, t8=1.7, mu=0.25, soc=0.0, soc_z=0.0, c=1.0),
    'roig-mnf2': _spec(3, 'tetragonal', 4, 'dxy', ROIG,
                       'One-orbital Eqs. (3), (5)–(6), Fig. S4 and Table 5 energy defaults in eV',
                       _roig_rutile, t1=0.0, t2=0.13, t3=0.0, t4=-0.02, t5=0.015,
                       t6=0.0, t7=0.03, t8=0.33, mu=-0.01, j=0.05, soc=0.0, soc_z=0.0, c=1.0),
    'roig-fesb2': _spec(3, 'orthorhombic', 4, 'dxy', ROIG,
                        'Eqs. (6), (30)–(31), Fig. 7 and Table 5 energy defaults in eV',
                        _roig_fesb2, t1x=-0.1, t1y=-0.05, t2=-0.05, t3=0.06, t4x=0.1,
                        t4y=0.05, t5=-0.05, t6=0.05, t7=-0.1, t8=0.15, mu=-0.12,
                        j=0.1, soc_x=0.0, soc_y=0.0, soc_z=0.0, b=1.0, c=1.0),
    'roig-kappa-cl': _spec(2, 'rectangular', 4, 'dxy', ROIG,
                           'Eqs. (37)–(38), Fig. S7 and Table 5 energy defaults in eV',
                           _roig_kappa_cl, t1x=0.08, t1y=-0.01, t3=-0.03, t6=0.05,
                           t8=0.3, mu=-0.1, j=0.05, b=1.0),
}

for _label, _offset, _form, _geometry in (
    ('l15-2ab', (0.5, 0), 'dxy', 'rectangular'),
    ('l16-2ab', (0, 0.5), 'dxy', 'rectangular'),
    ('l17-2ab', (0.5, 0.5), 'dxy', 'rectangular'),
    ('l40-2ab', (0.5, 0), 'dxy', 'rectangular'),
    ('l44-2ab', (0.5, 0.5), 'dxy', 'rectangular'),
    ('l51-2c', (0.5, 0.5), 'd-mixed', 'square'),
    ('l61-2c', (0.5, 0.5), 'dx2-y2', 'square'),
    ('l63-2a', (0.5, 0.5), 'g', 'square'),
    ('l63-2b', (0.5, 0.5), 'dxy', 'square'),
):
    _extra = {'b': 1.0} if _geometry == 'rectangular' else {}
    if _form == 'd-mixed':
        _extra['td_xy'] = 0.15
    _MODELS[f'gondolf-{_label}'] = _spec(
        2, _geometry, 4, _form, GONDOLF,
        'Appendix B Table II intersite and anisotropic harmonics, representative other terms',
        _layer_group(_offset, _form), t=1.0, tx=0.8, td=0.3, **_extra,
    )

MODEL_REGISTRY = MappingProxyType(_MODELS)


def make_model(name, *, neel=(0.0, 0.0, 1.0), **parameters):
    """Instantiate a registry model; `j` is the signed exchange energy."""
    if name not in MODEL_REGISTRY:
        raise ValueError(f'unknown model {name!r}; choose from {", ".join(MODEL_REGISTRY)}')
    spec = MODEL_REGISTRY[name]
    unknown = parameters.keys() - spec.defaults.keys()
    if unknown:
        raise ValueError(f'unknown parameter(s) for {name}: {", ".join(sorted(unknown))}')
    params = {key: float(value) for key, value in (spec.defaults | parameters).items()}
    if not all(np.isfinite(value) for value in params.values()):
        raise ValueError('model parameters must be finite')
    if any(params.get(axis, 1.0) <= 0 for axis in ('a', 'b', 'c')):
        raise ValueError('lattice lengths must be positive')
    return BlochModel(name, spec, MappingProxyType(params), _unit_vector(neel))
