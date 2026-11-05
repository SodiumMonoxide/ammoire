"""Plane-wave continuum model for a twisted bilayer with a common parent lattice.

Reciprocal vectors are rows in Cartesian coordinates, and angles are in degrees.
The two layers rotate actively by +theta/2 and -theta/2 about z. The basis is
ordered as (layer, reciprocal replica, monolayer basis), with the top layer
first. The Hamiltonian uses cell-periodic monolayer matrix conventions.
"""

import operator
from types import MappingProxyType

import numpy as np


def _rotation(angle):
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s], [s, c]])


def _vector2(value, name):
    vector = np.asarray(value, dtype=float)
    if vector.shape != (2,) or not np.all(np.isfinite(vector)):
        raise ValueError(f"{name} must be a finite Cartesian vector of length 2")
    return vector


def _lattice_star(reciprocal):
    """Return primitive integer first-star shifts and the geometric lattice type."""
    reciprocal = np.asarray(reciprocal, dtype=float)
    if reciprocal.shape != (2, 2) or not np.all(np.isfinite(reciprocal)):
        raise ValueError("reciprocal must contain two finite Cartesian row vectors")
    lengths = np.linalg.norm(reciprocal, axis=1)
    if np.any(lengths == 0):
        raise ValueError("reciprocal vectors must be nonzero")
    directions = reciprocal / lengths[:, None]
    cosine = directions[0] @ directions[1]
    star = [(1, 0), (0, 1)]
    if np.isclose(cosine, 0, atol=1e-10):
        geometry = "square" if np.isclose(*lengths, rtol=1e-10, atol=0) else "rectangular"
    elif np.isclose(abs(cosine), 0.5, atol=1e-10) and np.isclose(
        *lengths, rtol=1e-10, atol=0
    ):
        geometry = "hexagonal"
        star.append((1, -int(np.sign(cosine))))
    else:
        raise ValueError("use square, rectangular, or 60°/120° hexagonal primitive rows")
    return tuple(star + [(-n, -m) for n, m in star]), geometry


def tunneling_harmonics(reciprocal, *, constant=0.0, first_star=None):
    """Build a uniform Fourier ansatz with equal amplitudes on the first star.

    Set either amplitude to None to omit it. Amplitudes can be scalars or
    rectangular orbital matrices. Use an explicit dictionary for distinct
    amplitudes or phases, including three-transfer hexagonal valley models.
    The latter also need valley-specific momentum anchors and orbital gauge.
    For a rectangular lattice the star contains both primitive directions.
    """
    harmonics = {} if constant is None else {(0, 0): constant}
    if first_star is not None:
        shifts, _ = _lattice_star(reciprocal)
        harmonics.update((shift, first_star) for shift in shifts)
    return harmonics


def _shell_indices(star, shells):
    frontier = {(0, 0)}
    waves = [(0, 0)]
    seen = frontier.copy()
    for _ in range(shells):
        frontier = {(i + n, j + m) for i, j in frontier for n, m in star} - seen
        waves.extend(sorted(frontier))
        seen.update(frontier)
    return np.array(waves, dtype=int)


def _tunneling_matrix(harmonic, top_size, bottom_size):
    matrix = np.asarray(harmonic, dtype=complex)
    if matrix.ndim == 0:
        if top_size != bottom_size:
            raise ValueError("different layer sizes require explicit tunneling matrices")
        matrix = matrix * np.eye(top_size)
    if matrix.shape != (top_size, bottom_size) or not np.all(np.isfinite(matrix)):
        raise ValueError(f"tunneling must have shape ({top_size}, {bottom_size}) and be finite")
    matrix = matrix.copy()
    matrix.setflags(write=False)
    return matrix


def _integer_shift(shift):
    try:
        n, m = map(operator.index, shift)
    except (TypeError, ValueError) as exc:
        raise ValueError("harmonic keys must be integer pairs") from exc
    return n, m


def _intralayer_potential(waves, size, harmonics, uniform):
    """Return U(G-G') in one layer. One of ±g suffices for each pair."""
    if not harmonics:
        return None
    result = np.zeros((len(waves) * size,) * 2, complex)
    amplitudes = {}
    for shift, amplitude in harmonics.items():
        n, m = _integer_shift(shift)
        if uniform and (n, m) != (0, 0):
            raise ValueError("zero twist accepts only the G=0 harmonic")
        matrix = _tunneling_matrix(amplitude, size, size)
        reverse = amplitudes.get((-n, -m), matrix if (n, m) == (0, 0) else None)
        if reverse is not None and not np.allclose(matrix, reverse.conj().T):
            raise ValueError("intralayer amplitudes must satisfy U(-g)=U(g)†")
        if reverse is not None:
            matrix = (matrix + reverse.conj().T) / 2
        amplitudes[(n, m)] = matrix
        amplitudes[(-n, -m)] = matrix.conj().T
    indices = {tuple(wave): i for i, wave in enumerate(waves)}
    for (n, m), matrix in amplitudes.items():
        for (i, j), row in indices.items():
            col = indices.get((i - n, j - m))
            if col is not None:
                result[row * size:(row + 1) * size, col * size:(col + 1) * size] = matrix
    return result


class MoireBilayer:
    """Truncated continuum bilayer of two Bloch Hamiltonians.

    Each layer supplies dimension=2, size, reciprocal, and hamiltonian(k).
    Primitive reciprocal rows must match before rotating. Each shell adds
    plane waves one step further in the first-star graph, starting at G=0.

    harmonics[(n, m)] is the top-to-bottom Fourier amplitude T(n*g1 + m*g2).
    H[top,G; bottom,G'] = T(G-G'). An amplitude is a (top.size, bottom.size)
    complex matrix in the layer Bloch bases. A scalar acts as identity in
    equal-size bases. The reverse block is fixed by Hermiticity. Bias is
    E_top - E_bottom, implemented as shifts of +bias/2 and -bias/2.
    potential_top and potential_bottom map harmonics to intralayer U(G-G').
    Supply one or both of ±g per pair; the missing partner is inferred by
    Hermiticity. U(0) is Hermitian. All potentials use the final layer basis.

    anchor_top and anchor_bottom are independent local Cartesian momenta at
    k=G=0, defaulting to Γ. As a shortcut, valley specifies one parent-lattice
    Cartesian valley shared by both layers. In this gauge k=0 is the mean of
    the rotated valley centers, and both layers use a common absolute momentum.

    By default the layer spin basis is treated as a common lab basis, with only
    the momentum argument rotated. With local_spin_frame=True, both layers
    must use orbital⊗physical spin with spin innermost. Then rotate the spin
    matrices (including the local Néel vector and SOC) to the common lab basis.
    T is always supplied in the final basis. Additional orbital rotations or
    basis gauges, if needed, are encoded by the caller in the models and T.

    This is a long-wavelength model with prescribed Fourier amplitudes. The
    finite plane-wave truncation can break mini-zone periodicity at its edge.
    A zero-twist bilayer is represented with shells=0 and only G=0 tunneling.
    """

    dimension = 2

    def __init__(
        self,
        top,
        bottom=None,
        *,
        twist_deg=2.0,
        shells=1,
        harmonics=None,
        potential_top=None,
        potential_bottom=None,
        anchor_top=None,
        anchor_bottom=None,
        valley=None,
        bias=0.0,
        local_spin_frame=False,
    ):
        bottom = top if bottom is None else bottom
        if top.dimension != 2 or bottom.dimension != 2:
            raise ValueError("both layers must be two-dimensional")
        top_size, bottom_size = operator.index(top.size), operator.index(bottom.size)
        if min(top_size, bottom_size) < 1:
            raise ValueError("layer sizes must be positive")
        parent_reciprocal = np.asarray(top.reciprocal, dtype=float)
        star, geometry = _lattice_star(parent_reciprocal)
        bottom_reciprocal = np.asarray(bottom.reciprocal, dtype=float)
        scale = np.max(np.abs(parent_reciprocal))
        if bottom_reciprocal.shape != (2, 2) or not np.allclose(
            parent_reciprocal, bottom_reciprocal, rtol=1e-10, atol=scale * 1e-12
        ):
            raise ValueError("both layers must use the same parent reciprocal basis")
        shells = operator.index(shells)
        if shells < 0:
            raise ValueError("shells must be nonnegative")
        twist_deg, bias = float(twist_deg), float(bias)
        if not np.isfinite(twist_deg) or not -180 < twist_deg < 180:
            raise ValueError("twist_deg must lie between -180 and 180 degrees")
        if not np.isfinite(bias):
            raise ValueError("bias must be finite")
        if twist_deg == 0 and shells:
            raise ValueError("zero twist requires shells=0")
        if local_spin_frame and (top_size % 2 or bottom_size % 2):
            raise ValueError("local_spin_frame requires complete physical spin doublets")

        half_angle = np.deg2rad(twist_deg) / 2
        top_rotation, bottom_rotation = _rotation(half_angle), _rotation(-half_angle)
        reciprocal = 2 * np.sin(half_angle) * (parent_reciprocal @ _rotation(np.pi / 2).T)
        if valley is not None:
            if anchor_top is not None or anchor_bottom is not None:
                raise ValueError("choose valley or explicit anchors")
            origin = _vector2(valley, "valley") * np.cos(half_angle)
            anchor_top, anchor_bottom = origin @ top_rotation, origin @ bottom_rotation
        self.top = top
        self.bottom = bottom
        self.twist_deg = twist_deg
        self.shells = shells
        self.bias = bias
        self.geometry = geometry
        self.reciprocal = reciprocal
        self.plane_waves = _shell_indices(star, shells)
        self.offsets = self.plane_waves @ reciprocal
        self.reciprocal.setflags(write=False)
        self.plane_waves.setflags(write=False)
        self.offsets.setflags(write=False)
        self._top_dimension = len(self.plane_waves) * top_size
        self.size = len(self.plane_waves) * (top_size + bottom_size)
        self._potentials = (
            _intralayer_potential(self.plane_waves, top_size, potential_top, twist_deg == 0),
            _intralayer_potential(self.plane_waves, bottom_size, potential_bottom, twist_deg == 0),
        )

        self._layers = []
        for model, size, angle, rotation, anchor, shift, first in (
            (top, top_size, half_angle, top_rotation, anchor_top, bias / 2, 0),
            (bottom, bottom_size, -half_angle, bottom_rotation, anchor_bottom, -bias / 2, self._top_dimension),
        ):
            anchor = np.zeros(2) if anchor is None else _vector2(anchor, "anchor").copy()
            indices = first + np.arange(len(self.plane_waves) * size).reshape(-1, size)
            spin_phases = None
            if local_spin_frame:
                phase = np.tile(np.exp(-0.5j * angle * np.array([1, -1])), size // 2)
                spin_phases = phase[:, None] * phase.conj()[None, :]
            self._layers.append((model, rotation, anchor, shift * np.eye(size), indices, spin_phases))

        amplitudes = {}
        self._coupling = np.zeros((self._top_dimension, self.size - self._top_dimension), dtype=complex)
        wave_index = {tuple(wave): i for i, wave in enumerate(self.plane_waves)}
        for shift, amplitude in ({} if harmonics is None else harmonics).items():
            n, m = _integer_shift(shift)
            if twist_deg == 0 and (n, m) != (0, 0):
                raise ValueError("zero twist accepts only the G=0 harmonic")
            matrix = _tunneling_matrix(amplitude, top_size, bottom_size)
            amplitudes[(n, m)] = matrix
            for (i, j), top_index in wave_index.items():
                bottom_index = wave_index.get((i - n, j - m))
                if bottom_index is None:
                    continue
                rows = slice(top_index * top_size, (top_index + 1) * top_size)
                cols = slice(bottom_index * bottom_size, (bottom_index + 1) * bottom_size)
                self._coupling[rows, cols] = matrix
        self.harmonics = MappingProxyType(amplitudes)

    def hamiltonian(self, k):
        """Return H(k) for one Cartesian mini-momentum or an array (..., 2)."""
        k = np.asarray(k, dtype=float)
        if k.ndim == 0 or k.shape[-1] != 2 or not np.all(np.isfinite(k)):
            raise ValueError("k must be a finite array of shape (..., 2)")
        global_momenta = k[..., None, :] + self.offsets
        out = np.zeros(k.shape[:-1] + (self.size, self.size), dtype=complex)
        for model, rotation, anchor, bias, indices, spin_phases in self._layers:
            local_momenta = global_momenta @ rotation + anchor
            blocks = model.hamiltonian(local_momenta) + bias
            if spin_phases is not None:
                blocks = blocks * spin_phases
            out[..., indices[:, :, None], indices[:, None, :]] = blocks
        split = self._top_dimension
        if self._potentials[0] is not None:
            out[..., :split, :split] += self._potentials[0]
        if self._potentials[1] is not None:
            out[..., split:, split:] += self._potentials[1]
        out[..., :split, split:] = self._coupling
        out[..., split:, :split] = self._coupling.conj().T
        return out

    def high_symmetry_path(self):
        """Return geometric path labels and Cartesian nodes in the mini-zone."""
        if self.twist_deg == 0:
            raise ValueError("zero twist has no moiré mini-zone")
        g1, g2 = self.reciprocal
        gamma = np.zeros(2)
        if self.geometry == "hexagonal":
            corner = (g1 + g2) / 3 if g1 @ g2 > 0 else (2 * g1 + g2) / 3
            return ("Γ", "M", "K", "Γ"), np.array([gamma, g1 / 2, corner, gamma])
        if self.geometry == "rectangular":
            return ("Γ", "X", "M", "Y", "Γ"), np.array(
                [gamma, g1 / 2, (g1 + g2) / 2, g2 / 2, gamma]
            )
        return ("Γ", "X", "M", "Γ"), np.array([gamma, g1 / 2, (g1 + g2) / 2, gamma])
