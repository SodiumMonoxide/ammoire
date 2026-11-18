"""Finite-sample band tests for spin–ferroelectric domain transformations.

Momenta and the operation use a common Cartesian frame. The two no-SOC,
collinear domains must use the same physical energy units and zero. All
comparisons use ascending energy ranks independently within each spin
sector; this remains well defined at band crossings and degeneracies.
"""

import numpy as np

from .commensurate import CommensurateTriangularBilayer
from .compensated import CompensatedTriangularBilayer
from .models import MomentumSlice, PAULI, _momentum, _unit_vector
from .moire import MoireBilayer
from .observables import layer_operator


def _spin_component(model, axis):
    """Physical Pauli spin along one common laboratory axis."""
    if isinstance(model, MomentumSlice):
        return _spin_component(model.parent, axis)
    if isinstance(model, MoireBilayer):
        return layer_operator(model, _spin_component(model.top, axis),
                              _spin_component(model.bottom, axis))
    if any(value != 0 for name, value in getattr(model, 'params', {}).items()
           if name in ('soc', 'soc_x', 'soc_y', 'soc_z')):
        raise ValueError('the spin–ferroelectric audit requires a no-SOC Hamiltonian')
    if isinstance(model, (CommensurateTriangularBilayer, CompensatedTriangularBilayer)):
        return np.kron(np.eye(model.size // 2), np.einsum('i,ijk->jk', axis, PAULI))
    operators = getattr(model, 'spin_operators', None)
    if operators is None or np.shape(operators) != (3, model.size, model.size):
        raise ValueError('each domain must supply physical Pauli spin_operators')
    return np.einsum('i,ijk->jk', axis, operators)


def _sector_bases(model, axis):
    spin = _spin_component(model, axis)
    if not np.all(np.isfinite(spin)) or not np.allclose(spin, spin.conj().T, atol=1e-12, rtol=0):
        raise ValueError('the physical spin component must be finite and Hermitian')
    values, vectors = np.linalg.eigh(spin)
    half = model.size // 2
    if (half == 0 or model.size != 2 * half
            or not np.allclose(values[:half], -1, atol=1e-10, rtol=0)
            or not np.allclose(values[half:], 1, atol=1e-10, rtol=0)):
        raise ValueError('physical Pauli spin must have equally many +1 and -1 eigenvalues')
    return vectors[:, half:], vectors[:, :half]


def _sector_energies(model, momenta, bases):
    h = np.asarray(model.hamiltonian(momenta), dtype=complex)
    if h.shape != (len(momenta), model.size, model.size) or not np.all(np.isfinite(h)):
        raise ValueError('domain Hamiltonians must return finite square matrices for every k')
    matrix_tol = 1e-12 + 32 * np.finfo(float).eps * model.size * np.max(abs(h), axis=(-2, -1))
    if np.any(np.max(abs(h - h.conj().swapaxes(-2, -1)), axis=(-2, -1)) > matrix_tol):
        raise ValueError('domain Hamiltonians must be Hermitian')
    up, down = bases
    if np.any(np.max(abs(up.conj().T @ h @ down), axis=(-2, -1)) > matrix_tol):
        raise ValueError('Hamiltonian mixes physical spin sectors along the declared spin axis')
    return np.stack([np.linalg.eigvalsh(basis.conj().T @ h @ basis) for basis in bases], axis=1)


def _sample_inputs(kpoints, dimension, spin_axis, atol):
    points = _momentum(kpoints, dimension)
    if points.ndim != 2 or not len(points):
        raise ValueError('supply a nonempty array of k points with shape (n, dimension)')
    atol = float(atol)
    if not np.isfinite(atol) or atol < 0:
        raise ValueError('atol must be finite and nonnegative')
    return points, _unit_vector(spin_axis), atol


def audit_group_spectra(model, kpoints, group, *, spin_axis=(0, 0, 1), atol=1e-8):
    """Test every signed operation of a SpinLaueGroup on one sampled domain.

    Physical-spin, energy-rank, no-SOC and Cartesian conventions are the
    same as for audit_domain_spectra. A finite sample tests consistency
    with the supplied symmetry; agreement alone cannot determine the group.
    """
    if model.dimension != group.dimension:
        raise ValueError('group and model must have the same dimension')
    points, axis, atol = _sample_inputs(kpoints, model.dimension, spin_axis, atol)
    bases = _sector_bases(model, axis)
    reference = _sector_energies(model, points, bases)
    residual = 0.0
    operations = group.operations
    for operation, sign in operations:
        partner = _sector_energies(model, points @ operation.T, bases)
        if sign == -1:
            partner = partner[:, ::-1, :]
        residual = max(residual, float(np.max(np.abs(reference - partner))))
    return {'max_residual': residual, 'matches': residual <= atol, 'tested_operations': len(operations)}


def audit_domain_spectra(positive, negative, kpoints, operation, *, spin_axis=(0, 0, 1), atol=1e-8):
    """Compare independently supplied domain spectra under an active domain map.

    At row points k, test E_positive(k,s) = E_negative(k @ operation.T,s).
    Separately compare both domains at identical k with fixed or exchanged
    spin. Neither the magnetic order nor the ferroelectric operation is
    inferred from the bands. For a 2D section, supply the induced in-plane
    action and corresponding sections of both domains.

    Models implement dimension, size, batched hamiltonian(k), and three
    physical Pauli spin_operators in their own orthonormal basis. Native
    commensurate and moiré bilayers of physical-spin layers are supported.
    For custom models, the caller establishes the no-SOC assumption and
    aligns momentum frames, units, and energy references. Spin conservation
    is checked at the evaluated points. atol is absolute in energy units.
    Residuals are maximum absolute differences over all samples and ranks.
    """
    dimension = positive.dimension
    if dimension not in (2, 3) or negative.dimension != dimension or positive.size != negative.size:
        raise ValueError('domains must have matching dimension and band counts (2D or 3D)')
    for attr in ('energy_unit', 'length_unit'):
        units = [getattr(model, attr, None) for model in (positive, negative)]
        if None not in units and units[0] != units[1]:
            raise ValueError(f'domains must have a common {attr}')
    points, axis, atol = _sample_inputs(kpoints, dimension, spin_axis, atol)
    operation = np.asarray(operation)
    if (np.iscomplexobj(operation) or operation.shape != (dimension, dimension)
            or not np.all(np.isfinite(operation))
            or not np.allclose(operation.T @ operation, np.eye(dimension), atol=1e-8, rtol=0)):
        raise ValueError('operation must be a real orthogonal Cartesian matrix')
    plus = _sector_energies(positive, points, _sector_bases(positive, axis))
    minus_bases = _sector_bases(negative, axis)
    minus = _sector_energies(negative, points, minus_bases)
    mapped = _sector_energies(negative, points @ operation.T, minus_bases)

    def maximum(array):
        return float(np.max(np.abs(array)))

    mapped_residual = maximum(plus - mapped)
    preserved = maximum(plus - minus)
    reversed_spin = maximum(plus - minus[:, ::-1, :])
    positive_split = plus[:, 0] - plus[:, 1]
    negative_split = minus[:, 0] - minus[:, 1]
    max_plus, max_minus = maximum(positive_split), maximum(negative_split)
    if max(max_plus, max_minus) <= atol or max(preserved, reversed_spin) <= atol:
        relation = 'unresolved'
    elif preserved <= atol:
        relation = 'preserved'
    elif reversed_spin <= atol:
        relation = 'reversed'
    else:
        relation = 'neither'
    return {
        'sample_count': len(points),
        'bands_per_spin': positive.size // 2,
        'spin_axis': list(axis), 'atol': atol,
        'operation_cartesian': operation.tolist(),
        'mapped_same_spin_residual': mapped_residual,
        'mapped_same_spin_matches': mapped_residual <= atol,
        'same_k_preserved_residual': preserved,
        'same_k_reversed_residual': reversed_spin,
        'same_k_relation': relation,
        'max_splitting_positive': max_plus, 'max_splitting_negative': max_minus,
        'same_k_splitting_preserved_residual': maximum(positive_split - negative_split),
        'same_k_splitting_reversed_residual': maximum(positive_split + negative_split),
    }
