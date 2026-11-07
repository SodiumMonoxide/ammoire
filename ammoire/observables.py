"""Operators and degeneracy-aware diagnostics in an orthonormal basis."""

import numpy as np


def layer_operator(bilayer, top, bottom):
    """Lift k-independent layer operators to (layer, G, orbital, label)."""
    n = len(bilayer.plane_waves)
    top, bottom = (np.asarray(op, complex) for op in (top, bottom))
    dt, db = bilayer.top.size, bilayer.bottom.size
    if top.shape != (dt, dt) or bottom.shape != (db, db):
        raise ValueError("operator shape must match its layer's basis")
    out = np.zeros((bilayer.size, bilayer.size), complex)
    split = n * dt
    out[:split, :split] = np.kron(np.eye(n), top)
    out[split:, split:] = np.kron(np.eye(n), bottom)
    return out


def commutator_norm(hamiltonian, operator):
    """Frobenius norm of [H, Q], in H units for dimensionless Q."""
    return np.linalg.norm(hamiltonian @ operator - operator @ hamiltonian, axis=(-2, -1))


def expectation(vectors, operator):
    """Expectation per column of an orthonormal eigenvector matrix."""
    return np.sum(vectors.conj() * (operator @ vectors), axis=-2).real


def resolve_bands(hamiltonian, operator, *, tolerance=1e-10):
    """Return ascending E, basis-resolved <Q>, and eigenvectors at one k.

    Within numerically degenerate groups, diagonalize the projection of Q.
    If Q does not commute with H, the result is an expectation, not a sector.
    The tolerance is absolute in H units.

    H and Q must be finite and Hermitian to 1e-12 plus 16 machine epsilons
    times their largest entry. Accepted rounding errors are symmetrized.
    This input check is independent of the degeneracy grouping tolerance.
    """
    if tolerance < 0 or not np.isfinite(tolerance):
        raise ValueError("tolerance must be finite and nonnegative")
    hamiltonian, operator = (np.asarray(matrix, complex) for matrix in (hamiltonian, operator))
    if (hamiltonian.ndim != 2 or hamiltonian.shape != operator.shape
            or hamiltonian.shape[0] != hamiltonian.shape[1]
            or hamiltonian.shape[0] == 0):
        raise ValueError("H and Q must be equally sized nonempty square matrices")
    for name, matrix in (("H", hamiltonian), ("Q", operator)):
        if not np.all(np.isfinite(matrix)):
            raise ValueError(f"{name} must be finite")
        limit = 1e-12 + 16 * np.finfo(float).eps * np.max(np.abs(matrix))
        if np.max(np.abs(matrix - matrix.conj().T)) > limit:
            raise ValueError(f"{name} must be Hermitian")
    hamiltonian = 0.5 * hamiltonian + 0.5 * hamiltonian.conj().T
    operator = 0.5 * operator + 0.5 * operator.conj().T
    energies, vectors = np.linalg.eigh(hamiltonian)
    first = 0
    while first < len(energies):
        last = np.searchsorted(energies, energies[first] + tolerance, side="right")
        if last - first > 1:
            group = vectors[:, first:last]
            _, rotation = np.linalg.eigh(group.conj().T @ operator @ group)
            vectors[:, first:last] = group @ rotation
        first = last
    return energies, expectation(vectors, operator), vectors
