"""Spin Laue cosets for ferroelectric switching in collinear altermagnets.

The classification follows Eqs. (2)–(9) of Sun et al., Nat. Commun. 17,
3101 (2026), doi:10.1038/s41467-026-69635-2. Matrices act on Cartesian
momentum columns. For an array of momentum rows, use k @ operation.T.
"""

from collections import deque
from operator import index

import numpy as np


_ATOL = 1e-8


def _integer(value, name):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be an integer")
    try:
        return index(value)
    except TypeError as error:
        raise ValueError(f"{name} must be an integer") from error


def _matching_sign(operations, matrix):
    for member, sign in operations:
        if np.allclose(matrix, member, atol=_ATOL, rtol=0):
            return sign
    return None


class SpinLaueGroup:
    """Construct the signed spatial actions of a collinear spin Laue group.

    Supply generators of the full band symmetry group as (matrix, sign)
    pairs. A sign of +1 preserves the spin channel and −1 exchanges it.
    Band inversion (−I, +1) is implicit in the nonrelativistic theory.
    The signs must define an index-two subgroup of the spatial group.

    For a 2D layer in the xy plane, both 2×2 matrices and layer-preserving
    3×3 operations are accepted. Classification uses their action on the
    in-plane momentum. Bulk models use full 3×3 matrices.
    """

    def __init__(self, dimension, generators):
        self._dimension = _integer(dimension, "dimension")
        if self._dimension not in (2, 3):
            raise ValueError("dimension must be 2 or 3")
        identity = np.eye(self._dimension)
        signed_generators = [(-identity, 1)]
        for matrix, sign in generators:
            sign = _integer(sign, "spin sign")
            if sign not in (-1, 1):
                raise ValueError("spin sign must be +1 or -1")
            signed_generators.append((self.momentum_action(matrix), sign))
        if all(sign == 1 for _, sign in signed_generators):
            raise ValueError("at least one spin-exchanging generator is required")

        members = [(identity, 1)]
        pending = deque(members)
        max_order = 12 if self._dimension == 2 else 48
        while pending:
            member, member_sign = pending.popleft()
            for generator, generator_sign in signed_generators:
                product = member @ generator
                sign = member_sign * generator_sign
                known_sign = _matching_sign(members, product)
                if known_sign is not None:
                    if sign != known_sign:
                        raise ValueError("conflicting spin signs for the same momentum action")
                    continue
                # A lattice automorphism has integer trace in every basis.
                trace = np.trace(product)
                if abs(trace - round(trace)) > self._dimension * _ATOL or len(members) >= max_order:
                    raise ValueError("generators do not close as a crystallographic point group")
                members.append((product, sign))
                pending.append((product, sign))
        self._operations = tuple(members)

    @property
    def dimension(self):
        """Number of momentum components used for classification."""
        return self._dimension

    @property
    def operations(self):
        """All (momentum action, spin sign) pairs in the generated group."""
        return tuple((matrix.copy(), sign) for matrix, sign in self._operations)

    def momentum_action(self, operation):
        """Return the Cartesian orthogonal action in the momentum space used."""
        if np.iscomplexobj(operation):
            raise ValueError("Cartesian operation must be real")
        matrix = np.asarray(operation, dtype=float)
        if matrix.shape not in ((self.dimension, self.dimension), (3, 3)):
            raise ValueError(f"invalid operation shape for dimension {self.dimension}")
        if not np.isfinite(matrix).all():
            raise ValueError("Cartesian operation must be finite")
        if not np.allclose(matrix.T @ matrix, np.eye(len(matrix)), atol=_ATOL, rtol=0):
            raise ValueError("Cartesian operation must be orthogonal")
        if self.dimension == 2 and matrix.shape == (3, 3):
            if np.any(np.abs(matrix[2, :2]) > _ATOL) or np.any(np.abs(matrix[:2, 2]) > _ATOL):
                raise ValueError("a 2D operation must preserve the xy layer plane")
            matrix = matrix[:2, :2]
        return matrix.copy()

    def classify(self, operation):
        """Classify a specified domain operation that preserves magnetic order.

        Return 'I' for the preserving subgroup, 'II' for the exchanging
        coset, or 'III' for an action outside the supplied group. Membership
        does not establish that a polarization switch exists in a material.
        """
        sign = _matching_sign(self._operations, self.momentum_action(operation))
        return {1: "I", -1: "II", None: "III"}[sign]
