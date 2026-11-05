"""Import a cell-gauge Wannier90 Hamiltonian into the Bloch model interface.

The lattice contains Cartesian real-space primitive vectors as rows. With
the lattice in angstroms, k is in inverse angstroms and H(k) in eV. Spin and
orbital order follow the file. The Hamiltonian alone specifies no spin,
hyperspin, orbital-position, or interlayer tunneling matrix elements.

File specification: https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/files/#seedname_hrdat
"""

from dataclasses import dataclass
from operator import index
from pathlib import Path

import numpy as np


def _check_wsvec(path):
    if not path.name.endswith("_hr.dat"):
        return False
    sidecar = path.with_name(path.name.removesuffix("_hr.dat") + "_wsvec.dat")
    if not sidecar.is_file():
        return False
    with sidecar.open(encoding="utf-8-sig") as file:
        header = "".join(file.readline().lower().split())
    expanded = "write_ndegen_applied=.true." in header
    if not expanded and "use_ws_distance=.false." not in header:
        raise ValueError(
            f"{sidecar.name} indicates possible orbital-dependent Wigner–Seitz shifts. "
            "Export with write_ndegen_applied=true or use_ws_distance=false."
        )
    return expanded


def _read_hr(path):
    with path.open(encoding="utf-8-sig") as file:
        if not file.readline():
            raise ValueError("empty Wannier90 Hamiltonian file")
        records = (line.rstrip("\r\n") for line in file if line.strip())
        try:
            size, count = int(next(records)), int(next(records))
        except (StopIteration, ValueError) as exc:
            raise ValueError("expected the Wannier count and R-vector count") from exc
        if min(size, count) < 1:
            raise ValueError("the Wannier count and R-vector count must be positive")
        degeneracies = []
        while len(degeneracies) < count:
            try:
                degeneracies.extend(map(int, next(records).split()))
            except (StopIteration, ValueError) as exc:
                raise ValueError("invalid Wigner–Seitz degeneracy block") from exc
        if len(degeneracies) != count or min(degeneracies) < 1:
            raise ValueError("expected one positive Wigner–Seitz degeneracy per R vector")

        translations = []
        hopping = np.empty((count, size, size), dtype=complex)
        for block in hopping:
            seen = np.zeros((size, size), dtype=bool)
            translation = None
            for _ in range(size**2):
                line = next(records, None)
                if line is None:
                    raise ValueError("truncated Wannier90 hopping block")
                fields = line.split()
                if len(fields) != 7 and len(line) == 49:
                    # Standard Fortran (5I5,2F12.6), including touching fields.
                    fields = [line[i:i + 5] for i in range(0, 25, 5)] + [line[25:37], line[37:]]
                try:
                    rx, ry, rz, i, j = map(int, fields[:5])
                    re, im = (float(part.replace("D", "E").replace("d", "e")) for part in fields[5:])
                except ValueError as exc:
                    raise ValueError(f"invalid Wannier90 hopping record {line!r}") from exc
                r = (rx, ry, rz)
                if translation is None:
                    translation = r
                if translation != r or not (1 <= i <= size and 1 <= j <= size):
                    raise ValueError("each R block must contain one complete Wannier matrix")
                if seen[i - 1, j - 1]:
                    raise ValueError(f"duplicate matrix entry at R={r}, m={i}, n={j}")
                seen[i - 1, j - 1] = True
                block[i - 1, j - 1] = complex(re, im)
            translations.append(translation)
        if next(records, None) is not None:
            raise ValueError("extra records after the declared Wannier90 hopping blocks")
    if not np.all(np.isfinite(hopping)):
        raise ValueError("Wannier90 hoppings must be finite")
    return translations, np.array(degeneracies), hopping


def _hermitian_weights(translations, degeneracies, hopping, tolerance):
    indices = {r: i for i, r in enumerate(translations)}
    if len(indices) != len(translations):
        raise ValueError("duplicate Wannier90 R blocks")
    hopping /= degeneracies[:, None, None]
    for r, i in indices.items():
        j = indices.get(tuple(-x for x in r))
        if j is None or degeneracies[i] != degeneracies[j]:
            raise ValueError(f"expected a matching -R block with equal degeneracy for R={r}")
        if i > j:
            continue
        reverse = hopping[j].conj().T
        residual = np.max(np.abs(hopping[i] - reverse))
        if residual > tolerance:
            raise ValueError(f"H(R) != H(-R)† at R={r}, residual {residual:.3g} eV")
        # Remove residual accepted rounding error to preserve a Hermitian Bloch matrix.
        hopping[i] = 0.5 * hopping[i] + 0.5 * reverse
        hopping[j] = hopping[i].conj().T


@dataclass(frozen=True)
class WannierModel:
    """Bloch Hamiltonian in the orthonormal Wannier basis of the input file."""

    name: str
    dimension: int
    lattice: np.ndarray
    reciprocal: np.ndarray
    translations: np.ndarray
    degeneracies: np.ndarray
    weighted_hoppings: np.ndarray
    _cartesian_translations: np.ndarray

    @property
    def size(self):
        return self.weighted_hoppings.shape[1]

    @classmethod
    def from_hr(cls, path, lattice, *, dimension, tolerance=1e-6):
        """Load standard *_hr.dat with H(k) = Σ_R exp(i k·R) H(R) / ndegen(R).

        Use files produced with use_ws_distance=false or with distance shifts
        already included by write_ndegen_applied=true. The hr file alone cannot
        reveal omitted shifts. A standard adjacent *_wsvec.dat header is checked
        when present. For dimension=2, supply a 2×2 in-plane lattice and all
        third integer translation components must be zero. Tolerance (eV)
        applies to individual weighted hopping entries before symmetrization.
        """
        dimension = index(dimension)
        lattice = np.array(lattice, dtype=float, copy=True)
        if dimension not in (2, 3) or lattice.shape != (dimension, dimension):
            raise ValueError("dimension=2 or 3 requires a matching square row-vector lattice")
        if not np.all(np.isfinite(lattice)) or np.linalg.matrix_rank(lattice) != dimension:
            raise ValueError("lattice vectors must be finite and linearly independent")
        tolerance = float(tolerance)
        if not np.isfinite(tolerance) or tolerance < 0:
            raise ValueError("tolerance must be finite and nonnegative")
        path = Path(path)
        expanded = _check_wsvec(path)
        translations, degeneracies, hopping = _read_hr(path)
        if expanded and np.any(degeneracies != 1):
            raise ValueError("expanded Wannier90 exports must have unit degeneracies")
        if dimension == 2 and any(r[2] != 0 for r in translations):
            raise ValueError("dimension=2 requires Rz=0 in every Wannier90 hopping block")
        _hermitian_weights(translations, degeneracies, hopping, tolerance)
        translations = np.array(translations)[:, :dimension]
        reciprocal = 2 * np.pi * np.linalg.inv(lattice).T
        cartesian = translations @ lattice
        for array in (lattice, reciprocal, translations, degeneracies, hopping, cartesian):
            array.setflags(write=False)
        return cls(path.stem.removesuffix("_hr"), dimension, lattice, reciprocal,
                   translations, degeneracies, hopping, cartesian)

    def hamiltonian(self, k):
        """Evaluate at Cartesian k of shape (..., dimension) in the cell gauge."""
        k = np.asarray(k, dtype=float)
        if k.ndim == 0 or k.shape[-1] != self.dimension or not np.all(np.isfinite(k)):
            raise ValueError(f"k must have finite shape (..., {self.dimension})")
        phase = np.exp(1j * (k @ self._cartesian_translations.T))
        return np.tensordot(phase, self.weighted_hoppings, axes=1)


def _self_test():
    from tempfile import TemporaryDirectory

    onsite = np.array([[0.2, 0.15j], [-0.15j, -0.3]])
    forward = np.array([[0.4, 0.1 + 0.2j], [0.3 - 0.4j, 0.5]])
    translations = ((-1, 1, 0), (0, 0, 0), (1, -1, 0))
    matrices = (forward.conj().T, onsite, forward)
    lines = ["synthetic Wannier90 fixture", "2", "3", "2 1 2"]
    for r, matrix in zip(translations, matrices):
        for j in range(2):
            for i in range(2):
                re, im = matrix[i, j].real, matrix[i, j].imag
                lines.append(f"{r[0]} {r[1]} {r[2]} {i + 1} {j + 1} {re} {im}")
    with TemporaryDirectory(prefix="wannier-fixture-", dir=Path.cwd()) as directory:
        path = Path(directory) / "fixture_hr.dat"
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        lattice = np.array([[2, 0], [0.3, 1.5]])
        model = WannierModel.from_hr(path, lattice, dimension=2)
        k = np.array([[0.2, 0.4], [-0.1, 0.3]])
        phase = np.exp(1j * k @ (lattice[0] - lattice[1]))[:, None, None]
        expected = onsite + (phase * forward + phase.conj() * forward.conj().T) / 2
        np.testing.assert_allclose(model.hamiltonian(k), expected, atol=1e-14)
        np.testing.assert_allclose(model.hamiltonian(k + model.reciprocal[0]), expected, atol=1e-14)
        lattice3 = np.eye(3)
        lattice3[:2, :2] = lattice
        model3 = WannierModel.from_hr(path, lattice3, dimension=3)
        np.testing.assert_allclose(model3.hamiltonian(np.c_[k, [0.2, 1.3]]), expected, atol=1e-14)
        path.write_text("\n".join(lines[:-1] + ["1 -1 0 2 2 2.0 0.0"]), encoding="utf-8")
        with np.testing.assert_raises_regex(ValueError, "residual"):
            WannierModel.from_hr(path, lattice, dimension=2)
        path.with_name("fixture_wsvec.dat").write_text("use_ws_distance=.true.\n", encoding="utf-8")
        with np.testing.assert_raises_regex(ValueError, "Wigner"):
            WannierModel.from_hr(path, lattice, dimension=2)
    print("Wannier90 self-test passed")


if __name__ == "__main__":
    _self_test()
