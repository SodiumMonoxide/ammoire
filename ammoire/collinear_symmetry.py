"""Audit collinear spin spectra at specified Cartesian momentum partners.

Run ``python -m ammoire.collinear_symmetry MANIFEST.json [--out REPORT.json]``.
The manifest supplies the direct lattice as three rows in angstroms and
explicitly names each probe and symmetry operation. The EIGENVAL files
must come from static ISPIN=2 calculations on a common cell and energy
reference. The audit compares ascending energy ranks independently for
each physical spin. Orbital identities and spin textures require separate
wavefunction or orbital information.
"""

import argparse
from dataclasses import dataclass
import json
from operator import index
from pathlib import Path

import numpy as np

from .vasp_bridge import read_eigenval


@dataclass(frozen=True)
class CollinearProbe:
    """One Cartesian k point in inverse angstroms and (bands, up/down) energies in eV."""

    k_cartesian: np.ndarray
    energies: np.ndarray


@dataclass(frozen=True)
class SymmetryPair:
    """Map k_source to k_target with a Cartesian operation and optional G shift.

    For row vectors, ``k_target = k_source @ operation_cartesian.T + G``.
    The shift gives integer coefficients in the supplied reciprocal row
    basis. Set swap_spin for an operation that exchanges up and down.
    """

    source: str
    target: str
    operation_cartesian: np.ndarray
    swap_spin: bool
    reciprocal_shift: tuple[int, int, int] = (0, 0, 0)


def _real_array(value, shape, name):
    if np.iscomplexobj(value):
        raise ValueError(f"{name} must be real")
    array = np.asarray(value, dtype=float)
    if (array.ndim != len(shape) or any(n is not None and n != actual for n, actual in zip(shape, array.shape))
            or not np.all(np.isfinite(array))):
        raise ValueError(f"{name} must be a finite real array with shape {shape}")
    return array


def _integer(value, name):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be an integer")
    try:
        return index(value)
    except TypeError as exc:
        raise ValueError(f"{name} must be an integer") from exc


def _invert_lattice(value, name):
    matrix = _real_array(value, (3, 3), name)
    try:
        inverse = np.linalg.inv(matrix)
    except np.linalg.LinAlgError as exc:
        raise ValueError(f"{name} must be invertible") from exc
    if not np.all(np.isfinite(inverse)):
        raise ValueError(f"{name} has a non-finite inverse")
    return matrix, inverse


def _maximum(values):
    band, spin = np.unravel_index(np.argmax(np.abs(values)), values.shape)
    difference = float(values[band, spin])
    result = {"max_abs_eV": abs(difference), "signed_eV": difference, "band_rank_1based": int(band + 1)}
    if values.shape[1] == 2:
        result["source_spin"] = ("up", "down")[spin]
    return result


def audit_collinear_symmetry(probes, pairs, *, occupied_bands, nodes=(), bands=(),
                             reciprocal=None, k_tolerance_inverse_angstrom=1e-7):
    """Report spin symmetry residuals and same-k differences in eV.

    ``probes`` maps names to CollinearProbe objects. ``pairs`` contains
    SymmetryPair objects. Every pair must have matching band counts and
    Cartesian momenta within the specified tolerance. ``reciprocal`` gives
    3D reciprocal row vectors in inverse angstroms when shifts are used.

    ``occupied_bands`` is the supplied count of occupied ranks in each
    spin channel. It defines the reported low-energy prefix. ``bands``
    selects additional one-based ranks for reporting signed same-k splittings.
    ``nodes`` designates probes for an explicit same-k degeneracy check.
    Include a self-pair with spin swap to also check the node momentum's
    invariance under the specified operation.

    Pair residuals are source minus target after the requested spin map.
    Same-k differences are up minus down at identical ascending ranks.
    """
    occupied = _integer(occupied_bands, "occupied_bands")
    selected = tuple(_integer(rank, "selected band rank") for rank in bands)
    if occupied < 1 or any(rank < 1 for rank in selected) or len(set(selected)) != len(selected):
        raise ValueError("occupied_bands must be positive and selected band ranks positive and unique")
    tolerance = float(k_tolerance_inverse_angstrom)
    if not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError("momentum tolerance must be finite and nonnegative")
    if reciprocal is not None:
        reciprocal, _ = _invert_lattice(reciprocal, "reciprocal")

    momenta, spectra, probe_reports = {}, {}, {}
    for name, probe in probes.items():
        if not isinstance(name, str) or not name.strip():
            raise ValueError("probe names must be nonempty strings")
        momenta[name] = _real_array(probe.k_cartesian, (3,), f"k at {name}")
        energies = _real_array(probe.energies, (None, 2), f"energies at {name}")
        if len(energies) < max((occupied, *selected)):
            raise ValueError(f"{name}: requested occupied or selected rank exceeds {len(energies)} bands")
        spectra[name] = np.sort(energies, axis=0)
        difference = spectra[name][:, :1] - spectra[name][:, 1:]
        probe_reports[name] = {
            "k_cartesian_inverse_angstrom": momenta[name].tolist(),
            "band_count_per_spin": len(energies),
            "same_k_max_occupied": _maximum(difference[:occupied]),
            "same_k_max_all": _maximum(difference),
            "selected_rank_splitting_eV": {str(rank): float(difference[rank - 1, 0]) for rank in selected},
        }
    if not probe_reports:
        raise ValueError("provide at least one probe")

    pair_reports = []
    for pair in pairs:
        label = f"{pair.source} -> {pair.target}"
        if pair.source not in momenta or pair.target not in momenta:
            raise ValueError(f"{label}: unknown probe name")
        if not isinstance(pair.swap_spin, (bool, np.bool_)):
            raise ValueError(f"{label}: swap_spin must be a boolean")
        operation = _real_array(pair.operation_cartesian, (3, 3), f"Cartesian operation for {label}")
        if not np.allclose(operation @ operation.T, np.eye(3), atol=1e-8, rtol=0):
            raise ValueError(f"{label}: Cartesian operation must be orthogonal")
        shift = tuple(_integer(value, f"reciprocal shift for {label}") for value in pair.reciprocal_shift)
        if len(shift) != 3 or (any(shift) and reciprocal is None):
            raise ValueError(f"{label}: provide three integer shift coefficients and a basis for a nonzero shift")
        predicted = momenta[pair.source] @ operation.T
        if any(shift):
            predicted += np.asarray(shift) @ reciprocal
        k_error = float(np.linalg.norm(predicted - momenta[pair.target]))
        if k_error > tolerance:
            raise ValueError(f"{label}: Cartesian momentum mismatch {k_error:.6g} A^-1 exceeds {tolerance:.6g}")
        source, target = spectra[pair.source], spectra[pair.target]
        if source.shape != target.shape:
            raise ValueError(f"{label}: full-spectrum comparison requires equal band counts")
        difference = source - (target[:, ::-1] if pair.swap_spin else target)
        pair_reports.append({
            "source": pair.source,
            "target": pair.target,
            "spin_relation": "exchange" if pair.swap_spin else "same_spin",
            "operation_cartesian": operation.tolist(),
            "reciprocal_shift": list(shift),
            "k_error_inverse_angstrom": k_error,
            "band_count_per_spin": len(source),
            "max_occupied": _maximum(difference[:occupied]),
            "max_all": _maximum(difference),
        })

    node_reports = {}
    for name in nodes:
        if name not in probe_reports or name in node_reports:
            raise ValueError(f"unknown or repeated node probe {name!r}")
        node_reports[name] = {
            "max_occupied": probe_reports[name]["same_k_max_occupied"],
            "max_all": probe_reports[name]["same_k_max_all"],
        }
    return {
        "energy_convention": "ascending within-spin energy ranks in eV",
        "occupied_bands_per_spin": occupied,
        "k_tolerance_inverse_angstrom": tolerance,
        "probes": probe_reports,
        "pairs": pair_reports,
        "nodes": node_reports,
    }


def audit_manifest(path):
    """Read a JSON manifest and static collinear EIGENVAL probes.

    Each item under ``probes`` specifies ``eigenval`` and optionally a
    one-based ``k_index``. The index defaults to one for a one-point file.
    Paths are interpreted relative to the manifest unless absolute.
    ``pairs`` entries take the fields of SymmetryPair. Optional top-level
    keys are ``nodes``, ``bands``, and ``k_tolerance_inverse_angstrom``.
    EIGENVAL occupations are summarized to help assess the chosen prefix.
    """
    path = Path(path)
    manifest = json.loads(path.read_text(encoding="utf-8-sig"))
    lattice, inverse = _invert_lattice(manifest["lattice_angstrom"], "lattice_angstrom")
    reciprocal = 2 * np.pi * inverse.T
    occupied = _integer(manifest["occupied_bands"], "occupied_bands")
    if occupied < 1:
        raise ValueError("occupied_bands must be positive")
    cache, probes, occupation_reports = {}, {}, {}
    for name, probe in manifest["probes"].items():
        source = (path.parent / probe["eigenval"]).resolve()
        if source not in cache:
            cache[source] = read_eigenval(source)
        data = cache[source]
        point_count = len(data.fractional_k)
        if "k_index" not in probe and point_count != 1:
            raise ValueError(f"{name}: k_index is required for a multi-point EIGENVAL")
        k_index = _integer(probe.get("k_index", 1), f"k_index for {name}")
        if not 1 <= k_index <= point_count:
            raise ValueError(f"{name}: k_index must lie in 1..{point_count}")
        i = k_index - 1
        energies = data.energies[i]
        probes[name] = CollinearProbe(data.fractional_k[i] @ reciprocal, energies)
        order = np.argsort(energies, axis=0, kind="stable")
        fillings = np.take_along_axis(data.occupations[i], order, axis=0)
        occupation_reports[name] = {
            "prefix_min": float(np.min(fillings[:occupied])),
            "remaining_max": float(np.max(fillings[occupied:])) if occupied < len(fillings) else None,
        }
    report = audit_collinear_symmetry(
        probes, (SymmetryPair(**pair) for pair in manifest["pairs"]),
        occupied_bands=occupied, nodes=manifest.get("nodes", ()), bands=manifest.get("bands", ()),
        reciprocal=reciprocal,
        k_tolerance_inverse_angstrom=manifest.get("k_tolerance_inverse_angstrom", 1e-7),
    )
    report["lattice_angstrom"] = lattice.tolist()
    report["reciprocal_inverse_angstrom"] = reciprocal.tolist()
    for name, occupation in occupation_reports.items():
        report["probes"][name]["eigenval_occupations"] = occupation
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, help="JSON manifest for lattice, EIGENVAL probes, and Cartesian partners")
    parser.add_argument("--out", type=Path, help="write JSON to a new file instead of stdout")
    args = parser.parse_args(argv)
    try:
        report = json.dumps(audit_manifest(args.manifest), indent=2, ensure_ascii=False, allow_nan=False)
        if args.out is None:
            print(report)
        else:
            with args.out.open("x", encoding="utf-8") as stream:
                stream.write(report + "\n")
    except (OSError, ValueError, TypeError, KeyError, AttributeError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
