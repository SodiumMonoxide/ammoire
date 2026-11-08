"""Compare imported H(k) with the same seed's Wannier90 band interpolation.

Use the original gnuplot-format band files from standard interpolation. The
comparison uses the exported eV reference directly and measures agreement with
Wannier90's interpolant. Independent DFT bands require a separate comparison.

File layout: https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/files/#seedname_banddat
Wannier90 3.1.0 writer: https://github.com/wannier-developers/wannier90/blob/v3.1.0/src/plot.F90
"""

import argparse
import json
from pathlib import Path

import numpy as np

from ._wannier_io import seed_file
from .wannier import WannierModel


def _numeric_rows(file, path, column_counts, first_line):
    for number, line in enumerate(file, start=first_line):
        fields = line.split()
        if not fields:
            continue
        if len(fields) not in column_counts:
            raise ValueError(f"{path.name}:{number}: unexpected number of columns")
        try:
            values = [float(field.replace("D", "E").replace("d", "e")) for field in fields]
        except ValueError as exc:
            raise ValueError(f"{path.name}:{number}: invalid numeric record") from exc
        if not np.all(np.isfinite(values)):
            raise ValueError(f"{path.name}:{number}: non-finite value")
        yield values


def _read_band_kpoints(path):
    with path.open(encoding="utf-8-sig") as file:
        try:
            count = int(file.readline())
        except ValueError as exc:
            raise ValueError(f"{path.name}: expected a k-point count on the first line") from exc
        if count < 1:
            raise ValueError(f"{path.name}: expected a positive k-point count")
        points = [row[:3] for row in _numeric_rows(file, path, (3, 4), 2)]
    if len(points) != count:
        raise ValueError(f"{path.name}: declared {count} k-points but read {len(points)}")
    return np.asarray(points)


def _read_band_energies(path, num_kpoints, num_bands):
    with path.open(encoding="utf-8-sig") as file:
        rows = [row[:2] for row in _numeric_rows(file, path, (2, 3), 1)]
    if len(rows) != num_kpoints * num_bands:
        raise ValueError(f"{path.name}: expected {num_bands} bands of {num_kpoints} k-points")
    bands = np.asarray(rows).reshape(num_bands, num_kpoints, 2)
    distances = bands[..., 0]
    if not np.all(distances == distances[0]):
        raise ValueError(f"{path.name}: band blocks have different path coordinates")
    return bands[..., 1].T


def validate(seed, *, tolerance=1e-4):
    """Return eV errors for energy ranks and reciprocal-periodic H(k)."""
    seed = Path(seed)
    tolerance = float(tolerance)
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("tolerance must be positive and finite")
    model = WannierModel.from_seed(seed)
    fractional_k = _read_band_kpoints(seed_file(seed, "_band.kpt"))
    reference = _read_band_energies(seed_file(seed, "_band.dat"), len(fractional_k), model.size)
    k = fractional_k @ model.reciprocal
    residual = np.empty_like(reference)
    periodicity = 0.0
    # Limit temporary matrix arrays for long paths in large supercells.
    batch = max(1, 524288 // model.size**2)
    for first in range(0, len(k), batch):
        section = slice(first, first + batch)
        points = k[section]
        h = model.hamiltonian(points)
        residual[section] = np.linalg.eigvalsh(h) - reference[section]
        for g in model.reciprocal:
            periodicity = max(periodicity, float(np.max(np.abs(model.hamiltonian(points + g) - h))))
    max_abs = float(np.max(np.abs(residual)))
    return {
        "seed": seed.name,
        "comparison": "same_seed_wannier90_interpolation",
        "interpolation": model.interpolation,
        "num_kpoints": len(k),
        "num_bands": model.size,
        "tolerance_eV": tolerance,
        "rmse_eV": float(np.sqrt(np.mean(residual**2))),
        "max_abs_eV": max_abs,
        "periodicity_residual_eV": periodicity,
        "pass": max(max_abs, periodicity) <= tolerance,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("seed", type=Path, help="Wannier90 seed path, e.g. run/wannier90.1")
    parser.add_argument("--tolerance", type=float, default=1e-4, metavar="EV",
                        help="maximum allowed eV error (default: %(default)g)")
    parser.add_argument("--out", type=Path, metavar="JSON", help="also save the JSON report")
    args = parser.parse_args(argv)
    try:
        report = validate(args.seed, tolerance=args.tolerance)
        output = json.dumps(report, indent=2, allow_nan=False) + "\n"
        if args.out is not None:
            input_suffixes = (".win", "_hr.dat", "_wsvec.dat", "_band.kpt", "_band.dat")
            if args.out.resolve() in {
                seed_file(args.seed, suffix).resolve() for suffix in input_suffixes
            }:
                raise ValueError("report path is a Wannier90 input to this comparison")
            args.out.write_text(output, encoding="utf-8")
    except (OSError, ValueError, np.linalg.LinAlgError) as exc:
        parser.error(str(exc))
    print(output, end="")
    return 0 if report["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
