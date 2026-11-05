"""Compare collinear VASP spin splittings with a spin-conserving bulk TB model.

This reads the five-column band records of a static ISPIN=2 EIGENVAL. Supply
orbital-matched VASP band indices and use the same reciprocal cell and axis
conventions in the calculation and model. EIGENVAL alone does not establish
orbital character, magnetic order, or convergence.
"""

import argparse
import csv
from dataclasses import dataclass
import json
from operator import index
from pathlib import Path

import numpy as np

from .models import MODEL_REGISTRY, make_model


@dataclass(frozen=True)
class Eigenvalues:
    fractional_k: np.ndarray
    weights: np.ndarray
    energies: np.ndarray
    occupations: np.ndarray


def _read_floats(fields, count, context):
    if len(fields) != count:
        raise ValueError(f"expected {count} columns for {context}")
    try:
        values = [float(s.replace("D", "E").replace("d", "e")) for s in fields]
    except ValueError as exc:
        raise ValueError(f"invalid numeric value for {context}") from exc
    if not np.all(np.isfinite(values)):
        raise ValueError(f"non-finite value for {context}")
    return values


def read_eigenval(path):
    """Read fractional k and arrays shaped (nk, nband, 2), spin up then down."""
    with Path(path).open(encoding="utf-8-sig") as file:
        first = file.readline().split()
        if len(first) != 4 or first[-1] != "2":
            raise ValueError("expected collinear ISPIN=2 in the EIGENVAL header")
        for _ in range(4):
            if not file.readline():
                raise ValueError("truncated EIGENVAL header")
        try:
            _, nk, nb = file.readline().split()
            nk, nb = int(nk), int(nb)
        except ValueError as exc:
            raise ValueError("invalid EIGENVAL k-point and band counts") from exc
        if min(nk, nb) < 1:
            raise ValueError("EIGENVAL must contain k-points and bands")
        lines = (line.split() for line in file if line.strip())
        points, bands = [], []
        for ik in range(1, nk + 1):
            points.append(_read_floats(next(lines, []), 4, f"k-point {ik}"))
            block = []
            for ib in range(1, nb + 1):
                values = _read_floats(next(lines, []), 5, f"k-point {ik}, band {ib}")
                if values[0] != ib:
                    raise ValueError(f"unexpected band index {values[0]:g} at k-point {ik}")
                block.append(values[1:])
            bands.append(block)
        if next(lines, None) is not None:
            raise ValueError("extra records after declared EIGENVAL k-point blocks")
    points, bands = np.asarray(points), np.asarray(bands)
    return Eigenvalues(points[:, :3], points[:, 3], bands[..., :2], bands[..., 2:])


def compare_eigenvalues(data, model, *, band_up, band_down=None, branch="lower",
                        node_tolerance=1e-5, energy_scale=1.0):
    """Return CSV columns and unweighted errors for one pair of VASP bands.

    The model must have a three-dimensional primitive reciprocal basis and
    spin as its innermost index. The same ascending within-spin TB branch
    is chosen for both spins. Model parameters and domain are supplied by
    the caller. energy_scale expresses eV per model energy unit.
    """
    if model.dimension != 3 or model.size % 2:
        raise ValueError("comparison requires a 3D bulk model with spin doublets")
    if any(value != 0 for key, value in model.params.items() if key.startswith("soc")):
        raise ValueError("comparison requires a TB model without SOC")
    if branch not in {"lower", "upper"}:
        raise ValueError("branch must be lower or upper")
    if not np.isfinite(node_tolerance) or node_tolerance < 0:
        raise ValueError("node tolerance must be finite and nonnegative")
    if not np.isfinite(energy_scale) or energy_scale <= 0:
        raise ValueError("energy scale must be positive and finite")
    band_up = index(band_up)
    band_down = band_up if band_down is None else index(band_down)
    if not (1 <= band_up <= data.energies.shape[1] and 1 <= band_down <= data.energies.shape[1]):
        raise ValueError("one-based VASP band indices must lie in EIGENVAL's band range")
    k = data.fractional_k @ model.reciprocal
    h = model.hamiltonian(k)
    if np.max(np.abs(h[..., ::2, 1::2])) > 1e-10:
        raise ValueError("TB Hamiltonian mixes up and down spin channels")
    which = 0 if branch == "lower" else -1
    tb_up = energy_scale * np.linalg.eigvalsh(h[..., ::2, ::2])[..., which]
    tb_down = energy_scale * np.linalg.eigvalsh(h[..., 1::2, 1::2])[..., which]
    vasp_up = data.energies[:, band_up - 1, 0]
    vasp_down = data.energies[:, band_down - 1, 1]
    dft_delta = vasp_up - vasp_down
    tb_delta = tb_up - tb_down
    residual = dft_delta - tb_delta
    near_zero = np.abs(tb_delta) <= node_tolerance
    columns = {
        "k_index": np.arange(1, len(k) + 1),
        **{f"k_frac_{axis}": column for axis, column in zip("123", data.fractional_k.T)},
        **{f"k_cart_{axis}": column for axis, column in zip("xyz", k.T)},
        "k_weight": data.weights,
        "vasp_up_eV": vasp_up,
        "vasp_down_eV": vasp_down,
        "vasp_occ_up": data.occupations[:, band_up - 1, 0],
        "vasp_occ_down": data.occupations[:, band_down - 1, 1],
        "tb_up_eV": tb_up,
        "tb_down_eV": tb_down,
        "vasp_up_minus_down_eV": dft_delta,
        "tb_up_minus_down_eV": tb_delta,
        "vasp_minus_tb_splitting_eV": residual,
        "tb_near_zero": near_zero,
    }
    summary = {
        "tb_model": model.name,
        "tb_parameters": dict(model.params),
        "tb_neel": model.neel,
        "tb_branch": branch,
        "tb_eV_per_model_unit": energy_scale,
        "vasp_band_up": band_up,
        "vasp_band_down": band_down,
        "k_points": len(k),
        "node_tolerance_eV": node_tolerance,
        "tb_near_zero_count": int(np.count_nonzero(near_zero)),
        "max_abs_vasp_splitting_at_tb_near_zero_eV": (
            float(np.max(np.abs(dft_delta[near_zero]))) if near_zero.any() else None
        ),
        "unweighted_splitting_rmse_eV": float(np.sqrt(np.mean(residual**2))),
        "max_abs_splitting_residual_eV": float(np.max(np.abs(residual))),
    }
    return columns, summary


def _parameter(value):
    try:
        key, number = value.split("=", 1)
        return key, float(number)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("parameter overrides use KEY=VALUE") from exc


def _self_test():
    from tempfile import TemporaryDirectory

    fixture = (
        " 6 6 1 2\n header\n header\n CAR\n test\n 20 2 2\n"
        "\n 0 0 0 0.5\n 1 -1.50 -1.50 1 1\n 2 0.90 0.90 0 0\n"
        "\n 0.21 0.17 0.3 0.5\n 1 -1.20 -1.25 1 1\n 2 0.80 0.90 0 0\n"
    )
    with TemporaryDirectory(dir=Path.cwd()) as folder:
        path = Path(folder) / "EIGENVAL"
        path.write_text(fixture, encoding="utf-8")
        data = read_eigenval(path)
        model = make_model("roig-mnf2", a=4.8, c=3.3)
        columns, summary = compare_eigenvalues(data, model, band_up=1, band_down=1)
        np.testing.assert_allclose(columns["vasp_up_minus_down_eV"], [0, 0.05])
        np.testing.assert_allclose(columns["k_cart_z"], [0, 0.3 * 2 * np.pi / 3.3])
        assert columns["tb_near_zero"].tolist() == [True, False]
        assert summary["max_abs_vasp_splitting_at_tb_near_zero_eV"] == 0
        shifted, _ = compare_eigenvalues(data, model, band_up=1, band_down=2)
        np.testing.assert_allclose(shifted["vasp_up_minus_down_eV"], [-2.4, -2.1])
        flipped, _ = compare_eigenvalues(data, model.with_n((0, 0, -1)), band_up=1)
        np.testing.assert_allclose(flipped["tb_up_minus_down_eV"], -columns["tb_up_minus_down_eV"])
        path.write_text(fixture.replace(" 6 6 1 2", " 6 6 1 1"), encoding="utf-8")
        with np.testing.assert_raises(ValueError):
            read_eigenval(path)
        path.write_text(fixture.replace(" 2 0.80 0.90 0 0\n", ""), encoding="utf-8")
        with np.testing.assert_raises(ValueError):
            read_eigenval(path)
    print("VASP bridge synthetic tests passed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("eigenval", type=Path, nargs="?")
    parser.add_argument("--band-up", type=int, help="one-based, selected using orbital character")
    parser.add_argument("--band-down", type=int, help="default: same index as --band-up")
    parser.add_argument("--branch", choices=("lower", "upper"), default="lower")
    parser.add_argument("--model", choices=tuple(
        name for name, spec in MODEL_REGISTRY.items() if spec.dimension == 3
    ), default="roig-mnf2")
    parser.add_argument("--param", type=_parameter, action="append", default=[])
    parser.add_argument("--neel-sign", type=int, choices=(1, -1), default=1)
    parser.add_argument("--tb-ev-per-unit", type=float, help="required for dimensionless presets")
    parser.add_argument("--node-tol", type=float, default=1e-5, help="TB near-zero threshold (eV)")
    parser.add_argument("--out", type=Path, help="optional per-k CSV path, must not already exist")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        _self_test()
        return
    if args.eigenval is None or args.band_up is None:
        parser.error("provide EIGENVAL and --band-up, or use --self-test")
    scale = args.tb_ev_per_unit
    if scale is None:
        if args.model not in {"roig-mnf2", "roig-ruo2", "roig-fesb2"}:
            parser.error("provide --tb-ev-per-unit for this preset")
        scale = 1.0
    try:
        model = make_model(args.model, neel=(0, 0, args.neel_sign), **dict(args.param))
        columns, summary = compare_eigenvalues(
            read_eigenval(args.eigenval), model, band_up=args.band_up,
            band_down=args.band_down, branch=args.branch, node_tolerance=args.node_tol,
            energy_scale=scale,
        )
        if args.out is not None:
            with args.out.open("x", encoding="utf-8", newline="") as file:
                writer = csv.writer(file)
                writer.writerow(columns)
                writer.writerows(zip(*columns.values()))
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
