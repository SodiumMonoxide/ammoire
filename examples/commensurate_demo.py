"""Export the triangular bilayer demonstration and eight-angle controls.

Run from the repository root with ``python -m examples.commensurate_demo``
or, after installing AMMoire, ``python examples/commensurate_demo.py``.
Only NumPy is required; add ``--plot`` for Matplotlib figures. The model
uses representative parameters and fixed magnetic order, not a NiCl2 fit.
"""

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from ammoire.commensurate import COMMENSURATE_PAIRS, CommensurateTriangularBilayer


def diagnostics(model):
    """Numerical symmetry residuals and an isolated-branch i-wave check."""
    k = (np.random.default_rng(20261107).random((30, 2)) - 0.5) @ model.reciprocal
    h = model.hamiltonian(k)
    bands = model.spin_bands(k)
    radius = 0.4 * np.linalg.norm(model.reciprocal[0])
    angles = np.arange(12) * np.pi / 6
    nodes = radius * np.column_stack([np.cos(angles), np.sin(angles)])
    lobes = radius * np.column_stack([np.cos(angles + np.pi / 12), np.sin(angles + np.pi / 12)])
    lobe_bands = model.spin_bands(lobes)
    radii = radius * np.array([0.23, 0.30, 0.37, 0.44])
    ray = radii[:, None] * [np.cos(np.pi / 12), np.sin(np.pi / 12)]
    phi = 2 * np.pi / 3
    rotation = np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])
    residuals = {
        "hermiticity_eV": float(np.max(np.abs(h - h.conj().swapaxes(-1, -2)))),
        "reciprocal_periodicity_eV": float(max(np.max(np.abs(h - model.hamiltonian(k + g)))
                                                for g in model.reciprocal)),
        "C2_spin_exchange_eV": float(np.max(np.abs(bands - model.spin_bands(k * [1, -1])[:, ::-1]))),
        "C3_eV": float(np.max(np.abs(bands - model.spin_bands(k @ rotation.T)))),
        "inversion_eV": float(np.max(np.abs(bands - model.spin_bands(-k)))),
        "i_wave_nodes_eV": float(np.max(np.abs(model.splitting(nodes)))),
        "zero_twist_splitting_eV": float(np.max(np.abs(CommensurateTriangularBilayer(False).splitting(k)))),
        "zero_exchange_splitting_eV": float(np.max(np.abs(CommensurateTriangularBilayer(exchange=0).splitting(k)))),
    }
    if max(residuals.values()) > 1e-10:
        raise AssertionError(f"Symmetry validation failed: {residuals}")
    lobe_splitting = lobe_bands[:, 0, 0] - lobe_bands[:, 1, 0]
    exponent = float(np.polyfit(np.log(radii), np.log(np.abs(model.splitting(ray)[:, 0])), 1)[0])
    isolation = float(np.min(np.diff(lobe_bands[..., :2], axis=-1)))
    if not np.all(lobe_splitting * np.roll(lobe_splitting, 1) < 0) or abs(exponent - 6) > 0.03 or isolation < 0.1:
        raise AssertionError("Expected an isolated lowest branch with alternating i-wave lobes")
    return {
        "residuals": residuals,
        "lowest_branch_lobe_radius_inverse_a": radius,
        "lowest_branch_lobe_splitting_eV": lobe_splitting.tolist(),
        "lowest_branch_gap_to_next_at_lobes_eV": isolation,
        "lowest_branch_small_k_log_log_slope": exponent,
        "small_k_fit_radii_inverse_a": radii.tolist(),
        "generic_max_rank_paired_spin_splitting_eV": float(np.max(np.abs(model.splitting(k)))),
    }


def angle_scan():
    """Compare identical Cartesian radius separately from mini-zone fraction.

    The matching-Q table deliberately keeps primitive same-Q differences
    separate from folded lowest-rank differences: band folding can reverse
    their signs even when interlayer hopping is zero.
    """
    phi = np.arange(72) * 2 * np.pi / 72
    directions = np.column_stack([np.cos(phi), np.sin(phi)])
    basis6 = np.sin(6 * phi)
    fixed_radius = 0.45
    q = (9.14 / np.sqrt(7)) * np.array([0.2976694799071473, 0.10311573261235644])
    indices = np.array([(u, v) for u in range(-2, 3) for v in range(-2, 3)])
    rows, matched = [], []
    for m, n in COMMENSURATE_PAIRS:
        model = CommensurateTriangularBilayer(m=m, n=n)
        decoupled = CommensurateTriangularBilayer(m=m, n=n, interlayer=0)
        fraction_radius = 0.4 * np.linalg.norm(model.reciprocal[0])
        fixed_split = model.splitting(fixed_radius * directions)[:, 0]
        fixed_decoupled_split = decoupled.splitting(fixed_radius * directions)[:, 0]
        primitive_split = model.primitive_dispersion(fixed_radius * directions, 0) - model.primitive_dispersion(
            fixed_radius * directions, 1)
        np.testing.assert_allclose(fixed_decoupled_split, primitive_split, atol=1e-12)
        fraction_bands = model.spin_bands(fraction_radius * directions)
        fraction_split = fraction_bands[:, 0, 0] - fraction_bands[:, 1, 0]
        rows.append({
            "angle_deg": model.twist_deg, "m": m, "n": n, "sites_per_layer": model.sites_per_layer,
            "fixed_ka": fixed_radius,
            "fixed_full_sin6_meV": float(2000 * np.mean(fixed_split * basis6)),
            "fixed_decoupled_sin6_meV": float(2000 * np.mean(fixed_decoupled_split * basis6)),
            "fixed_leading_sin6_meV": float(1000 * model.hopping * np.sin(3 * np.deg2rad(model.twist_deg)) / 1920 * fixed_radius**6),
            "fraction_ka": float(fraction_radius),
            "fraction_full_max_meV": float(1000 * np.max(np.abs(fraction_split))),
            "fraction_decoupled_max_meV": float(1000 * np.max(np.abs(decoupled.splitting(fraction_radius * directions)[:, 0]))),
            "fraction_lowest_gap_min_eV": float(np.min(np.diff(fraction_bands[..., :2], axis=-1))),
        })
        candidates = q + indices @ model.reciprocal
        nearest_index = int(np.argmin(np.linalg.norm(candidates, axis=1)))
        nearest = candidates[nearest_index]
        np.testing.assert_allclose(model.spin_bands(q), model.spin_bands(nearest), atol=1e-12)
        matched.append({
            "angle_deg": model.twist_deg, "qx_a": float(q[0]), "qy_a": float(q[1]),
            "nearest_kx_a": float(nearest[0]), "nearest_ky_a": float(nearest[1]),
            "added_b1": int(indices[nearest_index, 0]), "added_b2": int(indices[nearest_index, 1]),
            "primitive_unfolded_same_q_delta_meV": float(1000 * (model.primitive_dispersion(q, 0) - model.primitive_dispersion(q, 1))),
            "primitive_nearest_representative_delta_meV": float(1000 * (model.primitive_dispersion(nearest, 0) - model.primitive_dispersion(nearest, 1))),
            "folded_ordered_lowest_tperp0_delta_meV": float(1000 * decoupled.splitting(q)[0]),
            "folded_ordered_lowest_full_delta_meV": float(1000 * model.splitting(q)[0]),
        })
    return rows, matched


def write_rows(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("results/commensurate"))
    parser.add_argument("--points", type=int, default=40, help="samples per path segment")
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args(argv)
    if args.points < 2:
        parser.error("--points must be at least two")
    if args.plot:
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            parser.error("--plot requires Matplotlib; install ammoire[plot]")
    model = CommensurateTriangularBilayer()
    control = CommensurateTriangularBilayer(twisted=False)
    b1, b2 = model.reciprocal
    qprime, q = np.array([[0.30, -0.15], [0.45, 0.15]]) @ model.reciprocal
    vertices = np.array([qprime, [0.0, 0.0], q, (2 * b1 + b2) / 3])
    labels = ["Q prime", "Gamma", "Q", "K"]
    points = np.concatenate([np.linspace(a, b, args.points, endpoint=False)
                             for a, b in zip(vertices[:-1], vertices[1:])] + [vertices[-1:]])
    distances = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]
    bands = {"AA control": control.spin_bands(points), "21.7868 degree twist": model.spin_bands(points)}
    summary = {
        "model": "CommensurateTriangularBilayer",
        "parameter_status": "illustrative geometry model; not fitted to NiCl2 or an exact paper reproduction",
        "twist_deg": model.twist_deg,
        "area_in_primitive_cells": 7,
        "states": model.size,
        "parameters_eV": {"hopping": model.hopping, "interlayer": model.interlayer, "exchange": model.exchange},
        "interlayer_width_a": model.width,
        "interlayer_cutoff_a": 1.75,
        "real_cell_a": model.lattice.tolist(),
        "reciprocal_cell_inverse_a": model.reciprocal.tolist(),
        "spin_difference_convention": "up minus down, same within-spin ascending energy rank; not orbital tracking",
        "Q_Qprime_fractional": [[0.45, 0.15], [0.30, -0.15]],
        "Q_Qprime_splitting_eV": model.splitting(np.array([[0.45, 0.15], [0.30, -0.15]]) @ model.reciprocal).tolist(),
        "checks": diagnostics(model),
        "angle_scan": {
            "validated_pairs": [list(pair) for pair in COMMENSURATE_PAIRS],
            "radius_conventions": "fixed_ka=0.45 is the same Cartesian radius; fraction_ka=0.8*Gamma-M varies with the supercell",
            "folding_control": "tb_matched_q.csv separates same primitive Q from folded lowest energy rank; a sign difference can persist at t_perp=0",
            "magnetic_order": "fixed A-type collinear AFM for every angle, not a prediction of magnetic stability",
        },
    }
    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / "bands.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["case", "k_index", "kx_inverse_a", "ky_inverse_a", "k_distance_inverse_a",
                         "spin", "energy_rank_zero_based", "energy_eV"])
        for label, spectrum in bands.items():
            for i, point in enumerate(points):
                for spin, spin_name in enumerate(("up", "down")):
                    for rank, energy in enumerate(spectrum[i, spin]):
                        writer.writerow([label, i, *point, distances[i], spin_name, rank, energy])
    angle_rows, matched_rows = angle_scan()
    write_rows(args.out / "angle_scan.csv", angle_rows)
    write_rows(args.out / "tb_matched_q.csv", matched_rows)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    if args.plot:
        figure, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True, layout="constrained")
        ticks = distances[::args.points]
        for axis, (label, spectrum) in zip(axes, bands.items()):
            up = axis.plot(distances, spectrum[:, 0], color="#a12643", lw=0.9)
            down = axis.plot(distances, spectrum[:, 1], color="#145a87", ls="--", lw=0.8)
            axis.legend([up[0], down[0]], ["spin up", "spin down"], frameon=False, fontsize=8)
            axis.set_xticks(ticks, labels)
            axis.set(title=label, xlabel="Cartesian k path")
            axis.grid(axis="x", alpha=0.25)
        axes[0].set_ylabel("Illustrative TB energy (eV)")
        figure.savefig(args.out / "bands.png", dpi=180)
        plt.close(figure)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
