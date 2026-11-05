"""Run AM/HAM moiré benchmark and export band data."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from .hyperspin import SquareHAM, SZ, charge_channels
from .models import MODEL_REGISTRY, make_model
from .moire import MoireBilayer
from .observables import commutator_norm, expectation, layer_operator, resolve_bands


def mixed_bilayer(angle, shells, mixing, coupling, bias):
    """Square HAM on top, d-wave AM below, common illustrative energy unit."""
    # Place the upper HAM branch near the AM band edge at Γ.
    top = SquareHAM(offset=1.17)
    bottom = make_model("square-dx2-y2")
    q_top, q_bottom = top.charge, SZ
    seed = np.vstack((np.eye(2), 0.65 * np.eye(2)))
    same, opposite = charge_channels(seed, q_top, q_bottom)
    tunneling = coupling * (same + mixing * opposite)
    harmonics = {(0, 0): tunneling}
    if angle != 0:
        harmonics.update({g: 0.4 * tunneling for g in ((1, 0), (-1, 0), (0, 1), (0, -1))})
    bilayer = MoireBilayer(top, bottom, twist_deg=angle, shells=shells,
                          bias=bias, harmonics=harmonics)
    return bilayer, layer_operator(bilayer, q_top, q_bottom)


def sample_path(bilayer, points_per_segment):
    if bilayer.twist_deg == 0:
        g1, g2 = bilayer.top.reciprocal
        labels = ("Γ", "X", "M", "Γ")
        vertices = np.array([np.zeros(2), g1 / 2, (g1 + g2) / 2, np.zeros(2)])
    else:
        labels, vertices = bilayer.high_symmetry_path()
    points = np.concatenate([
        np.linspace(start, end, points_per_segment, endpoint=False)
        for start, end in zip(vertices[:-1], vertices[1:])
    ] + [vertices[-1:]])
    distance = np.r_[0, np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]
    return points, distance, labels, distance[::points_per_segment]


def calculate(bilayer, charge, points):
    n = len(bilayer.plane_waves)
    top_weight = layer_operator(bilayer, np.eye(bilayer.top.size),
                                np.zeros((bilayer.bottom.size, bilayer.bottom.size)))
    central = np.flatnonzero(np.all(bilayer.plane_waves == 0, axis=1)).item()
    central_indices = np.r_[
        central * bilayer.top.size + np.arange(bilayer.top.size),
        n * bilayer.top.size + central * bilayer.bottom.size + np.arange(bilayer.bottom.size),
    ]
    energies, charges, top_weights, central_weights, residuals = [], [], [], [], []
    for k in points:
        h = bilayer.hamiltonian(k)
        e, q, v = resolve_bands(h, charge)
        energies.append(e)
        charges.append(q)
        top_weights.append(expectation(v, top_weight))
        central_weights.append(np.sum(abs(v[central_indices])**2, axis=0))
        residuals.append(commutator_norm(h, charge))
    return tuple(map(np.array, (energies, charges, top_weights, central_weights, residuals)))


def cutoff_difference(coarse, fine, broadening):
    """Compare G=0 projected Gaussian spectra on one common energy grid."""
    ec, _, _, _, _ = coarse
    ef, _, _, _, _ = fine
    lo, hi = min(ec.min(), ef.min()), max(ec.max(), ef.max())
    grid = np.linspace(lo - 6 * broadening, hi + 6 * broadening,
                       max(400, int((hi - lo) / broadening * 8)))

    def spectral_density(results):
        energies, _, _, weights, _ = results
        x = (grid[None, None, :] - energies[..., None]) / broadening
        return np.sum(weights[..., None] * np.exp(-0.5 * x**2), axis=1)

    reference = spectral_density(fine)
    discrepancy = np.abs(spectral_density(coarse) - reference).sum(axis=-1)
    return float(np.max(discrepancy / reference.sum(axis=-1)))


def plot_results(path, distance, labels, ticks, cases):
    import matplotlib.pyplot as plt

    fig, panels = plt.subplots(1, len(cases), figsize=(6 * len(cases), 4.5),
                               sharey=True, layout="constrained", squeeze=False)
    for axis, (name, results) in zip(panels[0], cases.items()):
        e, q, _, w, comm = results
        x = np.broadcast_to(distance[:, None], e.shape)
        selected = w > 0.01
        plot = axis.scatter(x[selected], e[selected], c=q[selected],
                            s=2 + 12 * w[selected], cmap="coolwarm", vmin=-1, vmax=1)
        axis.set_xticks(ticks, labels)
        axis.grid(axis="x", alpha=0.25)
        axis.set_title(f"{name}\nmax ||[H,Q]|| = {comm.max():.3g}")
    panels[0, 0].set_ylabel("Energy (model units)")
    fig.colorbar(plot, ax=panels.ravel().tolist(), label=r"$\langle Q\rangle$")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("models", help="list native models")
    demo = commands.add_parser("demo", help="compare same- and mixed-charge moiré tunneling")
    demo.add_argument("--angle", type=float, default=4.0, help="twist angle in degrees")
    demo.add_argument("--shells", type=int, default=2)
    demo.add_argument("--mixing", type=float, default=0.7, help="opposite-Q channel scale")
    demo.add_argument("--coupling", type=float, default=0.08, help="tunneling amplitude")
    demo.add_argument("--bias", type=float, default=0.05, help="top-minus-bottom energy shift")
    demo.add_argument("--points", type=int, default=24, help="samples per path segment")
    demo.add_argument("--out", type=Path, default=Path("results"))
    demo.add_argument("--plot", action="store_true")
    demo.add_argument("--check-cutoff", action="store_true",
                      help="compare G=0 spectra at path vertices with one additional shell")
    demo.add_argument("--broadening", type=float, default=0.04,
                      help="Gaussian width in model units for --check-cutoff")
    args = parser.parse_args()
    if args.command == "models":
        print(*MODEL_REGISTRY, "square-ham", sep="\n")
        return
    if args.points < 1 or not np.all(np.isfinite([args.mixing, args.coupling, args.bias])):
        parser.error("points must be positive and model parameters finite")
    if args.check_cutoff and (args.angle == 0 or not 0.001 <= args.broadening <= 10):
        parser.error("cutoff comparison requires nonzero twist and broadening in [0.001, 10]")
    if args.plot:
        try:
            import matplotlib
        except ImportError:
            parser.error("--plot needs matplotlib (install ammoire[plot])")
    cases = {}
    summary = {"angle_deg": args.angle, "shells": args.shells, "bias": args.bias,
               "coupling": args.coupling, "opposite_charge_scale": args.mixing}
    for label, scale in (("Same-Q tunneling", 0.0), ("With opposite-Q tunneling", args.mixing)):
        bilayer, charge = mixed_bilayer(args.angle, args.shells, scale, args.coupling, args.bias)
        points, distance, labels, ticks = sample_path(bilayer, args.points)
        cases[label] = calculate(bilayer, charge, points)
        _, q, _, w, residuals = cases[label]
        summary[label] = {
            "max_commutator_frobenius": float(residuals.max()),
            "mean_central_weighted_Q_variance": float(np.sum(w * np.maximum(0, 1 - q**2)) / w.sum()),
        }
        if args.check_cutoff:
            finer, finer_charge = mixed_bilayer(
                args.angle, args.shells + 1, scale, args.coupling, args.bias
            )
            vertices = points[:-1:args.points]
            coarse_at_vertices = tuple(x[:-1:args.points] for x in cases[label])
            summary[label]["max_relative_G0_spectral_L1_change"] = cutoff_difference(
                coarse_at_vertices, calculate(finer, finer_charge, vertices), args.broadening
            )
    if args.check_cutoff:
        summary["cutoff_comparison"] = {
            "finer_shells": args.shells + 1, "gaussian_sigma": args.broadening,
            "vertices": labels[:-1],
        }
    summary.update(plane_waves=len(bilayer.plane_waves), bands=bilayer.size,
                   k_points=len(points), ham_offset=bilayer.top.offset)
    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / "bands.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(("case", "k_index", "kx", "ky", "k_distance", "band", "energy",
                         "Q_expectation", "top_weight", "G0_weight", "commutator_norm"))
        for label, (e, q, t, w, residuals) in cases.items():
            for i, (kx, ky) in enumerate(points):
                for j in range(e.shape[1]):
                    writer.writerow((label, i, kx, ky, distance[i], j, e[i, j], q[i, j],
                                     t[i, j], w[i, j], residuals[i]))
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    if args.plot:
        plot_results(args.out / "bands.png", distance, labels, ticks, cases)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
