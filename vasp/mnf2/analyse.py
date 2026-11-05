"""Analyze a completed run_wsl.sh directory.

From the ammoire directory:
    python vasp/mnf2/analyse.py vasp/mnf2/results-smoke
"""

import argparse
import csv
import json
from pathlib import Path
import re
import sys

import numpy as np

# Allow invocation directly from the source tree without installing the package.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from ammoire.models import make_model
from ammoire.vasp_bridge import compare_eigenvalues, read_eigenval


def stage_status(folder):
    outcar = (folder / "OUTCAR").read_text(encoding="utf-8", errors="replace")
    if (folder / "exit_code.txt").read_text().strip() != "0":
        raise ValueError(f"{folder.name} VASP exited with an error")
    if "General timing and accounting informations" not in outcar:
        raise ValueError(f"{folder.name} OUTCAR lacks the normal completion footer")
    if "aborting loop because EDIFF is reached" not in outcar:
        raise ValueError(f"{folder.name} OUTCAR lacks the EDIFF convergence marker")
    ediff = re.findall(r"EDIFF\s*=\s*(\S+)\s+stopping-criterion", outcar)
    energy = re.findall(r"free\s+energy\s+TOTEN\s*=\s*(\S+)\s+eV", outcar)
    if not ediff or not energy:
        raise ValueError(f"missing EDIFF or TOTEN in {folder.name} OUTCAR")
    return {"ediff_reached": True, "ediff_eV": float(ediff[-1]),
            "final_free_energy_eV": float(energy[-1])}, outcar


def magnetic_moments(outcar, oszicar):
    _, separator, table = outcar.rpartition("magnetization (x)")
    if not separator:
        raise ValueError("missing collinear site moments in SCF OUTCAR")
    block = re.split(r"^\s*tot\s", table, maxsplit=1, flags=re.MULTILINE)[0]
    moments = [float(fields[-1]) for line in block.splitlines()
               if len(fields := line.split()) == 5 and fields[0].isdigit()]
    if len(moments) != 6:
        raise ValueError("expected six final site moments in MnF2 SCF OUTCAR")
    mag = re.findall(r"\bmag=\s*(\S+)", oszicar)
    if not mag:
        raise ValueError("missing net moment in SCF OSZICAR")
    return {"Mn_MuB": moments[:2], "F_MuB": moments[2:],
            "net_cell_MuB_oszicar": float(mag[-1])}


def read_probe_points(path):
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    if len(lines) < 4 or lines[2].strip().lower() != "reciprocal":
        raise ValueError("probe KPOINTS must use explicit reciprocal coordinates")
    rows = [line.partition("!") for line in lines[3:] if line.strip()]
    points = np.array([[float(x) for x in data.split()[:3]] for data, _, _ in rows])
    if points.shape != (int(lines[1]), 3):
        raise ValueError("probe KPOINTS count mismatch")
    return points, [label.strip() for _, _, label in rows]


def reciprocal_equivalent(points, target):
    delta = points - target
    return np.all(np.abs(delta - np.rint(delta)) <= 1e-7, axis=-1)


def occupied_spectra(data):
    occupied = data.occupations > 0.5
    counts = occupied.sum(axis=1)
    nocc = int(counts[0, 0])
    if (np.max(abs(data.occupations - occupied)) > 1e-3 or
            not np.all(counts == nocc) or not 0 < nocc < data.energies.shape[1]):
        raise ValueError("occupied-spectrum benchmark needs saturated occupations and equal spin counts")
    e_occ = np.sort(np.where(occupied, data.energies, np.inf), axis=1)[:, :nocc]
    cbm = np.min(data.energies[~occupied])
    vbm = np.max(e_occ)
    return e_occ, float(cbm - vbm)


def symmetry_pairs(points, spectra, labels):
    operations = {
        "mirror_x_spin_swap": ([-1, 1, 1], True),
        "mirror_y_spin_swap": ([1, -1, 1], True),
        "inversion_same_spin": ([-1, -1, -1], False),
    }
    rows = []
    for name, (sign, swap) in operations.items():
        for i, k in enumerate(points):
            for j in np.flatnonzero(reciprocal_equivalent(points, k * sign)):
                if i > j:
                    continue
                difference = spectra[i] - (spectra[j, :, ::-1] if swap else spectra[j])
                rows.append((name, i + 1, j + 1, labels[i], labels[j],
                             float(np.max(np.abs(difference)))))
    return rows


def read_mn_d_weights(path, band, data):
    """Collect Mn d projections of one fixed VASP band from both PROCAR blocks."""
    result = np.full((len(data.fractional_k), 2, 2), np.nan)
    spin, k, selected, d_columns = 0, -1, False, []
    with path.open(encoding="utf-8") as file:
        for line in file:
            k_match = re.match(r"\s*k-point\s+(\d+)\s*:", line)
            band_match = re.match(r"\s*band\s+(\d+)\s*# energy\s+(\S+)", line)
            if k_match:
                next_k = int(k_match[1]) - 1
                if next_k == 0 and k >= 0:
                    spin += 1
                if not 0 <= next_k < len(result) or spin >= 2:
                    raise ValueError("PROCAR point count or spin block count disagrees with EIGENVAL")
                k, selected = next_k, False
            elif band_match:
                selected = int(band_match[1]) == band
                if selected and not np.isclose(
                    float(band_match[2]), data.energies[k, band - 1, spin], atol=2e-6, rtol=0
                ):
                    raise ValueError("PROCAR and EIGENVAL disagree on the representative band")
            elif selected:
                fields = line.split()
                if fields and fields[0] == "ion":
                    d_columns = [i for i, value in enumerate(fields)
                                 if value.startswith("d") or value == "x2-y2"]
                    if len(d_columns) != 5:
                        raise ValueError("expected lm-decomposed PROCAR with five d columns")
                if fields and fields[0] in {"1", "2"} and d_columns:
                    result[k, spin, int(fields[0]) - 1] = sum(float(fields[i]) for i in d_columns)
    if not np.all(np.isfinite(result)):
        raise ValueError("representative Mn d projection was not found for all k and spins")
    return result


def write_csv(path, header, rows):
    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(header)
        writer.writerows(rows)


def analyse(run, output, band=None, branch="lower", neel_sign=1):
    scf, outcar = stage_status(run / "scf")
    probe, _ = stage_status(run / "probe")
    scf.update(magnetic_moments(outcar, (run / "scf" / "OSZICAR").read_text()))
    points, labels = read_probe_points(run / "probe" / "KPOINTS")
    data = read_eigenval(run / "probe" / "EIGENVAL")
    if data.fractional_k.shape != points.shape or not np.all(
        np.abs((data.fractional_k - points) - np.rint(data.fractional_k - points)) <= 1e-7
    ):
        raise ValueError("EIGENVAL k-points do not match the executed KPOINTS order")
    spectra, gap = occupied_spectra(data)
    nocc = spectra.shape[1]
    pairs = symmetry_pairs(points, spectra, labels)
    splitting = spectra[..., 0] - spectra[..., 1]
    nodal = np.any(np.isclose(2 * points[:, :2], np.rint(2 * points[:, :2]), atol=1e-8), axis=1)
    max_pos = np.unravel_index(np.argmax(np.abs(splitting)), splitting.shape)
    band = nocc if band is None else band
    if not 1 <= band <= data.energies.shape[1]:
        raise ValueError("--band is outside the EIGENVAL range")
    if np.any(data.occupations[:, band - 1] <= 0.5):
        raise ValueError("--band must be occupied at all probe points in both spin channels")

    poscar = (run / "probe" / "POSCAR").read_text().splitlines()
    if len(poscar) < 8 or poscar[5].split() != ["Mn", "F"] or poscar[6].split() != ["2", "4"]:
        raise ValueError("expected Mn F site order with two Mn and four F atoms")
    lattice = float(poscar[1]) * np.array([[float(x) for x in line.split()] for line in poscar[2:5]])
    a, b, c = lattice.diagonal()
    if not np.allclose(lattice, np.diag([a, b, c])) or not np.isclose(a, b) or min(a, b, c) <= 0:
        raise ValueError("expected the MnF2 tetragonal conventional-axis primitive cell")
    model = make_model("roig-mnf2", a=a, c=c, neel=(0, 0, neel_sign))
    columns, comparison = compare_eigenvalues(data, model, band_up=band, branch=branch)
    mn_d = read_mn_d_weights(run / "probe" / "PROCAR", band, data)
    columns = {"label": labels, **columns}
    for s, spin in enumerate(("up", "down")):
        for atom in (0, 1):
            columns[f"Mn{atom + 1}_d_projection_{spin}"] = mn_d[:, s, atom]
    comparison.update({
        "description": "Fixed VASP index and one TB branch for illustration. Mn d projections alone "
                       "do not identify a unique orbital-matched branch across k.",
        "min_max_Mn_d_projection_sum": [float(mn_d.sum(axis=2).min()),
                                       float(mn_d.sum(axis=2).max())],
    })
    symmetry = {}
    for name, _, _, _, _, residual in pairs:
        stats = symmetry.setdefault(name, {"n_pairs": 0, "max_abs_residual_eV": 0.0})
        stats["n_pairs"] += 1
        stats["max_abs_residual_eV"] = max(stats["max_abs_residual_eV"], residual)
    scf["gamma_mesh"] = [
        int(x) for x in (run / "scf" / "KPOINTS").read_text().splitlines()[3].split()
    ]
    summary = {
        "scf": scf,
        "fixed_density_probe": probe,
        "n_kpoints": len(points),
        "occupied_bands_per_spin": nocc,
        "sampled_occupied_unoccupied_gap_eV": gap,
        "nodal_k_count": int(nodal.sum()),
        "max_abs_occupied_splitting_at_nodal_k_eV": (
            float(np.max(abs(splitting[nodal]))) if nodal.any() else None
        ),
        "largest_occupied_splitting": {
            "eV_signed": float(splitting[max_pos]), "k_index": max_pos[0].item() + 1,
            "label": labels[max_pos[0]], "occupied_energy_rank": max_pos[1].item() + 1,
        },
        "occupied_spectrum_symmetry": symmetry,
        "fixed_index_tb_comparison": comparison,
    }
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output / "occupied_splitting.csv",
              ("k_index", "label", "k_frac_1", "k_frac_2", "k_frac_3", "nodal_k",
               "occupied_energy_rank", "up_eV", "down_eV", "up_minus_down_eV"),
              ((i + 1, labels[i], *k, nodal[i], j + 1, *spectra[i, j], splitting[i, j])
               for i, k in enumerate(points) for j in range(nocc)))
    write_csv(output / "symmetry_pairs.csv",
              ("operation", "k1_index", "k2_index", "k1_label", "k2_label", "max_residual_eV"), pairs)
    write_csv(output / "fixed_index_tb_comparison.csv", columns, zip(*columns.values()))
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path, help="completed run_wsl.sh output directory")
    parser.add_argument("--out", type=Path, help="defaults to RUN/analysis")
    parser.add_argument("--band", type=int, help="illustrative VASP band index, defaults to highest occupied")
    parser.add_argument("--branch", choices=("lower", "upper"), default="lower")
    parser.add_argument("--neel-sign", type=int, choices=(1, -1), default=1)
    args = parser.parse_args()
    try:
        summary = analyse(args.run, args.out or args.run / "analysis", args.band, args.branch, args.neel_sign)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
