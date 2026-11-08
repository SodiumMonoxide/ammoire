import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from ammoire.collinear_symmetry import (
    CollinearProbe, SymmetryPair, audit_collinear_symmetry, audit_manifest, main,
)


class CollinearSymmetryTests(unittest.TestCase):
    def setUp(self):
        self.lattice = np.array([[4.5, 0, 0], [1.2, 5.0, 0], [0.4, -0.7, 21]])
        self.reciprocal = 2 * np.pi * np.linalg.inv(self.lattice).T
        self.k = np.array([0.18, 0.22, 0.08])
        self.reflection = np.diag([-1, 1, -1])
        self.shift = (1, -1, 0)
        self.reflected = self.k @ self.reflection.T + np.array(self.shift) @ self.reciprocal
        self.energy = np.array([[-3, -2.95], [-1.23, -1.2], [1.2, 1.5]])
        reflected_energy = self.energy[:, ::-1].copy()
        reflected_energy[0, 0] += 0.0003
        reflected_energy[-1, 1] += 0.08
        inverse_energy = self.energy.copy()
        inverse_energy[-1, 0] += 0.03
        self.probes = {
            "Q": CollinearProbe(self.k, self.energy),
            "partner": CollinearProbe(self.reflected, reflected_energy),
            "minus_Q": CollinearProbe(-self.k, inverse_energy),
            "Gamma": CollinearProbe([0, 0, 0], [[-3, -3], [-1, -1.0001], [1.1, 1.2]]),
        }
        self.pairs = [
            SymmetryPair("Q", "partner", self.reflection, True, self.shift),
            SymmetryPair("Q", "minus_Q", -np.eye(3), False),
            SymmetryPair("Gamma", "Gamma", self.reflection, True),
        ]

    def test_audit_separates_occupied_and_all_ranks_for_pairs_and_nodes(self):
        result = audit_collinear_symmetry(self.probes, self.pairs, occupied_bands=2, bands=[2],
                                          nodes=["Gamma"], reciprocal=self.reciprocal)
        exchange, inversion, node = result["pairs"]
        self.assertEqual(exchange["spin_relation"], "exchange")
        self.assertEqual(exchange["max_all"]["source_spin"], "up")
        self.assertEqual(exchange["max_all"]["band_rank_1based"], 3)
        self.assertAlmostEqual(exchange["max_all"]["signed_eV"], -0.08)
        self.assertAlmostEqual(exchange["max_occupied"]["max_abs_eV"], 0.0003)
        self.assertEqual(inversion["spin_relation"], "same_spin")
        self.assertAlmostEqual(inversion["max_all"]["max_abs_eV"], 0.03)
        self.assertEqual(inversion["max_occupied"]["max_abs_eV"], 0)
        self.assertAlmostEqual(node["max_occupied"]["max_abs_eV"], 0.0001)
        self.assertAlmostEqual(result["nodes"]["Gamma"]["max_all"]["max_abs_eV"], 0.1)
        self.assertAlmostEqual(result["probes"]["Q"]["selected_rank_splitting_eV"]["2"], -0.03)

    def test_sorts_ranks_independently_without_inferring_orbitals(self):
        shuffled = self.energy.copy()
        shuffled[:, 1] = shuffled[[2, 0, 1], 1]
        probes = {
            "first": CollinearProbe(self.k, self.energy),
            "reordered": CollinearProbe(self.k, shuffled),
        }
        pair = SymmetryPair("first", "reordered", np.eye(3), False)
        result = audit_collinear_symmetry(probes, [pair], occupied_bands=2)
        self.assertEqual(result["pairs"][0]["max_all"]["max_abs_eV"], 0)

    def test_rejects_mismatched_coordinates_or_band_counts(self):
        for pair in (
            SymmetryPair("Q", "partner", self.reflection, True),
            SymmetryPair("Q", "partner", np.diag([1, 1, -1]), True, self.shift),
        ):
            with self.subTest(pair=pair), self.assertRaisesRegex(ValueError, "Cartesian momentum mismatch"):
                audit_collinear_symmetry(self.probes, [pair], occupied_bands=2, reciprocal=self.reciprocal)
        with self.assertRaisesRegex(ValueError, "basis for a nonzero shift"):
            audit_collinear_symmetry(self.probes, self.pairs, occupied_bands=2)
        bad_pair = SymmetryPair("Q", "partner", 2 * self.reflection, True, self.shift)
        with self.assertRaisesRegex(ValueError, "must be orthogonal"):
            audit_collinear_symmetry(self.probes, [bad_pair], occupied_bands=2, reciprocal=self.reciprocal)
        self.probes["partner"] = CollinearProbe(self.reflected, self.energy[:2])
        with self.assertRaisesRegex(ValueError, "requires equal band counts"):
            audit_collinear_symmetry(self.probes, self.pairs, occupied_bands=2, reciprocal=self.reciprocal)

    def test_static_eigenval_manifest_from_nonorthogonal_cell(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as folder:
            root = Path(folder)
            eigenval = root / "EIGENVAL"
            names = ("Q", "partner", "minus_Q")
            records = ["2 2 1 2", "header", "header", "CAR", "test", "4 3 3"]
            inverse = np.linalg.inv(self.reciprocal)
            for name in names:
                probe = self.probes[name]
                q = probe.k_cartesian @ inverse
                records += ["", " ".join([*(f"{x:.15f}" for x in q), "1.0"])]
                for i, (up, down) in enumerate(probe.energies, 1):
                    occupation = int(i <= 2)
                    records.append(f"{i} {up:.10f} {down:.10f} {occupation} {occupation}")
            eigenval.write_text("\n".join(records) + "\n", encoding="utf-8")
            manifest = {
                "lattice_angstrom": self.lattice.tolist(),
                "occupied_bands": 2,
                "probes": {name: {"eigenval": "EIGENVAL", "k_index": i} for i, name in enumerate(names, 1)},
                "pairs": [
                    {"source": pair.source, "target": pair.target, "swap_spin": pair.swap_spin,
                     "operation_cartesian": pair.operation_cartesian.tolist(), "reciprocal_shift": list(pair.reciprocal_shift)}
                    for pair in self.pairs[:2]
                ],
            }
            path = root / "manifest.json"
            path.write_text(json.dumps(manifest), encoding="utf-8")
            report_file = root / "report.json"
            main([str(path), "--out", str(report_file)])
            report = json.loads(report_file.read_text(encoding="utf-8"))
            self.assertAlmostEqual(report["pairs"][0]["max_occupied"]["max_abs_eV"], 0.0003)
            self.assertLess(report["pairs"][0]["k_error_inverse_angstrom"], 1e-13)
            self.assertEqual(report["probes"]["Q"]["eigenval_occupations"], {"prefix_min": 1.0, "remaining_max": 0.0})
            del manifest["probes"]["partner"]["k_index"]
            path.write_text(json.dumps(manifest), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "k_index is required"):
                audit_manifest(path)


if __name__ == "__main__":
    unittest.main()
