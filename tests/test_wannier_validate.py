import io
import json
from contextlib import redirect_stdout
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from ammoire.wannier_validate import main, validate


class WannierValidationTests(unittest.TestCase):
    def setUp(self):
        directory = TemporaryDirectory(prefix="wannier-validate-", dir=Path.cwd())
        self.addCleanup(directory.cleanup)
        self.seed = Path(directory.name) / "wannier90.2"
        self._file(".win").write_text(
            "num_wann = 2\nuse_ws_distance = false\n"
            "begin unit_cell_cart\nang\n2 0 0\n0.7 1.8 0\n0 0 4\nend unit_cell_cart\n",
            encoding="utf-8",
        )
        rows = ["synthetic dispersive two-orbital model", "2", "3", "1 1 1"]
        for translation, diagonal in [((-1, 1, 0), (0.17, -0.13)),
                                      ((0, 0, 0), (-0.7, 1.2)),
                                      ((1, -1, 0), (0.17, -0.13))]:
            r = " ".join(map(str, translation))
            rows.extend([f"{r} 1 1 {diagonal[0]} 0", f"{r} 2 1 0 0",
                         f"{r} 1 2 0 0", f"{r} 2 2 {diagonal[1]} 0"])
        self._file("_hr.dat").write_text("\n".join(rows) + "\n", encoding="utf-8")
        self._file("_band.kpt").write_text(
            "4\n0.0 0.0 0.0 1.0\n0.2 0.1 0.0 1.0\n"
            "0.3 -0.1 0.0 1.0\n0.4 -0.2 0.0 1.0\n", encoding="utf-8",
        )

    def _file(self, suffix):
        return self.seed.with_name(self.seed.name + suffix)

    def _write_reference(self, *, energy_shift=0.0, mismatched_path=False):
        # Distinct k points on either side of a path jump share an x coordinate.
        path = [0.0, 0.2, 0.2, 0.6]
        phase = 2 * np.pi * np.array([0.0, 0.1, 0.4, 0.6])
        energies = [-0.7 + 0.34 * np.cos(phase), 1.2 - 0.26 * np.cos(phase)]
        blocks = []
        for index, band in enumerate(energies):
            distances = [0.0, 0.2, 0.3, 0.6] if index and mismatched_path else path
            blocks.append("\n".join(
                f"{x:g} {energy + energy_shift:.8E} 0.25".replace("E", "D")
                for x, energy in zip(distances, band)
            ))
        self._file("_band.dat").write_text("\n\n".join(blocks) + "\n", encoding="utf-8")

    def test_band_major_reference_with_projection_weight(self):
        self._write_reference()
        result_path = self._file(".json")
        stdout = io.StringIO()
        with redirect_stdout(stdout):
            status = main([str(self.seed), "--tolerance", "1e-7", "--out", str(result_path)])
        report = json.loads(result_path.read_text(encoding="utf-8"))
        self.assertEqual(status, 0)
        self.assertEqual(json.loads(stdout.getvalue()), report)
        self.assertEqual(report["comparison"], "same_seed_wannier90_interpolation")
        self.assertEqual((report["num_kpoints"], report["num_bands"]), (4, 2))
        self.assertTrue(report["pass"])
        self.assertLess(report["max_abs_eV"], 1e-8)
        self.assertLess(report["periodicity_residual_eV"], 1e-12)

    def test_constant_energy_shift_fails_without_alignment(self):
        self._write_reference(energy_shift=0.05)
        stdout = io.StringIO()
        with redirect_stdout(stdout):
            status = main([str(self.seed)])
        report = json.loads(stdout.getvalue())
        self.assertEqual(status, 1)
        self.assertFalse(report["pass"])
        self.assertAlmostEqual(report["max_abs_eV"], 0.05, places=7)
        self.assertAlmostEqual(report["rmse_eV"], 0.05, places=7)

    def test_inconsistent_path_is_rejected(self):
        self._write_reference(mismatched_path=True)
        with self.assertRaisesRegex(ValueError, "different path coordinates"):
            validate(self.seed)

    def test_truncated_kpoint_file_is_rejected(self):
        self._file("_band.kpt").write_text("4\n0.0 0.0 0.0\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "declared 4 k-points but read 1"):
            validate(self.seed)


if __name__ == "__main__":
    unittest.main()
