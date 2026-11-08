import tempfile
import unittest
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose

from ammoire.wannier import WannierModel


class WannierSeedTests(unittest.TestCase):
    def setUp(self):
        scratch = tempfile.TemporaryDirectory(dir=Path.cwd())
        self.addCleanup(scratch.cleanup)
        self.seed = Path(scratch.name) / "case.1"
        self.win = self.seed.with_name("case.1.win")
        self.hr = self.seed.with_name("case.1_hr.dat")
        self.ws = self.seed.with_name("case.1_wsvec.dat")
        self.cell_text = "Begin Unit_Cell_Cart ! comment\nBohr\n4d0 0 0\n1, 3, 0\n0 0 12\nEnd Unit_Cell_Cart\n"
        self.win.write_text(self.cell_text, encoding="utf-8")
        self.raw_hr = ["synthetic complex fixture", "2", "3", "2 1 2"]
        self.raw_ws = ["## fixture with use_ws_distance=.true."]
        forward = np.array([[0.5, 4 + 2j], [0, -1.]])
        onsite = np.diag([-2., 1.])
        for rx, matrix in ((-1, forward.conj().T), (0, onsite), (1, forward)):
            for i, j in np.ndindex(2, 2):
                value = matrix[i, j]
                self.raw_hr.append(f"{rx} 0 0 {i+1} {j+1} {value.real} {value.imag}")
                shifts = [0]
                if (rx, i, j) in ((1, 0, 1), (-1, 1, 0)):
                    shifts.append(2 * rx)
                self.raw_ws.extend([f"{rx} 0 0 {i+1} {j+1}", str(len(shifts))])
                self.raw_ws.extend(f"{shift} 0 0" for shift in shifts)
        self.hr.write_text("\n".join(self.raw_hr) + "\n", encoding="utf-8")

    def write_sidecar(self, rows=None):
        self.ws.write_text("\n".join(self.raw_ws if rows is None else rows) + "\n", encoding="utf-8")

    def test_orbital_pair_shifts_change_off_grid_interpolation(self):
        self.write_sidecar()
        model = WannierModel.from_seed(self.seed)
        self.assertEqual(model.interpolation, "wsvec")
        assert_allclose(model.lattice, .529177210903 * np.array([[4, 0, 0], [1, 3, 0], [0, 0, 12]]))
        q = np.array([.2, -.11, .17])
        k = q @ model.reciprocal
        phase = np.exp(2j * np.pi * q[0])
        # One off-diagonal element receives both independent degeneracy
        # weights. The diagonal terms still use the unshifted R pair.
        off = (4 + 2j) * (phase + phase**3) / 4
        expected = [[-2 + .5 * phase.real, off], [off.conjugate(), 1 - phase.real]]
        assert_allclose(model.hamiltonian(k), expected, atol=1e-13)
        assert_allclose(model.hamiltonian(k + model.reciprocal), np.tile(expected, (3, 1, 1)), atol=1e-13)
        assert_allclose(WannierModel.from_hr(self.hr, model.lattice, dimension=3).hamiltonian(k), expected)

    def test_seed_requires_matching_sidecar_or_explicit_disabled_mode(self):
        with self.assertRaisesRegex(ValueError, "requires.*_wsvec.dat"):
            WannierModel.from_seed(self.seed)
        self.win.write_text(self.cell_text + "use_ws_distance : F # deliberately disabled\n", encoding="utf-8")
        self.assertEqual(WannierModel.from_seed(self.seed).interpolation, "degeneracy")
        self.write_sidecar()
        with self.assertRaisesRegex(ValueError, "differs between"):
            WannierModel.from_seed(self.seed)

    def test_rejects_bad_ws_entry_and_out_of_plane_effective_shift(self):
        self.write_sidecar(self.raw_ws[:-1])
        with self.assertRaisesRegex(ValueError, "truncated Wigner"):
            WannierModel.from_seed(self.seed)
        self.write_sidecar(self.raw_ws[:1] + self.raw_ws[1:4] + self.raw_ws[1:4] + self.raw_ws[4:])
        with self.assertRaisesRegex(ValueError, "duplicate Wigner"):
            WannierModel.from_seed(self.seed)
        self.write_sidecar(["0 0 1" if row == "2 0 0" else "0 0 -1" if row == "-2 0 0" else row
                            for row in self.raw_ws])
        with self.assertRaisesRegex(ValueError, "effective Rz=0"):
            WannierModel.from_hr(self.hr, np.eye(2), dimension=2)

    def test_expanded_export_does_not_apply_sidecar_twice(self):
        self.win.write_text(self.cell_text + "write_ndegen_applied=.TRUE.\n", encoding="utf-8")
        self.write_sidecar(["## with use_ws_distance=.true. and write_ndegen_applied=.true."])
        with self.assertRaisesRegex(ValueError, "unit degeneracies"):
            WannierModel.from_seed(self.seed)
        records = self.raw_hr.copy()
        records[3] = "1 1 1"
        self.hr.write_text("\n".join(records) + "\n", encoding="utf-8")
        model = WannierModel.from_seed(self.seed)
        self.assertEqual(model.interpolation, "expanded")
        self.assertEqual(model.hamiltonian([0, 0, 0])[0, 1], 4 + 2j)


if __name__ == "__main__":
    unittest.main()
