import unittest

import numpy as np

from ammoire import MoireBilayer, SquareHAM, make_model
from ammoire.fe_spectra import audit_domain_spectra, audit_group_spectra
from ammoire.fe_symmetry import SpinLaueGroup


POINTS = np.array([[0.36, -0.71], [0.51, 0.28], [-0.62, 0.32]])
MIRROR = np.diag([-1, 1])


class DomainSpectraTests(unittest.TestCase):
    def test_preserving_and_reversing_spin_along_a_common_tilted_axis(self):
        axis = [1, 1, 0]
        plus = make_model('square-dxy', neel=axis, j=0.4)
        preserved = audit_domain_spectra(plus, plus, POINTS, -np.eye(2), spin_axis=axis)
        self.assertTrue(preserved['mapped_same_spin_matches'])
        self.assertEqual(preserved['same_k_relation'], 'preserved')

        minus = make_model('square-dxy', neel=axis, j=-0.4)
        reversed_spin = audit_domain_spectra(plus, minus, POINTS, MIRROR, spin_axis=axis)
        self.assertEqual(reversed_spin['same_k_relation'], 'reversed')
        self.assertTrue(reversed_spin['mapped_same_spin_matches'])
        self.assertGreater(reversed_spin['same_k_preserved_residual'], 0.01)
        self.assertLess(reversed_spin['same_k_splitting_reversed_residual'], 1e-12)

    def test_no_splitting_on_a_nodal_line(self):
        model = make_model('square-dxy')
        report = audit_domain_spectra(model, model, [[0, 0.5], [0, 0.8]], MIRROR)
        self.assertTrue(report['mapped_same_spin_matches'])
        self.assertEqual(report['same_k_relation'], 'unresolved')
        self.assertLess(report['max_splitting_positive'], 1e-12)

    def test_active_rotation_distinguishes_the_two_momentum_maps(self):
        plus = make_model('gondolf-dxy')
        theta = 2 * np.pi / 3
        rotation = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])

        class RotatedDomain:
            spin_operators = plus.spin_operators
            dimension, size = plus.dimension, plus.size

            def hamiltonian(self, k):
                return plus.hamiltonian(k @ rotation)

        minus = RotatedDomain()
        report = audit_domain_spectra(plus, minus, POINTS, rotation)
        self.assertTrue(report['mapped_same_spin_matches'])
        self.assertEqual(report['same_k_relation'], 'neither')
        self.assertEqual(report['bands_per_spin'], 2)
        self.assertFalse(audit_domain_spectra(plus, minus, POINTS, rotation.T)['mapped_same_spin_matches'])

    def test_common_energy_reference(self):
        plus = make_model('square-dx2-y2', mu=0)
        shifted = make_model('square-dx2-y2', mu=1.7)
        report = audit_domain_spectra(plus, shifted, POINTS, np.eye(2))
        self.assertAlmostEqual(report['mapped_same_spin_residual'], 1.7)
        self.assertAlmostEqual(report['same_k_splitting_preserved_residual'], 0)
        self.assertEqual(report['same_k_relation'], 'neither')

    def test_reject_unphysical_spin_or_spin_mixing(self):
        model = make_model('square-dxy')
        with self.assertRaisesRegex(ValueError, 'spin_operators'):
            audit_domain_spectra(SquareHAM(), SquareHAM(), POINTS, MIRROR)
        with self.assertRaisesRegex(ValueError, 'mixes physical spin'):
            audit_domain_spectra(model, model, POINTS, MIRROR, spin_axis=[0, 1, 0])
        soc = make_model('lieb-dx2-y2', soc=0.2)
        with self.assertRaisesRegex(ValueError, 'no-SOC'):
            audit_domain_spectra(soc, soc, [[0, 0]], MIRROR)

    def test_band_counts_points_operation_and_tolerance(self):
        model = make_model('square-dxy')
        with self.assertRaisesRegex(ValueError, 'band counts'):
            audit_domain_spectra(model, make_model('brekke-lieb'), POINTS, MIRROR)
        for points, operation, tolerance in [
            (POINTS, [[2, 0], [0, 1]], 1e-8), (POINTS, MIRROR, np.nan),
            (POINTS[:0], MIRROR, 0), (POINTS * np.nan, MIRROR, 0)
        ]:
            with self.subTest(points=points, operation=operation, tolerance=tolerance):
                with self.assertRaises(ValueError):
                    audit_domain_spectra(model, model, points, operation, atol=tolerance)

    def test_bilayer_and_group_operations(self):
        layer = make_model('square-dxy')
        bilayer = MoireBilayer(layer, harmonics={(0, 0): 0.05}, twist_deg=0, shells=0)
        group = SpinLaueGroup(2, [(MIRROR, -1)])
        report = audit_group_spectra(bilayer, POINTS, group)
        self.assertTrue(report['matches'])
        self.assertEqual(report['tested_operations'], 4)
        incompatible = make_model('square-dx2-y2')
        self.assertFalse(audit_group_spectra(incompatible, POINTS, group)['matches'])

        mixed = MoireBilayer(SquareHAM(), layer, twist_deg=0, shells=0,
                             harmonics={(0, 0): np.ones((4, 2))})
        with self.assertRaisesRegex(ValueError, 'spin_operators'):
            audit_group_spectra(mixed, POINTS, group)


if __name__ == '__main__':
    unittest.main()
