import unittest

import numpy as np
from numpy.testing import assert_allclose

from ammoire.compensated import CompensatedTriangularBilayer


class CompensatedTriangularBilayerTests(unittest.TestCase):
    def test_periodic_spin_blocks_and_rotating_spectrum(self):
        for pair, sites in (((2, 1), 7), ((5, 2), 13)):
            with self.subTest(pair=pair):
                model = CompensatedTriangularBilayer(m=pair[0], n=pair[1])
                self.assertEqual(model.sites_per_sublayer, sites)
                self.assertEqual(model.size, 8 * sites)
                k = np.array([[0.21, 0.13], [-0.32, 0.29]]) @ model.reciprocal
                h = model.hamiltonian(k)
                assert_allclose(h, h.swapaxes(-1, -2).conj(), atol=1e-12)
                for g in model.reciprocal:
                    assert_allclose(h, model.hamiltonian(k + g), atol=1e-12)
                assert_allclose(h[:, ::2, 1::2], 0, atol=1e-12)
                bands = model.spin_bands(k)
                assert_allclose(bands, model.spin_bands(k * [1, -1])[:, ::-1], atol=1e-12)
                assert_allclose(bands, model.spin_bands(-k), atol=1e-12)
                assert_allclose(np.sort(bands.reshape(len(k), model.size)), np.linalg.eigvalsh(h), atol=1e-12)

    def test_isolated_parents_and_exchange_off_preserve_degeneracy(self):
        points = np.array([[0.13, 0.27], [0.29, -0.33]])
        for options in ({"interlayer": 0, "facing_potential": 0}, {"exchange": 0}, {"m": 2, "n": 2}, {"m": 1, "n": 0}):
            with self.subTest(options=options):
                model = CompensatedTriangularBilayer(**options)
                assert_allclose(model.splitting(points), 0, atol=1e-12)

    def test_interface_ablations_reproduce_the_reference_ring(self):
        model = CompensatedTriangularBilayer()
        phi = np.arange(120) * 2 * np.pi / 120
        radius = 0.4 * np.linalg.norm(model.reciprocal[0])
        ring = radius * np.column_stack([np.cos(phi), np.sin(phi)])
        expected_meV = (
            ({}, 0.10473144),
            ({"facing_potential": 0}, 0.01133083),
            ({"interlayer": 0}, 0.13543365),
            ({"interlayer": 0, "registry_scale": 0}, 0.13543540),
            ({"interlayer": 0, "uniform_scale": 0}, 0.00187638),
        )
        for options, expected in expected_meV:
            with self.subTest(options=options):
                model = CompensatedTriangularBilayer(**options)
                bands = model.spin_bands(ring)
                splitting = bands[:, 0, 0] - bands[:, 1, 0]
                assert_allclose(splitting[::10], 0, atol=1e-12)
                self.assertAlmostEqual(1000 * np.max(np.abs(splitting)), expected, delta=5e-8)
                if options.get("uniform_scale") == 0:
                    self.assertLess(np.min(np.diff(bands[..., :2], axis=-1)), 2e-6)

    def test_interface_potential_decomposition_is_additive_in_hamiltonian(self):
        k = [0.31, -0.38]
        isolated = CompensatedTriangularBilayer(interlayer=0, facing_potential=0).hamiltonian(k)
        potential = CompensatedTriangularBilayer(interlayer=0).hamiltonian(k) - isolated
        uniform = CompensatedTriangularBilayer(interlayer=0, registry_scale=0).hamiltonian(k) - isolated
        registry = CompensatedTriangularBilayer(interlayer=0, uniform_scale=0).hamiltonian(k) - isolated
        assert_allclose(potential, uniform + registry, atol=1e-13)
        assert_allclose(potential, np.diag(np.diag(potential)), atol=1e-13)
        per_plane = np.diag(registry).reshape(4, 14)
        assert_allclose(per_plane.mean(axis=1), 0, atol=1e-13)
        assert_allclose(per_plane[[0, 3]], 0, atol=1e-13)


if __name__ == "__main__":
    unittest.main()
