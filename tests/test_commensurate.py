import unittest

import numpy as np
from numpy.testing import assert_allclose

from ammoire.commensurate import COMMENSURATE_PAIRS, CommensurateTriangularBilayer


class CommensurateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = CommensurateTriangularBilayer()
        cls.points = (np.random.default_rng(20261107).random((20, 2)) - 0.5) @ cls.model.reciprocal

    def test_cell_and_hermitian_periodic_hamiltonian(self):
        model, k = self.model, self.points.reshape(4, 5, 2)
        primitive_area = np.sqrt(3) / 2
        self.assertAlmostEqual(np.linalg.det(model.lattice) / primitive_area, 7)
        assert_allclose(model.lattice @ model.reciprocal.T, 2 * np.pi * np.eye(2), atol=1e-13)
        h = model.hamiltonian(k)
        self.assertEqual(h.shape, (4, 5, 28, 28))
        assert_allclose(h, h.conj().swapaxes(-1, -2), atol=1e-12)
        for g in model.reciprocal:
            assert_allclose(model.hamiltonian(k + g), h, atol=1e-12)
        spin_z = np.kron(np.eye(14), np.diag([1, -1]))
        assert_allclose(h @ spin_z - spin_z @ h, 0, atol=1e-13)

    def test_spin_exchange_mirror_C3_and_inversion(self):
        model, k = self.model, self.points
        bands = model.spin_bands(k)
        assert_allclose(bands, model.spin_bands(k * [1, -1])[:, ::-1, :], atol=1e-12)
        phi = 2 * np.pi / 3
        rotation = np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])
        assert_allclose(bands, model.spin_bands(k @ rotation.T), atol=1e-12)
        assert_allclose(bands, model.spin_bands(-k), atol=1e-12)
        assert_allclose(np.sort(bands.reshape(len(k), 28)),
                        np.linalg.eigvalsh(model.hamiltonian(k)), atol=1e-12)
        self.assertGreater(np.max(np.abs(model.splitting(k))), 1e-3)

    def test_zero_exchange_and_zero_twist_spin_degeneracy(self):
        for model in (CommensurateTriangularBilayer(exchange=0),
                      CommensurateTriangularBilayer(twisted=False)):
            assert_allclose(model.splitting(self.points), 0, atol=1e-12)

    def test_zero_interlayer_and_atomic_exchange_limits(self):
        decoupled = CommensurateTriangularBilayer(interlayer=0)
        h = decoupled.hamiltonian(self.points)
        assert_allclose(h[..., :14, 14:], 0, atol=1e-13)
        assert_allclose(h[..., 14:, :14], 0, atol=1e-13)
        # Each physical layer's two spins differ by exactly its 2J field.
        for layer, sign in ((slice(0, 14), -1), (slice(14, 28), 1)):
            block = h[..., layer, layer]
            assert_allclose(block[..., ::2, ::2] - block[..., 1::2, 1::2],
                            np.broadcast_to(2 * sign * decoupled.exchange * np.eye(7), (len(h), 7, 7)),
                            atol=1e-13)
        atomic = CommensurateTriangularBilayer(hopping=0, interlayer=0)
        assert_allclose(atomic.hamiltonian([0.2, -0.4]),
                        atomic.exchange * atomic.exchange_operator, atol=1e-13)

    def test_i_wave_nodes_alternation_and_sixth_order_on_isolated_branch(self):
        model = self.model
        radius = 0.4 * np.linalg.norm(model.reciprocal[0])
        angles = np.arange(12) * np.pi / 6
        nodes = radius * np.column_stack([np.cos(angles), np.sin(angles)])
        assert_allclose(model.splitting(nodes), 0, atol=1e-12)
        lobes = radius * np.column_stack([np.cos(angles + np.pi / 12), np.sin(angles + np.pi / 12)])
        bands = model.spin_bands(lobes)
        split = bands[:, 0, 0] - bands[:, 1, 0]
        self.assertTrue(np.all(split * np.roll(split, 1) < 0))
        self.assertGreater(np.min(np.diff(bands[..., :2], axis=-1)), 0.1)
        radii = radius * np.array([0.23, 0.30, 0.37, 0.44])
        ray = radii[:, None] * [np.cos(np.pi / 12), np.sin(np.pi / 12)]
        exponent = np.polyfit(np.log(radii), np.log(np.abs(model.splitting(ray)[:, 0])), 1)[0]
        self.assertLess(abs(exponent - 6), 0.03)

    def test_existing_Q_pair_and_gamma_reduced_subspace(self):
        model = self.model
        points = np.array([[0.45, 0.15], [0.30, -0.15]]) @ model.reciprocal
        split = model.splitting(points)
        assert_allclose(split[0], -split[1], atol=1e-12)
        self.assertAlmostEqual(split[0, 0], -0.0001223, delta=1e-7)
        self.assertAlmostEqual(split[0, 4], -0.04917, delta=1e-5)
        control = CommensurateTriangularBilayer(twisted=False)
        primitive = control.primitive_layers[0]
        vectors = np.array([(u, v) for u in range(-2, 3) for v in range(-2, 3)]) @ primitive
        distances = np.linalg.norm(vectors, axis=1)
        t_gamma = control.interlayer * np.exp(-distances[distances < 1.75]**2 / (2 * control.width**2)).sum()
        # Uniform seven-site layer orbitals form an exact Gamma subspace.
        projector = np.kron(np.kron(np.eye(2), np.full((7, 1), 1 / np.sqrt(7))), np.eye(2))
        expected = np.kron([[6 * control.hopping, t_gamma], [t_gamma, 6 * control.hopping]], np.eye(2))
        expected += control.exchange * np.kron(np.diag([-1, 1]), np.diag([1, -1]))
        assert_allclose(control.hamiltonian([0, 0]) @ projector, projector @ expected, atol=1e-12)

    def test_reject_invalid_parameters_and_momenta(self):
        for parameters in ({"twisted": 4}, {"hopping": np.nan}, {"exchange": np.inf},
                           {"width": 0}, {"width": -1}, {"interlayer": np.nan},
                           {"m": 2.0}, {"m": 4, "n": 2}, {"m": True},
                           {"m": 3, "n": 2, "twisted": False}):
            with self.subTest(parameters=parameters), self.assertRaises(ValueError):
                CommensurateTriangularBilayer(**parameters)
        for momentum in (1, [1, 2, 3], [np.nan, 0], [[0, np.inf]]):
            with self.subTest(momentum=momentum), self.assertRaises(ValueError):
                self.model.hamiltonian(momentum)

    def test_eight_angle_cells_periodicity_and_spin_exchange_symmetry(self):
        expected_sites = (7, 19, 7, 13, 13, 7, 19, 7)
        expected_angles = (0, 13.173551107258918, 21.78678929826181, 27.795772496027972,
                           32.204227503972035, 38.21321070173819, 46.82644889274108, 60)
        for pair, sites, angle in zip(COMMENSURATE_PAIRS, expected_sites, expected_angles):
            with self.subTest(pair=pair):
                model = CommensurateTriangularBilayer(m=pair[0], n=pair[1])
                self.assertEqual(model.sites_per_layer, sites)
                self.assertEqual(model.size, 4 * sites)
                self.assertAlmostEqual(model.twist_deg, angle)
                k = np.array([[0.08, 0.17], [0.43, -0.21], [-0.19, -0.16]]) @ model.reciprocal
                h = model.hamiltonian(k)
                assert_allclose(h, h.conj().swapaxes(-1, -2), atol=1e-12)
                for g in model.reciprocal:
                    assert_allclose(model.hamiltonian(k + g), h, atol=1e-12)
                bands = model.spin_bands(k)
                # The aligned endpoints exchange layers at the same k; the
                # arbitrarily folded seven-site cell need not preserve k_y reflection.
                partner = k if pair in ((2, 2), (1, 0)) else k * [1, -1]
                assert_allclose(bands, model.spin_bands(partner)[:, ::-1], atol=1e-12)
                phi = 2 * np.pi / 3
                rotation = np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])
                assert_allclose(bands, model.spin_bands(k @ rotation.T), atol=1e-12)
                if pair in ((2, 2), (1, 0)):
                    assert_allclose(model.splitting(k), 0, atol=1e-12)

    def test_complementary_angles_rotate_and_swap_spin(self):
        phi = np.pi / 6
        rotation = np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])
        for first, second in ((0, 7), (1, 6), (2, 5), (3, 4)):
            left = CommensurateTriangularBilayer(m=COMMENSURATE_PAIRS[first][0], n=COMMENSURATE_PAIRS[first][1])
            right = CommensurateTriangularBilayer(m=COMMENSURATE_PAIRS[second][0], n=COMMENSURATE_PAIRS[second][1])
            with self.subTest(angles=(left.twist_deg, right.twist_deg)):
                k = (np.random.default_rng(9).random((12, 2)) - 0.5) @ left.reciprocal
                assert_allclose(left.spin_bands(k), right.spin_bands(k @ rotation.T)[:, ::-1], atol=1e-12)

    def test_decoupled_small_k_reference_and_folding_sign_controls(self):
        angles = np.arange(24) * 2 * np.pi / 24
        small_k = 0.45 * np.column_stack([np.cos(angles), np.sin(angles)])
        q = (9.14 / np.sqrt(7)) * np.array([0.2976694799071473, 0.10311573261235644])
        indices = np.array([(u, v) for u in range(-2, 3) for v in range(-2, 3)])
        for pair in COMMENSURATE_PAIRS:
            with self.subTest(pair=pair):
                model = CommensurateTriangularBilayer(m=pair[0], n=pair[1], interlayer=0)
                primitive = model.primitive_dispersion(small_k, 0) - model.primitive_dispersion(small_k, 1)
                assert_allclose(model.splitting(small_k)[:, 0], primitive, atol=1e-12)
                candidates = q + indices @ model.reciprocal
                nearest = candidates[np.argmin(np.linalg.norm(candidates, axis=1))]
                assert_allclose(model.spin_bands(q), model.spin_bands(nearest), atol=1e-12)
                direct_delta = model.primitive_dispersion(q, 0) - model.primitive_dispersion(q, 1)
                nearest_delta = model.primitive_dispersion(nearest, 0) - model.primitive_dispersion(nearest, 1)
                assert_allclose(model.splitting(q)[0], nearest_delta, atol=1e-12)
                if pair in ((3, 2), (5, 2), (3, 1), (7, 1)):
                    self.assertLess(direct_delta, 0)
                    self.assertGreater(model.splitting(q)[0], 0)
                if pair == (5, 2):
                    self.assertAlmostEqual(1000 * direct_delta, -0.13478314670734015, delta=1e-8)
                    self.assertAlmostEqual(1000 * model.splitting(q)[0], 0.07014948547201527, delta=1e-8)


if __name__ == "__main__":
    unittest.main()
