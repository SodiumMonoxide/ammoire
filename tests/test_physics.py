import unittest

import numpy as np
from numpy.testing import assert_allclose

from ammoire import MODEL_REGISTRY, MoireBilayer, SquareHAM, make_model
from ammoire.__main__ import calculate, cutoff_difference, mixed_bilayer, sample_path
from ammoire.observables import commutator_norm, resolve_bands


class PhysicsTests(unittest.TestCase):
    def test_registry_periodicity_and_exchange_limit(self):
        for name in MODEL_REGISTRY:
            with self.subTest(name=name):
                model = make_model(name)
                k = np.array([0.39, -0.53, 0.74])[:model.dimension]
                h = model.hamiltonian([k, -k])
                self.assertEqual(h.shape, (2, model.size, model.size))
                assert_allclose(h, h.conj().swapaxes(-1, -2), atol=1e-12)
                assert_allclose(model.hamiltonian(k + model.reciprocal), h[:1].repeat(model.dimension, axis=0),
                                atol=1e-12)
                zero = np.linalg.eigvalsh(make_model(name, j=0).hamiltonian(k))
                assert_allclose(zero[::2], zero[1::2], atol=1e-12)

    def test_ham_matrix_spectrum_and_odd_sector_exchange(self):
        model = SquareHAM()
        k = np.array([[0, 0], [0.4, -0.8], [1.1, 0.3]])
        h = model.hamiltonian(k)
        assert_allclose(h, model.hamiltonian(-k).conj(), atol=1e-13)
        assert_allclose(commutator_norm(h, model.charge), 0, atol=1e-13)
        bands = model.sector_bands(k)
        assert_allclose(bands, model.sector_bands(-k)[:, ::-1], atol=1e-13)
        assert_allclose(np.linalg.eigvalsh(h), np.sort(bands.reshape(-1, 4)), atol=1e-13)

    def test_degenerate_sector_resolution(self):
        model = SquareHAM()
        h = model.hamiltonian([0.0, 0.0])
        energies, charges, vectors = resolve_bands(h, model.charge)
        assert_allclose(h @ vectors, vectors * energies, atol=1e-13)
        assert_allclose(charges, [-1, 1, -1, 1], atol=1e-13)

    def test_complex_rectangular_harmonic_convention(self):
        top = SquareHAM()
        bottom = make_model("square-dxy")
        t = np.arange(8).reshape(4, 2) * (0.02 + 0.03j)
        bilayer = MoireBilayer(top, bottom, twist_deg=5, shells=1,
                              harmonics={(1, 0): t}, bias=0.2)
        h = bilayer.hamiltonian([0.0, 0.0])
        first = 4 * np.flatnonzero(np.all(bilayer.plane_waves == [1, 0], axis=1)).item()
        bottom_start = 4 * len(bilayer.plane_waves)
        assert_allclose(h[first:first+4, bottom_start:bottom_start+2], t)
        assert_allclose(h, h.conj().T)
        assert_allclose(h[:4, :4], top.hamiltonian([0, 0]) + 0.1 * np.eye(4))

    def test_moire_conservation_and_projected_weight(self):
        for mixing in (0.0, 0.7):
            bilayer, charge = mixed_bilayer(4, 1, mixing, 0.08, 0.05)
            points, _, _, _ = sample_path(bilayer, 3)
            results = calculate(bilayer, charge, points)
            _, q, w_top, w_central, residual = results
            assert_allclose(w_top.sum(axis=1), 4 * len(bilayer.plane_waves), atol=1e-11)
            assert_allclose(w_central.sum(axis=1), 6, atol=1e-11)
            self.assertEqual(cutoff_difference(results, results, 0.04), 0)
            if mixing == 0:
                assert_allclose(residual, 0, atol=1e-13)
                assert_allclose(abs(q), 1, atol=1e-11)
            else:
                self.assertGreater(residual.min(), 0)
                self.assertLess(abs(q).min(), 0.99)

    def test_intralayer_potential_adds_hermitian_partner(self):
        model = make_model("square-dxy")
        exchange = (0.06 + 0.02j) * model.spin_operators[0]
        modulated = MoireBilayer(model, twist_deg=4, shells=1, potential_top={(1, 0): exchange})
        bare = MoireBilayer(model, twist_deg=4, shells=1)
        dh = modulated.hamiltonian([0.02, 0.05]) - bare.hamiltonian([0.02, 0.05])
        pos = 2 * np.flatnonzero(np.all(modulated.plane_waves == [1, 0], axis=1)).item()
        assert_allclose(dh[pos:pos+2, :2], exchange)
        assert_allclose(dh[:2, pos:pos+2], exchange.conj().T)
        assert_allclose(dh, dh.conj().T)
        with self.assertRaises(ValueError):
            MoireBilayer(model, potential_top={(1, 0): exchange, (-1, 0): exchange})


if __name__ == "__main__":
    unittest.main()
