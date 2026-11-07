import unittest

import numpy as np
from numpy.testing import assert_allclose

from ammoire import MoireBilayer, make_model
from ammoire.cairo import cell_gauge, evaluate
from ammoire.models import PAULI


def spin_bands(model, k):
    h = model.hamiltonian(k)
    return np.stack([np.linalg.eigvalsh(h[..., s::2, s::2]) for s in (0, 1)], axis=-2)


def splitting(model, k):
    bands = spin_bands(model, k)
    return bands[..., 0, :] - bands[..., 1, :]


class CairoTests(unittest.TestCase):
    def test_formula_and_cell_gauge(self):
        model = make_model('cairo', delta=.1)
        k = np.array([[.31, -.79], [-1.1, .43]])
        raw = evaluate(k, model.params, PAULI[2], periodic=False)
        # Independently check three S7 bonds and signed exchange onsite.
        assert_allclose(raw[..., 0, 2], -.4 * (1.1 * np.cos(k.sum(axis=-1) / 2)
                                             + .9 * np.cos(np.diff(k, axis=-1)[..., 0] / 2)))
        assert_allclose(raw[..., 0, 4], -.9 * np.exp(.5j * k[:, 0]))
        assert_allclose(raw[..., 0, 8], -1.1 * np.exp(.5j * k[:, 1]))
        assert_allclose(raw[..., [0, 1, 2, 3], [0, 1, 2, 3]], [[-1, 1, 1, -1]] * 2)
        u = cell_gauge(k)
        h = model.hamiltonian(k)
        assert_allclose(h, u[..., :, None] * raw * u[..., None, :].conj())
        assert_allclose(np.linalg.eigvalsh(h), np.linalg.eigvalsh(raw), atol=1e-12)
        for g in model.reciprocal:
            assert_allclose(model.hamiltonian(k + g), h, atol=1e-12)

    def test_nodes_and_strain_power(self):
        x = np.linspace(-2.4, 2.4, 37)
        for strain in (0, .1):
            model = make_model('cairo', delta=strain)
            for sign in (-1, 1):
                assert_allclose(splitting(model, np.stack([x, sign * x], axis=-1)), 0, atol=1e-12)
        axis = np.stack([x, np.zeros_like(x)], axis=-1)
        assert_allclose(splitting(make_model('cairo'), axis), 0, atol=1e-12)
        self.assertGreater(np.max(np.abs(splitting(make_model('cairo', delta=.1), axis))), .1)
        r = np.array([.04, .02])
        ray = np.stack([2 * r, r], axis=-1)
        powers = []
        for strain in (0, .1):
            delta = splitting(make_model('cairo', delta=strain), ray)[:, 4]
            powers.append(np.log2(abs(delta[0] / delta[1])))
        assert_allclose(powers, [4, 2], atol=.06)

    def test_exchange_and_ligand_limits(self):
        k = np.random.default_rng(7).uniform(-2, 2, (15, 2))
        for parameter in ('j', 't', 'tn1', 'tn2'):
            assert_allclose(splitting(make_model('cairo', **{parameter: 0}), k), 0, atol=1e-12)
        # With ligand hopping off, the magnetic block has analytic 2-site spectra.
        model = make_model('cairo', t=0)
        raw = evaluate(k, model.params, PAULI[2], periodic=False)
        radius = np.sqrt(1 + (4 * .2 * np.cos(k[:, 0] / 2) * np.cos(k[:, 1] / 2)) ** 2)
        assert_allclose(np.linalg.eigvalsh(raw[..., :4, :4]),
                        np.stack([-radius, -radius, radius, radius], axis=-1), atol=1e-12)

    def test_spin_rotation_and_continuum_composition(self):
        model = make_model('cairo')
        k = [.2, .4]
        assert_allclose(np.linalg.eigvalsh(model.hamiltonian(k)),
                        np.linalg.eigvalsh(model.with_n((1, 2, 3)).hamiltonian(k)), atol=1e-12)
        bilayer = MoireBilayer(model, model.with_n((0, 0, -1)), twist_deg=4,
                              shells=0, harmonics={(0, 0): .03})
        h = bilayer.hamiltonian(k)
        self.assertEqual(h.shape, (24, 24))
        assert_allclose(h, h.conj().T, atol=1e-12)


if __name__ == '__main__':
    unittest.main()
