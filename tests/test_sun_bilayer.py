import unittest

import numpy as np
from numpy.testing import assert_allclose

from ammoire import make_model


class SunBilayerTests(unittest.TestCase):
    def test_real_space_neighbours(self):
        # Independent geometric enumeration in fractional direct coordinates.
        k = np.array([2.13, -1.07])
        cells = np.array([(i, j) for i in range(-2, 3) for j in range(-2, 3)])
        for stacking in range(6):
            model = make_model('sun-hex-bilayer', stacking=stacking)
            angle = stacking * np.pi / 3
            shift = np.array([np.cos(angle), np.sin(angle)]) / 3
            sites = np.array([[0, 0], [1/3, -2/3]]) @ model.lattice
            translations = cells @ model.lattice
            inter = np.zeros((2, 2), dtype=complex)
            for i, bottom in enumerate(sites):
                for j, top in enumerate(sites):
                    distances = np.linalg.norm(translations + top + shift - bottom, axis=1)
                    neighbours = translations[np.isclose(distances, 1/3)]
                    self.assertEqual(len(neighbours), 1)
                    inter[i, j] = 0.5 * np.exp(1j * (neighbours @ k)).sum()
            h = model.hamiltonian(k)
            assert_allclose(h[:4, 4:], np.kron(inter, np.eye(2)), atol=1e-13)
            assert_allclose(h.diagonal()[::2], [-1, 1, 1, -1])

    def test_fixed_frame_domain_spectra(self):
        reference = make_model('sun-hex-bilayer')
        points = np.array([[.38, 1.70], [2.1, -1.08], [1.3, .49], [4.1, 1.54]])
        mirror = np.diag([1, -1])
        c3 = np.array([[-0.5, -np.sqrt(3)/2], [np.sqrt(3)/2, -0.5]])

        def spectra(model, k):
            h = model.hamiltonian(k)
            return np.stack([np.linalg.eigvalsh(h[..., s::2, s::2]) for s in (0, 1)], axis=1)

        base = spectra(reference, points)
        assert_allclose(spectra(reference, -points), base, atol=1e-12)
        assert_allclose(spectra(reference, points @ mirror.T), base[:, ::-1], atol=1e-12)
        self.assertGreater(np.max(np.abs(base[:, 0] - base[:, 1])), 0.1)
        for stacking, sign, action in [(3, -1, np.eye(2)), (3, 1, mirror), (5, -1, c3)]:
            with self.subTest(stacking=stacking, sign=sign):
                switched = make_model('sun-hex-bilayer', stacking=stacking, layer_sign=sign)
                assert_allclose(spectra(switched, points @ action.T), base, atol=1e-12)
                same_k = spectra(switched, points)
                if stacking == 5:
                    self.assertGreater(np.max(abs(same_k - base)), 0.1)
                    self.assertGreater(np.max(abs(same_k - base[:, ::-1])), 0.1)
                else:
                    expected = base if sign == -1 else base[:, ::-1]
                    assert_allclose(same_k, expected, atol=1e-12)

    def test_disconnected_layers_and_stacking_parameters(self):
        model = make_model('sun-hex-bilayer', t_inter=0)
        h = model.hamiltonian([1.6, -1.2])
        assert_allclose(h[:4, 4:], 0)
        assert_allclose(np.linalg.eigvalsh(h[::2, ::2]), np.linalg.eigvalsh(h[1::2, 1::2]))
        for params in ({'stacking': 1.5}, {'stacking': 6}, {'layer_sign': 0}):
            with self.assertRaises(ValueError):
                make_model('sun-hex-bilayer', **params).hamiltonian([0, 0])
