import unittest

import numpy as np
from numpy.testing import assert_allclose

from ammoire.observables import resolve_bands


class ResolveBandsTests(unittest.TestCase):
    def test_rejects_asymmetry_hidden_by_default_allclose(self):
        malformed = np.array([[0., 1.], [1.000009, 0.]])
        with self.assertRaisesRegex(ValueError, "H must be Hermitian"):
            resolve_bands(malformed, np.diag([1., -1.]))
        with self.assertRaisesRegex(ValueError, "Q must be Hermitian"):
            resolve_bands(np.diag([0., 1.]), malformed)
        # Large entries must not restore NumPy's default relative tolerance.
        with self.assertRaisesRegex(ValueError, "H must be Hermitian"):
            resolve_bands(1e10 * malformed, np.diag([1., -1.]))

    def test_rejects_nonfinite_matrices_and_empty_basis(self):
        for value in (np.nan, np.inf, -np.inf):
            malformed = np.diag([value, 1.])
            with self.subTest(value=value, matrix="H"):
                with self.assertRaisesRegex(ValueError, "H must be finite"):
                    resolve_bands(malformed, np.eye(2))
            with self.subTest(value=value, matrix="Q"):
                with self.assertRaisesRegex(ValueError, "Q must be finite"):
                    resolve_bands(np.eye(2), malformed)
        with self.assertRaises(ValueError):
            resolve_bands(np.empty((0, 0)), np.empty((0, 0)))

    def test_hermitian_complex_inputs_at_multiple_energy_scales(self):
        h0 = np.array([[0.2, 0.3 + 0.4j], [0.3 - 0.4j, -0.7]])
        q = np.array([[1., 0.2j], [-0.2j, -1.]])
        for scale in (1., 1e12):
            with self.subTest(scale=scale):
                h = scale * h0
                e, charges, vectors = resolve_bands(h, q)
                residual = np.linalg.norm(h @ vectors - vectors * e)
                self.assertLess(residual / np.linalg.norm(h), 1e-14)
                assert_allclose(vectors.conj().T @ vectors, np.eye(2), atol=1e-14)
                assert_allclose(charges, np.diag(vectors.conj().T @ q @ vectors).real,
                                atol=1e-14)

    def test_rounding_noise_is_symmetrized_without_mutating_inputs(self):
        h = np.array([[0., 1.], [1. + 2e-13, 2.]])
        q = np.array([[1., 0.3j], [-0.3j + 2e-13, -1.]])
        saved_h, saved_q = h.copy(), q.copy()
        e, charges, vectors = resolve_bands(h, q)
        symmetric_h = (h + h.conj().T) / 2
        symmetric_q = (q + q.conj().T) / 2
        assert_allclose(symmetric_h @ vectors, vectors * e, atol=1e-14)
        assert_allclose(charges, np.diag(vectors.conj().T @ symmetric_q @ vectors).real,
                        atol=1e-14)
        np.testing.assert_array_equal(h, saved_h)
        np.testing.assert_array_equal(q, saved_q)


if __name__ == "__main__":
    unittest.main()
