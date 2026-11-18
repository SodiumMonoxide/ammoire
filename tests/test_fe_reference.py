import json
import unittest

import numpy as np

from ammoire.fe_reference import literature_cases, transform_conductivity


class DomainReferenceTests(unittest.TestCase):
    def test_all_published_domain_labels_and_response_types(self):
        cases = json.loads(json.dumps(literature_cases(), allow_nan=False))
        self.assertEqual(len(cases), 8)
        for order, responses in (("L1", ["I", "III", "III"]), ("L2", ["II", "III", "III"])):
            subcases = [c for c in cases if c["material"] == "MnPS3" and c["order"] == order]
            self.assertEqual([c["path"] for c in subcases], [1, 2, 3])
            self.assertEqual([c["response_type"] for c in subcases], responses)
            self.assertEqual([c["initial_state"] for c in subcases], ["P1+"] * 3)
            self.assertEqual([c["final_state"] for c in subcases], ["P1-", "P2-", "P3-"])
        bulk = [c for c in cases if c["material"] == "BiFeO3"]
        self.assertEqual([c["response_type"] for c in bulk], ["I", "II"])
        inversion, twofold = [np.array(c["momentum_action"]) for c in bulk]
        np.testing.assert_array_equal(inversion, -np.eye(3))
        np.testing.assert_array_equal(twofold @ [-1, 1, 0], [-1, 1, 0])
        np.testing.assert_array_equal(twofold @ [1, 1, 1], [-1, -1, -1])
        self.assertAlmostEqual(np.linalg.det(twofold), 1)

    def test_sense_and_factor_order_of_bilayer_path_rotations(self):
        ex = [1, 0, 0]
        ez = [0, 0, 1]
        for case in (c for c in literature_cases() if c["material"] == "MnPS3"):
            angle = 2 * np.pi * (case["path"] - 1) / 3
            if case["order"] == "L2":
                angle = -angle  # My acts after the rotation.
            operation = np.array(case["operation_cartesian"])
            np.testing.assert_allclose(operation @ ex, [np.cos(angle), np.sin(angle), 0], atol=1e-14)
            np.testing.assert_allclose(operation @ ez, [0, 0, -1])
            np.testing.assert_allclose(case["momentum_action"], operation[:2, :2])

    def test_conductivity_rotation_matches_rotated_current_and_field(self):
        tensors = np.array([[[1, 1j], [-2j, 4]], [[-3, -2j], [5j, 1]]])
        field = [1.3, -0.2]
        for case in (c for c in literature_cases() if c["material"] == "MnPS3"):
            operation = np.array(case["momentum_action"])
            result = transform_conductivity(tensors, operation)
            for source, target in zip(tensors, result):
                np.testing.assert_allclose(target @ (operation @ field), operation @ (source @ field))
        bulk = np.arange(18).reshape(2, 3, 3)
        np.testing.assert_array_equal(transform_conductivity(bulk, -np.eye(3)), bulk)

    def test_type_ii_exchange_follows_from_mirror_related_input_tensors(self):
        # For these channels My is a symmetry combined with spin exchange.
        up = np.array([[2, 0.3], [0.3, 1]])
        down = np.array([[2, -0.3], [-0.3, 1]])
        np.testing.assert_allclose(transform_conductivity([up, down], np.diag([1, -1])), [down, up])

    def test_invalid_conductivity_and_operations(self):
        for tensors, operation in ((np.zeros((2, 3)), np.eye(3)), (np.ones((2, 2, 2)), np.eye(3)),
                                   (np.ones((2, 2, 2)), 2 * np.eye(2)), (np.zeros((2, 2, 2)), 1j * np.eye(2)),
                                   (np.full((2, 2, 2), np.nan), np.eye(2))):
            with self.subTest(tensors=tensors, operation=operation), self.assertRaises(ValueError):
                transform_conductivity(tensors, operation)


if __name__ == "__main__":
    unittest.main()
