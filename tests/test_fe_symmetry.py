import unittest

import numpy as np

from ammoire.fe_symmetry import SpinLaueGroup


def rotation_z(degrees):
    angle = np.deg2rad(degrees)
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


class SpinLaueTests(unittest.TestCase):
    def test_layer_operations_in_sun_supplement(self):
        mirror_y = np.diag([1, -1, 1])
        mirror_z = np.diag([1, 1, -1])
        c3z = rotation_z(120)
        group = SpinLaueGroup(2, [(mirror_y, -1)])
        self.assertEqual(len(group.operations), 4)
        for operation, result in (
            (mirror_z, "I"),
            (mirror_z @ c3z, "III"),
            (mirror_z @ mirror_y, "II"),
            (mirror_z @ mirror_y @ c3z, "III"),
            (rotation_z(180), "I"),
            (np.diag([-1, 1, -1]), "II"),
        ):
            with self.subTest(operation=operation):
                self.assertEqual(group.classify(operation), result)
        np.testing.assert_allclose(group.momentum_action(c3z) @ [1, 0],
                                   [-0.5, np.sqrt(3) / 2])

    def test_bulk_actions_preserve_kz_distinctions(self):
        mirror_y = np.diag([1, -1, 1])
        group = SpinLaueGroup(3, [(mirror_y, -1)])
        self.assertEqual(group.classify(-mirror_y), "II")
        self.assertEqual(group.classify(np.diag([1, 1, -1])), "III")
        self.assertEqual(group.classify(-np.eye(3)), "I")
        self.assertEqual(group.dimension, 3)

    def test_closure_products_have_the_composed_spin_action(self):
        mirror = np.diag([1, -1])
        rotation = rotation_z(120)[:2, :2]
        group = SpinLaueGroup(2, [(rotation, 1), (mirror, -1)])
        self.assertEqual(sum(sign == 1 for _, sign in group.operations), 6)
        self.assertEqual(sum(sign == -1 for _, sign in group.operations), 6)
        self.assertEqual(group.classify(-rotation @ mirror), "II")
        self.assertEqual(group.classify(rotation.T), "I")
        self.assertEqual(group.classify(rotation_z(90)), "III")

    def test_full_cubic_group(self):
        cycle = np.eye(3)[[1, 2, 0]]
        group = SpinLaueGroup(3, [(cycle, 1), (rotation_z(90), -1)])
        self.assertEqual(len(group.operations), 48)
        self.assertEqual(group.classify(rotation_z(180)), "I")
        self.assertEqual(group.classify(np.eye(3)[[1, 0, 2]]), "II")

    def test_spin_labels_must_define_an_index_two_subgroup(self):
        mirror = np.diag([1, -1])
        with self.assertRaisesRegex(ValueError, "spin-exchanging"):
            SpinLaueGroup(2, [(mirror, 1)])
        for generators in (
            [(mirror, -1), (mirror, 1)],
            [(-np.eye(2), -1)],
            [(rotation_z(120), -1)],
            [(np.diag([1, 1, -1]), -1)],
        ):
            with self.subTest(generators=generators), self.assertRaisesRegex(ValueError, "conflicting"):
                SpinLaueGroup(2, generators)

    def test_rejects_noncrystallographic_generated_groups(self):
        mirror = np.diag([1, -1])
        rotation = rotation_z(36)[:2, :2]
        # The generators are two valid reflections. Their product has order 5.
        with self.assertRaisesRegex(ValueError, "crystallographic"):
            SpinLaueGroup(2, [(mirror, -1), (rotation @ mirror @ rotation.T, -1)])

    def test_invalid_dimensions_matrices_and_signs(self):
        mirror = np.diag([1, -1, 1])
        for dimension in (1, 4, True, 2.5):
            with self.subTest(dimension=dimension), self.assertRaises(ValueError):
                SpinLaueGroup(dimension, [(mirror, -1)])
        group = SpinLaueGroup(2, [(mirror, -1)])
        for operation in (np.eye(1), [[1, 2], [0, 1]], np.eye(2) * 1j, np.full((2, 2), np.nan),
                          np.eye(3)[[2, 1, 0]], np.diag([1, 1, 0])):
            with self.subTest(operation=operation), self.assertRaises(ValueError):
                group.classify(operation)
        for sign in (0, True, "exchange", -1.5):
            with self.subTest(sign=sign), self.assertRaisesRegex(ValueError, "spin sign"):
                SpinLaueGroup(2, [(mirror, sign)])


if __name__ == "__main__":
    unittest.main()
