import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

from ammoire import make_model
from ammoire.fe_response import bilayer_demo, domain_response
from ammoire.fe_symmetry import SpinLaueGroup


class FerroelectricResponseTests(unittest.TestCase):
    def test_cli_manifest_through_module_entry_point(self):
        root = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'report.json'
            result = subprocess.run(
                [sys.executable, '-m', 'ammoire', 'fe-response',
                 'examples/fe_square_dxy.json', '--out', str(output)],
                cwd=root, capture_output=True, text=True, check=True,
            )
            report = json.loads(output.read_text(encoding='utf-8'))
        self.assertEqual(json.loads(result.stdout), report)
        self.assertEqual(report['response_type'], 'II')
        self.assertEqual(report['sample_verdict'], 'consistent_on_sample')
        self.assertEqual(report['two_domain_spectra']['sample_count'], 4)
        self.assertEqual(report['positive_domain_group']['tested_operations'], 8)

    def test_bilayer_demo_distinguishes_all_three_predictions(self):
        paths = bilayer_demo()['paths']
        self.assertEqual([p['response_type'] for p in paths], ['I', 'II', 'III'])
        self.assertEqual([p['two_domain_spectra']['same_k_relation'] for p in paths],
                         ['preserved', 'reversed', 'neither'])
        for path in paths:
            self.assertEqual(path['sample_verdict'], 'consistent_on_sample')
            self.assertTrue(path['positive_domain_group']['matches'])
            self.assertTrue(path['two_domain_spectra']['mapped_same_spin_matches'])

    def test_group_incompatible_with_an_otherwise_consistent_path(self):
        model = make_model('square-dx2-y2')
        mirror = np.diag([1, -1])
        incorrect_group = SpinLaueGroup(2, [(mirror, -1)])
        report = domain_response(incorrect_group, np.eye(2), model, model, [[0.3, 0.9]])
        self.assertEqual(report['response_type'], 'I')
        self.assertTrue(report['two_domain_spectra']['mapped_same_spin_matches'])
        self.assertFalse(report['positive_domain_group']['matches'])
        self.assertEqual(report['sample_verdict'], 'inconsistent')

    def test_node_only_sample_does_not_resolve_spin_reversal(self):
        mirror = np.diag([1, -1])
        group = SpinLaueGroup(2, [(mirror, -1)])
        plus, minus = (make_model('square-dxy', j=j) for j in (0.2, -0.2))
        report = domain_response(group, mirror, plus, minus, [[0, 0], [0.7, 0]])
        self.assertEqual(report['response_type'], 'II')
        self.assertEqual(report['sample_verdict'], 'undistinguished_on_sample')

    def test_bulk_section_at_generic_kz_requires_full_3d_momenta(self):
        plus, minus = (make_model('cubic-i', j=j) for j in (0.2, -0.2))
        exchange_xy = np.eye(3)[[1, 0, 2]]
        group = SpinLaueGroup(3, [(exchange_xy, -1)])
        report = domain_response(group, exchange_xy, plus, minus, [[1.2, -0.9, 0.53]])
        self.assertEqual(report['sample_verdict'], 'consistent_on_sample')
        self.assertEqual(report['two_domain_spectra']['same_k_relation'], 'reversed')
        plane_group = SpinLaueGroup(2, [(exchange_xy, -1)])
        with self.assertRaisesRegex(ValueError, '3D bulk'):
            domain_response(plane_group, exchange_xy, plus.slice_kz(0.53), minus.slice_kz(0.53),
                            [[1.2, -0.9]])


if __name__ == '__main__':
    unittest.main()
