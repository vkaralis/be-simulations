# Author: Vangelis D. Karalis
"""Verify interactive scripts and their saved artifacts."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from PIL import Image


class CLITests(unittest.TestCase):
    def test_all_interactive_scripts(self):
        root = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as directory:
            env = dict(os.environ, MPLCONFIGDIR=directory, PYTHONDONTWRITEBYTECODE='1')
            for design in ('2x2', '2x3', '2x4', '3x3'):
                run = subprocess.run(
                    [sys.executable, str(root / f'BE{design}.py'),
                     '--seed', '42', '--output', directory],
                    input='24\n1000\n30\n0.8:1.3:0.05\n',
                    capture_output=True, text=True, env=env, check=True)
                self.assertIn('Monte Carlo', run.stdout)
                stem = Path(directory) / f'power_{design}'
                metadata = json.loads(stem.with_suffix('.json').read_text())
                self.assertEqual(metadata['author'], 'Vangelis D. Karalis')
                self.assertEqual(metadata['design'], design)
                self.assertEqual(len(stem.with_suffix('.csv').read_text().splitlines()), 12)
                with Image.open(stem.with_suffix('.png')) as image:
                    self.assertEqual(image.info['Author'], 'Vangelis D. Karalis')
                    self.assertEqual(image.info['Software'], 'BEsim')
                    image.verify()

    def test_invalid_input_has_clean_error(self):
        run = subprocess.run(
            [sys.executable, '-m', 'besim', '--design', '2x2',
             '--subjects', '5', '--trials', '100', '--cv', '30', '--gmr', '1'],
            capture_output=True, text=True)
        self.assertEqual(run.returncode, 2)
        self.assertIn('divisible', run.stderr)
        self.assertNotIn('Traceback', run.stderr)
