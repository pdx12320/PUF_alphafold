"""Boundary, missing-outcome and training-only transformation regression checks."""
import importlib.util
from pathlib import Path
import sys
import os
import unittest
import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location('binary_runner', Path(__file__).with_name('train.py'))
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


class Labels(unittest.TestCase):
    def test_boundary_and_increases(self):
        d = pd.DataFrame({'C388_delta_pp': [-12.50001, -12.5, -12.49999, 0., 20., np.nan]})
        np.testing.assert_array_equal(runner.binary_labels(d, 'C388', 12.5), [0, 1, 1, 1, 1, -1])

    def test_scope(self):
        with self.assertRaises(ValueError):
            runner.binary_labels(pd.DataFrame(), 'C295', 12.5)


@unittest.skipUnless(os.environ.get('C388_V4'), 'Set C388_V4 for archived-data leakage check')
class TrainOnlyTransform(unittest.TestCase):
    def test_heldout_values_do_not_change_selection_or_training_embedding(self):
        sys.path.insert(0, str(Path(os.environ['C388_V4']) / 'src'))
        import core
        ds = core.Dataset('C388', 'combined_PR')
        c = next(c for c in core.CFG['representations']
                 if c['family'] == 'PR_train_union_interface' and c['mode'] == 'raw'
                 and c['retain'] == 25 and c['ranking'] == 'variance'
                 and not c['pca'] and not c['cp_trm'])
        tr = np.arange(24)
        original = ds.transform(c, tr)
        ds.blocks['PR'][24:] = 0.99
        ds.basecache.clear()
        ds.statecache.clear()
        ds.reprcache.clear()
        ds.gramcache.clear()
        changed = ds.transform(c, tr)
        for key in ['selected', 'mean', 'scale', 'variance']:
            np.testing.assert_array_equal(original[key], changed[key])
        np.testing.assert_allclose(original['embedding'][tr], changed['embedding'][tr], atol=1e-10)


if __name__ == '__main__':
    unittest.main()
