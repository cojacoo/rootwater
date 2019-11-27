import unittest

import os
import numpy as np
import pandas as pd
from numpy.testing import assert_almost_equal

from rootwater import rw, sf

# get the basebath for test reference files
BASEPATH = os.path.abspath(os.path.dirname(__file__))

def _read_helper(fname):
    df = pd.read_csv(os.path.join(BASEPATH, fname),index_col=0)
    df.index = pd.to_datetime(df.index)
    return df


class TestIt(unittest.TestCase):
    def setUp(self):
        # read references
        self.SMtest = _read_helper('SM_test.csv')
        self.RWUtest = _read_helper('RWU_test.csv')
        self.SFtest = _read_helper('SF_test.csv')
        self.SVtest = _read_helper('SV_test.csv')

    def test_RWU(self):
        assert_almost_equal(
            pd.concat(rw.dfRWUc(self.SMtest)).dropna().values,
            self.RWUtest.values,
            decimal=2
            )
        

    def test_SF(self):
        assert_almost_equal(
            sf.sap_calc(self.SVtest,32.,0.95,'beech').values,
            self.SFtest.values,
            decimal=2
        )


if __name__ == '__main__':
    unittest.main()
