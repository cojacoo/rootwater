import unittest

import numpy as np
import pandas as pd
from numpy.testing import assert_almost_equal

import rootwater 

# read references
SMtest = pd.read_csv('SM_test.cvs',index_col=0)
SMtest.index = pd.to_datetime(SMtest.index)

RWUtest = pd.read_csv('RWU_test.cvs',index_col=0)
RWUtest.index = pd.to_datetime(RWUtest.index)

SFtest = pd.read_csv('SF_test.cvs',index_col=0)
SFtest.index = pd.to_datetime(SFtest.index)

SVtest = pd.read_csv('SV_test.cvs',index_col=0)
SVtest.index = pd.to_datetime(SVtest.index)



class TestIt(unittest.TestCase):
    def test_RWU(self):
        assert_almost_equal(
            pd.concat(rw.dfRWUc(SMtest)).dropna().values,
            RWUtest.values,
            decimal=2
            )
        

    def test_SF(self):
        assert_almost_equal(
            sf.sap_calc(SVtest,32.,'beech').values,
            SFtest.values,
            decimal=2
        )


if __name__ == '__main__':
    unittest.main()
