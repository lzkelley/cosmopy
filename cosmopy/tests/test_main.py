"""Test methods for `cosmopy.__main__.py`.

Can be run with:
    $ nosetests cosmopy/tests/test_main.py
    $ nosetests cosmopy/tests/test_main.py:TestMain.test_parse_input

"""
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np

from numpy.testing import run_module_suite
# import scipy as sp
# import scipy.stats
from nose.tools import assert_almost_equal  # , assert_true, assert_false, assert_raises


class TestMain(object):

    @classmethod
    def setup_class(cls):
        np.random.seed(9865)

    def test__parse_input__unitless(self):
        from cosmopy.__main__ import parse_input

        ins = ['-12', '+36', '2.2', '2.3e23', '-1.98e-1', 2.3]
        outs = [-12, 36, 2.2, 2.3e23, -1.98e-1, 2.3]

        for ii, oo in zip(ins, outs):
            print("In: '{}'".format(ii))
            rv = parse_input(ii)
            print("\tOut: '{}' (vs. '{}')".format(rv, oo))
            assert_almost_equal(rv, oo)

        return

    def test__parse_input__defunit(self):
        from cosmopy.__main__ import parse_input
        from astropy.units import Unit

        ins = ['2.3e12', '7.4e2']
        dunits = [Unit('s'), Unit('yr')]
        outs = [2.3e12, 23352624000.0]

        for ii, dd, oo in zip(ins, dunits, outs):
            print("In: '{}'".format(ii))
            rv = parse_input(ii, dd).cgs.value
            print("\tOut: '{}' (vs. '{}')".format(rv, oo))
            assert_almost_equal(rv, oo)

        return

    def test__parse_input__units(self):
        from cosmopy.__main__ import parse_input

        ins = ['-35 kg', '7.4e2 s', '2.2 cm', '2.3e23 yr', '-1.98e-1 pc']
        outs = [-35e3, 7.4e2, 2.2, 7.258248e+30, -6.10964161e+17]

        for ii, oo in zip(ins, outs):
            print("In: '{}'".format(ii))
            rv = parse_input(ii).cgs.value
            print("\tOut: '{}' (vs. '{}')".format(rv, oo))
            assert_almost_equal(rv, oo, delta=np.fabs(1e-6*oo))

        return


# Run all methods as if with `nosetests ...`
if __name__ == "__main__":
    run_module_suite()
