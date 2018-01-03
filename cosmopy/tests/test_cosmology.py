"""Test methods for `cosmopy.cosmology.py`.
"""
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np

from numpy.testing import run_module_suite
# import scipy as sp
# import scipy.stats
from nose.tools import assert_almost_equal  # , assert_true, assert_false, assert_raises


class TestCosmology(object):

    @classmethod
    def setup_class(cls):
        np.random.seed(986523)

    def test_init(self):
        from cosmopy.cosmology import Cosmology
        cosmo = Cosmology()
        print(cosmo)
        return

# Run all methods as if with `nosetests ...`
if __name__ == "__main__":
    run_module_suite()
