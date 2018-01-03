"""Test methods for `cosmopy.cosmology.py`.
"""
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import astropy as ap
import astropy.units

from numpy.testing import run_module_suite
from nose.tools import assert_almost_equal, assert_true, assert_raises, assert_equal
import unittest


class TestCosmology(unittest.TestCase):

    def setUp(self):
        np.random.seed(986523)

        self.funcs_forw = ["luminosity_distance", "comoving_distance", "age", "lookback_time"]
        self.funcs_back = ["dlum_to_z", "dcom_to_z", "tage_to_z", "tlbk_to_z"]
        self.z_low = 0.001
        self.vals_low = [1.3359e+25, 1.3346e+25, 4.3459e+17, 4.4494e+14]

    def test_init(self):
        from cosmopy.cosmology import Cosmology
        cosmo = Cosmology()
        print(cosmo)
        return

    def test_interp_grid(self):
        from cosmopy.cosmology import Cosmology
        cosmo = Cosmology()
        grid_z = cosmo._grid_z
        num_z = grid_z.size
        print("z-grid has size: {}".format(num_z))
        assert_almost_equal(grid_z[-1], 0.0)

        print("cosmo._grid_a.size, num_z = {}, {}".format(cosmo._grid_a.size, num_z))
        assert_true(cosmo._grid_a.size == num_z)
        print("cosmo._sort_z.size, num_z = {}, {}".format(cosmo._sort_z.size, num_z))
        assert_true(cosmo._sort_z.size == num_z)

        dz = np.diff(grid_z[cosmo._sort_z])
        assert_true(np.all(dz > 0.0))

        return

    def test_interp(self):
        from cosmopy.cosmology import Cosmology

        xx = np.array([1.0, 2.0])
        yy = np.array([10.0, 20.0])
        val = 15.0

        zz = Cosmology._interp(1.5, xx, yy)
        print("Interp = {}, should be = {}".format(zz, val))
        assert_almost_equal(zz, val)

        zz = Cosmology._interp(1.5, xx[::-1], yy[::-1])
        assert_almost_equal(zz, val)

        # Out of domain should be NaN
        bad = Cosmology._interp(0.5, xx, yy)
        assert_true(np.isnan(bad))

        return

    def test_a_z_conv(self):
        from cosmopy.cosmology import Cosmology
        # These are static-methods so we dont need an instance

        zz = [3.0, 1.0, 0.5, 0.0]
        aa = [0.25, 0.5, 0.666666666, 1.0]

        for z, a in zip(zz, aa):
            _a = Cosmology._z_to_a(z)
            _z = Cosmology._a_to_z(a)

            print("z = {} ===> a = {} ({})".format(z, _a, a))
            assert_almost_equal(z, _z)

            print("a = {} ===> z = {} ({})".format(a, _z, z))
            assert_almost_equal(a, _a)

        # Test array input
        assert_true(np.allclose(Cosmology._z_to_a(zz), aa))
        assert_true(np.allclose(Cosmology._a_to_z(aa), zz))

        # Test out-of-domain errors
        with assert_raises(ValueError) as cm:
            Cosmology._a_to_z(-0.1)
        with assert_raises(ValueError) as cm:
            Cosmology._a_to_z(1.1)
        with assert_raises(ValueError) as cm:
            Cosmology._z_to_a(-0.2)

        return

    def test_forward(self):
        """Test "forward" conversion from redshift to other parameters
        """
        from cosmopy.cosmology import Cosmology
        cosmo = Cosmology()

        z = self.z_low
        for vv, ff in zip(self.vals_low, self.funcs_forw):
            _v = getattr(cosmo, ff)(z).cgs.value
            exp = np.power(10.0, np.floor(np.log10(vv)))
            print("{}({}) = {} (should be: {})".format(ff, z, _v, vv))
            assert_almost_equal(_v/exp, vv/exp, 3)

        return

    def test_backward(self):
        """Test "backward" conversion from parameters to redshift
        """
        from cosmopy.cosmology import Cosmology
        cosmo = Cosmology()

        z = self.z_low
        for vv, ff in zip(self.vals_low, self.funcs_back):
            _z = getattr(cosmo, ff)(vv)
            print("{}({}) = {} (should be: {})".format(ff, vv, _z, z))
            assert_almost_equal(z, _z, 3)

        return

    def test__get_grid(self):
        from cosmopy.cosmology import Cosmology
        cosmo = Cosmology()

        grid, names, units = cosmo.get_grid()
        z_grid = cosmo._grid_z
        num = z_grid.size
        print(names)
        print(units)
        assert_equal(grid.shape, (num, 6,))
        assert_true(np.allclose(grid[:, 0], cosmo._grid_z))

        ind = names.index('dl')
        dl = cosmo._grid_dlum
        un = ap.units.Unit(units[ind])
        print("'dl' ind = {}, unit = {}".format(ind, un))
        grid_dl = grid[:, ind] * un
        grid_dl = grid_dl.cgs.value
        print("dl ind = {}, {}, {}".format(ind, dl[:4], grid_dl[:4]))
        assert_true(np.allclose(grid_dl, dl))


# Run all methods as if with `nosetests ...`
if __name__ == "__main__":
    run_module_suite()
