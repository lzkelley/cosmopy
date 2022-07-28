"""Test methods for `cosmopy.cosmology.py`.
"""

import numpy as np
import astropy as ap
import astropy.units

import pytest


class TestCosmology:

    funcs_forw = ["luminosity_distance", "comoving_distance", "age", "lookback_time"]
    funcs_back = ["dlum_to_z", "dcom_to_z", "tage_to_z", "tlbk_to_z"]
    z_low = 0.001
    vals_low = [1.3359e+25, 1.3346e+25, 4.3459e+17, 4.4494e+14]

    def test_init(self):
        from cosmopy.cosmology import Cosmology
        cosmo = Cosmology()
        print(cosmo)
        return

    def test__init_interp_grid(self):
        from cosmopy.cosmology import Cosmology

        zpnts = [10.0, 1.0]
        numz = 10

        grid = Cosmology._init_interp_grid(zpnts, numz)
        assert grid.size == numz * len(zpnts)

        with pytest.raises(ValueError):
            Cosmology._init_interp_grid(zpnts[::-1], numz)

        return

    def test_interp_grid(self):
        from cosmopy.cosmology import Cosmology
        cosmo = Cosmology()
        grid_z = cosmo._grid_z
        num_z = grid_z.size
        print("z-grid has size: {}".format(num_z))
        assert np.isclose(grid_z[-1], 0.0)

        print("cosmo._grid_a.size, num_z = {}, {}".format(cosmo._grid_a.size, num_z))
        assert (cosmo._grid_a.size == num_z)
        print("cosmo._sort_z.size, num_z = {}, {}".format(cosmo._sort_z.size, num_z))
        assert (cosmo._sort_z.size == num_z)

        dz = np.diff(grid_z[cosmo._sort_z])
        assert (np.all(dz > 0.0))

        return

    def test_interp(self):
        from cosmopy.cosmology import Cosmology

        xx = np.array([1.0, 2.0])
        yy = np.array([10.0, 20.0])
        val = 15.0

        zz = Cosmology._interp(1.5, xx, yy)
        print("Interp = {}, should be = {}".format(zz, val))
        assert np.isclose(zz, val)

        zz = Cosmology._interp(1.5, xx[::-1], yy[::-1])
        assert np.isclose(zz, val)

        # Out of domain should be NaN
        bad = Cosmology._interp(0.5, xx, yy)
        assert np.isnan(bad)

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
            assert np.isclose(z, _z)

            print("a = {} ===> z = {} ({})".format(a, _z, z))
            assert np.isclose(a, _a)

        # Test array input
        assert (np.allclose(Cosmology._z_to_a(zz), aa))
        assert (np.allclose(Cosmology._a_to_z(aa), zz))

        # Test out-of-domain errors
        with pytest.raises(ValueError):
            Cosmology._a_to_z(-0.1)
        with pytest.raises(ValueError):
            Cosmology._a_to_z(1.1)
        with pytest.raises(ValueError):
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
            assert np.isclose(_v/exp, vv/exp, 1e-3)

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
            assert np.isclose(z, _z, atol=1e-3)

        return

    def test__get_grid(self):
        from cosmopy.cosmology import Cosmology
        cosmo = Cosmology()

        grid, names, units = cosmo.get_grid()
        z_grid = cosmo._grid_z
        num = z_grid.size
        print(names)
        print(units)
        assert grid.shape == (num, 6,)
        assert np.allclose(grid[:, 0], cosmo._grid_z)

        ind = names.index('dl')
        dl = cosmo._grid_dlum
        un = ap.units.Unit(units[ind])
        print("'dl' ind = {}, unit = {}".format(ind, un))
        grid_dl = grid[:, ind] * un
        grid_dl = grid_dl.cgs.value
        print("dl ind = {}, {}, {}".format(ind, dl[:4], grid_dl[:4]))
        assert np.allclose(grid_dl, dl)
