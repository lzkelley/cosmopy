"""Test methods for `cosmopy.cosmology.py`.
"""

import numpy as np
import astropy as ap
import astropy.units

import pytest


def _stats(vv, percs=[0, 25, 50, 75, 100]):
    ss = np.percentile(vv, percs)
    msg = ", ".join(['{:.4e}'.format(_s) for _s in ss])
    return msg


class TestCosmology:

    funcs_forw = ["luminosity_distance", "comoving_distance", "age", "lookback_time"]
    funcs_back = ["dlum_to_z", "dcom_to_z", "tage_to_z", "tlbk_to_z"]
    z_low = 0.001
    vals_low = [1.3353e+25, 1.3340e+25, 4.3356e+17, 4.4475e+14]

    def test_init(self):
        from cosmopy.cosmology import Cosmology
        cosmo = Cosmology()
        print(cosmo)

        # Cannot provide both `h` and `H0`
        with pytest.raises(ValueError):
            Cosmology(h=0.697, H0=2.24e-18)

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
            _a = Cosmology.z_to_a(z)
            _z = Cosmology.a_to_z(a)

            print("z = {} ===> a = {} ({})".format(z, _a, a))
            assert np.isclose(z, _z)

            print("a = {} ===> z = {} ({})".format(a, _z, z))
            assert np.isclose(a, _a)

        # Test array input
        assert (np.allclose(Cosmology.z_to_a(zz), aa))
        assert (np.allclose(Cosmology.a_to_z(aa), zz))

        # Test out-of-domain errors
        with pytest.raises(ValueError):
            Cosmology.a_to_z(-0.1)
        with pytest.raises(ValueError):
            Cosmology.a_to_z(1.1)
        with pytest.raises(ValueError):
            Cosmology.z_to_a(-0.2)

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

        grid, names, units = cosmo._get_grid()
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

    def test_dvc_dz(self):
        from cosmopy.cosmology import Cosmology
        cosmo = Cosmology()

        NUM = 1000
        RTOL = 1e-6

        redz = 10**np.random.uniform(-2, 1, int(NUM))
        redz = np.sort(redz)

        vals = cosmo.dVcdz(redz, cgs=True)
        check = cosmo.differential_comoving_volume(redz).cgs.value * 4*np.pi
        print("redz  : ", _stats(redz))
        print("dVcdz : ", _stats(vals))
        print("check : ", _stats(check))
        assert np.allclose(vals, check, rtol=RTOL), "Error in dvcdz exceeds rtol = {}".format(RTOL)

        # use CGS units
        print("cgs = True")
        vals = cosmo.dVcdz(redz, cgs=True)
        check = cosmo.differential_comoving_volume(redz).cgs.value * 4*np.pi
        print("\ncheck = ", check[:4], check[-4:])

        print("redz  : ", _stats(redz))
        print("dVcdz : ", _stats(vals))
        print("check : ", _stats(check))
        assert np.allclose(vals, check, rtol=RTOL), "Error in dvcdz(cgs=True) exceeds rtol = {}".format(RTOL)

        # use standard units
        print("cgs = False")
        vals = cosmo.dVcdz(redz, cgs=False).value
        check = cosmo.differential_comoving_volume(redz).value * 4*np.pi
        print("\ncheck = ", check[:4], check[-4:])

        print("redz  : ", _stats(redz))
        print("dVcdz : ", _stats(vals))
        print("check : ", _stats(check))
        assert np.allclose(vals, check, rtol=RTOL), "Error in dvcdz(cgs=False) exceeds rtol = {}".format(RTOL)

        return

    def test_dtdz(self):
        from cosmopy.cosmology import Cosmology
        cosmo = Cosmology()

        NUM = 1000
        redz = 10 ** np.random.uniform(-2, 1, int(NUM))
        redz = np.sort(redz)

        def check(zz):
            zz = np.asarray(zz)
            efac = cosmo.efunc(zz)
            th = cosmo.hubble_time.to('s').value
            rv = th / (1.0 + zz) / efac
            return rv

        vals = cosmo.dtdz(redz)
        check = check(redz)
        assert np.allclose(vals, check)

        # precomputed values for [Omega0 = 0.2880, OmegaBaryon = 0.0472, HubbleParam = 0.6933]
        test = [0.001, 0.01, 0.1, 0.5, 1.0, 5.0, 10.0]
        check = [4.4443e+17, 4.3875e+17, 3.8660e+17, 2.2865e+17, 1.2814e+17, 9.3516e+15, 2.0647e+15]
        vv = cosmo.dtdz(test)
        assert np.allclose(vv, check, rtol=1e-4)

        return


def _get_err(aa, bb):
    err = (aa - bb) / np.minimum(aa, bb)
    err = np.fabs(err)
    return err.mean(), err.max()


def test_z_to__methods():
    from cosmopy.cosmology import Cosmology
    cosmo = Cosmology()

    NUM = 1000
    AVE_THRESH = 1e-6
    MAX_THRESH = 1e-4

    redz = 10**np.random.uniform(-2, 1, int(NUM))
    redz = np.sort(redz)

    funcs = ['tage', 'tlbk', 'dlum', 'dcom']
    funcs_check = ['age', 'lookback_time', 'luminosity_distance', 'comoving_distance']
    funcs_for = ["z_to_" + ff for ff in funcs]
    funcs_back = [ff + "_to_z" for ff in funcs]

    for ii, (func_for, check, func_back) in enumerate(zip(funcs_for, funcs_check, funcs_back)):
        _err = "{} :: Failed - ".format(funcs[ii])

        val_for = getattr(cosmo, func_for)(redz)
        chk = getattr(cosmo, check)(redz).cgs.value
        val_back = getattr(cosmo, func_back)(val_for)

        chk_ave, chk_max = _get_err(val_for, chk)
        print("forward errors: ave={:.4e}, max={:.4e}".format(chk_ave, chk_max))

        # --- Convert forwards, and check against astropy default calculation
        err = _err + "forward error exceeded average threshold {} > {}!".format(chk_ave, AVE_THRESH)
        assert chk_ave < AVE_THRESH, err

        err = _err + "forward error exceeded maximum threshold {} > {}!".format(chk_max, MAX_THRESH)
        assert chk_max < MAX_THRESH, err

        # --- Convert backwards, and check against input redshifts
        bck_ave, bck_max = _get_err(val_back, redz)
        print("backward errors: ave={:.4e}, max={:.4e}".format(bck_ave, bck_max))

        err = _err + "backward error exceeded average threshold {} > {}!".format(bck_ave, AVE_THRESH)
        assert bck_ave < AVE_THRESH, err

        err = _err + "backward error exceeded maximum threshold {} > {}!".format(bck_max, MAX_THRESH)
        assert bck_max < MAX_THRESH, err

    return


def test_a_to__methods():
    from cosmopy.cosmology import Cosmology
    cosmo = Cosmology()

    NUM = 1000
    AVE_THRESH = 1e-6
    MAX_THRESH = 1e-4

    _redz = 10**np.random.uniform(-2, 1, int(NUM))
    _redz = np.sort(_redz)
    scafa = cosmo.z_to_a(_redz)

    funcs = ['tage', 'tlbk', 'dlum', 'dcom']
    funcs_check = ['age', 'lookback_time', 'luminosity_distance', 'comoving_distance']
    funcs_for = ["a_to_" + ff for ff in funcs]
    funcs_back = [ff + "_to_a" for ff in funcs]

    for ii, (func_for, check, func_back) in enumerate(zip(funcs_for, funcs_check, funcs_back)):
        _err = "{} :: Failed - ".format(funcs[ii])

        val_for = getattr(cosmo, func_for)(scafa)
        chk = getattr(cosmo, check)(_redz).cgs.value
        val_back = getattr(cosmo, func_back)(val_for)
        print(funcs[ii], "scafa    : ", _stats(scafa))
        print(funcs[ii], "redz     : ", _stats(_redz))
        print(funcs[ii], "forward  : ", _stats(val_for))
        print(funcs[ii], "check    : ", _stats(chk))
        print(funcs[ii], "backward : ", _stats(val_back))

        chk_ave, chk_max = _get_err(val_for, chk)
        print("{} forward errors: ave={:.4e}, max={:.4e}".format(funcs[ii], chk_ave, chk_max))

        # --- Convert forwards, and check against astropy default calculation
        err = _err + "forward error exceeded average threshold {} > {}!".format(chk_ave, AVE_THRESH)
        assert chk_ave < AVE_THRESH, err

        err = _err + "forward error exceeded maximum threshold {} > {}!".format(chk_max, MAX_THRESH)
        assert chk_max < MAX_THRESH, err

        # --- Convert backwards, and check against input redshifts
        bck_ave, bck_max = _get_err(val_back, scafa)
        print("{} backward errors: ave={:.4e}, max={:.4e}".format(funcs[ii], bck_ave, bck_max))

        err = _err + "backward error exceeded average threshold {} > {}!".format(bck_ave, AVE_THRESH)
        assert bck_ave < AVE_THRESH, err

        err = _err + "backward error exceeded maximum threshold {} > {}!".format(bck_max, MAX_THRESH)
        assert bck_max < MAX_THRESH, err

    return
