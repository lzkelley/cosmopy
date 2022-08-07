"""Test methods for `cosmopy.__main__.py`.
"""

import runpy

import numpy as np
import astropy as ap
import astropy.units

import pytest


class TestMain:

    np.random.seed(986523)
    z_low = 0.001
    vals_low = [1.3353e+25, 1.3340e+25, 4.3356e+17, 4.4475e+14]
    arg_keys = ['dl', 'dc', 'ta', 'tl']

    def test__module(self):
        with pytest.raises(SystemExit):
            runpy.run_module('cosmopy', run_name='__main__')

        # globs = {'sys.argv': ['cosmopy', '-z', '0.1']}
        # print("\n\n\n=================================================\n\n\n")
        # runpy.run_module('cosmopy', run_name='__main__', init_globals=globs)
        return

    def test__parse_args(self):
        from cosmopy.__main__ import parse_args

        with pytest.raises(SystemExit):
            parse_args([])

        with pytest.raises(SystemExit):
            parse_args(['-v'])

        vals = self.vals_low + [0.1, 0.5]
        keys = self.arg_keys + ['z', 'a']

        for vv, kk in zip(vals, keys):
            args = ["-" + kk, "{:.4e}".format(vv)]
            print("args = ", args)
            parse_args(args)

        return

    def test__main(self):
        from cosmopy.__main__ import main

        with pytest.raises(SystemExit):
            main([])

        with pytest.raises(SystemExit):
            main(['-v'])

        keys = self.arg_keys + ['z', 'a']
        vals = self.vals_low + [0.1, 0.5]

        for kk, vv in zip(keys, vals):
            args = ["-" + kk, "{:.4e}".format(vv)]
            print("args = ", args)
            main(args)

        return

    def test__parse_input(self):
        from cosmopy.__main__ import parse_input
        from astropy.units import Unit

        ins = ['-12', '+36', '2.2', '2.3e23', '-1.98e-1', 2.3]
        outs = [-12, 36, 2.2, 2.3e23, -1.98e-1, 2.3]

        for ii, oo in zip(ins, outs):
            print("In: '{}'".format(ii))
            rv = parse_input(ii)
            print("\tOut: '{}' (vs. '{}')".format(rv, oo))
            assert np.isclose(rv, oo)

        ins = ['2.3e12', '7.4e2']
        dunits = [Unit('s'), Unit('yr')]
        outs = [2.3e12, 23352624000.0]

        for ii, dd, oo in zip(ins, dunits, outs):
            print("In: '{}'".format(ii))
            rv = parse_input(ii, dd).cgs.value
            print("\tOut: '{}' (vs. '{}')".format(rv, oo))
            assert np.isclose(rv, oo)

        ins = ['-35 kg', '7.4e2 s', '2.2 cm', '2.3e23 yr', '-1.98e-1 pc']
        outs = [-35e3, 7.4e2, 2.2, 7.258248e+30, -6.10964161e+17]

        for ii, oo in zip(ins, outs):
            print("In: '{}'".format(ii))
            rv = parse_input(ii).cgs.value
            print("\tOut: '{}' (vs. '{}')".format(rv, oo))
            assert np.isclose(rv, oo, atol=np.fabs(1e-6*oo))

        # Should fail with unrecognized units
        with pytest.raises(ValueError):
            parse_input('12 Monkeys')

        return

    def test__get_cosmology(self):
        from cosmopy.__main__ import get_cosmology

        cosmo = get_cosmology()
        print(cosmo)

        dl = cosmo.luminosity_distance(self.z_low).cgs.value
        print(dl, self.vals_low[0])
        assert np.isclose(dl, self.vals_low[0], rtol=1e-3)

        z = cosmo.dlum_to_z(dl)
        print(z, self.z_low)
        assert np.isclose(z, self.z_low)
        return

    def test__api(self):
        from cosmopy.__main__ import api

        # Use a key that shouldn't be recognized
        with pytest.raises(KeyError):
            api('m', 0.2)

        # These values were calculated using built-in astropy cosmology instance, for comparison
        outs = {
            'z': '0.1000', 'a': '0.9091', 'ta': '12.4377 Gyr', 'dl': '465.1654 Mpc',
            'da': '384.4342 Mpc', 'tl': '1.3150 Gyr', 'dc': '422.8776 Mpc',
            'arc': '1863.7897 pc', 'dm': '38.3380'
        }

        def compare(aa):
            rv = api(*aa)
            tests = []
            print("input args = ", args)
            for kk, v1 in outs.items():
                v2 = rv[kk]
                q1 = ap.units.Quantity(v1).cgs.value
                q2 = ap.units.Quantity(v2).cgs.value
                tt = np.isclose(q1, q2, rtol=1e-3)
                print("{} : {} vs {} :: {}".format(kk, v1, v2, tt))
                tests.append(tt)

            assert np.all(tests)

        args = ['z', '0.1']
        compare(args)

        args = ['dc', '423.1241Mpc']
        compare(args)

        args = ['dc', '0.4231241Gpc']
        compare(args)

        return

    def test__calc_basic(self):
        from cosmopy.__main__ import calc_basic, _RESULTS_PARS, get_cosmology
        import argparse

        cosmo = get_cosmology()

        # Should raise error with no settings set
        sets = {kk: None for kk in _RESULTS_PARS}
        args = argparse.Namespace(**sets)
        with pytest.raises(ValueError):
            calc_basic(cosmo, args)

        return
