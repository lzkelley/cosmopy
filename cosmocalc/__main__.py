"""
"""
from collections import namedtuple
from argparse import ArgumentParser

import numpy as np
import astropy as ap

from . cosmology import Cosmology
from . import PC, ARCSEC, U_SEC, U_CM

funcs = ['luminosity_distance', 'comoving_distance', 'lookback_time', 'age']

_RES_PARS = ["redz", "scale", "dlum", "dcom", "tlbk", "tage"]
Results = namedtuple("Results", _RES_PARS)


def calc(cosmo, args):
    scale = None
    _vals = [None, None, None, None]
    # Redshift
    if args.z is not None:
        redz = args.z
    # Scalefactor
    elif args.a is not None:
        scale = args.a
        redz = cosmo._a_to_z(scale)
    # Age of the Universe
    elif args.ta is not None:
        tage = parse_input(args.ta, U_SEC)
        redz = cosmo.tage_to_z(tage.cgs.value)
        _vals[3] = tage
    # Loockback Time
    elif args.tl is not None:
        tlbk = parse_input(args.tl, U_SEC)
        tage = cosmo.hubble_time - tlbk
        redz = cosmo.tage_to_z(tage.cgs.value)
        _vals[2] = tlbk
        _vals[3] = tage
    # Comoving Distance
    elif args.dc is not None:
        dcom = parse_input(args.dc, U_CM)
        redz = cosmo.dcom_to_z(dcom.cgs.value)
        _vals[1] = dcom
        # Calculate luminosity-distance manually
        _vals[0] = dcom * (1.0 + redz)
    # Luminosity Distance
    elif args.dl is not None:
        dlum = parse_input(args.dl, U_CM)
        redz = cosmo.dcom_to_z(dlum.cgs.value)
        _vals[1] = dlum / (1.0 + redz)
        # Calculate luminosity-distance manually
        _vals[0] = dlum
    else:
        print("__main__.calc()")
        print("\t`args` = '{}'".format(args))
        raise ValueError("No input given!")

    scale = cosmo._z_to_a(redz) if scale is None else scale

    _vals = [getattr(cosmo, ff)(redz) if vv is None else vv
             for ff, vv in zip(funcs, _vals)]
    results = Results(redz, scale, *_vals)
    return results


def output(results):
    """Given a redshift and fractional comoving distance (D_C/D_H), return cosmological parameters.
    """

    dang = results.dcom/(1.0 + results.redz)
    darc = dang*ARCSEC
    dist_mod = 5.0*np.log10(results.dlum.cgs.value/(10.0*PC))

    vals = [results.redz, results.scale, results.dcom, results.dlum,
            dang, darc, results.tlbk,
            results.tage, dist_mod]

    symbs = ['z', 'a', 'D_c', 'D_L',
             'D_A', 'Arcsec', 'T_lb',
             'T_a', 'DM']
    types = [None, None, 'Mpc', 'Mpc',
             'Mpc', 'pc', 'Gyr',
             'Gyr', None]
    names = ['Redshift', 'Scale-factor', 'Comoving Distance', 'Luminosity Distance',
             'Angular Diameter Distance', 'Arcsecond Scale', 'Lookback Time',
             'Age of the Universe', 'Distance Modulus']

    rets = []
    for vv, ss, tt, nn in zip(vals, symbs, types, names):
        try:
            vv = vv.item()
        except AttributeError:
            pass
        v_std = vv if tt is None else vv.to(tt)
        try:
            base = "{:>10s} = {:.4f}".format(ss, v_std)
        except Exception:
            print("Failed constructing `base` string on '{}' ({})".format(nn, ss))
            print("  Value = '{}', type = '{}'".format(v_std, type(v_std)))
            raise

        try:
            v_cgs = vv.cgs
            conv = " ~ {:.4e}".format(v_cgs)
        except Exception:
            conv = ""

        rstr = "{base:30s}{conv:20s} : {name}".format(base=base, conv=conv, name=nn)
        rets.append(rstr)

    print("\n" + "\n".join(rets) + "\n")
    return


def parse_input(inval, unit=None):
    """Convert an input string value into a number, possibly with units.

    e.g. "2.3e12 cm" will be parsed into an `astropy.units.quantity.Quantity`
    """
    try:
        val = np.float(inval)
    except ValueError:
        try:
            val = ap.units.quantity.Quantity(inval)
        except Exception:
            print("__main__.parse_input()")
            print("Failed to convert '{}'".format(inval))
            raise
    else:
        if unit is not None:
            val = val * unit

    # val = val.cgs.value

    return val


def parse_args():
    """Initialize desired command line options, retrieve their values.
    """

    parser = ArgumentParser(description="CosmoCalc: cosmological calculator.")

    # Target Parameters
    parser.add_argument('-z', type=float, default=None,
                        help='Target redshift z')
    parser.add_argument('-a', type=float, default=None,
                        help='Target scale factor a')
    parser.add_argument('-dc', '-cd', default=None,
                        help='Target coming distance D_C')
    parser.add_argument('-dl', '-ld', default=None,
                        help='Target luminosity distance D_L')
    parser.add_argument('-tl', '-lt', default=None,
                        help='Target look-back time T_L')
    parser.add_argument('-ta', '-at', default=None,
                        help='Target universe age T_A')

    '''
    # Modify Cosmological Parameters
    parser.add_argument('--H0', '--h0', type=float, metavar='', default=Parameters.H0, help='Hubble Constant (H0) [km/s/Mpc]')
    parser.add_argument('--T0', '--t0', type=float, metavar='', default=Parameters.T0, help='Hubble Time (T0) [1/s]')
    parser.add_argument('--ODM', '--darkmatter', type=float, metavar='OmegaDM', default=Parameters.OmegaDM, help='Dark Matter fraction (OmegaDM) [1/critical-density]')
    parser.add_argument('--OB', '--baryon', type=float, metavar='OmegaB', default=Parameters.OmegaB, help='Baryon      fraction (OmegaB)  [1/critical-density]')
    parser.add_argument('--OL', '--darkenergy', type=float, metavar='OmegaL', default=Parameters.OmegaL, help='Dark Energy fraction (OmegaL)  [1/critical-density]')
    parser.add_argument('--OR', '--radiation', type=float, metavar='OmegaR', default=Parameters.OmegaR, help='Radtion     fraction (OmegaR)  [1/critical-density]')
    '''
    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    cosmo = Cosmology(args)

    results = calc(cosmo, args)

    output(results)

    return


if __name__ == "__main__":
    main()
