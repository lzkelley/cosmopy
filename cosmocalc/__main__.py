"""
"""
import argparse

import numpy as np
import astropy as ap

from . cosmology import Cosmology
from . import PC, ARCSEC, U_SEC, U_CM

_RESULTS_FUNCS = ['luminosity_distance', 'comoving_distance', 'lookback_time', 'age']

_RESULTS_FUNC_PARS = ['dl', 'dc', 'tl', 'ta']
_RESULTS_PARS = ['z', 'a'] + _RESULTS_FUNC_PARS


def calc_basic(cosmo, args):
    """Given a `Cosmology` instance and input arguments, calculate basic cosmological values.

    Parameters
    ----------
    cosmo : `Cosmology`
        Class for calculating cosmological parameters, based on `astropy.cosmology` classes.
    args : `argparse.Namespace`
        Namespace containing input parameters.
        One of its parameters, {`a`, `z`, `dc`, `dl`, `tl`, `ta`} must be non-None.

    Returns
    -------
    results : `Results` namedtuple
        Namedtuple of basic cosmological parameter results which can be processed into outputs.

    """
    scale = None

    results = {kk: None for kk in _RESULTS_PARS}

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
        results['ta'] = tage

    # Loockback Time
    elif args.tl is not None:
        tlbk = parse_input(args.tl, U_SEC)
        tage = cosmo.hubble_time - tlbk
        redz = cosmo.tage_to_z(tage.cgs.value)
        results['ta'] = tage
        results['tl'] = tlbk

    # Comoving Distance
    elif args.dc is not None:
        dcom = parse_input(args.dc, U_CM)
        redz = cosmo.dcom_to_z(dcom.cgs.value)
        # Calculate luminosity-distance manually
        dlum = dcom * (1.0 + redz)
        results['dc'] = dcom
        results['dl'] = dlum

    # Luminosity Distance
    elif args.dl is not None:
        dlum = parse_input(args.dl, U_CM)
        redz = cosmo.dlum_to_z(dlum.cgs.value)
        # Calculate comoving-distance manually
        dcom = dlum / (1.0 + redz)
        results['dc'] = dcom
        results['dl'] = dlum

    # No arguments set, raise error
    else:
        print("__main__.calc()")
        print("\t`args` = '{}'".format(args))
        raise ValueError("No input given!")

    # Calculate scale-factor if needed
    scale = cosmo._z_to_a(redz) if scale is None else scale
    results['a'] = scale
    results['z'] = redz

    # Calculate the values not already set
    for pp, ff in zip(_RESULTS_FUNC_PARS, _RESULTS_FUNCS):
        if results[pp] is None:
            func = getattr(cosmo, ff)
            results[pp] = func(redz)

    return results


def calc_derived(results):
    dang = results['dc']/(1.0 + results['z'])
    darc = dang*ARCSEC
    dist_mod = 5.0*np.log10(results['dl'].cgs.value/(10.0*PC))

    results['da'] = dang
    results['arc'] = darc
    results['dm'] = dist_mod

    return results


def output(results, print_output=True):
    """Format and print the given cosmological parameters to stdout.

    Arguments
    ---------
    results : `Results` namedtuple
        Cosmological parameters to output.

    """

    keys = ['z', 'a', 'dc', 'dl',
            'da', 'arc', 'tl',
            'ta', 'dm']
    symbs = ['z', 'a', 'D_c', 'D_L',
             'D_A', 'Arcsec', 'T_lb',
             'T_a', 'DM']
    types = [None, None, 'Mpc', 'Mpc',
             'Mpc', 'pc', 'Gyr',
             'Gyr', None]
    names = ['Redshift', 'Scale-factor', 'Comoving Distance', 'Luminosity Distance',
             'Angular Diameter Distance', 'Arcsecond Scale', 'Lookback Time',
             'Age of the Universe', 'Distance Modulus']

    outs = []
    retvals = {}
    for kk, ss, tt, nn in zip(keys, symbs, types, names):
        vv = results[kk]
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
        retvals[kk] = rstr
        outs.append(rstr)

    if print_output:
        print("\n" + "\n".join(outs) + "\n")

    return retvals


def parse_input(inval, unit=None):
    """Convert an input string value into a number, possibly with units.

    If `unit` is given and the input is a unitless float, then the given unit is added.
    e.g. "2.3e12 cm" will be parsed into an `astropy.units.quantity.Quantity`

    Arguments
    ---------
    inval : str
    unit : `astropy.unit.Unit` or `None`
        Unit to be added to the given input (if it doesn't already specify a unit value).

    Returns
    -------
    val : `astropy.unit.quantity.Quantity`
        Parsed value.  This value will have units attached if either the input string specified
        them, or if a `unit` value was given.

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

    return val


def parse_args():
    """Initialize desired command line options, retrieve their values.
    """

    parser = argparse.ArgumentParser(description="CosmoCalc: cosmological calculator.")

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


def api(key, val, cosmo=None):
    """Primary method for external API access to this package.
    """
    args = {kk: None for kk in _RESULTS_PARS}
    if key not in args:
        print("Valid keys: '{}'".format(_RESULTS_PARS))
        raise ValueError("Unrecognized key '{}'".format(key))

    # Set the given value
    args[key] = val
    # Convert to `Namespace` object (expected by the `calc_basic` method)
    args = argparse.Namespace(**args)

    # Load cosmology object if needed
    if cosmo is None:
        cosmo = get_cosmology()

    # Calculate
    # ----------------
    # Calculate basic cosmological values
    results = calc_basic(cosmo, args)
    # Calculate additional derived values
    results = calc_derived(results)
    # Format output string and return as dictionary
    retvals = output(results, print_output=False)

    return retvals


def get_cosmology():
    cosmo = Cosmology()
    return cosmo


def main():
    """Primary computational method for command-line access to calculations.

    Loads command line arguments and a `Cosmology` instance.
    Calculates cosmological distance measures.
    Formats and prints the output to stdout.

    """
    # Initialize
    # ----------------
    # Load arguments
    args = parse_args()
    # Load cosmology object
    cosmo = get_cosmology()

    # Calculate
    # ----------------
    # Calculate basic cosmological values
    results = calc_basic(cosmo, args)
    # Calculate additional derived values
    results = calc_derived(results)
    # Format and print output
    output(results)

    return


if __name__ == "__main__":
    main()
