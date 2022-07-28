"""Interface to the cosmopy package via either command-line or python API.
"""
import argparse
import os
import sys

import click

import numpy as np
import astropy as ap

from . cosmology import Cosmology
from . import PC, ARCSEC, U_SEC, U_CM, COLORED_OUTPUT

_RESULTS_FUNCS = ['luminosity_distance', 'comoving_distance', 'lookback_time', 'age']

_RESULTS_FUNC_PARS = ['dl', 'dc', 'tl', 'ta']
_RESULTS_PARS = ['z', 'a'] + _RESULTS_FUNC_PARS

_COLOR_BASE = "magenta"  # "cyan"
_COLOR_CONV = "green"  # "yellow"
_COLOR_NAME = "red"

# Enable imperial units (e.g. 'miles')
ap.units.imperial.enable()


def calc_basic(cosmo, sets):
    """Given a `Cosmology` instance and input arguments, calculate basic cosmological values.

    Parameters
    ----------
    cosmo : `Cosmology`
        Class for calculating cosmological parameters, based on `astropy.cosmology` classes.
    sets : `argparse.Namespace`
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
    if sets.z is not None:
        redz = parse_input(sets.z)

    # Scalefactor
    elif sets.a is not None:
        scale = parse_input(sets.a)
        redz = cosmo._a_to_z(scale)

    # Age of the Universe
    elif sets.ta is not None:
        tage = parse_input(sets.ta, U_SEC)
        redz = cosmo.tage_to_z(tage.cgs.value)
        results['ta'] = tage

    # Loockback Time
    elif sets.tl is not None:
        tlbk = parse_input(sets.tl, U_SEC)
        redz = cosmo.tlbk_to_z(tlbk.cgs.value)
        results['tl'] = tlbk

    # Comoving Distance
    elif sets.dc is not None:
        dcom = parse_input(sets.dc, U_CM)
        redz = cosmo.dcom_to_z(dcom.cgs.value)
        # Calculate luminosity-distance manually
        dlum = dcom * (1.0 + redz)
        results['dc'] = dcom
        results['dl'] = dlum

    # Luminosity Distance
    elif sets.dl is not None:
        dlum = parse_input(sets.dl, U_CM)
        redz = cosmo.dlum_to_z(dlum.cgs.value)
        # Calculate comoving-distance manually
        dcom = dlum / (1.0 + redz)
        results['dc'] = dcom
        results['dl'] = dlum

    # No arguments set, raise error
    else:
        print("__main__.calc()")
        print("\t`sets` = '{}'".format(sets))
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
    """Calculate derived cosmological parameters based on the standard ones.

    Adds:
        - 'da': the Angular diameter distance
        - 'arc': the physical size corresponding to one arcsecond
        - 'dm': the distance modulus

    Arguments
    ---------
    results : dict
        Dictionary of standard cosmological parameters (e.g. dl, z, tl, etc)

    Returns
    -------
    results : dict
        Dictionary of parameters with derived quantities added.

    """
    dang = results['dc']/(1.0 + results['z'])
    darc = dang*ARCSEC
    dist_mod = 5.0*np.log10(results['dl'].cgs.value/(10.0*PC))

    results['da'] = dang
    results['arc'] = darc
    results['dm'] = dist_mod

    return results


def output_print(results, print_output=True):
    """Format and print the given cosmological parameters to stdout.

    If the global `COLORED_OUTPUT` parameter is true, `click` is used to print colored output.

    Arguments
    ---------
    results : `Results` namedtuple
        Cosmological parameters to output.
    print_output : bool
        Print output to stdout in addition to returning dictionary of values

    Returns
    -------
    retvals : dict
        Dictionary of cosmological parameter outputs in formatted str form.

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

    retvals = {}
    printvals = []
    outs = []
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
        printvals.append({"base": base, "conv": conv, "name": nn})

    if print_output:
        print("")

        if COLORED_OUTPUT:
            for val in printvals:
                click.secho('{base:30s}'.format(base=val['base']),
                            bold=True, nl=False, fg=_COLOR_BASE)
                click.secho('{conv:20s}'.format(conv=val['conv']),
                            bold=False, nl=False, fg=_COLOR_CONV)
                click.secho(' : ', bold=True, nl=False)  # , fg="green")
                click.secho('{name}'.format(name=val["name"]),
                            bold=False, nl=True, fg=_COLOR_NAME)
        else:
            print("\n".join(outs))

        print("")

    return retvals


def output_api(results):
    """Convert cosmological parameters from numeric values to formatted-strings.

    Arguments
    ---------
    results : `dict`
        Cosmological parameters in numerical form.

    Returns
    -------
    retvals : `dict`
        Cosmological parameters as nicely formatted strings.

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

    retvals = {}
    for kk, ss, tt, nn in zip(keys, symbs, types, names):
        vv = results[kk]
        try:
            vv = vv.item()
        except AttributeError:
            pass
        v_std = vv if tt is None else vv.to(tt)
        try:
            base = "{:.4f}".format(v_std)
        except Exception:
            print("Failed constructing `base` string on '{}' ({})".format(nn, ss))
            print("  Value = '{}', type = '{}'".format(v_std, type(v_std)))
            raise

        retvals[kk] = base

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
        val = float(inval)
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


def parse_args(args):
    """Initialize desired command line options, retrieve their values.
    """

    parser = argparse.ArgumentParser(description="cosmopy: cosmological calculator.")

    # Target Parameters
    parser.add_argument('-z', type=float, default=None,
                        help='target redshift z')
    parser.add_argument('-a', type=float, default=None,
                        help='target scale factor a')
    parser.add_argument('-dc', '-cd', default=None,
                        help='target coming distance D_C')
    parser.add_argument('-dl', '-ld', default=None,
                        help='target luminosity distance D_L')
    parser.add_argument('-tl', '-lt', default=None,
                        help='target look-back time T_L')
    parser.add_argument('-ta', '-at', default=None,
                        help='target universe age T_A')

    parser.add_argument('-v', '--version', action="store_true", default=False,
                        help='print version information.')

    # print("args = ", args)
    sets = parser.parse_args(args)

    if sets.version:
        from . import __version__
        print("cosmopy :: https://github.com/lzkelley/cosmopy/")
        print(os.path.abspath(__file__))
        print(__version__)
        sys.exit(1)

    arg_given = False
    arg_vars = vars(sets)
    arg_vars.pop('version')
    # print(arg_vars)
    for kk, vv in arg_vars.items():
        if vv is not None:
            arg_given = True

    if not arg_given:
        print("No argument given.  Must provide an input argument.\n")
        parser.print_help()
        print("")
        sys.exit(1)

    return sets


def api(key, val, cosmo=None):
    """Primary method for external API access to this package.

    Arguments
    ---------
    key : str
        Specification of target cosmological parameter.  Must be one of `_RESULTS_PARS`.
    val : str or scalar
        The value of the target cosmological parameter.
        Can be a scalar value in CGS units, or a string including unit specification.

    Returns
    -------
    retvals : dict
        Dictionary of cosmological parameters as formatted strings.

    """
    sets = {kk: None for kk in _RESULTS_PARS}
    if key not in sets:
        print("Valid keys: '{}'".format(_RESULTS_PARS))
        raise KeyError("Unrecognized key '{}'".format(key))

    # Set the given value
    sets[key] = val
    # Convert to `Namespace` object (expected by the `calc_basic` method)
    sets = argparse.Namespace(**sets)

    # Load cosmology object if needed
    if cosmo is None:
        cosmo = get_cosmology()

    # Calculate
    # ----------------
    # Calculate basic cosmological values
    results = calc_basic(cosmo, sets)
    # Calculate additional derived values
    results = calc_derived(results)
    # Format output string and return as dictionary
    retvals = output_api(results)

    return retvals


def get_cosmology():
    """Construct an instance of the `Cosmology` class and return.
    """
    cosmo = Cosmology()
    return cosmo


def main(args=None):
    """Primary computational method for command-line access to calculations.

    Loads command line arguments and a `Cosmology` instance.
    Calculates cosmological distance measures.
    Formats and prints the output to stdout.

    """
    # Initialize
    # ----------------
    # Load arguments
    if args is None:
        args = sys.argv[1:]
    sets = parse_args(args)
    # Load cosmology object
    cosmo = get_cosmology()

    # Calculate
    # ----------------
    # Calculate basic cosmological values
    results = calc_basic(cosmo, sets)
    # Calculate additional derived values
    results = calc_derived(results)
    # Format and print outputs
    output_print(results)

    return


if __name__ == "__main__":
    main()
