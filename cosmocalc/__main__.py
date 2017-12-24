"""
"""

from argparse import ArgumentParser

from . import cosmocalc
from . settings import Settings
from . parameters import Parameters


def ParseArgs():
    """
    Initialize the argparse object with desired command line options, then retrieve their values.
    """

    parser = ArgumentParser()

    # Target Parameters
    parser.add_argument('-z', type=float, default=Settings.z, help='Target redshift z')
    parser.add_argument('-a', type=float, default=Settings.a, help='Target scale factor a')
    parser.add_argument('-cd', '-dc', type=float, default=Settings.cd, help='Target coming distance D_C')
    parser.add_argument('-ld', '-dl', type=float, default=Settings.ld, help='Target luminosity distance D_L')
    parser.add_argument('-lt', '-tl', type=float, default=Settings.lt, help='Target look-back time T_L')
    parser.add_argument('-ta', '-at', type=float, default=Settings.t, help='Target universe age T_A')

    # Modifiers
    parser.add_argument('--pc', action='store_true', default=Settings.use_pc, help="Use provided distance ('-cd' or '-ld') in parsecs")
    parser.add_argument('--mpc', action='store_true', default=Settings.use_mpc, help="Use provided distance ('-cd' or '-ld') in megaparsecs")
    parser.add_argument('--ly', action='store_true', default=Settings.use_ly, help="Use provided distance ('-cd' or '-ld') in lightyears")
    parser.add_argument('--yr', action='store_true', default=Settings.use_yr, help="Use provided time ('-lt' or '-at') in years")
    parser.add_argument('--myr', action='store_true', default=Settings.use_myr, help="Use provided time ('-lt' or '-at') in megayears")

    # Behavior
    parser.add_argument('--print', dest='prt', action='store_true', default=Settings.print_flag, help="Print defaul cosmological parameters")

    # Modify Cosmological Parameters
    parser.add_argument('--H0', '--h0', type=float, metavar='', default=Parameters.H0, help='Hubble Constant (H0) [km/s/Mpc]')
    parser.add_argument('--T0', '--t0', type=float, metavar='', default=Parameters.T0, help='Hubble Time (T0) [1/s]')
    parser.add_argument('--ODM', '--darkmatter', type=float, metavar='OmegaDM', default=Parameters.OmegaDM, help='Dark Matter fraction (OmegaDM) [1/critical-density]')
    parser.add_argument('--OB', '--baryon', type=float, metavar='OmegaB', default=Parameters.OmegaB, help='Baryon      fraction (OmegaB)  [1/critical-density]')
    parser.add_argument('--OL', '--darkenergy', type=float, metavar='OmegaL', default=Parameters.OmegaL, help='Dark Energy fraction (OmegaL)  [1/critical-density]')
    parser.add_argument('--OR', '--radiation', type=float, metavar='OmegaR', default=Parameters.OmegaR, help='Radtion     fraction (OmegaR)  [1/critical-density]')

    args = parser.parse_args()
    Settings.print_flag = args.prt

    if args.z is not None:
        Settings.z = float(args.z)
    if args.a is not None:
        Settings.a = float(args.a)
    if args.cd is not None:
        Settings.cd = float(args.cd)
    if args.ld is not None:
        Settings.ld = float(args.ld)
    if args.lt is not None:
        Settings.lt = float(args.lt)
    if args.ta is not None:
        Settings.t = float(args.ta)

    Settings.use_pc = args.pc
    Settings.use_mpc = args.mpc
    Settings.use_ly = args.ly
    Settings.use_yr = args.yr
    Settings.use_myr = args.myr

    Parameters.H0 = args.H0
    Parameters.T0 = args.T0
    Parameters.OmegaDM = args.ODM
    Parameters.OmegaB = args.OB
    Parameters.OmegaL = args.OL
    Parameters.OmegaR = args.OR

    return sets, pars


def main():

    print("\nCosmoCalc")

    global sets, pars
    # Create a Settings object
    sets = Settings()
    # Create a Cosmological Parameters Object
    pars = Parameters()
    # Modify settings and parameters based on user arguments
    sets, pars = ParseArgs()

    # Print Cosmological parameters
    if Settings.print_flag:
        print(Parameters.ParameterString())

    pars, errs = cosmocalc.solve(
        z=sets.z, a=sets.a, cd=sets.cd, ld=sets.ld, lt=sets.lt, t=sets.t, pout=True)


if __name__ == "__main__":
    main()
