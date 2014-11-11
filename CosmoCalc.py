# ==================================================================================================
# CosmoCalc.py
# ------------
#
#
# ------------------
# Luke Zoltan Kelley
# LKelley@cfa.harvard.edu
# ==================================================================================================

import os
import sys

import numpy as np
import scipy as sp
import scipy.optimize
import scipy.integrate

from argparse import ArgumentParser


# Global Variables (initialized in 'main')
#sets = None
#pars = None



class Settings(object):
    """
    Object to store the current configuration of CosmoCalc.
    """

    verbose    = True                                                                               # Print excess output to stdout
    print_flag = False                                                                              # Print the cosmological parameters
    build_flag = False                                                                              # Default whether to rebuild integration table

    use_pc     = False                                                                              # Use input distance in parsecs
    use_mpc    = False                                                                              # Use input distance in Megaparsecs
    use_ly     = False                                                                              # Use input distance in lightyears

    use_yr     = False                                                                              # Use input time     in years
    use_myr    = False                                                                              # Use input time     in megayears

    # Target epochs (calculate epoch based on given parameter below)
    z   = None                                                                               # Redshift (z)
    a   = None                                                                               # Scale Factor (a)
    cd  = None                                                                               # Comoving Distance
    ld  = None                                                                               # Luminosity Distance
    lt  = None                                                                               # Lookback time
    t   = None                                                                               # Age of the universe (time)

    root_tol   = 1.0e-10                                                                            # Root finding tolerance
    quad_iter  = 1000                                                                               # Number of iterations for quadrature





class Parameters(object):
    """
    Standard Cosmological Parameters.

    Cosmological Parameters taken from
    [1] '1212.5226v3 - Nine-Year WMAP Observations: Cosmological Parameter Results'
    [2] '1002.3488   - The Cosmological Parameters 2010'
    """

    # Default Cosmological Parameters
    H0         = 69.7                      # Hubble's Constant           [km/s/Mpc]
    T0         = 13.76                     # Age of the Universe         [Gyr]
    OmegaDM    = 0.235                     # Dark Matter Density         [1]
    OmegaB     = 0.0464                    # Baryonic Density            [1]
    OmegaL     = 0.7185                    # Dark Energy Density         [1]
    OmegaR     = 4.984e-4                  # Radiation Density           [2] omegar*h^2 = 2.47e-5; h = 0.704

    # Physical Constants
    c          = 2.99792458e10             # Speed of Light              [cm/s]
    G          = 6.67259e-8                # Gravitational Constant      [cm^3/g/s^2]
    mProton    = 1.6726231e-24             # Proton Mass                 [g]
    mElectron  = 9.1093829e-28             # Electron Mass               [g]
    planckH    = 6.6260755e-27             # Planck's Constant           [erg s]
    parsec     = 3.08567758e+18            # Parsec                      [cm]
    arcsec     = 4.84813681e-06            # Arcsecond                   [radians]
    year       = 3.15569520e+07            # Year                        [s]

    @staticmethod
    def HubbleTime():      return 1.0*(1.0e6*Parameters.parsec/1.0e5)/Parameters.H0

    @staticmethod
    def HubbleDistance():  return Parameters.c*Parameters.HubbleTime()

    @staticmethod
    def CriticalDensity(): return 3*np.square(Parameters.H0)/(8.*np.pi*Parameters.G)

    @staticmethod
    def OmegaMatter():     return Parameters.OmegaB+Parameters.OmegaDM

    @staticmethod
    def ParameterString():
        parstr  = "\n"
        parstr += "Hubble Constant,           H0  = %f\n" % (Parameters.H0)
        parstr += "Age of the Universe,       T0  = %f\n" % (Parameters.T0)
        parstr += "Baryon      Fraction,  OmegaB  = %f\n" % (Parameters.OmegaB )
        parstr += "Dark Matter Fraction,  OmegaDM = %f\n" % (Parameters.OmegaDM)
        parstr += "Dark Energy Fraction,  OmegaL  = %f\n" % (Parameters.OmegaL )
        parstr += "Radiation   Fraction,  OmegaR  = %f\n" % (Parameters.OmegaR )
        return parstr

    @staticmethod
    def HubbleFunction(zz): 
        """
        The E(z) function from Hogg1999.

        Incorporates matter (dark and light), radiation, and dark energy.  This
        version of the equation assumes that OmegaK (curvature) is zero.
        """

        # Determine the constituents at the given redshifts
        t_matter     = Parameters.OmegaMatter()*np.power((1.0+zz),3)
        t_radiation  = Parameters.OmegaR*np.power((1.0+zz),4)
        t_darkenergy = Parameters.OmegaL
        return np.sqrt( t_matter + t_radiation + t_darkenergy )


    @staticmethod
    def HubbleDistanceFunction(zz): return 1.0/Parameters.HubbleFunction(zz)
    @staticmethod
    def HubbleTimeFunction(zz): return 1.0/( (1.0+zz)*Parameters.HubbleFunction(zz) )








def DistToCM(indist):
    """
    Based on unit settings, convert distance to centimeters.

    Setup to accept [default] centimeters (cm), parsecs (pc), megaparsecs (mpc),
    or lightyears (ly).  Stored in Settings object (tst).
    """
    if(   Settings.use_pc  ): indist *= Parameters.parsec
    elif( Settings.use_mpc ): indist *= Parameters.parsec*(1.0e6)
    elif( Settings.use_ly  ): indist *= Parameters.year*Parameters.c
    return indist


def TimeToS(intime):
    """
    Based on unit settings, convert time to seconds.

    Setup to accept [default] seconds (s), years (yr), megayears (myr).  Stored
    in Settings object (tst).
    """
    if(   Settings.use_yr  ): intime *= Parameters.year
    elif( Settings.use_myr ): intime *= Parameters.year*(1.0e6)
    return intime






def Integrate(func, z0, z1):
    """Integrate a general function (func) from 'z0' to 'z1'."""
    val, err = sp.integrate.quadrature(func, z0, z1, maxiter=Settings.quad_iter)
    return val, err

def IntegrateDistance(z0, z1):
    """Integrate distance between given redshifts."""
    dist, derr = Integrate(Parameters.HubbleDistanceFunction, z0, z1)
    return dist, derr

def IntegrateTime(z0, z1):
    """Integrate time between given redshifts"""
    time, terr = Integrate(Parameters.HubbleTimeFunction, z0, z1)
    return time, terr



def RootFind(func): return sp.optimize.newton( func, 1.0, tol=Settings.root_tol )

def InvertScaleFactor(ascale):
    """Find redshift from Scale-Factor"""
    return (1.0/ascale) - 1.0

def InvertComDist(cdist):
    """Given a comoving-distance, find redshift."""
    rhs = lambda x: IntegrateDistance(0.0, x)[0] - cdist
    return RootFind(rhs)

def InvertLumDist(ldist):
    """Given a luminosity-distance, find redshift."""
    rhs = lambda x: (1.0+x)*IntegrateDistance(0.0, x)[0] - ldist
    return RootFind(rhs)

def InvertLookTime(ltime):
    """Given a look-back time, find redshift."""
    rhs = lambda x: IntegrateTime(0.0, x)[0] - ltime
    return RootFind(rhs)

def InvertAgeTime(atime):
    """Given an age of the universe, find redshift."""
    rhs = lambda x: 1.0 - IntegrateTime(0.0, x)[0] - atime
    return RootFind(rhs)


    

def CosmologicalParameters(redz, pout=False):
    """
    Given a redshift and fractional comoving distance (D_C/D_H), return cosmological parameters.
    """

    pars = []
    errs = []
    retstr = ""

    # Store Redshift
    pars.append(redz)
    errs.append(0.0)
    rstr = "z = %e" % (redz)
    retstr += rstr + "\n"
    #strs.append(rstr)

    # Store scale-factor
    scaleFactor = 1.0/(1.0+redz)
    pars.append(scaleFactor)
    errs.append(0.0)
    rstr = "a = %e" % (scaleFactor)
    retstr += rstr + "\n"
    #strs.append(rstr)

    # Integrate over redshift to find distance and time with errors
    dist, dist_err = IntegrateDistance(0.0, redz)
    time, time_err = IntegrateTime(0.0, redz)

    # Comoving Distance
    comDist       = Parameters.HubbleDistance()*dist
    comDist_err   = Parameters.HubbleDistance()*dist_err
    pars.append(comDist)
    errs.append(comDist_err)
    rstr          = "D_C  =  %+e +- %e  [cm]  ~  %+e [Mpc] : Comoving Distance" % \
                    (comDist, comDist_err, comDist/(1.0e6*Parameters.parsec) )
    retstr       += rstr + "\n"
    #strs.append(rstr)

    # Luminosity Distance
    lumDist       = (1.0 + redz)*comDist
    lumDist_err   = (1.0 + redz)*comDist_err
    pars.append(lumDist)
    errs.append(lumDist_err)
    rstr          = "D_L  =  %+e +- %e  [cm]  ~  %+e [Mpc] : Luminosity Distance" % \
                    (lumDist, lumDist_err, lumDist/(1.0e6*Parameters.parsec) )
    retstr       += rstr + "\n"
    #strs.append(rstr)

    # Angular Diameter Distance
    angDist       = comDist/(1.0 + redz)
    angDist_err   = comDist_err/(1.0 + redz)
    pars.append(angDist)
    errs.append(angDist_err)
    rstr          = "D_A  =  %+e +- %e  [cm]  ~  %+e [Mpc] : Angular Diameter Distance" % \
                    (angDist, angDist_err, lumDist/(1.0e6*Parameters.parsec)) 
    retstr       += rstr + "\n"
    #strs.append(rstr)

    # Arcsecond size
    arcSize       = Parameters.arcsec*angDist
    arcSize_err   = Parameters.arcsec*angDist_err
    pars.append(arcSize)
    errs.append(arcSize_err)
    rstr          = "Arc  =  %+e +- %e  [cm]  ~  %+e [pc]  : Arcsecond Transverse Distance" % \
                    (arcSize, arcSize_err, arcSize/(Parameters.parsec) )
    retstr       += rstr + "\n"
    #strs.append(rstr)

    # Lookback Time
    lookTime      = Parameters.HubbleTime()*time
    lookTime_err  = Parameters.HubbleTime()*time_err
    pars.append(lookTime)
    errs.append(lookTime_err)
    rstr          = "T_L  =  %+e +- %e  [s]   ~  %+e [Myr] : Lookback Time" % \
                    (lookTime, lookTime_err, lookTime/(1.0e6*Parameters.year) )
    retstr       += rstr + "\n"
    #strs.append(rstr)

    # Age
    ageTime       = Parameters.HubbleTime()*(1.0-time)
    ageTime_err   = lookTime_err
    pars.append(ageTime)
    errs.append(ageTime_err)
    rstr          = "T_A  =  %+e +- %e  [s]   ~  %+e [Myr] : Age of the Universe" % \
                    (ageTime, ageTime_err, ageTime/(1.0e6*Parameters.year) )
    retstr       += rstr + "\n"
    #strs.append(rstr)

    # Distance Modulus
    distMod       = 5.0*np.log10(lumDist/(Parameters.parsec*10.0))
    distMod_err   = 5.0*np.log10(lumDist_err/(Parameters.parsec*10.0)) if lumDist_err > 0.0 else 0.0
    pars.append(distMod)
    errs.append(distMod_err)
    rstr          = "DM   =  %+e +- %e  []    Distance Modulus" % (distMod, np.abs(distMod_err) )
    retstr       += rstr + "\n"
    #strs.append(rstr)

    if( pout ): print retstr

    return pars, errs

    


def solve(z=None, a=None, cd=None, ld=None, lt=None, t=None, verbose=False, pout=False):

    # Redshift
    if( z != None ):
        if( verbose ): print " - Redshift"
        redz = z

    # Scale-factor
    elif( a != None ):
        if( verbose ): print " - Scale-factor"
        redz = InvertScaleFactor(a)

    # Comoving Distance
    elif( cd != None ):
        t_cd = DistToCM(cd)
        if( verbose ): print " - Comoving Distance, %.2e [cm]" % (t_cd)
        redz = InvertComDist( t_cd/Parameters.HubbleDistance() )

    # Luminosity Distance
    elif( ld != None ):
        t_ld = DistToCM(ld)
        if( verbose ): print " - Luminosity Distance, %.2e [cm]" % (t_ld)
        redz = InvertLumDist( t_ld/Parameters.HubbleDistance() )

    # Look-back Time
    elif( lt != None ):
        t_tl = TimeToS(tl)
        if( verbose ): print " - Look-back Time, %.2e [s]" % (t_tl)
        redz = InvertLookTime( t_tl/Parameters.HubbleTime() )

    # Universe Age Time
    elif( t != None ):
        t_ta = TimeToS(t)
        if( verbose ): print " - Age of the Universe, %.2e [s]" % (t_ta)
        redz = InvertAgeTime( t_ta/Parameters.HubbleTime() )

    else:
        print Parameters.ParameterString()
        raise RuntimeError("No parameters given!")


    # Make sure redz is an iterable
    if( not np.iterable(redz) ): redz = [ redz ]


    pars = []
    errs = []
    for zz in redz:
        pp, ee = CosmologicalParameters(zz, pout=pout)    
        pars.append(pp)
        errs.append(ee)

    pars = np.array(pars)
    errs = np.array(errs)

    return pars.T, errs.T



def ParseArgs():
    """
    Initialize the argparse object with desired command line options, then retrieve their values.
    """

    parser          = ArgumentParser()

    # Target Parameters
    parser.add_argument('-z',         type=float, default=Settings.z,       help='Target redshift z'                                        )
    parser.add_argument('-a',         type=float, default=Settings.a,       help='Target scale factor a'                                    )
    parser.add_argument('-cd','-dc',  type=float, default=Settings.cd,      help='Target coming distance D_C'                               )
    parser.add_argument('-ld','-dl',  type=float, default=Settings.ld,      help='Target luminosity distance D_L'                           )
    parser.add_argument('-lt','-tl',  type=float, default=Settings.lt,      help='Target look-back time T_L'                                )
    parser.add_argument('-ta','-at',  type=float, default=Settings.t,       help='Target universe age T_A'                                  )

    # Modifiers
    parser.add_argument('--pc',                 action='store_true', default=Settings.use_pc , help="Use provided distance ('-cd' or '-ld') in parsecs"     )
    parser.add_argument('--mpc',                action='store_true', default=Settings.use_mpc, help="Use provided distance ('-cd' or '-ld') in megaparsecs" )
    parser.add_argument('--ly',                 action='store_true', default=Settings.use_ly , help="Use provided distance ('-cd' or '-ld') in lightyears"  )
    parser.add_argument('--yr',                 action='store_true', default=Settings.use_yr , help="Use provided time ('-lt' or '-at') in years"           )
    parser.add_argument('--myr',                action='store_true', default=Settings.use_myr, help="Use provided time ('-lt' or '-at') in megayears"       )

    # Behavior
    parser.add_argument('--print', dest='prt',  action='store_true', default=Settings.print_flag, help="Print defaul cosmological parameters"                  )

    # Modify Cosmological Parameters
    parser.add_argument('--H0','--h0',          type=float, metavar='',        default=Parameters.H0,           help='Hubble Constant (H0) [km/s/Mpc]'                         )
    parser.add_argument('--T0','--t0',          type=float, metavar='',        default=Parameters.T0,           help='Hubble Time (T0) [1/s]'                                  )
    parser.add_argument('--ODM','--darkmatter', type=float, metavar='OmegaDM', default=Parameters.OmegaDM,      help='Dark Matter fraction (OmegaDM) [1/critical-density]'     )
    parser.add_argument('--OB', '--baryon',     type=float, metavar='OmegaB',  default=Parameters.OmegaB,       help='Baryon      fraction (OmegaB)  [1/critical-density]'     )
    parser.add_argument('--OL','--darkenergy',  type=float, metavar='OmegaL',  default=Parameters.OmegaL,       help='Dark Energy fraction (OmegaL)  [1/critical-density]'     )
    parser.add_argument('--OR','--radiation',   type=float, metavar='OmegaR',  default=Parameters.OmegaR,       help='Radtion     fraction (OmegaR)  [1/critical-density]'     )


    args            = parser.parse_args()
    Settings.print_flag = args.prt

    if( args.z  != None ): Settings.z   = float(args.z)
    if( args.a  != None ): Settings.a   = float(args.a)
    if( args.cd != None ): Settings.cd  = float(args.cd)
    if( args.ld != None ): Settings.ld  = float(args.ld)
    if( args.lt != None ): Settings.lt  = float(args.lt)
    if( args.ta != None ): Settings.t   = float(args.ta)

    Settings.use_pc     = args.pc
    Settings.use_mpc    = args.mpc
    Settings.use_ly     = args.ly
    Settings.use_yr     = args.yr
    Settings.use_myr    = args.myr

    Parameters.H0         = args.H0
    Parameters.T0         = args.T0
    Parameters.OmegaDM    = args.ODM
    Parameters.OmegaB     = args.OB
    Parameters.OmegaL     = args.OL
    Parameters.OmegaR     = args.OR

    return sets, pars



def main():

    print "\nCosmoCalc"


    global sets, pars
    sets         = Settings()                                                                       # Create a Settings object
    pars         = Parameters()                                                                     # Create a Cosmological Parameters Object
    sets, pars   = ParseArgs()                                                                      # Modify settings and parameters based on user arguments

    # Print Cosmological parameters
    if( Settings.print_flag ): print Parameters.ParameterString()

    pars, errs = solve(z=sets.z, a=sets.a, cd=sets.cd, ld=sets.ld, lt=sets.lt, t=sets.t, pout=True)




if __name__ == "__main__": main()

