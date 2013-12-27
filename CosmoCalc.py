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
sets = None
pars = None


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
    z_target   = -1.0                                                                               # Redshift
    cd_target  = -1.0                                                                               # Comoving Distance
    ld_target  = -1.0                                                                               # Luminosity Distance
    tl_target  = -1.0                                                                               # Lookback time
    ta_target  = -1.0                                                                               # Age of the universe (time)

    root_tol   = 1.0e-10                                                                            # Root finding tolerance
    quad_iter  = 1000                                                                               # Number of iterations for quadrature

    #table_filename      = ".cc-table.dat"                                                           # Name of integration table file
    #script_path         = os.path.dirname( os.path.abspath(__file__) )                              # Path to this python script
    #full_table_filename = script_path + "/" + table_filename




class Parameters(object):
    """
    Standard Cosmological Parameters.

    Cosmological Parameters taken from
    [1] '1212.5226v3 - Nine-Year WMAP Observations: Cosmological Parameter Results'
    [2] '1002.3488   - The Cosmological Parameters 2010'
    """

    H0         = 69.7                      # Hubble's Constant           [km/s/Mpc]
    T0         = 13.76                     # Age of the Universe         [Gyr]
    
    OmegaDM    = 0.235                     # Dark Matter Density         [1]
    OmegaB     = 0.0464                    # Baryonic Density            [1]
    OmegaL     = 0.7185                    # Dark Energy Density         [1]
    OmegaR     = 4.984e-4                  # Radiation Density           [2] omegar*h^2 = 2.47e-5; h = 0.704

    c          = 2.99792458e10             # Speed of Light              [cm/s]
    G          = 6.67259e-8                # Gravitational Constant      [cm^3/g/s^2]
    mProton    = 1.6726231e-24             # Proton Mass                 [g]
    mElectron  = 9.1093829e-28             # Electron Mass               [g]
    planckH    = 6.6260755e-27             # Planck's Constant           [erg s]
    parsec     = 3.08567758e+18            # Parsec                      [cm]
    arcsec     = 4.84813681e-06            # Arcsecond                   [radians]
    year       = 3.15569520e+07            # Year                        [s]

    @staticmethod
    def HubbleTime():      return 1.0*(1.0e6*pars.parsec/1.0e5)/pars.H0
    @staticmethod
    def HubbleDistance():  return pars.c*pars.HubbleTime()
    @staticmethod
    def CriticalDensity(): return 3*np.square(pars.H0)/(8.*np.pi*pars.G)
    @staticmethod
    def OmegaMatter():     return pars.OmegaB+pars.OmegaDM

    @staticmethod
    def ParameterString():
        parstr  = "\n"
        parstr += "Hubble Constant,           H0  = %f\n" % (pars.H0)
        parstr += "Age of the Universe,       T0  = %f\n" % (pars.T0)
        parstr += "Baryon      Fraction,  OmegaB  = %f\n" % (pars.OmegaB )
        parstr += "Dark Matter Fraction,  OmegaDM = %f\n" % (pars.OmegaDM)
        parstr += "Dark Energy Fraction,  OmegaL  = %f\n" % (pars.OmegaL )
        parstr += "Radiation   Fraction,  OmegaR  = %f\n" % (pars.OmegaR )

        return parstr

    @staticmethod
    def HubbleFunction(zz): 
        """
        The E(z) function from Hogg1999.

        Incorporates matter (dark and light), radiation, and dark energy.  This
        version of the equation assumes that OmegaK (curvature) is zero.
        """
        return np.sqrt( pars.OmegaMatter()*np.power((1.0+zz),3) + pars.OmegaR*np.power((1.0+zz),4) + pars.OmegaL )

    @staticmethod
    def HubbleDistanceFunction(zz):
        return 1.0/pars.HubbleFunction(zz)

    @staticmethod
    def HubbleTimeFunction(zz):
        return 1.0/( (1.0+zz)*pars.HubbleFunction(zz) )






def DistToCM(tst, indist):
    """
    Based on unit settings, convert distance to centimeters.

    Setup to accept [default] centimeters (cm), parsecs (pc), megaparsecs (mpc),
    or lightyears (ly).  Stored in Settings object (tst).
    """
    if(   tst.use_pc  ): indist *= pars.parsec
    elif( tst.use_mpc ): indist *= pars.parsec*(1.0e6)
    elif( tst.use_ly  ): indist *= pars.year*pars.c
    return indist


def TimeToS(tst, intime):
    """
    Based on unit settings, convert time to seconds.

    Setup to accept [default] seconds (s), years (yr), megayears (myr).  Stored
    in Settings object (tst).
    """
    if(   tst.use_yr  ): intime *= pars.year
    elif( tst.use_myr ): intime *= pars.year*(1.0e6)
    return intime






def Integrate(func, z0, z1):
    """Integrate a general function (func) from 'z0' to 'z1'."""
    val, err = sp.integrate.quadrature(func, z0, z1, maxiter=sets.quad_iter)
    return val, err

def IntegrateDistance(z0, z1):
    """Integrate distance between given redshifts."""
    dist, derr = Integrate(pars.HubbleDistanceFunction, z0, z1)
    return dist, derr

def IntegrateTime(z0, z1):
    """Integrate time between given redshifts"""
    time, terr = Integrate(pars.HubbleTimeFunction, z0, z1)
    return time, terr



def RootFind(func): return sp.optimize.newton( func, 1.0, tol=sets.root_tol )

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
    
    dist, dist_err = IntegrateDistance(0.0, redz)
    time, time_err = IntegrateTime(0.0, redz)

    retstr = "z = %e\n" % (redz)

    # Comoving Distance
    comDist       = pars.HubbleDistance()*dist
    comDist_err   = pars.HubbleDistance()*dist_err
    retstr       += "D_C  =  %+e +- %e  [cm]  ~  %+e [Mpc] : Comoving Distance\n" % (comDist, comDist_err, comDist/(1.0e6*pars.parsec) )

    # Luminosity Distance
    lumDist       = (1.0 + redz)*comDist
    lumDist_err   = (1.0 + redz)*comDist_err
    retstr       += "D_L  =  %+e +- %e  [cm]  ~  %+e [Mpc] : Luminosity Distance\n" % (lumDist, lumDist_err, lumDist/(1.0e6*pars.parsec) )

    # Angular Diameter Distance
    angDist       = comDist/(1.0 + redz)
    angDist_err   = comDist_err/(1.0 + redz)
    retstr       += "D_A  =  %+e +- %e  [cm]  ~  %+e [Mpc] : Angular Diameter Distance\n" % (angDist, angDist_err, lumDist/(1.0e6*pars.parsec)) 

    # Arcsecond size
    arcSize       = pars.arcsec*angDist
    arcSize_err   = pars.arcsec*angDist_err
    retstr       += "Arc  =  %+e +- %e  [cm]  ~  %+e [pc]  : Arcsecond Transverse Distance\n" % (arcSize, arcSize_err, arcSize/(pars.parsec) )

    # Lookback Time
    lookTime      = pars.HubbleTime()*time
    lookTime_err  = pars.HubbleTime()*time_err
    retstr       += "T_L  =  %+e +- %e  [s]   ~  %+e [Myr] : Lookback Time\n" % (lookTime, lookTime_err, lookTime/(1.0e6*pars.year) )

    # Age
    ageTime       = pars.HubbleTime()*(1.0-time)
    ageTime_err   = lookTime_err
    retstr       += "T_A  =  %+e +- %e  [s]   ~  %+e [Myr] : Age of the Universe\n" % (ageTime, ageTime_err, ageTime/(1.0e6*pars.year) )

    # Distance Modulus
    distMod       = 5.0*np.log10(lumDist/(pars.parsec*10.0))
    distMod_err   = 5.0*np.log10(lumDist_err/(pars.parsec*10.0))
    retstr       += "DM   =  %+e +- %e  []    Distance Modulus\n" % (distMod, np.abs(distMod_err) )

    if( pout ): print retstr

    return [ [comDist , comDist_err , "Comoving Distance [cm]"        ] , \
             [lumDist , lumDist_err , "Luminosity Distance [cm]"      ] , \
             [angDist , angDist_err , "Angular Diameter Distance [cm]"] , \
             [arcSize , arcSize_err , "Arcsecond Transverse Size [cm]"] , \
             [lookTime, lookTime_err, "Lookback Time [s]"             ] , \
             [ageTime , ageTime_err , "Age of the Universe [s]"       ] , \
             [distMod , distMod_err , "Distance Modulus []"           ] ]




def RedshiftFromTarget(t_sets):
    """
    Given a settings object, determine which parameter is being targeted.

    The settings object contains the default, and command-line options.  Use this information
    to determine which cosmological parameter ('target') is being provided (e.g. redshift),
    and return that target parameter, along with the appropriate function to obtain that
    parameter from a line from the table (e.g. 'RedshiftFromLine').
    """

    all_targets = [ t_sets.z_target, t_sets.cd_target, t_sets.ld_target, t_sets.tl_target, t_sets.ta_target ]

    # Make sure only one target is specified (>= 0.0)
    if( sum(tarz >= 0.0 for tarz in all_targets) != 1 ): 
        raise RuntimeError('Must specify a single target parameter.')


    # Redshift
    if( t_sets.z_target >= 0.0 ): 
        if( t_sets.verbose ): print " - Redshift"
        return t_sets.z_target


    # Comoving Distance
    if( t_sets.cd_target >= 0.0 ):
        t_cd = DistToCM(t_sets, t_sets.cd_target)
        if( t_sets.verbose ): print " - Comoving Distance, %.2e [cm]" % (t_cd)
        return InvertComDist( t_cd/pars.HubbleDistance() )


    # Luminosity Distance
    if( t_sets.ld_target >= 0.0 ):
        t_ld = DistToCM(t_sets, t_sets.ld_target)
        if( t_sets.verbose ): print " - Luminosity Distance, %.2e [cm]" % (t_ld)
        return InvertLumDist( t_ld/pars.HubbleDistance() )


    # Look-back Time
    if( t_sets.tl_target >= 0.0 ):
        t_tl = TimeToS(t_sets, t_sets.tl_target)
        if( t_sets.verbose ): print " - Look-back Time, %.2e [s]" % (t_tl)
        return InvertLookTime( t_tl/pars.HubbleTime() )


    # Universe Age Time
    if( t_sets.ta_target >= 0.0 ):
        t_ta = TimeToS(t_sets, t_sets.ta_target)
        if( t_sets.verbose ): print " - Age of the Universe, %.2e [s]" % (t_ta)
        return InvertAgeTime( t_ta/pars.HubbleTime() )


    raise RuntimeError('Must specify a single target parameter.')

    


def InitArgParse():
    """
    Initialize the argparse object with desired command line options, then retrieve their values.
    """

    parser          = ArgumentParser()

    parser.add_argument('-z',         type=float, default=sets.z_target,       help='Target redshift z'                                        )
    parser.add_argument('-cd','-dc',  type=float, default=sets.cd_target,      help='Target coming distance D_C'                               )
    parser.add_argument('-ld','-dl',  type=float, default=sets.ld_target,      help='Target luminosity distance D_C'                           )
    parser.add_argument('-tl','-lt',  type=float, default=sets.tl_target,      help='Target look-back time T_L'                                )
    parser.add_argument('-ta','-at',  type=float, default=sets.ta_target,      help='Target universe age T_A'                                  )

    parser.add_argument('--print', dest='prt',  action='store_true', default=sets.print_flag, help="Print defaul cosmological parameters"                  )
    #parser.add_argument('--plot',               action='store_true', default=False, help="Plot cosmological parameters"                          )

    parser.add_argument('--pc',                 action='store_true', default=sets.use_pc , help="Use provided distance ('-cd' or '-ld') in parsecs"     )
    parser.add_argument('--mpc',                action='store_true', default=sets.use_mpc, help="Use provided distance ('-cd' or '-ld') in megaparsecs" )
    parser.add_argument('--ly',                 action='store_true', default=sets.use_ly , help="Use provided distance ('-cd' or '-ld') in lightyears"  )

    parser.add_argument('--yr',                 action='store_true', default=sets.use_yr , help="Use provided time ('-lt' or '-at') in years"           )
    parser.add_argument('--myr',                action='store_true', default=sets.use_myr, help="Use provided time ('-lt' or '-at') in megayears"       )


    args            = parser.parse_args()
    sets.print_flag = args.prt
    #sets.plot_flag  = args.plot

    sets.z_target   = float(args.z)
    sets.cd_target  = float(args.cd)
    sets.ld_target  = float(args.ld)
    sets.tl_target  = float(args.tl)
    sets.ta_target  = float(args.ta)

    sets.use_pc     = args.pc
    sets.use_mpc    = args.mpc
    sets.use_ly     = args.ly
    sets.use_yr     = args.yr
    sets.use_myr    = args.myr


def main():

    print "\nCosmoCalc"

    global sets, pars
    sets         = Settings()
    pars         = Parameters()
    args         = InitArgParse()

    # Print Cosmological parameters
    if( sets.print_flag ):
        print pars.ParameterString()
        exit(1)

    redshift     = RedshiftFromTarget(sets)
    cosmo_params = CosmologicalParameters(redshift, pout=True)




if __name__ == "__main__": main()





# 
