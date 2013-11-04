# ==================================================================================================
# CosmoCalc.py
# ------------
#
#
# ------------------
# Luke Zoltan Kelley
# LKelley@cfa.harvard.edu
# ==================================================================================================


#from astropy.io import ascii as at
#from astropy.utils.console import ProgressBar
#from linecache import getline as lc_getline

import numpy as np
#from scipy.interpolate import interp1d
import scipy as sp
import scipy.optimize
import scipy.integrate

import os
import sys
from argparse import ArgumentParser


class Sets(object):
    """
    Object to store the current configuration of CosmoCalc.
    """

    verbose    = True                                                                               # Print excess output to stdout
    print_flag = False                                                                              # Print the cosmological parameters
    build_flag = False                                                                              # Default whether to rebuild integration table
    use_pc     = False                                                                              # Use input distance in parsecs
    use_mpc    = False                                                                              # Use input distance in Megaparsecs
    use_ly     = False                                                                              # Use input distance in lightyears

    z_target   = -1.0
    cd_target  = -1.0
    ld_target  = -1.0

    root_tol   = 1.0e-10
    quad_iter  = 1000

    table_filename      = ".cc-table.dat"                                                           # Name of integration table file
    script_path         = os.path.dirname( os.path.abspath(__file__) )                              # Path to this python script
    full_table_filename = script_path + "/" + table_filename




class Pars(object):
    """
    Standard Cosmological Parameters.

    Cosmological Parameters taken from
    '1212.5226v3 - Nine-Year WMAP Observations: Cosmological Parameter Results'
    """

    H0         = 69.7                      # Hubble's Constant           [km/s/Mpc]
    T0         = 13.76                     # Age of the Universe         [Gyr]
    
    OmegaDM    = 0.235                     # Dark Matter Density         [1/critical-density]
    OmegaB     = 0.0464                    # Baryonic Density            [1/critical-density]
    OmegaL     = 0.7185                    # Dark Energy Density         [1/critical-density]

    c          = 2.99792458e10             # Speed of Light              [cm/s]
    G          = 6.67259e-8                # Gravitational Constant      [cm^3/g/s^2]
    mProton    = 1.6726231e-24             # Proton Mass                 [g]
    mElectron  = 9.1093829e-28             # Electron Mass               [g]
    planckH    = 6.6260755e-27             # Planck's Constant           [erg s]
    parsec     = 3.08567758e+18            # Parsec                      [cm]
    arcsec     = 4.84813681e-06            # Arcsecond                   [radians]
    year       = 3.15569520e+07            # Year                        [s]

    @staticmethod
    def HubbleTime():      return 1.0*(1.0e6*Pars.parsec/1.0e5)/Pars.H0
    @staticmethod
    def HubbleDistance():  return Pars.c*Pars.HubbleTime()
    @staticmethod
    def CriticalDensity(): return 3*np.square(Pars.H0)/(8.*np.pi*Pars.G)
    @staticmethod
    def OmegaMatter():     return Pars.OmegaB+Pars.OmegaDM

    @staticmethod
    def ParameterString():
        parstr  = "\n"
        parstr += "Hubble Constant,           H0  = %e\n" % (Pars.H0)
        parstr += "Age of the Universe,       T0  = %e\n" % (Pars.T0)
        parstr += "Baryon Fraction,       OmegaB  = %e\n" % (Pars.OmegaB )
        parstr += "Dark Matter Fraction,  OmegaDM = %e\n" % (Pars.OmegaDM)
        parstr += "Dark Energy Fraction,  OmegaL  = %e\n" % (Pars.OmegaL )

        return parstr

    @staticmethod
    def HubbleFunction(zz): 
        """
        The E(z) function from Hogg1999.

        This version of the equation assumes that OmegaK (curvature) is zero.
        """
        return np.sqrt( Pars.OmegaMatter()*np.power((1+zz),3) + Pars.OmegaL )

    @staticmethod
    def HubbleDistanceFunction(zz):
        return 1.0/Pars.HubbleFunction(zz)

    @staticmethod
    def HubbleTimeFunction(zz):
        return 1.0/( (1.0+zz)*Pars.HubbleFunction(zz) )




def ERROR(func=None, errstr=None, errnum=None, ex=False):
    """Print error message."""
    where_str = "CosmoCalc.py"
    error_str = "ERROR"
    if( func   != None ): where_str += " - %s" % (func  )
    if( errnum != None ): error_str += " %d -" % (errnum)
    if( errstr != None ): error_str += " %s"   % (errstr)

    print '\n[%s] %s !!\n\n' % (where_str, error_str)

    if( ex ):
        sys.exit(errnum)
        exit(errnum)

    return



def Integrate(func, z0, z1):
    val, err = sp.integrate.quadrature(func, z0, z1, maxiter=Sets.quad_iter)
    return val, err

def IntegrateDistance(z0, z1):
    dist, derr = Integrate(Pars.HubbleDistanceFunction, z0, z1)
    return dist, derr

def IntegrateTime(z0, z1):
    time, terr = Integrate(Pars.HubbleTimeFunction, z0, z1)
    return time, terr



def InvertComDist(cdist):
    rhs = lambda x: IntegrateDistance(0.0, x)[0] - cdist
    return sp.optimize.newton( rhs, 1.0, tol=Sets.root_tol )

def InvertLumDist(ldist):
    rhs = lambda x: (1.0+x)*IntegrateDistance(0.0, x)[0] - ldist
    return sp.optimize.newton( rhs, 1.0, tol=Sets.root_tol )


    

def CosmologicalParameters(redz, pout=False):
    """
    Given a redshift and fractional comoving distance (D_C/D_H), return cosmological parameters.
    """
    
    dist, dist_err = IntegrateDistance(0.0, redz)
    time, time_err = IntegrateTime(0.0, redz)

    retstr = "z = %e\n" % (redz)

    # Comoving Distance
    comDist       = Pars.HubbleDistance()*dist
    comDist_err   = Pars.HubbleDistance()*dist_err
    retstr       += "D_C  =  %+e +- %e  [cm]  ~  %+e [Mpc] : Comoving Distance\n" % (comDist, comDist_err, comDist/(1.0e6*Pars.parsec) )

    # Luminosity Distance
    lumDist       = (1.0 + redz)*comDist
    lumDist_err   = (1.0 + redz)*comDist_err
    retstr       += "D_L  =  %+e +- %e  [cm]  ~  %+e [Mpc] : Luminosity Distance\n" % (lumDist, lumDist_err, lumDist/(1.0e6*Pars.parsec) )

    # Angular Diameter Distance
    angDist       = comDist/(1.0 + redz)
    angDist_err   = comDist_err/(1.0 + redz)
    retstr       += "D_A  =  %+e +- %e  [cm]  ~  %+e [Mpc] : Angular Diameter Distance\n" % (angDist, angDist_err, lumDist/(1.0e6*Pars.parsec)) 

    # Arcsecond size
    arcSize       = Pars.arcsec*angDist
    arcSize_err   = Pars.arcsec*angDist_err
    retstr       += "Arc  =  %+e +- %e  [cm]  ~  %+e [pc]  : Arcsecond Transverse Distance\n" % (arcSize, arcSize_err, arcSize/(Pars.parsec) )

    # Lookback Time
    lookTime      = Pars.HubbleTime()*time
    lookTime_err  = Pars.HubbleTime()*time_err
    retstr       += "T_L  =  %+e +- %e  [s]   ~  %+e [Myr] : Lookback Time\n" % (lookTime, lookTime_err, lookTime/(1.0e6*Pars.year) )

    # Age
    ageTime       = Pars.HubbleTime()*(1.0-time)
    ageTime_err   = lookTime_err
    retstr       += "t    =  %+e +- %e  [s]   ~  %+e [Myr] : Age of the Universe\n" % (ageTime, ageTime_err, ageTime/(1.0e6*Pars.year) )

    # Distance Modulus
    distMod       = 5.0*np.log10(lumDist/(Pars.parsec*10.0))
    distMod_err   = 5.0*np.log10(lumDist_err/(Pars.parsec*10.0))
    retstr       += "DM   =  %+e +- %e  []    Distance Modulus\n" % (distMod, np.abs(distMod_err) )

    if( pout ): print retstr

    return [ [comDist, comDist_err, "Comoving Distance [cm]"] , \
             [lumDist, lumDist_err, "Luminosity Distance [cm]"] , \
             [angDist, angDist_err, "Angular Diameter Distance [cm]"] , \
             [arcSize, arcSize_err, "Arcsecond Transverse Size [cm]"] , \
             [lookTime, lookTime_err, "Lookback Time [s]"] , \
             [ageTime, ageTime_err, "Age of the Universe [s]"] , \
             [distMod, distMod_err, "Distance Modulus []"] ]




def RedshiftFromTarget(t_sets):
    """
    Given a settings object, determine which parameter is being targeted.

    The settings object contains the default, and command-line options.  Use this information
    to determine which cosmological parameter ('target') is being provided (e.g. redshift),
    and return that target parameter, along with the appropriate function to obtain that
    parameter from a line from the table (e.g. 'RedshiftFromLine').
    """

    all_targets = [ t_sets.z_target, t_sets.cd_target, t_sets.ld_target ]

    # Make sure only one target is specified (>= 0.0)
    if( sum(tarz >= 0.0 for tarz in all_targets) != 1 ): 
        ERROR(func='DetermineTarget', errstr='Must specify a single target parameter.', ex=True)


    # Redshift
    if( t_sets.z_target >= 0.0 ): 
        if( t_sets.verbose ): print " - Redshift"
        return t_sets.z_target #, RedshiftFromLine


    # Comoving Distance
    if( t_sets.cd_target >= 0.0 ):
        if( t_sets.verbose ): print " - Comoving Distance"
        t_cd = t_sets.cd_target
        if(   t_sets.use_pc  ): t_cd *= Pars.parsec
        elif( t_sets.use_mpc ): t_cd *= Pars.parsec*(1.0e6)
        elif( t_sets.use_ly  ): t_cd *= Pars.year*Pars.c

        return InvertComDist( t_cd/Pars.HubbleDistance() )


    # Luminosity Distance
    if( t_sets.ld_target >= 0.0 ):
        if( t_sets.verbose ): print " - Luminosity Distance"
        t_ld = t_sets.ld_target
        if(   t_sets.use_pc  ): t_ld *= Pars.parsec
        elif( t_sets.use_mpc ): t_ld *= Pars.parsec*(1.0e6)
        elif( t_sets.use_ly  ): t_ld *= Pars.year*Pars.c

        return InvertLumDist( t_ld/Pars.HubbleDistance() )


    ERROR(func='DetermineTarget', errstr='Must specify a single target parameter.', ex=True)

    


def InitArgParse():
    """
    Initialize the argparse object with desired command line options, then retrieve their values.
    """

    parser          = ArgumentParser()

    parser.add_argument('-z',      type=float, default=Sets.z_target,       help='Target redshift z'                                        )
    parser.add_argument('-cd',     type=float, default=Sets.cd_target,      help='Target coming distance D_C'                               )
    parser.add_argument('-ld',     type=float, default=Sets.ld_target,      help='Target luminosity distance D_C'                           )

    parser.add_argument('--print', dest='prt',  action='store_true', default=False, help="Print defaul cosmological parameters"                  )
    parser.add_argument('--plot',               action='store_true', default=False, help="Plot cosmological parameters"                          )
    parser.add_argument('--mpc',                action='store_true', default=False, help="Use provided distance ('-cd' or '-ld') in megaparsecs" )
    parser.add_argument('--ly',                 action='store_true', default=False, help="Use provided distance ('-cd' or '-ld') in lightyears"  )
    parser.add_argument('--pc',                 action='store_true', default=False, help="Use provided distance ('-cd' or '-ld') in parsecs"     )

    args            = parser.parse_args()
    Sets.print_flag = args.prt
    Sets.plot_flag  = args.plot
    Sets.z_target   = float(args.z)
    Sets.cd_target  = float(args.cd)
    Sets.ld_target  = float(args.ld)
    Sets.use_pc     = args.pc
    Sets.use_mpc    = args.mpc
    Sets.use_ly     = args.ly
 



def main():

    args         = InitArgParse()

    # Print Cosmological parameters
    if( Sets.print_flag ):
        print Pars.ParameterString()


    redshift     = RedshiftFromTarget(Sets)
    cosmo_params = CosmologicalParameters(redshift, pout=True)




if __name__ == "__main__": main()





# 
