# ==================================================================================================
# CosmoCalc.py
# ------------
#
#
# ------------------
# Luke Zoltan Kelley
# LKelley@cfa.harvard.edu
# ==================================================================================================


from astropy.io import ascii as at
from astropy.utils.console import ProgressBar
from linecache import getline as lc_getline

import numpy as np
#import matplotlib as mpl
#from matplotlib import pyplot as plt

import sys
from argparse import ArgumentParser


class Sets(object):
    z_steps  = 1e5
    z_max    = 50.0

    build_flag = False
    table_filename = "cc-table.dat"

    z_target = -1.0
    header   = 4

    # Indices of parameters for output table file
    IND_NUM  = 0                                                                                    # Step number 
    IND_Z    = 1                                                                                    # Redshift
    IND_DF   = 2                                                                                    # Distance fraction
    IND_DFE  = 3                                                                                    # Distance fraction error
    IND_TF   = 4                                                                                    # Time     fraction
    IND_TFE  = 5                                                                                    # Time     fraction error



class Pars(object):
    """
    Standard Cosmological Parameters.

    Cosmological Parameters taken from
    '1212.5226v3 - Nine-Year WMAP Observations: Cosmological Parameter Results'
    """

    H0      = 69.7                      # Hubble's Constant           [km/s/Mpc]
    T0      = 13.76                     # Age of the Universe         [Gyr]
    OmegaDM = 0.235                     # Dark Matter Density         [1/critical-density]
    OmegaB  = 0.0464                    # Baryonic Density            [1/critical-density]
    OmegaL  = 0.7185                    # Dark Energy Density         [1/critical-density]
    c       = 2.99792458e10             # Speed of Light              [cm/s]
    G       = 6.67259e-8                # Gravitational Constant      [cm^3/g/s^2]
    mProton = 1.6726231e-24             # Proton Mass                 [g]
    mElectron = 1.6726231e-24           # Electron Mass               [g]
    planckH = 6.6260755e-27             # Planck's Constant           [erg s]
    parsec  = 3.08567758e+18            # Parsec                      [cm]
    arcsec  = 4.84813681e-06            # Arcsecond                   [radians]
    year    = 3.15569520e+07            # Year                        [s]

    @staticmethod
    def HubbleTime():      return 1.0*(1.0e6*Pars.parsec/1.0e5)/Pars.H0
    @staticmethod
    def HubbleDistance():  return Pars.c*Pars.HubbleTime()
    @staticmethod
    def CriticalDensity(): return 3*np.square(Pars.H0)/(8.*np.pi*Pars.G)
    @staticmethod
    def OmegaMatter():     return Pars.OmegaB+Pars.OmegaDM



def ERROR(errstr='', errnum=1, ex=False):
    """Print error message."""
    print '\n[CosmoCalc.py] ERROR: %d - %s !!\n\n' % (errnum, errstr)
    if( ex ):
        sys.exit(errnum)
        exit(errnum)


def HubbleFunction(zz): 
    """
    The E(z) function from Hogg1999.

    This version of the equation assumes that OmegaK (curvature) is zero.
    """
    return np.sqrt( Pars.OmegaMatter()*np.power((1+zz),3) + Pars.OmegaL )
    



def CreateTable(zsteps, zmax, filename=Sets.table_filename):
    """
    Perform a riemann sum to build a table of unscaled comoving distances as a function of redshift.

    The coming distance is given by the expression,
    astro-ph/9905116v4 : Hogg - Distance Measures in Cosmology - 1999
      D_C = D_H Int[0,z]( 1.0/E(z') dz' )
      s.t.
      D_H is the Hubble Distance (c/H0), and
      E(z) = sqrt( Omega_M*(1+z)^3 - OmegaK*(1+z)^2 + OmegaL )
    The table calculated here is D_C/D_H which can then be normalized and used as desired.

    The difference between left and right riemann sums are used as an absolute upper-limit to the
    error on the integral for each redshift value.
    """

    count    = 0                                                                                    # Count the step number
    redz     = 0.0                                                                                  # Current Redshift 
    dz       = 1.0/zsteps                                                                           # Calculate z stepsize
    numSteps = 1.0*zmax/dz                                                                          # Total Number of steps

    # Distance integration Parameters
    dist_leftsum  = 0.0
    dist_rightsum = 0.0
    dist_trapsum  = 0.0
    dist_err      = 0.0
    dist_larg     = 1.0/HubbleFunction(redz)
    dist_rarg     = dist_larg

    # Time integration Parameters
    time_leftsum  = 0.0
    time_rightsum = 0.0
    time_trapsum  = 0.0
    time_err      = 0.0
    time_larg     = 1.0/( (1.0+redz)*HubbleFunction(redz) )
    time_rarg     = time_larg

    # Initialize output file and print header
    outfile = open(filename, 'w')
    outfile.write("#  zmax = %e\n#  dz = %e\n" % (zmax, dz) )
    outfile.write("#  Num       z        dist-frac    dist-err    time-frac     time-err\n")
    outfile.write("# ---------------------------------------------------------------------\n")


    # Integrate
    with ProgressBar(numSteps) as bar:
        while( redz+dz*0.99 <= zmax ):
            # Integrate left Riemann Sums
            dist_leftsum  += dist_larg*dz
            time_leftsum  += time_larg*dz                                                           

            # Move Right terms to left
            dist_larg      = dist_rarg
            time_larg      = time_rarg

            # Step forward in redshift
            redz          += dz

            # Calculate new right terms
            dist_rarg      = 1.0/HubbleFunction(redz)
            time_rarg      = dist_rarg/(1.0+redz)

            # Integrate right Riemann sums
            dist_rightsum += dist_rarg*dz
            time_rightsum += time_rarg*dz

            # Integrate trapezoid rule
            dist_trapsum   = 0.5*(dist_leftsum + dist_rightsum)
            time_trapsum   = 0.5*(time_leftsum + time_rightsum)

            # Calculate erros as left-right differences
            dist_err       = np.abs(dist_leftsum - dist_rightsum)
            time_err       = np.abs(time_leftsum - time_rightsum)

            # Print to data file   (Make sure these match the ordering in Sets.IND_*
            outfile.write("%6d %e %e %e %e %e\n" % (count, redz, dist_trapsum, dist_err, time_trapsum, time_err) ) 

            # Increment counter, and update progress bar
            count    += 1
            bar.update()


    return count, redz, dist_trapsum, dist_err, time_trapsum, time_err




def ParseSettingsLine(inline):
    """Parse a line of format '... <VAR> = <VALUE>' to retrieve the value."""
    inline = inline.split('=')
    return inline[-1].strip()


def ParseTableLine(inline):
    return [ float(word) for word in inline.split() ]


def RedshiftsFromTable(tarz, filename=Sets.table_filename):
    """
    Retrieve the table values nearest in reshift to the one provided.
    """

    # Open file with table
    try: infile = open(filename)
    except:
        ERROR("RedshiftsFromTable() - could not open '%s', %s" % (filename, sys.exc_info()[0]), 1)
        raise

    # Retrieve the parameters that made this table
    Sets.z_max = float( ParseSettingsLine(infile.readline()) )
    Sets.dz    = float( ParseSettingsLine(infile.readline()) )

    # Estimate lines around target redshift
    num_lines  = 1.0*Sets.z_max/Sets.dz
    line_frac  = tarz/Sets.z_max
    # Don't allow calculations beyond table
    if( line_frac > 1.0 ):
        ERROR("RedshiftsFromTable() - tarz %e seems larger than maximum in file %e !\n" % \
              (tarz, Sets.z_max), 2, ex=True)
        
    guess_low  = Sets.header+int(np.floor(num_lines*line_frac))-1
    guess_low  = np.max([Sets.header, guess_low])+1
    line_low   = ParseTableLine( lc_getline(filename, guess_low) )                                  # Note: getline seems to 1-index lines

    if( line_low[1] > tarz and guess_low > Sets.header+1):
        line_low   = ParseTableLine( lc_getline(filename, guess_low-1) )
        guess_low -= 1
        
    guess_high = guess_low+1
    line_high  = ParseTableLine( lc_getline(filename, guess_high) )
    
    if( line_high[1] < tarz ):
        line_high   = ParseTableLine( lc_getline(filename, guess_high+1) )
        line_low    = ParseTableLine( lc_getline(filename, guess_high  ) )
        guess_high += 1

    return [ line_low, line_high ]


def InterpolateRedshifts(tarz, bef, aft):
    """
    Interpolate to the target redshift from the two given sets of values.

    Recall the parameter order is given by Sets.IND_*
    Target z sould never occur after the second set, but may occur before the first.
    """

    # Calculate interpolation distance
    deltaz     = tarz - bef[1]

    # Calculate Interpolant slopes
    dist_slope = (aft[Sets.IND_DF]-bef[Sets.IND_DF])/(aft[Sets.IND_Z]-bef[Sets.IND_Z])
    time_slope = (aft[Sets.IND_TF]-bef[Sets.IND_TF])/(aft[Sets.IND_Z]-bef[Sets.IND_Z])

    # interpolate
    dist_new   = bef[Sets.IND_DF] + deltaz*dist_slope
    time_new   = bef[Sets.IND_TF] + deltaz*dist_slope

    # Retrieve larger of possible errors (as upper limit)
    dist_err   = np.max( [bef[Sets.IND_DFE],aft[Sets.IND_DFE]] )
    time_err   = np.max( [bef[Sets.IND_TFE],aft[Sets.IND_TFE]] )
    
    # Give a fractional integration number (between interpolants)
    frac = bef[Sets.IND_NUM] + (deltaz/(aft[Sets.IND_Z]-bef[Sets.IND_Z]))

    # Make sure this parameter order matches Sets.IND_*
    return [ frac, tarz, dist_new, dist_err, time_new, time_err ]

    

def CosmologicalParameters(zscales, pout=False):
    """
    Given a redshift and fractional comoving distance (D_C/D_H), return cosmological parameters.
    """
    
    zz        = zscales[Sets.IND_Z]
    dist      = zscales[Sets.IND_DF]
    dist_err  = zscales[Sets.IND_DFE]
    time      = zscales[Sets.IND_TF]
    time_err  = zscales[Sets.IND_TFE]

    retstr = "z = %e\n" % (zz)

    # Comoving Distance
    comDist       = Pars.HubbleDistance()*dist
    comDist_err   = Pars.HubbleDistance()*dist_err
    retstr       += "D_C  =  %+e +- %e  [cm]  ~  %+e [Mpc] : Comoving Distance\n" % (comDist, comDist_err, comDist/(1.0e6*Pars.parsec) )

    # Luminosity Distance
    lumDist       = (1.0 + zz)*comDist
    lumDist_err   = (1.0 + zz)*comDist_err
    retstr       += "D_L  =  %+e +- %e  [cm]  ~  %+e [Mpc] : Luminosity Distance\n" % (lumDist, lumDist_err, lumDist/(1.0e6*Pars.parsec) )

    # Angular Diameter Distance
    angDist       = comDist/(1.0 + zz)
    angDist_err   = comDist_err/(1.0 + zz)
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


def CalculateScales(tarz, pout=False):
    before, after = RedshiftsFromTable(tarz)
    zmatch        = InterpolateRedshifts(tarz, before, after)
    zscales       = CosmologicalParameters(zmatch, pout=pout)
    return zscales


def InitArgParse():
    """
    Initialize the argparse object with desired command line options, then retrieve their values.
    """

    # parser = argparse.ArgumentParser()
    parser = ArgumentParser()

    # parser.add_argument('-o',                   metavar='OUTDIR',          type=str, default=None, required=True,        help='outpur directory' )
    # parser.add_argument('--flux',                      type=str, default=None,       choices=RIEMANN_CHOICES,     help='Riemann solver')
    # parser.add_argument('-d',              action='store_true', default=False,               help='display density'                                     )
    parser.add_argument('--build',  nargs=2, metavar=('zstep','zmax'), type=float, default=None,  help='Rebuild CosmoCalc table (cc-table.dat)')
    parser.add_argument('-z', type=float, default=-1.0, help='Target redshift z')

    args          = parser.parse_args()
    temp_build    = args.build
    Sets.z_target = float(args.z)

    if( temp_build != None ):
        Sets.build_flag          = True
        Sets.z_steps, Sets.z_max = temp_build




def main():
    
    args = InitArgParse()

    # Rebuild the integration table
    if( Sets.build_flag ): 
        res = CreateTable(Sets.z_steps, Sets.z_max)
        print "Steps %d, Final:  z = %e, dist = %e +- %e, time = %e +- %e\n" % \
            ( res[0], res[1], res[2]*Pars.HubbleDistance(), \
              np.abs(res[3]*Pars.HubbleDistance), res[4]*Pars.HubbleTime(), np.abs(res[5]*Pars.HubbleTime()) )

    if( Sets.z_target >= 0.0 ): 
        scales = CalculateScales(Sets.z_target, pout=True)


if __name__ == "__main__": main()





# 
