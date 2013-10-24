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

import os
import sys
from argparse import ArgumentParser


class Sets(object):
    """
    Object to store the current configuration of CosmoCalc.
    """

    z_steps  = 1e3                                                                                  # Default number of integration steps per 1z
    z_max    = 50.0                                                                                 # Default maximum redshift for integration

    build_flag = False                                                                              # Default whether to rebuild integration table
    use_pc     = False                                                                              # Use input distance in parsecs
    use_mpc    = False                                                                              # Use input distance in Megaparsecs
    use_ly     = False                                                                              # Use input distance in lightyears

    table_filename      = ".cc-table.dat"                                                           # Name of integration table file
    script_path         = os.path.dirname( os.path.abspath(__file__) )                              # Path to this python script
    full_table_filename = script_path + "/" + table_filename

    z_target  = -1.0
    cd_target = -1.0
    ld_target = -1.0
    header    = 4                                                                                   # Number of header lines in integration table

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



def CreateTable(zsteps, zmax, filename=Sets.full_table_filename):
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
            outfile.write("%8d %.12e %.12e %.12e %.12e %.12e\n" % \
                          (count, redz, dist_trapsum, dist_err, time_trapsum, time_err) ) 

            # Increment counter, and update progress bar
            count    += 1
            bar.update()


    return count, redz, dist_trapsum, dist_err, time_trapsum, time_err



def CheckTable(filename=Sets.full_table_filename):
    return os.path.isfile(filename)


def ParseSettingsLine(inline):
    """Parse a line of format '... <VAR> = <VALUE>' to retrieve the value."""
    inline = inline.split('=')
    return inline[-1].strip()


def ParseTableLine(inline):
    return [ float(word) for word in inline.split() ]


def RedshiftFromLine(inline):
    return np.array([ inline[Sets.IND_Z ] , 0.0 ])


def ComDistFromLine(inline):
    return np.array([ Pars.HubbleDistance()*inline[Sets.IND_DF ] , 
                      Pars.HubbleDistance()*inline[Sets.IND_DFE] ]) 


def LumDistFromLine(inline):
    return np.array([ (1.+inline[Sets.IND_Z])*Pars.HubbleDistance()*inline[Sets.IND_DF ] , 
                      (1.+inline[Sets.IND_Z])*Pars.HubbleDistance()*inline[Sets.IND_DFE] ]) 



def BoundsFromTable(targ, func, filename=Sets.full_table_filename):
    """
    Retrieve the table values nearest in reshift to the one provided.

    Note: getline seems to 1-index lines
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
    int_num_lines = int(np.floor(num_lines))
    par_max, dummy = func(  ParseTableLine( lc_getline(filename, int_num_lines) )  )
    # line_frac  = tarz/Sets.z_max
    targ_frac   = targ/par_max

    # Don't allow calculations beyond table
    if( targ_frac > 1.0 or targ > par_max ):
        ERROR("BoundsFromTable() - targ %e seems larger than maximum in file %e !\n" % \
              (targ, par_max), 2, ex=True)

        
    guess_low  = Sets.header+int(np.floor(int_num_lines*targ_frac))-1
    # guess_low  = np.max([Sets.header, guess_low])+1
    line_low, dummy = func(  ParseTableLine( lc_getline(filename, guess_low) )  )

    it_count = 0
    while( line_low > targ and guess_low > Sets.header+1):
        guess_low -= 1
        line_low, dummy = func(  ParseTableLine( lc_getline(filename, guess_low) )  )
        it_count += 1
        if( it_count > 100 ): 
            ERROR("BoundsFromTable() - @low it_count = %d, guess_low = %d, line_low = %e!\n" % \
                  it_count, guess_low, line_low, ex=True )
        
    guess_high = guess_low+1
    line_high, dummy = func(  ParseTableLine( lc_getline(filename, guess_high) )  )
    
    while( line_high < targ and guess_high < int_num_lines ):
        guess_low   = guess_high
        guess_high += 1
        line_high   = func(  ParseTableLine( lc_getline(filename, guess_high ) ) )
        line_low    = func(  ParseTableLine( lc_getline(filename, guess_low  ) ) )
        it_count += 1
        if( it_count > 100 ): 
            ERROR("BoundsFromTable() - @high it_count = %d, guess_high = %d, line_high = %e!\n" % \
                  it_count, guess_high, line_high, ex=True )


    if( it_count > 10 or line_high < targ or line_low > targ ):
        ERROR("BoundsFromTable() - it_count = %d  low,high = %d,%d  low,high = %e, %e!\n" % \
              it_count, guess_low, guess_high, line_low, line_high, ex=True )


    return [ ParseTableLine( lc_getline(filename, guess_low  ) ), 
             ParseTableLine( lc_getline(filename, guess_high ) ) ]




def InterpolateParameter(targ, func, bef, aft):
    """
    Interpolate to the target redshift from the two given sets of values.

    Recall the parameter order is given by Sets.IND_*
    Target z sould never occur after the second set, but may occur before the first.
    """

    # Calculate interpolation distance
    deltat     = targ - bef[1]

    # Calculate Interpolant slopes
    dist_slope = (aft[Sets.IND_DF]-bef[Sets.IND_DF])/(aft[Sets.IND_Z]-bef[Sets.IND_Z])
    time_slope = (aft[Sets.IND_TF]-bef[Sets.IND_TF])/(aft[Sets.IND_Z]-bef[Sets.IND_Z])

    # interpolate
    dist_new   = bef[Sets.IND_DF] + deltat*dist_slope
    time_new   = bef[Sets.IND_TF] + deltat*dist_slope

    # Retrieve larger of possible errors (as upper limit)
    dist_err   = np.max( [bef[Sets.IND_DFE],aft[Sets.IND_DFE]] )
    time_err   = np.max( [bef[Sets.IND_TFE],aft[Sets.IND_TFE]] )
    
    # Give a fractional integration number (between interpolants)
    frac = bef[Sets.IND_NUM] + (deltat/(aft[Sets.IND_Z]-bef[Sets.IND_Z]))

    # Make sure this parameter order matches Sets.IND_*
    return [ frac, tarz, dist_new, dist_err, time_new, time_err ]

    

def CosmologicalParameters(zscales, pout=False):
    """
    Given a redshift and fractional comoving distance (D_C/D_H), return cosmological parameters.
    """
    
    zz        = zscales[Sets.IND_Z  ]
    dist      = zscales[Sets.IND_DF ]
    dist_err  = zscales[Sets.IND_DFE]
    time      = zscales[Sets.IND_TF ]
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
    before, after = BoundsFromTable(tarz, RedshiftFromLine)

    # before, after = RedshiftsFromTable(tarz)
    zmatch        = InterpolateRedshifts(tarz, before, after)
    zscales       = CosmologicalParameters(zmatch, pout=pout)
    return zscales



def ConvertTargetsToRedshift(t_sets):
    all_targets = [ Sets.z_target, Sets.cd_target, Sets.ld_target ]

    # Make sure only one target is specified (>= 0.0)
    if( sum(tarz >= 0.0 for tarz in all_targets) != 1 ): 
        ERROR(errstr='Must specify a single target parameter.', ex=True)


    # Redshift
    if( Sets.z_target >= 0.0 ): 
        return Sets.z_target


    # Comoving Distance
    if( Sets.cd_target >= 0.0 ):

        t_cd = Sets.cd_target
        if(   Sets.use_pc  ): t_cd *= Pars.parsec
        elif( Sets.use_mpc ): t_cd *= Pars.parsec*(1.0e6)
        elif( Sets.use_ly  ): t_cd *= Pars.year*Pars.c

        return ComDistToRedshift(t_cd)


    # Luminosity Distance
    if( Sets.ld_target >= 0.0 ):

        t_ld = Sets.cd_target
        if(   Sets.use_pc  ): t_ld *= Pars.parsec
        elif( Sets.use_mpc ): t_ld *= Pars.parsec*(1.0e6)
        elif( Sets.use_ly  ): t_ld *= Pars.year*Pars.c

        return LumDistToRedshift(t_ld)



def InitArgParse():
    """
    Initialize the argparse object with desired command line options, then retrieve their values.
    """

    parser = ArgumentParser()
    parser.add_argument('--build',  nargs=2, metavar=('zstep','zmax'), type=float, default=None, help='Rebuild CosmoCalc table (cc-table.dat)' )
    parser.add_argument('-z', type=float, default=Sets.z_target, help='Target redshift z')
    parser.add_argument('-cd', type=float, default=Sets.cd_target, help='Target coming distance D_C'    )
    parser.add_argument('-ld', type=float, default=Sets.ld_target, help='Target luminosity distance D_C')
    parser.add_argument('--mpc',  action='store_true', default=False, help="Use provided distance ('-cd' or '-ld') in megaparsecs")
    parser.add_argument('--ly',   action='store_true', default=False, help="Use provided distance ('-cd' or '-ld') in lightyears" )
    parser.add_argument('--pc',   action='store_true', default=False, help="Use provided distance ('-cd' or '-ld') in parsecs"    )

    args           = parser.parse_args()
    build_args     = args.build
    Sets.z_target  = float(args.z)
    Sets.cd_target = float(args.cd)
    Sets.ld_target = float(args.ld)
    Sets.use_pc    = args.pc
    Sets.use_mpc   = args.mpc
    Sets.use_ly    = args.ly
 
    if( build_args != None ):
        Sets.build_flag          = True
        Sets.z_steps, Sets.z_max = build_args




def main():

    args = InitArgParse()

    # Rebuild the integration table
    if( Sets.build_flag ): 
        res = CreateTable(Sets.z_steps, Sets.z_max)
        print "Steps %d, Final:  z = %e, dist = %e +- %e, time = %e +- %e\n" % \
            ( res[0], res[1], res[2]*Pars.HubbleDistance(), \
              np.abs(res[3]*Pars.HubbleDistance), res[4]*Pars.HubbleTime(), \
              np.abs(res[5]*Pars.HubbleTime()) )
        return

    # Check that integration table exists
    if( not CheckTable() ):
        print "[CosmoCalc.py] Integration table (%s) does not exist." % ( Sets.full_table_filename )
        # Prompt user for action
        usr_inp = raw_input('Create default integration table?  y/[n]:')
        if( usr_inp != 'y' ): exit(0)

        # Create new table with default parameters
        res = CreateTable(Sets.z_steps, Sets.z_max)
        print "Steps %d, Final:  z = %e, dist = %e +- %e, time = %e +- %e\n" % \
            ( res[0], res[1], res[2]*Pars.HubbleDistance(), \
              np.abs(res[3]*Pars.HubbleDistance()), res[4]*Pars.HubbleTime(), \
              np.abs(res[5]*Pars.HubbleTime()) )
        
    """
    all_targets = [ Sets.z_target, Sets.cd_target, Sets.ld_target ]
    if( Sets.z_target >= 0.0 ): ERROR(errstr='Can only target one parameter.', ex=True)

    # Comoving Distance
    if( Sets.cd_target >= 0.0 ):

        use_dist = Sets.cd_target
        if(   Sets.use_pc  ): use_dist *= Pars.parsec
        elif( Sets.use_mpc ): use_dist *= Pars.parsec*(1.0e6)
        elif( Sets.use_ly  ): use_dist *= Pars.year*Pars.c
    """

    # Redshift
    if( Sets.z_target >= 0.0 ): 
        scales = CalculateScales(Sets.z_target, pout=True)
        return




if __name__ == "__main__": main()





# 
