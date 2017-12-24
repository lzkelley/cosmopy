"""
"""

import numpy as np
import scipy as sp
import scipy.optimize
import scipy.integrate  # noqa


from . settings import Settings
from . parameters import Parameters


def DistToCM(indist):
    """
    Based on unit settings, convert distance to centimeters.

    Setup to accept [default] centimeters (cm), parsecs (pc), megaparsecs (mpc),
    or lightyears (ly).  Stored in Settings object (tst).
    """
    if Settings.use_pc:
        indist *= Parameters.parsec
    elif Settings.use_mpc:
        indist *= Parameters.parsec*(1.0e6)
    elif Settings.use_ly:
        indist *= Parameters.year*Parameters.c
    return indist


def TimeToS(intime):
    """
    Based on unit settings, convert time to seconds.

    Setup to accept [default] seconds (s), years (yr), megayears (myr).  Stored
    in Settings object (tst).
    """
    if Settings.use_yr: intime *= Parameters.year
    elif Settings.use_myr: intime *= Parameters.year*(1.0e6)
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


def RootFind(func):
    return sp.optimize.newton(func, 1.0, tol=Settings.root_tol)


def InvertScaleFactor(ascale):
    """Find redshift from Scale-Factor"""
    return (1.0/ascale) - 1.0


def InvertComDist(cdist):
    """Given a comoving-distance, find redshift."""
    return RootFind(lambda x: IntegrateDistance(0.0, x)[0] - cdist)


def InvertLumDist(ldist):
    """Given a luminosity-distance, find redshift."""
    return RootFind(lambda x: (1.0+x)*IntegrateDistance(0.0, x)[0] - ldist)


def InvertLookTime(ltime):
    """Given a look-back time, find redshift."""
    return RootFind(lambda x: IntegrateTime(0.0, x)[0] - ltime)


def InvertAgeTime(atime):
    """Given an age of the universe, find redshift."""
    return RootFind(lambda x: 1.0 - IntegrateTime(0.0, x)[0] - atime)


def CosmologicalParameters(redz, pout=False):
    """
    Given a redshift and fractional comoving distance (D_C/D_H), return cosmological parameters.
    """

    pars = []
    errs = []
    retstr = ""

    # Store Redshift [0]
    pars.append(redz)
    errs.append(0.0)
    rstr = "z = {:e}".format(redz)
    retstr += rstr + "\n"

    # Store scale-factor [1]
    scaleFactor = 1.0/(1.0+redz)
    pars.append(scaleFactor)
    errs.append(0.0)
    rstr = "a = {:e}".format(scaleFactor)
    retstr += rstr + "\n"

    # Integrate over redshift to find distance and time with errors
    dist, dist_err = IntegrateDistance(0.0, redz)
    time, time_err = IntegrateTime(0.0, redz)

    # Comoving Distance [2]
    comDist = Parameters.HubbleDistance()*dist
    comDist_err = Parameters.HubbleDistance()*dist_err
    pars.append(comDist)
    errs.append(comDist_err)
    rstr = "D_C =  {:+e} +- {:e}  [cm]  ~  {:+e} [Mpc] : Comoving Distance".format(
        comDist, comDist_err, comDist/(1.0e6*Parameters.parsec))
    retstr += rstr + "\n"

    # Luminosity Distance [3]
    lumDist = (1.0 + redz)*comDist
    lumDist_err = (1.0 + redz)*comDist_err
    pars.append(lumDist)
    errs.append(lumDist_err)
    rstr = "D_L =  {:+e} +- {:e}  [cm]  ~  {:+e} [Mpc] : Luminosity Distance".format(
        lumDist, lumDist_err, lumDist/(1.0e6*Parameters.parsec))
    retstr += rstr + "\n"

    # Angular Diameter Distance [4]
    angDist = comDist/(1.0 + redz)
    angDist_err = comDist_err/(1.0 + redz)
    pars.append(angDist)
    errs.append(angDist_err)
    rstr = "D_A =  {:+e} +- {:e}  [cm]  ~  {:+e} [Mpc] : Angular Diameter Distance".format(
        angDist, angDist_err, lumDist/(1.0e6*Parameters.parsec))
    retstr += rstr + "\n"

    # Arcsecond size [5]
    arcSize = Parameters.arcsec*angDist
    arcSize_err = Parameters.arcsec*angDist_err
    pars.append(arcSize)
    errs.append(arcSize_err)
    rstr = "Arc =  {:+e} +- {:e}  [cm]  ~  {:+e} [pc]  : Arcsecond Transverse Distance".format(
        arcSize, arcSize_err, arcSize/(Parameters.parsec))
    retstr += rstr + "\n"

    # Lookback Time [6]
    lookTime = Parameters.HubbleTime()*time
    lookTime_err = Parameters.HubbleTime()*time_err
    pars.append(lookTime)
    errs.append(lookTime_err)
    rstr = "T_L =  {:+e} +- {:e}  [s]   ~  {:+e} [Myr] : Lookback Time".format(
        lookTime, lookTime_err, lookTime/(1.0e6*Parameters.year))
    retstr += rstr + "\n"

    # Age [7]
    ageTime = Parameters.HubbleTime()*(1.0-time)
    ageTime_err = lookTime_err
    pars.append(ageTime)
    errs.append(ageTime_err)
    rstr = "T_A =  {:+e} +- {:e}  [s]   ~  {:+e} [Myr] : Age of the Universe".format(
        ageTime, ageTime_err, ageTime/(1.0e6*Parameters.year))
    retstr += rstr + "\n"

    # Distance Modulus [8]
    distMod = 5.0*np.log10(lumDist/(Parameters.parsec*10.0)) if lumDist > 0.0 else 0.0
    distMod_err = 5.0*np.log10(lumDist_err/(Parameters.parsec*10.0)) if lumDist_err > 0.0 else 0.0
    pars.append(distMod)
    errs.append(distMod_err)
    rstr = "DM =  {:+e} +- {:e}  []    Distance Modulus".format(distMod, np.abs(distMod_err))
    retstr += rstr + "\n"

    if pout:
        print(retstr)

    return pars, errs


def solve(z=None, a=None, cd=None, ld=None, lt=None, t=None, verbose=False, pout=False):

    # Redshift
    if z is not None:
        if verbose:
            print(" - Redshift")
        redz = z

    # Scale-factor
    elif a is not None:
        if verbose:
            print(" - Scale-factor")
        redz = InvertScaleFactor(a)

    # Comoving Distance
    elif cd is not None:
        t_cd = DistToCM(cd)
        if verbose:
            print(" - Comoving Distance, {:.2e} [cm]".format(t_cd))
        redz = InvertComDist(t_cd/Parameters.HubbleDistance())

    # Luminosity Distance
    elif ld is not None:
        t_ld = DistToCM(ld)
        if verbose:
            print(" - Luminosity Distance, {:.2e} [cm]".format(t_ld))
        redz = InvertLumDist(t_ld/Parameters.HubbleDistance())

    # Look-back Time
    elif lt is not None:
        t_tl = TimeToS(lt)
        if verbose:
            print(" - Look-back Time, {:.2e} [s]".format(t_tl))
        redz = InvertLookTime(t_tl/Parameters.HubbleTime())

    # Universe Age Time
    elif t is not None:
        t_ta = TimeToS(t)
        if verbose:
            print(" - Age of the Universe, {:.2e} [s]".format(t_ta))
        redz = InvertAgeTime(t_ta/Parameters.HubbleTime())

    else:
        print(Parameters.ParameterString())
        raise RuntimeError("No parameters given!")

    # Make sure redz is an iterable
    if not np.iterable(redz):
        redz = [redz]

    pars = []
    errs = []
    for zz in redz:
        pp, ee = CosmologicalParameters(zz, pout=pout)
        pars.append(pp)
        errs.append(ee)

    pars = np.array(pars)
    errs = np.array(errs)

    return pars.T, errs.T
