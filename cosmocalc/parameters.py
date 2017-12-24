"""
"""

import numpy as np


class Parameters(object):
    """Standard Cosmological Parameters.

    Cosmological Parameters taken from
    [1] '1212.5226v3 - Nine-Year WMAP Observations: Cosmological Parameter Results'
    [2] '1002.3488   - The Cosmological Parameters 2010'
    """

    # Default Cosmological Parameters
    # -------------------------------
    # Hubble's Constant      [km/s/Mpc]
    H0 = 69.7
    # Age of the Universe    [Gyr]
    T0 = 13.76
    # Dark Matter Density    [1]
    OmegaDM = 0.235
    # Baryonic Density       [1]
    OmegaB = 0.0464
    # Dark Energy Density    [1]
    OmegaL = 0.7185
    # Radiation Density      [2] omegar*h^2 = 2.47e-5; h = 0.704
    OmegaR = 4.984e-4

    # Physical Constants
    # ------------------
    # Speed of Light         [cm/s]
    c = 2.99792458e10
    # Gravitational Constant [cm^3/g/s^2]
    G = 6.67259e-8
    # Proton Mass            [g]
    mProton = 1.6726231e-24
    # Electron Mass          [g]
    mElectron = 9.1093829e-28
    # Planck's Constant      [erg s]
    planckH = 6.6260755e-27
    # Parsec                 [cm]
    parsec = 3.08567758e+18
    # Arcsecond              [radians]
    arcsec = 4.84813681e-06
    # Year                   [s]
    year = 3.15569520e+07

    @staticmethod
    def HubbleTime():
        return 1.0*(1.0e6*Parameters.parsec/1.0e5)/Parameters.H0

    @staticmethod
    def HubbleDistance():
        return Parameters.c*Parameters.HubbleTime()

    @staticmethod
    def CriticalDensity():
        return 3*np.square(Parameters.H0)/(8.*np.pi*Parameters.G)

    @staticmethod
    def OmegaMatter():
        return Parameters.OmegaB+Parameters.OmegaDM

    @staticmethod
    def ParameterString():
        parstr = "\n"
        parstr += "Hubble Constant, H0 = {:f}\n".format(Parameters.H0)
        parstr += "Age of the Universe, T0 = {:f}\n".format(Parameters.T0)
        parstr += "Baryon      Fraction, OmegaB = {:f}\n".format(Parameters.OmegaB)
        parstr += "Dark Matter Fraction, OmegaDM = {:f}\n".format(Parameters.OmegaDM)
        parstr += "Dark Energy Fraction, OmegaL = {:f}\n".format(Parameters.OmegaL)
        parstr += "Radiation   Fraction, OmegaR = {:f}\n".format(Parameters.OmegaR)
        return parstr

    @staticmethod
    def HubbleFunction(zz):
        """The E(z) function from Hogg1999.

        Incorporates matter (dark and light), radiation, and dark energy.  This
        version of the equation assumes that OmegaK (curvature) is zero.
        """

        # Determine the constituents at the given redshifts
        t_matter = Parameters.OmegaMatter()*np.power((1.0+zz), 3)
        t_radiation = Parameters.OmegaR*np.power((1.0+zz), 4)
        t_darkenergy = Parameters.OmegaL
        return np.sqrt(t_matter + t_radiation + t_darkenergy)

    @staticmethod
    def HubbleDistanceFunction(zz):
        return 1.0/Parameters.HubbleFunction(zz)

    @staticmethod
    def HubbleTimeFunction(zz):
        return 1.0/((1.0+zz)*Parameters.HubbleFunction(zz))
