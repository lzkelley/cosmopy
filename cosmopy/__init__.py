"""
"""

__version__ = "3.2.1"
__author__ = "Luke Zoltan Kelley <lzkelley@northwestern.edu>"
__license__ = "MIT"

import astropy as ap
import astropy.constants  # noqa
import astropy.units  # noqa


NWTG = ap.constants.G.cgs.value
SPLC = ap.constants.c.cgs.value
MSOL = ap.constants.M_sun.cgs.value
LSOL = ap.constants.L_sun.cgs.value
RSOL = ap.constants.R_sun.cgs.value
PC = ap.constants.pc.cgs.value
AU = ap.constants.au.cgs.value

ARCSEC = ap.units.arcsec.cgs.scale       # arcsecond in radians
YR = ap.units.year.to(ap.units.s)       # year in seconds

U_SEC = ap.units.Unit('s')       # second unit
U_CM = ap.units.Unit('cm')       # centimeter unit

COLORED_OUTPUT = True


from . __main__ import get_cosmology, api, Cosmology  # noqa
