"""
"""

import os

import numpy as np
import scipy as sp
import scipy.interpolate  # noqa
import astropy as ap
import astropy.cosmology  # noqa

INTERP = "quadratic"                     # Type of interpolation for scipy
FLT_TYPE = np.float32
IMPOSE_FLOOR = True                # Enforce a minimum for certain parameters (e.g. comDist)
MAXIMUM_SCALE_FACTOR = 0.9999      # Scalefactor for minimum of certain parameters (e.g. comDist)


class Cosmology(ap.cosmology.FlatLambdaCDM):
    """

    Methods
    -------
    -   scale_to_age         - Convert from scale-factor to age of the universe [seconds].
    -   age_to_scale         - Convert from age of the universe [seconds] to scale-factor.

    """
    Omega0 = 0.2726
    OmegaLambda = 0.7274
    OmegaBaryon = 0.0456
    HubbleParam = 0.704
    H0 = HubbleParam * 100.0

    # z=0.0 is added automatically
    _Z_GRID = [1000.0, 10.0, 4.0, 2.0, 1.0, 0.5, 0.1, 0.01]
    _INTERP_POINTS = 20

    def __init__(self, args):
        # Initialize parent class
        super().__init__(H0=self.H0, Om0=self.Omega0, Ob0=self.OmegaBaryon)

        # Create grids for interpolations
        # -------------------------------
        #    Create a grid in redshift
        zgrid = self._init_interp_grid(self._Z_GRID, self._INTERP_POINTS)
        self._grid_z = zgrid
        self._sort_z = np.argsort(zgrid)
        # Calculate corresponding values in desired parameters
        #    Ages in seconds
        self._grid_age = self.age(zgrid).cgs.value
        self._sort_age = np.argsort(self._grid_age)
        #    Comoving distances in centimeters
        # self._grid_comdist = self.comoving_distance(zgrid).cgs.value
        # self._sort_comdist = np.argsort(self._grid_comdist)
        #    Comoving volume of the universe
        # self._grid_comvol = self.comoving_volume(zgrid).cgs.value
        return

    def _init_interp_grid(self, z_pnts, num_pnts):
        """Create a grid in redshift for interpolation.

        Arguments
        ---------
        z_pnts : (N,) array_like of scalar
            Monotonically decreasing redshift values greater than zero.  Inbetween each pair of
            values, `num_pnts` points are evenly added in log-space.  Also, `num_pnts` are added
            from the minimum given value down to redshift `0.0`.
        num_pnts : int
            Number of grid points to insert between each pair of values in `z_pnts`.

        Returns
        -------
        zgrid : (M,) array of scalar
            The total length of the final array will be ``N * num_pnts + 1``; that is, `num_pnts`
            for each value given in `z_pnts`, plus redshift `0.0` at the end.

        """
        # Make sure grid is monotonically decreasing
        if not np.all(np.diff(z_pnts) < 0.0):
            err_str = "Non-monotonic z_grid = {}".format(z_pnts)
            raise ValueError(err_str)

        z0 = z_pnts[0]
        zgrid = None
        # Create a log-spacing of points between each pair of values in the given points
        for z1 in z_pnts[1:]:
            temp = np.logspace(*np.log10([z0, z1]), num=num_pnts, endpoint=False)
            if zgrid is None:
                zgrid = temp
            else:
                zgrid = np.append(zgrid, temp)
            z0 = z1

        # Linearly space from the last point up to `0.0` (cant reach `0.0` with log-spacing)
        zgrid = np.append(zgrid, np.linspace(z0, 0.0, num=num_pnts))
        return zgrid

    @staticmethod
    def _a_to_z(sf):
        """Convert from scale-factor to redshift.
        """
        return (1.0/sf) - 1.0

    @staticmethod
    def _z_to_a(redz):
        """Convert from redshift to scale-factor.
        """
        return 1.0/(redz + 1.0)

    @staticmethod
    def _interp(vals, xx, yy, inds=None):
        if inds is None:
            inds = np.argsort(xx)
        # PchipInterpolate guarantees monotonicity with higher order.  Gives more accurate results.
        arrs = sp.interpolate.PchipInterpolator(xx[inds], yy[inds], extrapolate=False)(vals)
        return arrs

    '''
    def scale_to_age(self, sf):
        """Convert from scale-factor to age of the universe [seconds].
        """
        zz = self._scale_to_z(sf)
        age = self.age(zz).cgs.value
        return age

    def scale_to_comdist(self, sf):
        """Convert from scale-factor to comoving distance [cm].
        """
        zz = self._scale_to_z(sf)
        comdist = self.comoving_distance(zz).cgs.value
        return comdist
    '''

    def tage_to_z(self, age):
        """Convert from age of the universe [seconds] to redshift.
        """
        zz = self._interp(age, self._grid_age, self._grid_z, self._sort_age)
        return zz
