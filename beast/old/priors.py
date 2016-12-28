"""
Prior estimation
================

Computing posterior values on a grid are affected by the intrisinc sampling of
that given grid. If the sampling is uniform in all dimension, then the prior is uniform.
Otherwise the prior value is the local density at each grid point (normalized to the entire space)

In this package we explore the local density estimation using the nearest neighbors.
"""

import numpy as np
from scipy.spatial import cKDTree

try:
    from math import gamma as fn_gamma
except:
    from scipy.special import gamma as fn_gamma


class KDTreeDensityEstimator(object):
    """
    KDTreeInterpolator(points, values)

    Nearest-neighbours (barycentric) interpolation in N dimensions.

    This interpolator uses a KDTree to find the closest neighbours of a ND point
    and returns a barycentric interpolation (uses .bary.RbfInterpolator)

    if the number of neighbours is 1 or if the distance to the closest match is
    smaller than `eps`, the value of this point is returned instead.

    Parameters
    ----------
    points : (Npoints, Ndims) ndarray of floats
        Data point coordinates.

    normalize: bool
        Transform the data space to [0, 1] x ndim
        Makes each dimension weighted equally

    Notes
    -----
    Uses ``scipy.spatial.cKDTree``
    """

    def __init__(self, x, normalize=True):
        self.points = np.asarray(x)
        npoints, ndim = x.shape
        self.npoints = npoints
        self.ndim = ndim
        self._norm = False
        if normalize is True:
            self.normalize()
        self.tree = cKDTree(self.points)

    def normalize(self):
        """ Normalize the data space to [0, 1] x ndim """
        self._norm = True
        self.shift = self.points.min(axis=0)
        r = (self.points - self.shift)
        self.scale = r.max(axis=0)
        self.scale[self.scale == 0.] = 1.
        self.points = (self.points - self.shift) / self.scale

    def denormalize(self, x):
        """ Normalize the data space to [0, 1] x ndim """
        if self._norm is True:
            return x * self.scale + self.shift
        else:
            return x

    def __call__(self, *args, **kwargs):
        """
        Evaluate interpolator at given points.

        Parameters
        ----------
        xi : ndarray of float, shape (..., ndim)
            Points where to interpolate data at.

        k : integer
            The number of nearest neighbors to use.

        eps : non-negative float

            Return approximate nearest neighbors; the kth returned value
            is guaranteed to be no further than (1+eps) times the
            distance to the real k-th nearest neighbor.

        p : float, 1<=p<=infinity
            Which Minkowski p-norm to use.
            1 is the sum-of-absolute-values "Manhattan" distance
            2 is the usual Euclidean distance
            infinity is the maximum-coordinate-difference distance

        """
        xi = np.squeeze(np.asarray(args))
        s = xi.shape
        if s[1] != self.ndim:
            raise AttributeError('Points must have {0:d} dimensions, found {1:d}.'.format(self.ndim, s[1]))

        k = kwargs.get('k',10 )
        eps = kwargs.get('eps', 0.)
        return self.simple_density(xi, k=k, eps=eps)
    def simple_density(self, x, k=2, eps=0.):
        """ estimate density at x
            The density at a point x is estimated by n(x) ~ k / r_k^n

        Parameters
        ----------
        x: ndarray of float, shape (..., ndim)
            Points where to interpolate data at.

        k: int
            number of neighbors to use

        eps : non-negative float

            Return approximate nearest neighbors; the kth returned value
            is guaranteed to be no further than (1+eps) times the
            distance to the real k-th nearest neighbor.
        """
        xi, shape = self._prepare_x(x)
        d, ind = self.tree.query(xi, k=k, eps=eps)
        d = np.mean(d[:, 1:], axis=1)
        return float(k) / (hypersphere_volume(d, ndim=xi.shape[1]))

    def _prepare_x(self, x):
        """ make sure the data is correctly presented """
        xi = np.asarray(x)
        shape = xi.shape

        if self._norm is True:
            xi = (xi - self.shift) / self.scale

        if len(shape) < 1:
            raise ValueError('single value passed')

        if len(shape) == 1:
            return xi[..., np.newaxis], shape
        else:
            return xi, shape


def hypersphere_volume(r, ndim):
    """
    compute the hypersphere volume (or n-volume of a sphere) of radius r in n
    dimensions

    keywords
    --------
    r: float or ndarray
        radius or radii of the sphere(s)

    ndim: int
        dimensionality of the sphere
    """
    n = 0.5 * ndim
    return (np.pi) ** (n) / fn_gamma( n + 1 ) * (r ** ndim)


def testing():

    grid_def = dict(
        logt=[6.0, 10.13, 0.05],
        avs=[0., 8., 0.5],
        rvs=[2.1, 6.1, 1.],
        fbumps=[0., 1., 0.2],
        Z=0.019
    )

    def to_samples(d):
        r = np.atleast_1d(d)
        N = len(r)
        if N == 1:
            return r
        elif N == 3:
            return np.arange(d[0], d[1], d[2], dtype=float)
        elif N > 3:
            # assume already a sample
            return r
        else:
            raise ValueError('Do not understand defintion format')

    def make_fake_grid():
        names = grid_def.keys()
        q = [ to_samples(dk) for dk in grid_def.values() ]
        pts = np.nditer(np.ix_(*q))
        return names, pts

    n, q = make_fake_grid()
    qq = np.asarray([ k for k in q])
    t = KDTreeDensityEstimator(qq, normalize=True)

    print t[qq]
