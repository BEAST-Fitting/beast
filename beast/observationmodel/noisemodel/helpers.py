from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from scipy.spatial import cKDTree


def convert_dict_to_structured_ndarray(data):
    """convert_dict_to_structured_ndarray

    Parameters
    ----------

    data: dictionary like object
        data structure which provides iteritems and itervalues

    returns
    -------
    tab: structured ndarray
        structured numpy array
    """
    newdtype = []
    for key, dk in data.items():
        _dk = np.asarray(dk)
        dtype = _dk.dtype
        # unknown type is converted to text
        if dtype.type == np.object_:
            if len(data) == 0:
                longest = 0
            else:
                longest = len(max(_dk, key=len))
                _dk = _dk.astype('|%iS' % longest)
        if _dk.ndim > 1:
            newdtype.append((str(key), _dk.dtype, (_dk.shape[1],)))
        else:
            newdtype.append((str(key), _dk.dtype))
    tab = np.rec.fromarrays(iter(data.values()), dtype=newdtype)
    return tab


def _prepare_x(x):
        """ make sure the data is correctly presented """
        xi = np.asarray(x)
        shape = xi.shape

        if len(shape) < 1:
            raise ValueError('single value passed')

        if len(shape) == 1:
            return xi[..., np.newaxis]
        else:
            return xi


def nearest_neighbors(x, k=10,eps=0.):
    """ kd-tree for quick nearest-neighbor lookup

    ..note::
        using the default leafsize=10

    Parameters
    ----------
    x : array_like, last dimension self.m
        An array of points to query.

    k : integer
        The number of nearest neighbors to return.

    eps : non-negative float
        Return approximate nearest neighbors; the kth returned value
        is guaranteed to be no further than (1+eps) times the
        distance to the real k-th nearest neighbor.

    p : float, 1<=p<=infinity, default=2
        Which Minkowski p-norm to use.
        1 is the sum-of-absolute-values "Manhattan" distance
        2 is the usual Euclidean distance
        infinity is the maximum-coordinate-difference distance

    distance_upper_bound : nonnegative float, default=np.inf
        Return only neighbors within this distance.  This is used to prune
        tree searches, so if you are doing a series of nearest-neighbor
        queries, it may help to supply the distance to the nearest neighbor
        of the most recent point.

    Returns
    -------
    i : ndarray of ints
        The locations of the neighbors in self.data.
        If `x` has shape tuple+(self.m,), then `i` has shape tuple+(k,).
        Missing neighbors are indicated with self.n.
    """
    tree = cKDTree(x)
    d, ind = tree.query(x, k=k, eps=eps)
    return ind


#lazy functions
def toFlux(mag):
    return 10 ** (-0.4 * mag)
