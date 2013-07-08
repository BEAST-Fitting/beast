"""
Core of percentiles and statistics computations

This package implements common statistics often used on pdfs

weighted_percentile     Compute weighted percentiles.
expectation             Return the expectation value sum(p[i] * q[i] ) / sum(p[i])

Notes:
    Interestingly the cython implementation of expectation is almost as fast as numpy.mean (13% slower)

:author:        MF
:last update:   Tue Jun 11 12:55:16 PDT 2013
"""

from __future__ import division
import cython
cimport cython

import numpy as np
cimport numpy as np

DTYPE = np.float64
ctypedef np.float_t DTYPE_t

cdef extern from "numpy/npy_math.h":
    bint npy_isnan(DTYPE_t x)


@cython.boundscheck(False)
@cython.wraparound(False)
def weighted_percentile(np.ndarray[DTYPE_t, ndim=1] data,
                        np.ndarray[DTYPE_t, ndim=1] wt,
                        np.ndarray[DTYPE_t, ndim=1] percentiles):
    """Compute weighted percentiles.
    If the weights are equal, this is the same as normal percentiles.
    Elements of the data and wt arrays correspond to each other and must have
    equal length.

    INPUTS:
    -------
    data: ndarray[float, ndim=1]
        data points
    wt: ndarray[float, ndim=1]
        Weights of each point in data
        All the weights must be non-negative and the sum must be
        greater than zero.
    percentiles: ndarray[float, ndim=1]
        percentiles to use.
        (expecting values in the range of 0-1 rather than 0-100.)
    OUTPUTS:
    -------
    the weighted percentiles of the data.
    """
    cdef int n = len(data)
    cdef int pp = len(percentiles)
    cdef int ik, k
    cdef DTYPE_t norm = 1.0
    cdef DTYPE_t tmp = 0.0
    cdef DTYPE_t sw = 0.0

    cdef np.ndarray[DTYPE_t, ndim=1] sd = np.empty(n, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] aw = np.empty(n, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] _w = np.empty(n, dtype=DTYPE)
    cdef np.ndarray[long, ndim=1] i = np.argsort(data)

    cdef np.ndarray[DTYPE_t, ndim=1] o = np.empty(pp, dtype=DTYPE)

    for k in range(n):
        # make cumulative sum, sorted data arrays and normalized weights
        ## normalized weights are cumsum(sw) - 0.5 * sw, sw: data sorted weights
        ik = i[k]
        sw = wt[ik]
        sd[k] = data[ik]
        aw[k] = tmp
        tmp += sw
        _w[k] = aw[k] - 0.5 * sw

    cdef np.ndarray[long, ndim=1] spots = np.searchsorted(_w, percentiles * aw[n - 1])

    norm = 1. / aw[n - 1]
    for (pk, s, p) in zip(range(pp), spots, percentiles):
        if s == 0:
            o[pk] = sd[0]
        elif s == n:
            o[pk] = sd[n - 1]
        else:
            tmp = 1. / ((_w[s] - _w[s - 1]) * norm)
            f1 = (_w[s] * norm - p) * tmp
            f2 = (p - _w[s - 1] * norm) * tmp
            o[pk] = (sd[s - 1] * f1 + sd[s] * f2)
    return o

@cython.boundscheck(False)
@cython.wraparound(False)
def expectation(np.ndarray[DTYPE_t, ndim=1] q, np.ndarray[DTYPE_t, ndim=1] wt):
    """ return sum(p[i] * q[i] ) / sum(p[i])
    INPUTS
    ------
    q: ndarray[float, ndim=1]
        data from which we compute the expectation value
    wt: ndarray[float, ndim=1]
        weights associated to each data point

    OUTPUTS
    -------
    e: float
        expectation value
    """
    cdef DTYPE_t e = 0.0
    cdef DTYPE_t c = 0.0
    cdef int n=len(q)
    cdef int i

    for i in range(n):
        e += wt[i] * q[i]
        c += wt[i]

    return e / c

@cython.boundscheck(False)
@cython.wraparound(False)
def ptp(np.ndarray[DTYPE_t, ndim=1] q):
    """
    Range of values (maximum - minimum) along an axis.
    The name of the function comes from the acronym for 'peak to peak'.
    Equivalent to numpy.ptp

    Parameters
    ----------
    q : ndarray[float, ndim=1]
        Input values.

    Returns
    -------
    ptp : float
        max - min value of the array
    """
    cdef DTYPE_t m = q[0]
    cdef DTYPE_t M = q[0]
    cdef int n=len(q)
    cdef int i

    for i in range(n):
        if q[i] < m:
            m = q[i]
        elif q[i] > M:
            M = q[i]

    return M - m
