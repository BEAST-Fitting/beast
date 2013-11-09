""" This package implements many likelihoods based on the common chi2 statistics
It is meant to be fast even on large datasets, this is why it is in cython

Tests have shown that speeds of weave and cython codes are similar if
boundaries and negative indices checks are disabled.

c_N_chi2            Computes a classic (non-reduced) chi2 with normal errors
c_SN_chi2           Computes a chi2 with split errors pr non-symmetric errors
c_N_logLikelihood   Computes a normal likelihood (default, symmetric errors)
c_SN_logLikelihood  Computes a Split Normal likelihood (asymmetric errors)
c_getNorm_lnP       Compute the norm of a log-likelihood (overflow robust)

:author:        MF
:last update:   Mon May 20 11:43:27 PDT 2013
"""
from __future__ import division
import numpy as np
cimport numpy as np  # NOQA

from libc.math cimport log, exp

DTYPE = np.float64
ctypedef np.float_t DTYPE_t

cdef extern from "numpy/npy_math.h":
    bint npy_isnan(DTYPE_t x)


cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
def c_N_chi2(np.ndarray[DTYPE_t, ndim=1] flux,
           np.ndarray[DTYPE_t, ndim=1] fluxerr,
           np.ndarray[DTYPE_t, ndim=2] fluxmod,
           np.ndarray[int, ndim=1] mask):
    """ compute the non-reduced chi2 between data with uncertainties and
        perfectly known models

    INPUTS:
        flux:    np.ndarray[float, ndim=1]
            array of fluxes
        fluxerr: np.ndarray[float, ndim=1]
            array of flux errors
        fluxmod: np.ndarray[float, ndim=2]
            array of modeled fluxes (nfilters , nmodels)
        mask:    np.ndarray[bool, ndim=1]
            mask array to apply during the calculations mask.shape = flux.shape

    OUTPUTS:
        chi2:    np.ndarray[float, ndim=1]
            array of chi2 values (nmodels)
    """
    cdef int ni
    cdef int nj
    cdef int i
    cdef int j
    ni, nj = np.shape(fluxmod)
    cdef np.ndarray[DTYPE_t, ndim=1] out = np.empty(ni, dtype=DTYPE)

    cdef DTYPE_t temp = 0.0
    cdef DTYPE_t err = 0.0
    for i in range(ni):
        for j in range(nj):
            if mask[j] == 0:
                temp = ( flux[j] - fluxmod[i, j] )
                if fluxerr[j] > 0:
                    # if err>0 then divide else considere err=1.
                    temp /= fluxerr[j]
                temp *= temp
                out[i] += temp
    return out

# try:
#     from cython.parallel import prange
#
#     @cython.boundscheck(False)
#     @cython.wraparound(False)
#     def c_N_chi2_omp(np.ndarray[DTYPE_t, ndim=1] flux,
#                     np.ndarray[DTYPE_t, ndim=1] fluxerr,
#                     np.ndarray[DTYPE_t, ndim=2] fluxmod,
#                     np.ndarray[int, ndim=1] mask,
#                     int num_threads=0):
#         """ compute the non-reduced chi2 between data with uncertainties and
#             perfectly known models
#
#         INPUTS:
#             flux:    np.ndarray[float, ndim=1]
#                 array of fluxes
#             fluxerr: np.ndarray[float, ndim=1]
#                 array of flux errors
#             fluxmod: np.ndarray[float, ndim=2]
#                 array of modeled fluxes (nfilters , nmodels)
#             mask:    np.ndarray[bool, ndim=1]
#                 mask array to apply during the calculations mask.shape = flux.shape
#
#         OUTPUTS:
#             chi2:    np.ndarray[float, ndim=1]
#                 array of chi2 values (nmodels)
#         """
#         cdef int ni
#         cdef int nj
#         cdef int i
#         cdef int j
#         ni, nj = np.shape(fluxmod)
#         cdef np.ndarray[DTYPE_t, ndim=1] out = np.empty(ni, dtype=DTYPE)
#
#         cdef DTYPE_t temp = 0.0
#         cdef DTYPE_t err = 0.0
#
#         #normal behavior
#         if num_threads == 0:
#             for i in range(ni):
#                 for j in range(nj):
#                     if mask[j] == 0:
#                         temp = ( flux[j] - fluxmod[i, j] )
#                         if fluxerr[j] > 0:
#                             # if err>0 then divide else considere err=1.
#                             temp /= fluxerr[j]
#                         temp *= temp
#                         out[i] += temp
#         elif num_threads < 0:
#             for i in prange(ni, schedule='guided', nogil=True):
#                 for j in range(nj):
#                     if mask[j] == 0:
#                         temp = ( flux[j] - fluxmod[i, j] )
#                         if fluxerr[j] > 0:
#                             # if err>0 then divide else considere err=1.
#                             temp /= fluxerr[j]
#                         temp = temp * temp
#                         out[i] += temp
#         else:
#             for i in prange(ni, schedule='guided', nogil=True, num_threads=num_threads):
#                 for j in range(nj):
#                     if mask[j] == 0:
#                         temp = ( flux[j] - fluxmod[i, j] )
#                         if fluxerr[j] > 0:
#                             # if err>0 then divide else considere err=1.
#                             temp /= fluxerr[j]
#                         temp = temp * temp
#                         out[i] += temp
#         return out
# except ImportError:
#     pass
#

@cython.boundscheck(False)
@cython.wraparound(False)
def c_SN_chi2(np.ndarray[DTYPE_t, ndim=1] flux,
              np.ndarray[DTYPE_t, ndim=1] fluxerr_m,
              np.ndarray[DTYPE_t, ndim=1] fluxerr_p,
              np.ndarray[DTYPE_t, ndim=2] fluxmod,
              np.ndarray[int, ndim=1] mask):
    """ compute the non-reduced chi2 between data with uncertainties and
        perfectly known models

    INPUTS:
        flux:    np.ndarray[float, ndim=1]
            array of fluxes
        fluxerr_m: np.ndarray[float, ndim=1]
            array of flux errors on the left side (<= flux)
        fluxerr_p: np.ndarray[float, ndim=1]
            array of flux errors on the right side (>= flux)
        fluxmod: np.ndarray[float, ndim=2]
            array of modeled fluxes (nfilters , nmodels)
        mask:    np.ndarray[bool, ndim=1]
            mask array to apply during the calculations mask.shape = flux.shape

    OUTPUTS:
        chi2:    np.ndarray[float, ndim=1]
            array of chi2 values (nmodels)
    """
    cdef int ni
    cdef int nj
    cdef int i
    cdef int j
    ni, nj = np.shape(fluxmod)
    cdef np.ndarray[DTYPE_t, ndim=1] out = np.empty(ni, dtype=DTYPE)

    cdef DTYPE_t temp = 0.0
    cdef DTYPE_t err = 0.0
    for i in range(ni):
        for j in range(nj):
            if mask[j] == 0:
                temp = ( flux[j] - fluxmod[i, j] )
                # temp > 0 == left side
                if (temp > 0.)  & (fluxerr_m[j] > 0.):
                    # if err>0 then divide else considere err=1.
                    temp /= fluxerr_m[j]
                elif (temp < 0.)  & (fluxerr_p[j] > 0.):
                    temp /= fluxerr_p[j]
                temp *= temp
                out[i] += temp
    return out


@cython.boundscheck(False)
@cython.wraparound(False)
def c_SN_logLikelihood(np.ndarray[DTYPE_t, ndim=1] flux,
                       np.ndarray[DTYPE_t, ndim=1] fluxerr_m,
                       np.ndarray[DTYPE_t, ndim=1] fluxerr_p,
                       np.ndarray[DTYPE_t, ndim=2] fluxmod,
                       np.ndarray[int, ndim=1] mask,
                       DTYPE_t lnp_threshold=1000.):
    """ Compute the log of the chi2 likelihood between data with uncertainties and perfectly known models
    with split errors (or non symmetric errors)

    INPUTS:
        flux:    np.ndarray[float, ndim=1]
            array of fluxes
        fluxerr_m: np.ndarray[float, ndim=1]
            array of flux errors on the left side (<= flux)
        fluxerr_p: np.ndarray[float, ndim=1]
            array of flux errors on the right side (>= flux)
        fluxmod: np.ndarray[float, ndim=2]
            array of modeled fluxes (Nfilters , Nmodels)
        mask:    np.ndarray[bool, ndim=1]
            mask array to apply during the calculations mask.shape = flux.shape

    KEYWORDS:
        lnp_threshold:  float
            cut the values outside -x, x in lnp

    OUTPUTS:
        lnP:    np.ndarray[float, ndim=1]
            array of ln(P) values (Nmodels)

        with P = 1/[sqrt(pi/2) * (sig_p + sig_m)**2 ] * exp ( - 0.5 * chi2 )
            and chi2 uses sig_p or sig_m if fluxmod > flux (resp. <)
    """
    cdef DTYPE_t temp = 0.
    cdef DTYPE_t temp1 = 0.
    cdef DTYPE_t lnQ = 0.
    cdef int i
    cdef int j
    cdef int ni
    cdef int nj
    ni, nj = np.shape(fluxmod)
    cdef np.ndarray[DTYPE_t, ndim=1] lnP = np.zeros(ni, dtype=DTYPE)

    #compute the quality factor
    # lnQ = -0.5 * nj *  log( pi/2 ) - sum_j {log( err[j] ) }
    temp = 0.5 * log( 0.5 * np.pi )
    for j in range(nj):
        if mask[j] == 0:
            lnQ += temp
            temp1 = fluxerr_m[j] + fluxerr_p[j]
            if temp1 > 0.:
                lnQ += log(temp1)
    #lnQ is to be used * -1

    #compute the lnp = -lnQ - 0.5 * chi2
    # (Almost the c_chi2 code)
    for i in range(ni):
        for j in range(nj):
            if mask[j] == 0:
                temp = ( flux[j] - fluxmod[i, j] )
                # temp > 0 == left side
                if (temp > 0.)  & (fluxerr_m[j] > 0.):
                    # if err>0 then divide else considere err=1.
                    temp /= fluxerr_m[j]
                elif (temp < 0.)  & (fluxerr_p[j] > 0.):
                    temp /= fluxerr_p[j]
                temp *= temp
                lnP[i] += temp
        lnP[i] = -lnQ - 0.5 * lnP[i]
        if lnP[i] < -lnp_threshold:
            lnP[i] = -lnp_threshold

    return lnP


@cython.boundscheck(False)
@cython.wraparound(False)
def c_N_logLikelihood(np.ndarray[DTYPE_t, ndim=1] flux,
                      np.ndarray[DTYPE_t, ndim=1] fluxerr,
                      np.ndarray[DTYPE_t, ndim=2] fluxmod,
                      np.ndarray[int, ndim=1] mask,
                      DTYPE_t lnp_threshold=1000.):
    """ Compute the log of the chi2 likelihood between data with uncertainties and perfectly known models

    INPUTS:
        flux:    np.ndarray[float, ndim=1]
            array of fluxes
        fluxerr: np.ndarray[float, ndim=1]
            array of flux errors
        fluxmod: np.ndarray[float, ndim=2]
            array of modeled fluxes (Nfilters , Nmodels)
        mask:    np.ndarray[bool, ndim=1]
            mask array to apply during the calculations mask.shape = flux.shape

    KEYWORDS:
        lnp_threshold:  float
            cut the values outside -x, x in lnp

    OUTPUTS:
        lnP:    np.ndarray[float, ndim=1]
            array of ln(P) values (Nmodels)

        with P = 1/[sqrt(2pi) * sig**2 ] * exp ( - 0.5 * chi2 )
    """
    cdef DTYPE_t temp = 0.
    cdef DTYPE_t lnQ = 0.
    cdef int i
    cdef int j
    cdef int ni
    cdef int nj
    ni, nj = np.shape(fluxmod)
    cdef np.ndarray[DTYPE_t, ndim=1] lnP = np.zeros(ni, dtype=DTYPE)

    #compute the quality factor
    # lnQ = -0.5 * nj *  ln( 2 * pi) - sum_j {ln( err[j] ) }
    temp = 0.5 * log( 2. * np.pi )
    for j in range(nj):
        if mask[j] == 0:
            lnQ += temp
            if fluxerr[j] > 0.:
                lnQ += log(fluxerr[j])
    #lnQ is to be used * -1

    #compute the lnp = -lnQ - 0.5 * chi2
    # (Almost the c_chi2 code)
    for i in range(ni):
        for j in range(nj):
            if mask[j] == 0:
                temp = ( flux[j] - fluxmod[i, j] )
                if fluxerr[j] > 0:
                    # if err>0 then divide else considere err=1.
                    temp /= fluxerr[j]
                temp *= temp
                lnP[i] += temp
        lnP[i] = -lnQ - 0.5 * lnP[i]
        if lnP[i] < -lnp_threshold:
            lnP[i] = -lnp_threshold

    return lnP

@cython.boundscheck(False)
@cython.wraparound(False)
def c_N_covar_logLikelihood(np.ndarray[DTYPE_t, ndim=1] flux, 
                            np.ndarray[DTYPE_t, ndim=3] inv_cholesky_covar, 
                            np.ndarray[DTYPE_t, ndim=1] lnQ, 
                            np.ndarray[DTYPE_t, ndim=2] bias, 
                            np.ndarray[DTYPE_t, ndim=2] fluxmod):
    """
    Compute the log-likelihood given data, a covariance matrix, 
    and a bias term.
    INPUTS:
        flux:    np.ndarray([float, ndim=1])
             Measured fluxes
        inv_cholesky_covar:   np.ndarray([float, ndim=3])
             The inverses of the (lower triangular) Cholesky matrices of the covariance matrices
        lnQ:     np.ndarray([float, ndim=1])
             Logarithm of the determinants of the covariance matrices
        
    """
    cdef np.intp_t n_seds, len_seds, sed_i, len_i, len_j
    n_seds, len_seds = fluxmod.shape[0], fluxmod.shape[1]
    cdef np.ndarray[DTYPE_t, ndim=1] lnP = np.zeros([n_seds], dtype=np.double)
    cdef np.ndarray[DTYPE_t, ndim=1] off = np.zeros([len_seds], dtype=np.double)
    cdef DTYPE_t temp = 0

    
    for sed_i in range(n_seds):
        for len_i in range(len_seds):
            off[len_i] = flux[len_i] - (fluxmod[sed_i, len_i] + bias[sed_i, len_i]) 
        for len_i in range(len_seds):
            temp = 0
            for len_j in range(len_i+1):
                temp += inv_cholesky_covar[sed_i, len_i, len_j]*off[len_j]
            lnP[sed_i] += temp*temp
        lnP[sed_i] *= -0.5
        lnP[sed_i] -= lnQ[sed_i]
    return lnP

@cython.boundscheck(False)
@cython.wraparound(False)
def c_getNorm_lnP(np.ndarray[DTYPE_t, ndim=1] lnP):
    """ Compute the norm of a log-likelihood

    INPUTS:
        lnP:    np.ndarray[float, ndim=1]
            array of ln(P) values (Nmodels)

    OUTPUTS:
        norm:   float
            sum_i { exp( lnP[i] ) }
    """
    cdef DTYPE_t norm = 0.0
    cdef DTYPE_t K = 0.0
    cdef int ni = len(lnP)
    cdef int i

    # To make sure we don't have overflows, we normalize the sum by its max
    K = max(lnP)
    for i in range(ni):
        norm += exp(lnP[i] - K)

    norm *= exp(K)

    return norm

