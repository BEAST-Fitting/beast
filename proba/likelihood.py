"""
Core of likelihood computations

This package implements many likelihoods based on the common chi2 statistics

python/numpy version.

N_chi2            Computes a classic (non-reduced) chi2 with normal errors
SN_chi2           Computes a chi2 with split errors pr non-symmetric errors
N_logLikelihood   Computes a normal likelihood (default, symmetric errors)
SN_logLikelihood  Computes a Split Normal likelihood (asymmetric errors)
getNorm_lnP       Compute the norm of a log-likelihood (overflow robust)


:author:        MF
:last update:   Mon May 20 11:43:27 PDT 2013
"""
import numpy as np

def N_chi2(flux, fluxerr, fluxmod, mask=None):
    """ compute the non-reduced chi2 between data with uncertainties and
        perfectly known models

    INPUTS:
        flux:    np.ndarray[float, ndim=1]
            array of fluxes
        fluxerr: np.ndarray[float, ndim=1]
            array of flux errors
        fluxmod: np.ndarray[float, ndim=2]
            array of modeled fluxes (nfilters , nmodels)

    KEYWORDS:
        mask:    np.ndarray[bool, ndim=1]
            mask array to apply during the calculations mask.shape = flux.shape

    OUTPUTS:
        chi2:    np.ndarray[float, ndim=1]
            array of chi2 values (nmodels)
    """
    if mask is None:
        temp = flux[None, :] - fluxmod
        _e = fluxerr
    else:
        _m = ~mask.astype(bool)
        temp = flux[_m]
        temp = temp[None, :] - fluxmod[:, _m]
        _e = fluxerr[_m]

    for j in range(len(_e)):
        if _e[j] > 0:
            temp[:,j] /= _e[j]

    return (temp ** 2).sum(axis=1)


def SN_chi2(flux, fluxerr_m, fluxerr_p, fluxmod, mask=None):
    """ compute the non-reduced chi2 between data with asymmetric uncertainties and
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

    KEYWORDS:
        mask:    np.ndarray[bool, ndim=1]
            mask array to apply during the calculations mask.shape = flux.shape

    OUTPUTS:
        chi2:    np.ndarray[float, ndim=1]
            array of chi2 values (nmodels)
    """
    if mask is None:
        temp = flux[None, :] - fluxmod
        _em = fluxerr_m
        _ep = fluxerr_p
    else:
        _m = ~mask.astype(bool)
        _em = fluxerr_m[_m]
        _ep = fluxerr_p[_m]
        temp = flux[_m][None, :] - fluxmod[:, _m]

    for j in range(len(_em)):
        if _em[j] > 0:
            ind0 = np.where(temp[:, j] > 0.)
            temp[ind0, j] /= _em[j]
        if _ep[j] > 0:
            ind0 = np.where(temp[:, j] < 0.)
            temp[ind0, j] /= _ep[j]
    return (temp ** 2).sum(axis=1)


def SN_logLikelihood(flux, fluxerr_m, fluxerr_p, fluxmod, mask=None, lnp_threshold=1000.):
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

    KEYWORDS:
        mask:    np.ndarray[bool, ndim=1]
            mask array to apply during the calculations mask.shape = flux.shape
        lnp_threshold:  float
            cut the values outside -x, x in lnp

    OUTPUTS:
        lnP:    np.ndarray[float, ndim=1]
            array of ln(P) values (Nmodels)

        with P = 1/[sqrt(pi/2) * (sig_p + sig_m)**2 ] * exp ( - 0.5 * chi2 )
            and chi2 uses sig_p or sig_m if fluxmod > flux (resp. <)
    """
    ni, nj = np.shape(fluxmod)

    #compute the quality factor
    # lnQ = -0.5 * nj *  ln( pi/2 ) - sum_j {ln( err[j] ) }
    temp = 0.5 * np.log( 0.5 * np.pi )
    if mask is None:
        temp1 = fluxerr_m + fluxerr_p
    else:
        _m = ~mask.astype(bool)
        temp1 = fluxerr_m[_m] + fluxerr_p[_m]
    n = len(np.where(temp1 > 0)[0])
    lnQ = n * temp + np.sum(np.log(temp1))
    #lnQ is to be used * -1

    #compute the lnp = -lnQ - 0.5 * chi2
    _chi2 = SN_chi2(flux, fluxerr_m, fluxerr_p, fluxmod, mask=mask)

    lnP = -lnQ - 0.5 * _chi2

    return lnP


def N_logLikelihood(flux, fluxerr, fluxmod, mask=None, lnp_threshold=1000.):
    """ Compute the log of the chi2 likelihood between data with uncertainties and perfectly known models

    INPUTS:
        flux:    np.ndarray[float, ndim=1]
            array of fluxes
        fluxerr: np.ndarray[float, ndim=1]
            array of flux errors
        fluxmod: np.ndarray[float, ndim=2]
            array of modeled fluxes (Nfilters , Nmodels)

    KEYWORDS:
        mask:    np.ndarray[bool, ndim=1]
            mask array to apply during the calculations mask.shape = flux.shape
        lnp_threshold:  float
            cut the values outside -x, x in lnp

    OUTPUTS:
        lnP:    np.ndarray[float, ndim=1]
            array of ln(P) values (Nmodels)

        with P = 1/[sqrt(2pi) * sig**2 ] * exp ( - 0.5 * chi2 )
    """
    ni, nj = np.shape(fluxmod)

    #compute the quality factor
    # lnQ = -0.5 * nj *  ln( 2 * pi) - sum_j {ln( err[j] ) }
    temp = 0.5 * np.log( 2. * np.pi )
    if mask is None:
        temp1 = fluxerr
    else:
        _m = ~mask.astype(bool)
        temp1 = fluxerr[_m]
    n = len(np.where(temp1 > 0)[0])
    lnQ = n * temp + np.sum(np.log(temp1))
    #lnQ is to be used * -1

    #compute the lnp = -lnQ - 0.5 * chi2
    _chi2 = N_chi2(flux, fluxerr, fluxmod, mask=mask)

    lnP = -lnQ - 0.5 * _chi2
    #Removing Q factor for comparison with IDL SEDfitter
    #lnP = -0.5 * _chi2

    return lnP

def N_covar_logLikelihood(flux, inv_cholesky_covar, lnQ, bias, fluxmod):
    """
    Compute the log-likelihood given data, a covariance matrix,
    and a bias term. Very slow.
    INPUTS:
        flux:    np.ndarray([float, ndim=1])
             Measured fluxes
        inv_cholesky_covar:   np.ndarray([float, ndim=3])
             The inverses of the Cholesky decompositions of the covariance matrices
        lnQ:     np.ndarray([float, ndim=1])
             Logarithm of the determinants of the covariance matrices

    """
    lnP = np.zeros(fluxmod.shape)
    off = flux - (fluxmod + bias)
    for i in range(fluxmod.shape[0]):
        lnP[i] = np.dot(inv_cholesky_covar[i], off[i])
    lnP *= lnP
    lnP = -0.5*np.sum(lnP, axis=1)
    lnP -= lnQ
    return lnP


def getNorm_lnP(lnP):
    """ Compute the norm of a log-likelihood
    To make sure we don't have overflows, we normalize the sum by its max

    INPUTS:
        lnP:    np.ndarray[float, ndim=1]
            array of ln(P) values (Nmodels)

    OUTPUTS:
        norm:   float
            sum_i { exp( lnP[i] ) }
    """
    K = max(lnP)
    norm = np.exp(K) * np.sum(np.exp(lnP - K))
    return norm
