"""
Core of likelihood computations

This package implements many likelihoods based on the common chi2 statistics

python/numpy version.

N_chi2            Computes a classic (non-reduced) chi2 with normal errors
SN_chi2           Computes a chi2 with split errors pr non-symmetric errors
N_logLikelihood   Computes a normal likelihood (default, symmetric errors)
SN_logLikelihood  Computes a Split Normal likelihood (asymmetric errors)
getNorm_lnP       Compute the norm of a log-likelihood (overflow robust)
"""
import numpy as np
#import numexpr

def N_chi2(flux, fluxerr, fluxmod, mask=None):
    """ compute the non-reduced chi2 between data with uncertainties and
        perfectly known models

    Parameters
    ----------
    flux:    np.ndarray[float, ndim=1]
        array of fluxes

    fluxerr: np.ndarray[float, ndim=1]
        array of flux errors

    fluxmod: np.ndarray[float, ndim=2]
        array of modeled fluxes (nfilters , nmodels)

    mask:    np.ndarray[bool, ndim=1]
        mask array to apply during the calculations mask.shape = flux.shape

    Returns
    -------
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

    Parameters
    ----------
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

    Returns
    -------
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


def N_chi2_NM(flux, fluxmod, fluxerr, fluxbias, mask=None):
    """ compute the non-reduced chi2 between data and model taking into account
    the noise model computed from ASTs.

    Parameters
    ----------
    flux:    np.ndarray[float, ndim=1]
        array of fluxes

    fluxmod: np.ndarray[float, ndim=2]
        array of modeled fluxes (nfilters , nmodels)

    fluxerr: np.ndarray[float, ndim=2]
        array of dispersions per model (nfilters , nmodels)

    fluxbias: np.ndarray[float, ndim=2]
        array of fluxes biases per model (nfilters , nmodels)

    mask:    np.ndarray[bool, ndim=1]
        mask array to apply during the calculations mask.shape = flux.shape

    Returns
    -------
    chi2:    np.ndarray[float, ndim=1]
        array of chi2 values (nmodels)
    """
    if mask is None:
        temp = flux[None, :] - (fluxmod + fluxbias)
        #tflux = flux[None, :]
        #temp = numexpr.evaluate('tflux - fluxmod + fluxbias')
        _e = fluxerr
    else:
        _m = ~mask.astype(bool)
        temp = flux[_m]
        temp = temp[None, :] - (fluxmod[:, _m] + fluxbias[:,_m])
        _e = fluxerr[:,_m]

    for j in range(np.shape(_e)[1]):
        temp[:,j] /= _e[:,j]

    return (temp ** 2).sum(axis=1)


def N_covar_chi2(flux, fluxmod, fluxbias, icov_diag, icov_offdiag):
    """ compute the non-reduced chi2 between data and model using
    the full covariance matrix information computed from ASTs.

    Parameters
    ----------
    flux:    np.ndarray[float, ndim=1]
        array of fluxes

    fluxmod: np.ndarray[float, ndim=2]
        array of modeled fluxes (nfilters , nmodels)

    fluxbias: np.ndarray[float, ndim=2]
        array of fluxes biases per model (nfilters , nmodels)

    icov_diag: np.ndarray[float, ndim=2]
        array giving the diagnonal terms of the covariance matrix inverse

    icov_offdiag: np.ndarray[float, ndim=2]
        array giving the off diagnonal terms of the covariance matrix inverse

    Returns
    -------
    chi2:    np.ndarray[float, ndim=1]
        array of chi2 values (nmodels)

    Note
    ----
    Mask removed as it cannot be used with a precomputed inverse
    covariance matrix.  (KDG 29 Jan 2016)
    """
    # get the number of models and filters
    n_models, n_filters = fluxmod.shape
    
    # compute the difference in fluxes
    #    take into account the bias term from the AST results
    fluxdiff = flux[None, :] - (fluxmod + fluxbias)

    chisqr = icov_diag[:,n_filters-1]*np.square(fluxdiff[:,n_filters-1])
            
    m = 0
    for k in range(n_filters-1):
        tchisqr = icov_diag[:,k]*fluxdiff[:,k]
        for l in range(k+1,n_filters):
            tchisqr += (2.0*icov_offdiag[:,m])*fluxdiff[:,l]
            m += 1
        chisqr += tchisqr*fluxdiff[:,k]
        
    return chisqr


def SN_logLikelihood(flux, fluxerr_m, fluxerr_p, fluxmod, mask=None, 
                     lnp_threshold=1000.):
    """ Compute the log of the chi2 likelihood between data with 
    uncertainties and perfectly known models
    with split errors (or non symmetric errors)

    Parameters
    ----------
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

    lnp_threshold:  float
        cut the values outside -x, x in lnp

    Returns
    -------
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
    """ Compute the log of the chi2 likelihood between data with 
    uncertainties and perfectly known models

    Parameters
    ----------
    flux:    np.ndarray[float, ndim=1]
        array of fluxes

    fluxerr: np.ndarray[float, ndim=1]
        array of flux errors

    fluxmod: np.ndarray[float, ndim=2]
        array of modeled fluxes (Nfilters , Nmodels)

    mask:    np.ndarray[bool, ndim=1]
        mask array to apply during the calculations mask.shape = flux.shape

    lnp_threshold:  float
        cut the values outside -x, x in lnp

    Returns
    -------
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


def N_logLikelihood_NM(flux, fluxmod, fluxerr, fluxbias, mask=None,
                       lnp_threshold=1000.):
    """ Computes the log of the chi2 likelihood between data and model taking
    into account the noise model.

    Parameters
    ----------
    flux: np.ndarray[float, ndim=1]
        array of fluxes

    fluxmod: np.ndarray[float, ndim=2]
        array of modeled fluxes (nfilters , nmodels)

    fluxerr: np.ndarray[float, ndim=2]
        array of modeled fluxes (nfilters , nmodels)

    fluxbias: np.ndarray[float, ndim=2]
        array of modeled fluxes (nfilters , nmodels)

    mask:    np.ndarray[bool, ndim=1]
        mask array to apply during the calculations mask.shape = flux.shape

    lnp_threshold:  float
        cut the values outside -x, x in lnp

    Returns
    -------
    (lnp, chi2)
    lnP:    np.ndarray[float, ndim=1]
            array of ln(P) values (Nmodels)
    chi2:    np.ndarray[float, ndim=1]
            array of chi-squared values (Nmodels)

    .. math::

        P = \\frac{1}{\\sqrt{2pi} * \\sigma} \\times exp ( - \\chi^2 / 2 )
        and
        \\chi ^ 2 = \\sum_{k} (flux_{obs,k} - flux_{pred,k} - \mu_k) ^ 2 /
                    \sigma^2_{pred,k}
    """
    ni, nj = np.shape(fluxmod)

    #compute the quality factor
    # lnQ = -0.5 * nj *  ln( 2 * pi) - sum_j {ln( err[j] ) }
    temp = 0.5 * np.log( 2. * np.pi )
    if mask is None:
        temp1 = fluxerr
    else:
        _m = ~mask.astype(bool)
        temp1 = fluxerr[:,_m]

    # By definition errors computed from ASTs are positive.
    n = np.shape(temp1)[1]  
    # lnQ different for each model
    lnQ = n * temp + np.sum(np.log(temp1),axis=1)
    #lnQ is to be used * -1

    #compute the lnp = -lnQ - 0.5 * chi2
    _chi2 = N_chi2_NM(flux, fluxmod, fluxerr, fluxbias, mask=mask)

    lnP = -lnQ - 0.5 * _chi2

    return (lnP, _chi2)


def N_covar_logLikelihood(flux, fluxmod, fluxbias,
                          q_norm, icov_diag, icov_offdiag,
                          lnp_threshold=1000.):
    """ Computes the log of the chi2 likelihood between data and model taking
    into account the noise model.

    Parameters
    ----------
    flux: np.ndarray[float, ndim=1]
        array of fluxes

    fluxmod: np.ndarray[float, ndim=2]
        array of modeled fluxes (nfilters , nmodels)

    fluxbias: np.ndarray[float, ndim=2]
        array of modeled fluxes (nfilters , nmodels)

    q_norm: np.ndarray[float, ndim=2]
        array givign the q normalization of the likelihood
        q_norm = ln(1./Q) where Q = det(cov matrix)

    icov_diag: np.ndarray[float, ndim=2]
        array giving the diagnonal terms of the covariance matrix inverse

    icov_offdiag: np.ndarray[float, ndim=2]
        array giving the off diagnonal terms of the covariance matrix inverse

    lnp_threshold:  float
        cut the values outside -x, x in lnp

    Returns
    -------
    (lnp, chi2)
    lnP:    np.ndarray[float, ndim=1]
            array of ln(P) values (Nmodels)
    chi2:    np.ndarray[float, ndim=1]
            array of chi-squared values (Nmodels)
    """
    n_models, n_filters = np.shape(fluxmod)

    #compute the pi normalization term
    n_good_filters = n_filters
        
    pi_term = -0.5*n_good_filters*np.log(2.0*np.pi)

    # get the chi2 value
    _chi2 = N_covar_chi2(flux, fluxmod, fluxbias,
                         icov_diag, icov_offdiag)

    # compute the lnp = pi_term + q_norm - 0.5*chi2
    lnP = pi_term + q_norm - (0.5*_chi2)

    return (lnP, _chi2)


def N_covar_logLikelihood_cholesky(flux, inv_cholesky_covar, lnQ,
                                   bias, fluxmod):
    """
    Compute the log-likelihood given data, a covariance matrix,
    and a bias term. Very slow.

    Parameters
    ----------
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

    Parameters
    ----------
    lnP:    np.ndarray[float, ndim=1]
        array of ln(P) values (Nmodels)

    Returns
    -------
    norm:   float
        sum_i { exp( lnP[i] ) }
    """
    K = max(lnP)
    norm = np.exp(K) * np.sum(np.exp(lnP - K))
    return norm
