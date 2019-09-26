"""
Core of likelihood computations

This package implements many likelihoods based on the common chi2 statistics

python/numpy version.

N_logLikelihood   Computes a normal likelihood (default, symmetric errors)
SN_logLikelihood  Computes a Split Normal likelihood (asymmetric errors)
getNorm_lnP       Compute the norm of a log-likelihood (overflow robust)
"""
import numpy as np

__all__ = [
    "N_chi2_NM",
    "N_covar_chi2",
    "N_logLikelihood_NM",
    "N_covar_logLikelihood",
    "N_covar_logLikelihood_cholesky",
    "getNorm_lnP",
]


def N_chi2_NM(flux, fluxmod_wbias, ivar, mask=None):
    """ compute the non-reduced chi2 between data and model taking into account
    the noise model computed from ASTs.

    Parameters
    ----------
    flux:    np.ndarray[float, ndim=1]
        array of fluxes

    fluxmod_wbias: np.ndarray[float, ndim=2]
        array of modeled fluxes + ast-derived biases (nfilters , nmodels)

    ivar: np.ndarray[float, ndim=2]
        array of ast-derived inverse variances (nfilters , nmodels)

    mask:    np.ndarray[bool, ndim=1]
        mask array to apply during the calculations mask.shape = flux.shape

    Returns
    -------
    chi2:    np.ndarray[float, ndim=1]
        array of chi2 values (nmodels)
    """
    if (mask is None) or np.all(mask is False):
        temp = flux - fluxmod_wbias
        _ie = ivar
    else:
        _m = ~mask.astype(bool)
        temp = flux[_m] - fluxmod_wbias[:, _m]
        _ie = ivar[:, _m]

    return np.einsum("ij,ij,ij->i", temp, temp, _ie)


def N_covar_chi2(flux, fluxmod_wbias, icov_diag, two_icov_offdiag):
    """ compute the non-reduced chi2 between data and model using
    the full covariance matrix information computed from ASTs.

    Parameters
    ----------
    flux:    np.ndarray[float, ndim=1]
        array of fluxes

    fluxmod_wbias: np.ndarray[float, ndim=2]
        array of modeled fluxes (nfilters , nmodels)

    icov_diag: np.ndarray[float, ndim=2]
        array giving the diagnonal terms of the covariance matrix inverse

    two_icov_offdiag: np.ndarray[float, ndim=2]
        array giving 2x the off diagonal terms of the covariance matrix inverse

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
    n_models, n_filters = fluxmod_wbias.shape

    # compute the difference in fluxes
    fluxdiff = flux[None, :] - fluxmod_wbias

    # diagonal terms
    chisqr = np.einsum("ij,ij,ij->i", fluxdiff, fluxdiff, icov_diag)

    # off-diagonal terms
    m_start = 0
    for k in range(n_filters - 1):
        m_end = m_start + n_filters - k - 1
        tchisqr = np.einsum(
            "ij,ij->i", two_icov_offdiag[:, m_start:m_end], fluxdiff[:, k + 1 :]
        )
        tchisqr *= fluxdiff[:, k]
        chisqr += tchisqr
        m_start = m_end

    return chisqr


def N_logLikelihood_NM(flux, fluxmod_wbias, ivar, mask=None, lnp_threshold=1000.0):
    r""" Computes the log of the chi2 likelihood between data and model taking
    into account the noise model.

    Parameters
    ----------
    flux: np.ndarray[float, ndim=1]
        array of fluxes

    fluxmod_wbias: np.ndarray[float, ndim=2]
        array of modeled fluxes + ast-derived biases (nfilters, nmodels)

    ivar: np.ndarray[float, ndim=2]
        array of ast-derived inverse variances (nfilters , nmodels)

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
    ni, nj = np.shape(fluxmod_wbias)

    # compute the quality factor
    # lnQ = -0.5 * nj *  ln( 2 * pi) - sum_j {ln( err[j] ) }
    temp = 0.5 * np.log(2.0 * np.pi)
    if mask is None:
        temp1 = ivar  # fluxerr
    else:
        _m = ~mask.astype(bool)
        temp1 = ivar[:, _m]  # fluxerr[:,_m]

    # By definition errors computed from ASTs are positive.
    n = np.shape(temp1)[1]
    # lnQ different for each model
    lnQ = n * temp - 0.5 * np.sum(np.log(temp1), axis=1)
    # lnQ is to be used * -1

    # compute the lnp = -lnQ - 0.5 * chi2
    _chi2 = N_chi2_NM(flux, fluxmod_wbias, ivar, mask=mask)

    lnP = -lnQ - 0.5 * _chi2

    return (lnP, _chi2)


def N_covar_logLikelihood(
    flux, fluxmod_wbias, q_norm, icov_diag, two_icov_offdiag, lnp_threshold=1000.0
):
    """ Computes the log of the chi2 likelihood between data and model taking
    into account the noise model.

    Parameters
    ----------
    flux: np.ndarray[float, ndim=1]
        array of fluxes

    fluxmod_wbias: np.ndarray[float, ndim=2]
        array of modeled fluxes + ast-derived biases (nfilters , nmodels)

    q_norm: np.ndarray[float, ndim=2]
        array givign the q normalization of the likelihood
        q_norm = ln(1./Q) where Q = det(cov matrix)

    icov_diag: np.ndarray[float, ndim=2]
        array giving the diagnonal terms of the covariance matrix inverse

    two_icov_offdiag: np.ndarray[float, ndim=2]
        array giving 2x the off diagonal terms of the covariance matrix inverse

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
    n_models, n_filters = np.shape(fluxmod_wbias)

    # compute the pi normalization term
    n_good_filters = n_filters

    pi_term = -0.5 * n_good_filters * np.log(2.0 * np.pi)

    # get the chi2 value
    _chi2 = N_covar_chi2(flux, fluxmod_wbias, icov_diag, two_icov_offdiag)

    # compute the lnp = pi_term + q_norm - 0.5*chi2
    lnP = pi_term + q_norm - (0.5 * _chi2)

    return (lnP, _chi2)


def N_covar_logLikelihood_cholesky(flux, inv_cholesky_covar, lnQ, bias, fluxmod):
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
    lnP = -0.5 * np.sum(lnP, axis=1)
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
