import numpy as np

__all__ = [
    "_lognorm",
    "_two_lognorm",
    "_exponential",
    "_imf_salpeter",
    "_imf_kroupa",
    "_imf_flat",
]


def _lognorm(in_x, max_pos, sigma=0.5, N=1.0):
    """
    Lognormal distribution

    Parameters
    ----------
    in_x : vector
       x values

    max_pos : float
       Position of the lognormal function's maximum

    sigma : float
       Sigma of the lognormal function

    N : floats
       Multiplicative factor

    Returns
    -------
    lognormal computed on the x grid
    """
    sqrt_2pi = 1.0 / np.sqrt(2 * np.pi)
    mu = np.log(max_pos) + sigma ** 2

    # make it work for a single value or an array
    x = np.atleast_1d(in_x)

    # avoid zero or negative due to log
    pmask = x > 0

    lnorm = np.zeros(len(x))

    log_x = np.log(x[pmask])
    normalization = sqrt_2pi / (x[pmask] * sigma)

    lnorm[pmask] = N * normalization * np.exp(-0.5 * ((log_x - mu) / sigma) ** 2)

    if hasattr(x, "shape"):
        return lnorm
    else:
        return lnorm[0]


def _two_lognorm(xs, max_pos1, max_pos2, sigma1=0.5, sigma2=0.5, N1=1.0, N2=1.0):
    """
    Mixture of 2 lognormal functions

    Parameters
    ----------
    xs : vector
       x values

    max_pos1 : float
       Position of the lognormal function's maximum for component 1
    max_pos2 : float
       Position of the lognormal function's maximum for component 2

    sigma1 : float
       Sigma of the lognormal function for component 1
    sigma2 : float
       Sigma of the lognormal function for component 2

    N1 : floats
       Multiplicative factor for component 1
    N2 : floats
       Multiplicative factor for component 2

    Returns
    -------
    Mixture model: (LOGNORM1 + LOGNORM2) / INTEGRAL(LOGNORM1 + LOGNORM2)
    """
    pointwise = _lognorm(xs, max_pos1, sigma=sigma1, N=N1) + _lognorm(
        xs, max_pos2, sigma=sigma2, N=N2
    )
    normalization = np.trapz(pointwise, x=xs)
    return pointwise / normalization


def _exponential(x, tau=2.0, amp=1.0):
    """
    Exponential distribution

    Parameters
    ----------
    x : vector
       x values
    tau : float
       Decay Rate parameter in exp: N*e^-(x/tau)
    amp : float
       Amplitude for x = 0

    Returns
    -------
    exponential computed on the x grid
    """
    return amp * np.exp(-1.0 * x / tau)


def _imf_flat(x):
    """
    Compute a flat IMF (useful for simulations, not for normal BEAST runs)

    Parameters
    ----------
    x : numpy vector
      masses

    Returns
    -------
    imf : numpy vector
      unformalized IMF
    """
    return 1.0


def _imf_salpeter(x, slope=2.35):
    """
    Compute a Salpeter IMF

    Parameters
    ----------
    x : numpy vector
      masses

    slope : float
        powerlaw slope [default=2.35]

    Returns
    -------
    imf : numpy vector
      unformalized IMF
    """
    return x ** (-1.0 * slope)


def _imf_kroupa(in_x, alpha0=0.3, alpha1=1.3, alpha2=2.3, alpha3=2.3):
    """
    Compute a Kroupa IMF.

    Using equation 2 of Kroupa 2001, MNRAS, 322, 231

    Parameters
    ----------
    in_x : numpy vector
      masses

    alpha0,1,2,3 : float
        slopes between <0.08, 0.08-0.5, 0.5-1.0, >1.0 solar masses
        default = 0.3, 1.3, 2.3, 2.3

    Returns
    -------
    imf : numpy vector
      unformalized IMF
    """
    # allows for single float or an array
    x = np.atleast_1d(in_x)

    # mass break points
    m1 = 0.08
    m2 = 0.5
    m3 = 1.0

    ialpha0 = -1. * alpha0
    ialpha1 = -1. * alpha1
    ialpha2 = -1. * alpha2
    ialpha3 = -1. * alpha3

    imf = np.full((len(x)), 0.0)

    (indxs,) = np.where(x >= m3)
    if len(indxs) > 0:
        imf[indxs] = x[indxs] ** ialpha3

    (indxs,) = np.where((x >= m2) & (x < m3))
    fac1 = (m3 ** ialpha3) / (m3 ** ialpha2)
    if len(indxs) > 0:
        imf[indxs] = (x[indxs] ** ialpha2) * fac1

    (indxs,) = np.where((x >= m1) & (x < m2))
    fac2 = fac1 * (m2 ** ialpha2) / (m2 ** ialpha1)
    if len(indxs) > 0:
        imf[indxs] = (x[indxs] ** ialpha1) * fac2

    (indxs,) = np.where(x < m1)
    fac2 = fac2 * ((m1 ** ialpha1) / (m1 ** ialpha0))
    if len(indxs) > 0:
        imf[indxs] = (x[indxs] ** ialpha0) * fac2

    return imf
