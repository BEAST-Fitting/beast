import numpy as np

__all__ = [
    "_lognorm",
    "_two_lognorm",
    "_exponential",
    "_imf_salpeter",
    "_imf_kroupa",
    "_imf_flat",
]


def _lognorm(x, max_pos, sigma=0.5, N=1.0):
    """
    Lognormal distribution

    Parameters
    ----------
    xs: vector
       x values

    max_pos: float
       Position of the lognormal function's maximum

    sigma: float
       Sigma of the lognormal function

    N: floats
       Multiplicative factor

    Returns
    -------
    lognormal computed on the x grid
    """
    sqrt_2pi = 1.0 / np.sqrt(2 * np.pi)
    mu = np.log(max_pos) + sigma ** 2

    # avoid zero or negative due to log
    (indxs,) = np.where(x > 0)

    lnorm = np.zeros(len(x))

    log_x = np.log(x[indxs])
    normalization = sqrt_2pi / (x[indxs] * sigma)

    lnorm[indxs] = N * normalization * np.exp(-0.5 * ((log_x - mu) / sigma) ** 2)
    return lnorm


def _two_lognorm(xs, max_pos1, max_pos2, sigma1=0.5, sigma2=0.5, N1=1.0, N2=1.0):
    """
    Mixture of 2 lognormal functions

    Parameters
    ----------
    xs: vector
       x values

    max_pos1: float
       Position of the lognormal function's maximum for component 1
    max_pos2: float
       Position of the lognormal function's maximum for component 2

    sigma1: float
       Sigma of the lognormal function for component 1
    sigma2: float
       Sigma of the lognormal function for component 2

    N1: floats
       Multiplicative factor for component 1
    N2: floats
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


def _exponential(x, a=2.0, N=1.0):
    """
    Exponential distribution
    Parameters
    ----------
    x: vector
       x values
    a: float
       Decay Rate parameter in exp: N*e^-ax
       Distribution Mean = 1/a
    N: float
       Multiplicative factor
    Returns
    -------
    exponential computed on the x grid
    """
    return N * np.exp(-1.0 * a * x)


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


def _imf_salpeter(x):
    """
    Compute a Salpeter IMF

    Parameters
    ----------
    x : numpy vector
      masses

    Returns
    -------
    imf : numpy vector
      unformalized IMF
    """
    return x ** (-2.35)


def _imf_kroupa(in_x):
    """
    Compute a Kroupa IMF

    Parameters
    ----------
    in_x : numpy vector
      masses

    Returns
    -------
    imf : numpy vector
      unformalized IMF
    """
    # allows for single float or an array
    x = np.atleast_1d(in_x)

    m1 = 0.08
    m2 = 0.5
    alpha0 = -0.3
    alpha1 = -1.3
    alpha2 = -2.3
    imf = np.full((len(x)), 0.0)

    (indxs,) = np.where(x >= m2)
    if len(indxs) > 0:
        imf[indxs] = x[indxs] ** alpha2

    (indxs,) = np.where((x >= m1) & (x < m2))
    fac1 = (m2 ** alpha2) / (m2 ** alpha1)
    if len(indxs) > 0:
        imf[indxs] = (x[indxs] ** alpha1) * fac1

    (indxs,) = np.where(x < m1)
    fac2 = fac1 * ((m1 ** alpha1) / (m1 ** alpha0))
    if len(indxs) > 0:
        imf[indxs] = (x[indxs] ** alpha0) * fac2

    return imf
