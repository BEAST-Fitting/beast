import numpy as np
from scipy.stats import scoreatpercentile as sap
from . import kernels

kernel_switch = dict(gau=kernels.Gaussian, epa=kernels.Epanechnikov,
                    uni=kernels.Uniform, tri=kernels.Triangular,
                    biw=kernels.Biweight, triw=kernels.Triweight,
                    cos=kernels.Cosine)


#### Convenience Functions to be moved to kerneltools ####
def forrt(X, m=None):
    """
    RFFT with order like Munro (1976) FORTT routine.
    """
    if m is None:
        m = len(X)
    y = np.fft.rfft(X, m) / m
    return np.r_[y.real, y[1:-1].imag]


def revrt(X, m=None):
    """
    Inverse of forrt. Equivalent to Munro (1976) REVRT routine.
    """
    if m is None:
        m = len(X)
    y = X[:m / 2 + 1] + np.r_[0, X[m / 2 + 1:], 0] * 1j
    return np.fft.irfft(y) * m


def silverman_transform(bw, M, RANGE):
    """
    FFT of Gaussian kernel following to Silverman AS 176.

    Notes
    -----
    Underflow is intentional as a dampener.
    """
    J = np.arange(M / 2 + 1)
    FAC1 = 2 * (np.pi * bw / RANGE) ** 2
    JFAC = J ** 2 * FAC1
    BC = 1 - 1. / 3 * (J * 1. / M * np.pi) ** 2
    FAC = np.exp(-JFAC) / BC
    kern_est = np.r_[FAC, FAC[1:-1]]
    return kern_est


def linbin(X, a, b, M, trunc=1):
    """
    Linear Binning as described in Fan and Marron (1994)
    """
    gcnts = np.zeros(M)
    delta = (b - a) / (M - 1)

    for x in X:
        lxi = ((x - a) / delta)  # +1
        li = int(lxi)
        rem = lxi - li
        if (li > 1) and (li < M):
            gcnts[li] = gcnts[li] + 1 - rem
            gcnts[li + 1] = gcnts[li + 1] + rem
        if (li > M) and (trunc == 0):
            gcnts[M] = gcnts[M] + 1

    return gcnts


def _select_sigma(X):
    """
    Returns the smaller of std(X, ddof=1) or normalized IQR(X) over axis 0.

    References
    ----------
    Silverman (1986) p.47
    """
#    normalize = norm.ppf(.75) - norm.ppf(.25)
    normalize = 1.349
#    IQR = np.subtract.reduce(percentile(X, [75,25],
#                             axis=axis), axis=axis)/normalize
    IQR = (sap(X, 75) - sap(X, 25)) / normalize
    return np.minimum(np.std(X, axis=0, ddof=1), IQR)


## Univariate Rule of Thumb Bandwidths ##
def bw_scott(x):
    """
    Scott's Rule of Thumb

    Parameter
    ---------
    x : array-like
        Array for which to get the bandwidth

    Returns
    -------
    bw : float
        The estimate of the bandwidth

    Notes
    -----
    Returns 1.059 * A * n ** (-1/5.)

    A = min(std(x, ddof=1), IQR/1.349)
    IQR = np.subtract.reduce(np.percentile(x, [75,25]))

    References
    ---------- ::

    Scott, D.W. (1992) `Multivariate Density Estimation: Theory, Practice, and
        Visualization.`
    """
    A = _select_sigma(x)
    n = len(x)
    return 1.059 * A * n ** -.2


def bw_silverman(x):
    """f
    Silverman's Rule of Thumb

    Parameter
    ---------
    x : array-like
        Array for which to get the bandwidth

    Returns
    -------
    bw : float
        The estimate of the bandwidth

    Notes
    -----
    Returns .9 * A * n ** (-1/5.)

    A = min(std(x, ddof=1), IQR/1.349)
    IQR = np.subtract.reduce(np.percentile(x, [75,25]))

    References
    ---------- ::

    Silverman, B.W. (1986) `Density Estimation.`
    """
    A = _select_sigma(x)
    n = len(x)
    return .9 * A * n ** -.2


## Helper Functions ##
bandwidth_funcs = {'scott': bw_scott, 'silverman': bw_silverman}


def select_bandwidth(X, bw, kernel):
    """
    Selects bandwidth
    """
    bw = bw.lower()
    if bw not in ["scott", "silverman"]:
        raise ValueError("Bandwidth {} not understood".format(bw))
    return bandwidth_funcs[bw](X)


def counts(x, v):
    """
    Counts the number of elements of x that fall within the grid points v

    Notes
    -----
    Using np.digitize and np.bincount
    """
    idx = np.digitize(x, v)
    try:  # numpy 1.6
        return np.bincount(idx, minlength=len(v))
    except:
        bc = np.bincount(idx)
        return np.r_[bc, np.zeros(len(v) - len(bc))]


def kdesum(x, axis=0):
    return np.asarray([np.sum(x[i] - x, axis) for i in range(len(x))])


#### Kernel Density Estimator Functions ####
def kdensity(X, kernel="gau", bw="scott", weights=None, gridsize=None,
             adjust=1, clip=(-np.inf, np.inf), cut=3, retgrid=True, evalpts=None):
    """
    Rosenblatz-Parzen univariate kernel density estimator

    Parameters
    ----------
    X : array-like
        The variable for which the density estimate is desired.
    kernel : str
        The Kernel to be used. Choices are
        - "biw" for biweight
        - "cos" for cosine
        - "epa" for Epanechnikov
        - "gauss" for Gaussian.
        - "tri" for triangular
        - "triw" for triweight
        - "uni" for uniform
    bw : str, float
        "scott" - 1.059 * A * nobs ** (-1/5.), where A is min(std(X),IQR/1.34)
        "silverman" - .9 * A * nobs ** (-1/5.), where A is min(std(X),IQR/1.34)
        If a float is given, it is the bandwidth.
    weights : array or None
        Optional  weights. If the X value is clipped, then this weight is
        also dropped.
    gridsize : int
        If gridsize is None, max(len(X), 50) is used.
    adjust : float
        An adjustment factor for the bw. Bandwidth becomes bw * adjust.
    clip : tuple
        Observations in X that are outside of the range given by clip are
        dropped. The number of observations in X is then shortened.
    cut : float
        Defines the length of the grid past the lowest and highest values of X
        so that the kernel goes to zero. The end points are
        -/+ cut*bw*{min(X) or max(X)}
    retgrid : bool
        Whether or not to return the grid over which the density is estimated.

    Returns
    -------
    density : array
        The densities estimated at the grid points.
    grid : array, optional
        The grid points at which the density is estimated.

    Notes
    -----
    Creates an intermediate (`gridsize` x `nobs`) array. Use FFT for a more
    computationally efficient version.
    """
    X = np.asarray(X)
    if X.ndim == 1:
        X = X[:, None]
    clip_x = np.logical_and(X > clip[0], X < clip[1])
    X = X[clip_x]

    nobs = float(len(X))  # after trim

    if gridsize is None:
        gridsize = max(nobs, 50)  # don't need to resize if no FFT

        # handle weights
    if weights is None:
        weights = np.ones(nobs)
        q = nobs
    else:
        if len(weights) != len(clip_x):
            msg = "The length of the weights must be the same as the given X."
            raise ValueError(msg)
        weights = weights[clip_x.squeeze()]
        q = weights.sum()

    # if bw is None, select optimal bandwidth for kernel
    try:
        bw = float(bw)
    except:
        bw = select_bandwidth(X, bw, kernel)
    bw *= adjust

    if evalpts is None:
        a = np.min(X, axis=0) - cut * bw
        b = np.max(X, axis=0) + cut * bw
        grid = np.linspace(a, b, gridsize)
    else:
        grid = np.asarray(evalpts)

    k = (X.T - grid[:, None]) / bw  # uses broadcasting to make a gridsize x nobs

    # instantiate kernel class
    kern = kernel_switch[kernel](h=bw)
    # truncate to domain
    if kern.domain is not None:  # won't work for piecewise kernels like parzen
        z_lo, z_high = kern.domain
        domain_mask = (k < z_lo) | (k > z_high)
        k = kern(k)  # estimate density
        k[domain_mask] = 0
    else:
        k = kern(k)  # estimate density

    k[k < 0] = 0  # get rid of any negative values, do we need this?

    dens = np.dot(k, weights) / (q * bw)

    if retgrid:
        return dens, grid, bw
    else:
        return dens, bw
