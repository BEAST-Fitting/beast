from scipy.interpolate import interp1d
import numpy as np


def rebin_spectrum(wave, flux, resolution, wave_range):
    """
    Rebin spectrum to input resolution over input wavelength range
    High to lower resolution only.

    Parameters
    ----------
    wave : vector
        wavelengths of spectrum

    flux : vector
        spectrum in flux units

    resolution : float
        resolution of output spectrum

    wave_range : [float, float]
        wavelength range of output spectrum

    Outputs
    ------
    wave, flux, npts : tuple of vectors
        the model wavelength, flux, and npts at the requested wavelength

    """
    npts = int(
        np.log10(wave_range[1] / wave_range[0])
        / np.log10((1.0 + 2.0 * resolution) / (2.0 * resolution - 1.0))
    )

    twave = np.logspace(
        np.log10(wave_range[0]), np.log10(wave_range[1]), num=npts + 1, endpoint=True
    )
    full_wave_min = twave[0:-1]
    full_wave_max = twave[1:]
    full_wave = 0.5 * (full_wave_min + full_wave_max)

    full_flux = np.zeros((npts))
    full_npts = np.zeros((npts), dtype=int)

    for k in range(npts):
        (indxs,) = np.where((wave >= full_wave_min[k]) & (wave < full_wave_max[k]))
        n_indxs = len(indxs)
        if n_indxs > 0:
            full_flux[k] = np.sum(flux[indxs])
            full_npts[k] = n_indxs

    # divide by the # of points to create the final rebinned spectrum
    (indxs,) = np.where(full_npts > 0)
    if len(indxs):
        full_flux[indxs] = full_flux[indxs] / full_npts[indxs]

    # interpolate to fill in missing points in the rebinned spectrum
    #  e.g., the model spectrum is not computed at a high enough resolution
    #        at all the needed wavelengths
    (zindxs,) = np.where(full_npts <= 0)
    if len(zindxs):
        ifunc = interp1d(
            full_wave[indxs], full_flux[indxs], kind="linear", bounds_error=False
        )
        full_flux = ifunc(full_wave)
        full_npts[zindxs] = 1
        (nanindxs,) = np.where(~np.isfinite(full_flux))
        if len(nanindxs):
            full_flux[nanindxs] = 0.0
            full_npts[nanindxs] = 0

    return (full_wave, full_flux, full_npts)


def get_stellib_boundaries(s, dlogT=0.1, dlogg=0.3, closed=True):
    """
    Returns the closed boundary polygon around the stellar library with
    given margins.

    Parameters
    ----------
    s :  Stellib
        Stellar library object
    dlogT :  float
        margin in logT
    dlogg  : float
        margin in logg
    closed : bool
        if set, close the polygon

    Returns
    -------
    b : ndarray[float, ndim=2]
        (closed) boundary points: [logg, Teff]

    """
    leftb = [(k, np.max(s.logT[s.logg == k]) + dlogT) for k in np.unique(s.logg)]
    leftb += [(leftb[-1][0] + dlogg, leftb[-1][1])]
    leftb = [(leftb[0][0] - dlogg, leftb[0][1])] + leftb
    rightb = [(k, np.min(s.logT[s.logg == k]) - dlogT) for k in np.unique(s.logg)[::-1]]
    rightb += [(rightb[-1][0] - dlogg, rightb[-1][1])]
    rightb = [(rightb[0][0] + dlogg, rightb[0][1])] + rightb
    b = leftb + rightb
    if closed:
        b += [b[0]]
    b = np.array(b)
    outb = [b[:, 1], b[:, 0]]

    return np.transpose(np.array(outb))