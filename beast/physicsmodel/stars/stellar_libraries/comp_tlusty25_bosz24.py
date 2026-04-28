from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import QTable
from astropy.io import fits
import astropy.units as u
from astropy.io import ascii


def rebin_spectrum(wave, flux, resolution, wave_range):
    """
    Rebin spectrum to input resolution over input wavelength range
    High to lower resolution only

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


if __name__ == "__main__":  # pragma: no cover

    solar_z = 0.02

    # read the kurucz spectrum
    # wavelength first
    wave_filename = "/home/kgordon/Python/extstar_data/Models/BOSZ2024/r2000/bosz2024_wave_r2000.txt"
    mspec_lte_wave = ascii.read(
        wave_filename,
        format="no_header",
        #fast_reader={"exponent_style": "D"},
        names=["Wave"],
    )

    model_filename = "/home/kgordon/Python/extstar_data/Models/BOSZ2024/r2000/m+0.00/bosz2024_ap_t15000_g+4.0_m+0.00_a+0.00_c+0.00_v2_r2000_resam.txt.gz"
    mspec_lte = ascii.read(
        model_filename,
        format="no_header",
        #fast_reader={"exponent_style": "D"},
        names=["SFlux", "SCont"],
    )
    mspec_lte["Wave"] = mspec_lte_wave["Wave"] * u.angstrom

    # read the tlusty spectrum
    model_filename = "/home/kgordon/Python/extstar_data/Models/Tlusty_2023/z100t15000g400v2.spec.gz"
    # read in the model spectrum
    mspec = ascii.read(
        model_filename,
        format="no_header",
        fast_reader={"exponent_style": "D"},
        names=["Wave", "SFlux"],
    )

    # convert the type to float
    mspec["SFlux"] = mspec["SFlux"].astype(float)

    # set the units
    fluxunit = u.erg / (u.s * u.cm * u.cm * u.angstrom)
    mspec["Wave"].unit = u.angstrom
    mspec["SFlux"].unit = fluxunit

    # now extract the wave and flux colums
    mwave = mspec["Wave"]
    mflux = mspec["SFlux"]

    # rebin to R=4000 for speed
    rbres = 4000.0
    wave_rebin, flux_rebin, npts_rebin = rebin_spectrum(
        mwave.value, mflux.value, rbres, [200.0, 310000.0]
    )
    wave_rebin *= u.angstrom
    print(len(wave_rebin))

    # plot
    plt.plot(wave_rebin.to(u.micron), flux_rebin * 4 * np.pi)
    plt.plot(mspec_lte["Wave"].to(u.micron), mspec_lte["SFlux"] * 4 * np.pi)
    plt.yscale("log")
    plt.xscale("log")
    plt.show()
