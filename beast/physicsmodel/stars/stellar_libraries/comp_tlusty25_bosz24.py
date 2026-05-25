import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.io import ascii

from helpers import rebin_spectrum

if __name__ == "__main__":  # pragma: no cover

    solar_z = 0.02

    # read the kurucz spectrum
    # wavelength first
    wave_filename = "/home/kgordon/Python/extstar_data/Models/BOSZ2024/r2000/bosz2024_wave_r2000.txt"
    mspec_lte_wave = ascii.read(
        wave_filename,
        format="no_header",
        # fast_reader={"exponent_style": "D"},
        names=["Wave"],
    )

    model_filename = "/home/kgordon/Python/extstar_data/Models/BOSZ2024/r2000/m+0.00/bosz2024_ap_t15000_g+4.0_m+0.00_a+0.00_c+0.00_v2_r2000_resam.txt.gz"
    mspec_lte = ascii.read(
        model_filename,
        format="no_header",
        # fast_reader={"exponent_style": "D"},
        names=["SFlux", "SCont"],
    )
    mspec_lte["Wave"] = mspec_lte_wave["Wave"] * u.angstrom

    # read the tlusty spectrum
    model_filename = (
        "/home/kgordon/Python/extstar_data/Models/Tlusty_2023/z100t15000g400v2.spec.gz"
    )
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
