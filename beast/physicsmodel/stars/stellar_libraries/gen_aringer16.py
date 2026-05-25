import matplotlib.pyplot as plt

import numpy as np

from astropy.table import QTable
from astropy.io import fits
from astropy.io import ascii
from astropy import units as u
from astropy.modeling import models

from helpers import rebin_spectrum


def decode_params(filename):
    """
    Decode the tlusty filenames for the model parameters.

    Parameters
    ----------
    filename : str
        name of file to decode

    """
    model_params = {}

    slashpos = filename.rfind("/")

    tpos = filename.find("t", slashpos)
    gpos = filename.find("g", slashpos)
    masspos = filename.find("_m", slashpos)
    carbonpos = filename.find("_c", slashpos)
    oxygenpos = filename.find("_o", slashpos)
    ironpos = filename.find("_fe", slashpos)
    vpos = filename.find("_xi", slashpos)
    cpos = filename.find("_c", vpos)

    model_params["Teff"] = float(filename[tpos + 1 : gpos - 1])
    model_params["logg"] = float(filename[gpos + 1 : masspos]) * 0.01
    model_params["vturb"] = float(filename[vpos + 3 : cpos]) * 0.1
    model_params["Mass"] = float(filename[masspos + 2 : oxygenpos]) * 0.01
    model_params["logO"] = np.log10(float(f"0.{filename[oxygenpos + 2 : carbonpos]}"))
    model_params["logC"] = np.log10(float(f"0.{filename[carbonpos + 2 : ironpos]}"))
    model_params["C_to_O"] = (10 ** (model_params["logC"] - 12)) / (
        10 ** (model_params["logO"] - 12)
    )

    # metallicity from Fe (comments from Martha's code)

    # Aringer+2019 says "Since we use here the common definition, where log(eps(H)) is set
    # constant to 12, and no COMARCS models with special deviations of the Fe content exist,
    # the quantities [Z/H], [Fe/H], and log(Z/Zsun) are always equal."  So then all that
    # matters is what solar value I use.  The latest one that the CMD interface uses
    # is 0.0207... maybe I should use that?

    # SEE UPDATE!  KEEPING THIS TEXT FOR POSTERITY. "For the solar composition we adopted the
    # values from Anders & Grevesse (1989) except for C, N and O where we took the data from
    # Grevesse & Sauval (1994). This is in agreement with our previous work
    # (e.g. Aringer et al. 1999) and results in a Z⊙ of approximately
    # 0.02." --> So I will use 0.01922

    # UPDATE !!!!!! Asplund et al. 2021 say Zsun/Xsun = 0.0187.  Lodders et al. (2019)
    # say 0.0200.  Leo's CMD page says 0.0207.  I will use 0.02 here b/c Bernhard mentioned
    # 0.02 in his 2016 text AND because if you do the math on the Kurucz models, they use
    # Zsun=0.02 exactly.

    Fenum = float(f"0.{filename[ironpos + 3 : vpos]}") * 10.0
    # solar value is listed as 7.56 in Aringer+2019 (MNRAS, 487,2133)
    Fesun = 10.0 ** (7.56 - 12)
    Fe = 10.0 ** (Fenum - 12)
    FeH = np.round(np.log10(Fe) - np.log10(Fesun), 1)
    Z = (10.0**FeH) * 0.02  # see comments below
    model_params["Z"] = (10.0**FeH) * 0.02  # see comments below

    return model_params


if __name__ == "__main__":  # pragma: no cover

    solar_z = 0.02  # see comments in decode_params

    # get all vtrub=2 models
    path = "/home/kgordon/Python/extstar_data/Models/Aringer16/"
    mspec = ascii.read(f"{path}/mstar_spec_list.dat", data_start=0, format="no_header")
    files = [f"{path}{cfile}" for cfile in mspec["col1"].data]

    n_files = len(files)

    # fmt: off
    outtab = QTable(
        names=["Z", "Teff", "logg", "vtrub", "logT", "logz", "Mass", "logO", "logC", "R_Rs"],
        dtype=["f", "f", "f", "f", "f", "f", "f", "f", "f", "f"]
    )
    # fmt: on

    for k, cfile in enumerate(files):
        print(cfile)

        model_params = decode_params(cfile)
        Z = model_params["Z"] * solar_z

        # will also need the radius assumed when creating the spectrum
        # calculated from log g = log M − 2 log R + 4.437

        # radius in solar units
        R_Rs = 10.0 ** (
            (model_params["logg"] - 4.437 - np.log10(model_params["Mass"])) / (-2.0)
        )
        # radius in cm
        Reff = R_Rs * 6.957e10  # that's sun's radius in cm

        row = (
            Z,
            model_params["Teff"],
            model_params["logg"],
            model_params["vturb"],
            np.log10(model_params["Teff"]),
            model_params["Z"],
            model_params["Mass"],
            model_params["logO"],
            model_params["logC"],
            R_Rs,
        )
        print(row)
        outtab.add_row(row)

        # read in the model spectrum
        # wave units are angstrom
        # flux is in nu*L_nu [erg/s]
        # (radius has been used, unlike other stellar atmosphere results)
        mspec = ascii.read(
            cfile,
            names=["Wave", "d1", "SFlux", "d2"],
        )

        # convert the type to float
        mspec["SFlux"] = mspec["SFlux"].astype(float)

        # now extract the wave and flux columns
        mwave = mspec["Wave"]
        mflux = mspec["SFlux"]

        # spectra was in nu*L_nu [erg/s].  This conversion puts in erg/s/cm^2/AA
        mflux /= (4 * np.pi * (Reff**2.0)) / mwave

        # rebin to R=2000 for speed, common res and wave range with BOSZ LTE models
        rbres = 2000.0
        wave_rebin, flux_rebin, npts_rebin = rebin_spectrum(
            mwave.value, mflux.value, rbres, [500.0, 320000.0]
        )

        # add a scaled blackbody for the shorter and longer wavelengths that were not calculated
        bb = models.BlackBody(temperature = model_params["Teff"]*u.K,
                              scale=1.0 * u.erg / (u.cm * u.cm * u.s * u.AA * u.sr)
                              )
        modbb = bb(wave_rebin)

        # shorter wavelengths
        gvals = np.absolute(wave_rebin - 3370.0) < 10.0
        modbb *= np.nanmean(flux_rebin[gvals]) / np.nanmean(modbb[gvals])
        gvals = (wave_rebin < 5000.0) & (flux_rebin == 0.0)
        flux_rebin[gvals] = modbb[gvals]

        # longer wavelengths
        gvals = np.absolute(wave_rebin - (250000.0 - 150.0)) < 100.0
        modbb *= np.nanmean(flux_rebin[gvals]) / np.nanmean(modbb[gvals])
        gvals = (wave_rebin > 5000.0) & (flux_rebin == 0.0)
        flux_rebin[gvals] = modbb[gvals]

        if k == 0:
            outspec = np.zeros((len(files) + 1, len(wave_rebin)))
            outspec[-1, :] = wave_rebin
        outspec[k, :] = flux_rebin

        plt.plot(wave_rebin, flux_rebin)
    plt.xscale("log")
    plt.yscale("log")
    plt.show()

    # output the stellar library in the beast format
    spec_hdu = fits.PrimaryHDU(data=outspec)
    table_hdu = fits.BinTableHDU(data=outtab)
    table_hdu.name = "Aringer16"
    hdul = fits.HDUList([spec_hdu, table_hdu])
    hdul.writeto("aringer2016.grid.fits", overwrite=True)
