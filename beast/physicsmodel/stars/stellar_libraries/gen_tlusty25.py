import glob

import numpy as np

from astropy.table import QTable
from astropy.io import fits
from astropy.io import ascii

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
    periodpos = filename.rfind(".spec")

    zpos = filename.find("z", slashpos)
    tpos = filename.find("t", slashpos)
    gpos = filename.find("g", slashpos)
    vpos = filename.find("v", slashpos)

    if tpos - zpos > 4:
        model_params["Z"] = float(filename[zpos + 1 : tpos]) * 0.001
    else:
        model_params["Z"] = float(filename[zpos + 1 : tpos]) * 0.01
    model_params["Teff"] = float(filename[tpos + 1 : gpos])
    model_params["logg"] = float(filename[gpos + 1 : vpos]) * 0.01
    model_params["vturb"] = float(filename[vpos + 1 : periodpos])

    return model_params


if __name__ == "__main__":  # pragma: no cover

    solar_z = 0.02

    # get all vtrub=2 models
    path = "/home/kgordon/Python/extstar_data/Models/Tlusty_2025/"
    files = glob.glob(f"{path}/*v2.spec.gz")
    n_files = len(files)

    # fmt: off
    outtab = QTable(
        names=["Z", "Teff", "logg", "vturb", "logT", "logz"],
        dtype=["f", "f", "f", "f", "f", "f"]
    )
    # fmt: on

    for k, cfile in enumerate(files):
        print(cfile)

        # read in the model spectrum
        # wave units are angstrom
        # flux units are H, F = 4*pi*H in ergs/(cm^2 s A)
        mspec = ascii.read(
            cfile,
            format="no_header",
            fast_reader={"exponent_style": "D"},
            names=["Wave", "SFlux"],
        )

        # convert the type to float
        mspec["SFlux"] = mspec["SFlux"].astype(float)

        # now extract the wave and flux colums
        mwave = mspec["Wave"]
        mflux = mspec["SFlux"]

        # rebin to R=4000 for speed, common res and wave range with BOSZ LTE models
        rbres = 4000.0
        wave_rebin, flux_rebin, npts_rebin = rebin_spectrum(
            mwave.value, mflux.value, rbres, [500.0, 320000.0]
        )

        if k == 0:
            outspec = np.zeros((len(files) + 1, len(wave_rebin)))
            outspec[-1, :] = wave_rebin
        outspec[k, :] = flux_rebin * 4.0 * np.pi

        model_params = decode_params(cfile)
        Z = model_params["Z"] * solar_z
        row = (Z, model_params["Teff"], model_params["logg"], model_params["vturb"], np.log10(model_params["Teff"]),
               np.log10(Z))
        outtab.add_row(row)

    # output the stellar library in the beast format
    spec_hdu = fits.PrimaryHDU(data=outspec)
    table_hdu = fits.BinTableHDU(data=outtab)
    table_hdu.name = "TLUSTY"
    hdul = fits.HDUList([spec_hdu, table_hdu])
    hdul.writeto("tlusty2025.grid.fits", overwrite=True)
