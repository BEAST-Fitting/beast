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

    gpos = filename.find("_g", slashpos)
    zpos = filename.find("_m", gpos)
    tpos = filename.find("_t", slashpos)
    vpos = filename.find("_v", slashpos)

    # print(filename[zpos + 2 : zpos + 7])
    # print(filename[tpos + 2 : gpos])
    # print(filename[gpos + 2 : gpos + 6])
    # print(filename[vpos + 2 : vpos + 3])
    # exit()

    model_params["Z"] = float(filename[zpos + 2 : zpos + 7])
    model_params["Teff"] = float(filename[tpos + 2 : gpos])
    model_params["logg"] = float(filename[gpos + 2 : gpos + 6])
    model_params["vturb"] = float(filename[vpos + 2 : vpos + 3])

    return model_params


if __name__ == "__main__":  # pragma: no cover

    solar_z = 0.02

    # get all vtrub=2 models
    path = "/home/kgordon/Python/extstar_data/Models/BOSZ2024/r5000/r5000/"
    files = glob.glob(f"{path}/*/*a+0.00_c+0.00_v2_r5000_resam.txt.gz")
    n_files = len(files)

    # fmt: off
    outtab = QTable(
        names=["Z", "Teff", "logg", "vturb", "logT", "logz"],
        dtype=["f", "f", "f", "f", "f", "f"]
    )
    # fmt: on

    for k, cfile in enumerate(files):
        print(cfile)

        # get model parameters
        model_params = decode_params(cfile)

        # skip if ATLAS ("ap") and 8000 K or lower
        # MARCS models cover these temps
        usemod = True
        if ("ap" in cfile) and (model_params["Teff"] <= 8000.0):
            usemod = False
            print("skipped")

        if usemod:
            Z = (10 ** model_params["Z"]) * solar_z
            row = (Z, model_params["Teff"], model_params["logg"], model_params["vturb"], np.log10(model_params["Teff"]),
                np.log10(Z))
            outtab.add_row(row)
            print(row)

            # read in the model spectrum
            # wave units are angstrom
            # flux units are H, F = 4*pi*H in ergs/(cm^2 s A)
            # wavelength first
            wave_filename = f"{path}/bosz2024_wave_r5000.txt"
            mspec_lte_wave = ascii.read(
                wave_filename,
                format="no_header",
                names=["Wave"],
            )

            mspec = ascii.read(
                cfile,
                format="no_header",
                names=["SFlux", "SCont"],
            )
            mspec["Wave"] = mspec_lte_wave["Wave"]

            # convert the type to float
            mspec["SFlux"] = mspec["SFlux"].astype(float)

            # now extract the wave and flux columns
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

    # output the stellar library in the beast format
    # spec_hdu = fits.PrimaryHDU(data=outspec)
    # table_hdu = fits.BinTableHDU(data=outtab)
    # table_hdu.name = "BOSZ"
    # hdul = fits.HDUList([spec_hdu, table_hdu])
    # hdul.writeto("bosz2024.grid.fits", overwrite=True)
