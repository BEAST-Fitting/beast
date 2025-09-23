"""
Convert the MIST ascii files to FITS tables
"""
from glob import glob

import numpy as np
from tqdm import tqdm
from astropy.table import Table, Column

from beast.config import __ROOT__


def convert_ascii_to_fits_one_met(infiles, FeH, outfilename):
    """
    Convert the ascii files into FITS tables to speed up repeated reading

    Parameters
    ----------
    infiles : list
        ascii files names (one per intial mass)
    FeH : float
        Fe/H value for this set of tracks
    outfilename : str
        output filename
    """
    if isinstance(infiles, list):
        files = infiles
    else:
        files = [infiles]

    mass_act = np.array([])
    mass_ini = np.array([])
    logA = np.array([])
    logL = np.array([])
    logT = np.array([])
    logg = np.array([])
    phase = np.array([])
    eep = np.array([])
    for cfile in tqdm(files, desc=outfilename):
        a = Table.read(cfile, format="ascii", header_start=11)
        tmass = a["star_mass"].data
        mass_act = np.concatenate((mass_act, tmass))
        mass_ini = np.concatenate((mass_ini, np.full((len(tmass)), max(tmass))))
        logA = np.concatenate((logA, np.log10(a["star_age"].data)))
        logL = np.concatenate((logL, a["log_L"].data))
        logT = np.concatenate((logT, a["log_Teff"].data))
        logg = np.concatenate((logg, a["log_g"].data))
        phase = np.concatenate((phase, a["phase"].data))
        eep = np.concatenate((eep, range(len(a))))

    data = Table()
    data["log(M_act)"] = Column(np.log10(mass_act))
    data["log(M_ini)"] = Column(np.log10(mass_ini))
    data["logA"] = Column(logA)
    data["logL"] = Column(logL)
    data["logT"] = Column(logT)
    data["logg"] = Column(logg)
    data["phase"] = Column(phase)
    data["eep"] = Column(eep)

    met = np.full((len(eep)), FeH)
    data["met"] = Column(met, description="Fe/H values")

    data.write(outfilename, overwrite=True)


if __name__ == "__main__":
    # fmt: off
    FeH = [-4.00, -3.50, -3.00, -2.50, -2.00, -1.75,
           -1.50, -1.25, -1.00, -0.75, -0.25, 0.0,
           0.25, 0.5]
    FeH_str = ["m4.00", "m3.50", "m3.00", "m2.50", "m2.00", "m1.75",
               "m1.50", "m1.25", "m1.00", "m0.75", "m0.25", "p0.00",
               "p0.25", "p0.50"]
    # fmt: on

    for cFeH, cFeH_str in zip(FeH, FeH_str):
        cfiles = glob(
            f"{__ROOT__}/MIST/MIST_v1.2_feh_{cFeH_str}_afe_p0.0_vvcrit0.4_EEPS/*.eep"
        )
        outfile = f"MIST_FeH{cFeH:.2f}_vvcrit0.4.fits"
        convert_ascii_to_fits_one_met(cfiles, cFeH, outfile)
