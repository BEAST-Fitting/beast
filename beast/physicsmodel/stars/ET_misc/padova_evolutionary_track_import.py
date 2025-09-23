import glob
import os
import numpy as np
from astropy.table import Table, vstack, Column


# code from sometime ago, modified to run only from the command line
# updates needed
if __name__ == "__main__":
    data = Table()
    phil_phase_data = Table()
    directories = []

    data.meta["comments"] = [
        "Interpolated Padova evolutionary tracks",
        "For more information on data see http://philrosenfield.github.io/padova_tracks/",
        "Created for use in the BEAST code (https://github.com/BEAST-Fitting)",
        "Bolometric magnitudes have been calculated using: 4.77 - 2.5 * log L",
        "Surface gravities have been calculated using:  -10.616 + log mass + 4.0 * log Teff - log L",
        "0 : Beginning of Pre Main Sequence",
        "1 : Beginning of Main Sequence (MS)*",
        "2 : log Teff minimum on the Main Sequence*",
        "3 : Main Sequence turn off (MSTO)*",
        "4 : RG tip: Tip of the RGB (TRGB)*",
        "5 : Beginning Core He Fusion: Start central He-burning phase*",
        "6 : End of Core He Fusion*",
        "7 : Horizonzal Branch",
        "8 : Beginning of TP-AGB (if the star has the a phase)",
        "",
        "* EEP definitions followed Dotter et al. 2016,ApJ,222,8D as close as possible",
    ]

    # equivalent evolutionary points
    eeps = np.array(
        ([0] * 200)
        + [1] * 200
        + [2] * 200
        + [3] * 500
        + [4] * 30
        + [5] * 500
        + [6] * 100
        + [8] * 200
    )
    eeps_hb = np.array([7] * 600 + [8] * 200)

    # gets data from directories in current directory
    for dir in os.listdir("."):
        directories.append(dir)

    # combines data
    for dir in directories:
        print(dir)
        z = dir.split("Y")[0][1:]
        files = glob.glob(dir + "/*")
        for file in files:
            new_line = Table.read(file, format="ascii")
            z_col = Column(
                np.array([float(z)] * len(new_line)).tolist(),
                name="Z",
                description="Metallicity in solar metallicity",
            ).T
            index = Column(np.arange(0, len(z_col)), name="index")
            if "HB" in file:
                phase_col = eeps_hb[0 : len(new_line)]
            else:
                phase_col = eeps[0 : len(new_line)]
            # print(phase_col)
            phil_phase_array = Column(
                phase_col, name="evolutionary_point", description="rosenfield eep"
            ).T
            new_line.add_column(z_col)
            new_line.add_column(phil_phase_array)
            new_line.add_column(index)
            data = vstack([data, new_line])

    # specifies header information
    data.rename_column("logAge", "logA")
    data.rename_column("Mass", "M_ini")
    data.rename_column("logTe", "logT")
    data.rename_column("Mbol", "mbolmag")

    data["logA"] = Column(data["logA"], unit="yr", description="Age")
    data["M_ini"] = Column(data["M_ini"], unit="Msun", description="Initial mass")
    data["logT"] = Column(data["logT"], unit="K", description="Effective temperature")
    data["mbolmag"] = Column(
        data["mbolmag"],
        description="Bolomertric luminosities (magnitudes),\
                                        calculated using: 4.77 - 2.5 * log L",
    )
    data["logg"] = Column(
        data["logg"],
        unit="cm/s**2",
        description="surface gravity, calculated using:\
                                    -10.616 + log mass + 4.0 * log Teff - log L",
    )
    data["C/O"] = Column(data["C/O"], description="Carbon to oxygen ratio")

    data.write("Padova_ev_tracks.fits", overwrite=True)
