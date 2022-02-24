# script to generate BEAST library files
import numpy as np

import stsynphot as stsyn

import h5py

from beast.config import __ROOT__
from beast.observationmodel import phot


def make_libfile():
    """
    Extract filters from STSYNPHOT and save to main library file.
    """
    # wfc3_obsmodes_uvis
    wfc3_uvis = [
        "f200lp",
        "f218w",
        "f225w",
        "f275w",
        "f280n",
        "f300x",
        "f336w",
        "f343n",
        "f350lp",
        "f373n",
        "f390m",
        "f390w",
        "f395n",
        "f410m",
        "f438w",
        "f467m",
        "f469n",
        "f475w",
        "f475x",
        "f487n",
        "f502n",
        "f547m",
        "f555w",
        "f600lp",
        "f606w",
        "f621m",
        "f625w",
        "f631n",
        "f645n",
        "f656n",
        "f657n",
        "f658n",
        "f665n",
        "f673n",
        "f680n",
        "f689m",
        "f763m",
        "f775w",
        "f814w",
        "f845m",
        "f850lp",
        "f953n",
        "fq232n",
        "fq243n",
        "fq378n",
        "fq387n",
        "fq422m",
        "fq436n",
        "fq437n",
        "fq492n",
        "fq508n",
        "fq575n",
        "fq619n",
        "fq634n",
        "fq672n",
        "fq674n",
        "fq727n",
        "fq750n",
        "fq889n",
        "fq906n",
        "fq924n",
        "fq937n",
    ]

    wfc3_ir = [
        "f098m",
        "f105w",
        "f110w",
        "f125w",
        "f126n",
        "f127m",
        "f128n",
        "f130n",
        "f132n",
        "f139m",
        "f140w",
        "f153m",
        "f160w",
        "f164n",
        "f167n",
    ]

    # Open hd5 file for writing
    hf = h5py.File(__ROOT__ + "filters.hd5", "w")

    # Create group for nice hierarchical structure
    f = hf.create_group("filters")

    # Define arrays for "contents" / descriptive information
    tablenames = []
    observatories = []
    filternames = []
    names = []
    norms = []
    cwaves = []
    pwaves = []
    comments = []

    # Loop through WFC3_UVIS filters
    for mode in wfc3_uvis:

        # define uvis 1 and uvis2 modes
        mode_1 = "wfc3, uvis1, " + mode
        mode_2 = "wfc3, uvis2, " + mode

        # pull bandpasses from stsynphot for the two uvis modes
        bp_1 = stsyn.band(mode_1)
        bp_2 = stsyn.band(mode_2)

        # extract the wavelength array
        wave = bp_1.waveset

        # compute the average bandpass between uvis1 and uvis2
        bp_avg = np.average([bp_1(wave), bp_2(wave)], axis=0)

        # define the filter name
        filter_name = "HST_" + mode_1.split(",")[0].upper() + "_" + mode.upper()

        # build array of wavelength and throughput
        arr = np.array(
            list(zip(wave.value.astype(np.float64), bp_avg.astype(np.float64))),
            dtype=[("WAVELENGTH", "float64"), ("THROUGHPUT", "float64")],
        )

        # append dataset to the hdf5 filters group
        f.create_dataset(filter_name, data=arr)

        # generate filter instance to compute relevant info
        newfilt = phot.Filter(wave, bp_avg, name=filter_name)

        # populate contents lists with relevant information
        tablenames.append(filter_name)
        observatories.append("HST")
        filternames.append(mode.upper)
        names.append(newfilt.name)
        norms.append(newfilt.norm.value)
        cwaves.append(newfilt.cl.value)
        pwaves.append(newfilt.lpivot.value)
        comments.append("avg of uvis1 and uvis2")

    # Loop through WFC3_IR filters
    for mode in wfc3_ir:

        # define uvis 1 and uvis2 modes
        mode = "wfc3, ir, " + mode

        # pull bandpasses from stsynphot for the two uvis modes
        bp = stsyn.band(mode)

        # extract the wavelength array
        wave = bp.waveset

        # define the filter name
        filter_name = "HST_" + mode_1.split(",")[0].upper() + "_" + mode.upper()

        # build array of wavelength and throughput
        arr = np.array(
            list(zip(wave.value.astype(np.float64), bp(wave).astype(np.float64))),
            dtype=[("WAVELENGTH", "float64"), ("THROUGHPUT", "float64")],
        )

        # append dataset to the hdf5 filters group
        f.create_dataset(filter_name, data=arr)

        # generate filter instance to compute relevant info
        newfilt = phot.Filter(wave, bp(wave), name=filter_name)

        # populate contents lists with relevant information
        tablenames.append(filter_name)
        observatories.append("HST")
        filternames.append(mode.upper)
        names.append(newfilt.name)
        norms.append(newfilt.norm.value)
        cwaves.append(newfilt.cl.value)
        pwaves.append(newfilt.lpivot.value)
        comments.append("")

    # smash the contents arrays together
    contents = np.array(
        list(
            zip(
                tablenames,
                observatories,
                filternames,
                names,
                norms,
                cwaves,
                pwaves,
                comments,
            )
        ),
        dtype=[
            ("TABLENAME", "S40"),
            ("OBSERVATORY", "S30"),
            ("INSTRUMENT", "S30"),
            ("NAME", "S10"),
            ("NORM", "<f8"),
            ("CWAVE", "<f8"),
            ("PWAVE", "<f8"),
            ("COMMENT", "S100"),
        ],
    )

    # add the contents array as an hd5 dataset
    hf.create_dataset("content", data=contents)

    # close the file
    hf.close()


if __name__ == "__main__":  # pragma: no cover
    make_libfile()
