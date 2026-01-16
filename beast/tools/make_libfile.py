# script to generate BEAST library files
import os
import numpy as np

# suppress Vega file not found warning
import warnings
from astropy.utils.exceptions import AstropyUserWarning
warnings.simplefilter("ignore", category=AstropyUserWarning)

import stsynphot as stsyn

from pandeia.engine.instrument_factory import InstrumentFactory
import astropy.units as u
from astropy.table import QTable

import h5py

from beast.config import __ROOT__
from beast.observationmodel import phot


def make_filters_libfile():
    """
    Extract filters from STSYNPHOT and save to main library file.
    """

    # fmt: off
    wfc3_uvis = ["f218w", "f225w", "f275w", "f336w", "f390w", "f438w",
                 "f475w", "f555w", "f606w", "f625w", "f775w", "f814w",
                 "f390m", "f410m", "f467m", "f547m", "f621m", "f689m", "f763m", "f845m",
                 "f280n", "f343n", "f373n", "f395n", "f469n", "f487n", "f502n", "f631n",
                 "f645n", "f656n", "f657n", "f658n", "f665n", "f673n", "f680n", "f953n"]

    wfc3_ir = ["f105w", "f110w", "f125w", "f140w", "f160w",
               "f098m", "f127m", "f139m", "f153m"]

    wfpc2 = ["f157w", "f170w", "f185w", "f218w", "f255w",
             "f300w", "f336w", "f380w", "f439w", "f450w", "f555w",
             "f569w", "f606w", "f622w", "f675w", "f702w", "f791w", "f814w",
             "f122m", "f410m", "f467m", "f547m"]

    acs_wfc = ["f435w", "f475w", "f555w", "f606w", "f625w", "f775w", "f814w", "f850lp",
               "f502n", "f550m", "f658n"]

    acs_sbc = ["f115lp", "f125lp", "f140lp", "f150lp", "f165lp", "f122m"]

    jwst_nircam_sw = ["f150w2", "f070w", "f090w", "f115w", "f150w", "f200w",
                      "f140m", "f162m", "f182m", "f210m",
                      "f164n", "f187n", "f212n"]

    jwst_nircam_lw = ["f332w2", "f277w", "f356w", "f444w",
                      "f250m", "f300m", "f335m", "f360m", "f410m", "f430m", "f460m", "f480m",
                      "f323n", "f405n", "f466n", "f470n"]

    jwst_niriss = ["f090w", "f115w", "f150w", "f200w", "f277w", "f356w", "f444w",
                   "f140m", "f158m", "f380m", "f430m", "f480m"]

    jwst_miri = ["f560w", "f770w", "f1000w", "f1130w", "f1280w", "f1500w", "f1800w", "f2100w", "f2550w",
                 "f1065c", "f1140c", "f1550c", "f2300c"]

    roman_wfi = ["f062", "f087", "f106", "f129", "f158", "f184", "f213"]
    # fmt: on

    # galex
    galex = ["fuv", "nuv"]

    # Open hd5 file for writing
    hf = h5py.File(__ROOT__ + "filters.hd5", "w")

    # Create group for nice hierarchical structure
    f = hf.create_group("filters")

    # Define arrays for "contents" / descriptive information
    tablenames = []
    observatories = []
    instruments = []
    names = []
    norms = []
    cwaves = []
    pwaves = []
    comments = []

    # Loop through WFC3_UVIS filters
    for filt in wfc3_uvis:

        # define uvis 1 and uvis2 modes
        mode_1 = "wfc3, uvis1, " + filt
        mode_2 = "wfc3, uvis2, " + filt

        # pull bandpasses from stsynphot for the two uvis modes
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=AstropyUserWarning)
            bp_1 = stsyn.band(mode_1)
        bp_2 = stsyn.band(mode_2)

        # extract the wavelength array
        wave = bp_1.waveset

        # compute the average bandpass between uvis1 and uvis2
        bp_avg = np.average([bp_1(wave), bp_2(wave)], axis=0)

        # define the filter name
        filter_name = "HST_WFC3_" + filt.upper()

        # build array of wavelength and throughput
        arr = np.array(
            list(zip(wave.value.astype(np.float64), bp_avg.astype(np.float64))),
            dtype=[("WAVELENGTH", "float64"), ("THROUGHPUT", "float64")],
        )

        # append dataset to the hdf5 filters group
        f.create_dataset(filter_name, data=arr)

        # generate filter instance to compute relevant info
        newfilt = phot.Filter(wave, bp_avg, name=filt.upper())

        # populate contents lists with relevant information
        tablenames.append(filter_name)
        observatories.append("HST")
        instruments.append("WFC3")
        names.append(newfilt.name)
        norms.append(newfilt.norm.value)
        cwaves.append(newfilt.cl.value)
        pwaves.append(newfilt.lpivot.value)
        comments.append("avg of uvis1 and uvis2")

    # Loop through WFC3_IR filters
    for filt in wfc3_ir:

        # define ir mode
        mode = "wfc3, ir, " + filt

        # pull bandpasses from stsynphot for the two uvis modes
        bp = stsyn.band(mode)

        # extract the wavelength array
        wave = bp.waveset

        # define the filter name
        filter_name = "HST_WFC3_" + filt.upper()

        # build array of wavelength and throughput
        arr = np.array(
            list(zip(wave.value.astype(np.float64), bp(wave).astype(np.float64))),
            dtype=[("WAVELENGTH", "float64"), ("THROUGHPUT", "float64")],
        )

        # append dataset to the hdf5 filters group
        f.create_dataset(filter_name, data=arr)

        # generate filter instance to compute relevant info
        newfilt = phot.Filter(wave, bp(wave), name=filt.upper())

        # populate contents lists with relevant information
        tablenames.append(filter_name)
        observatories.append("HST")
        instruments.append("WFC3")
        names.append(newfilt.name)
        norms.append(newfilt.norm.value)
        cwaves.append(newfilt.cl.value)
        pwaves.append(newfilt.lpivot.value)
        comments.append("")

    # Loop through WFPC2 filters
    for filt in wfpc2:

        # define chips 1, 2, 3, 4 modes
        mode_1 = "wfpc2, 1, " + filt
        mode_2 = "wfpc2, 2, " + filt
        mode_3 = "wfpc2, 3, " + filt
        mode_4 = "wfpc2, 4, " + filt

        # pull bandpasses from stsynphot for the two uvis modes
        bp_1 = stsyn.band(mode_1)
        bp_2 = stsyn.band(mode_2)
        bp_3 = stsyn.band(mode_3)
        bp_4 = stsyn.band(mode_4)

        # extract the wavelength array
        wave = bp_1.waveset

        # compute the average bandpass between uvis1 and uvis2
        bp_avg = np.average([bp_1(wave), bp_2(wave), bp_3(wave), bp_4(wave)], axis=0)

        # define the filter name
        filter_name = "HST_WFPC2_" + filt.upper()

        # build array of wavelength and throughput
        arr = np.array(
            list(zip(wave.value.astype(np.float64), bp_avg.astype(np.float64))),
            dtype=[("WAVELENGTH", "float64"), ("THROUGHPUT", "float64")],
        )

        # append dataset to the hdf5 filters group
        f.create_dataset(filter_name, data=arr)

        # generate filter instance to compute relevant info
        newfilt = phot.Filter(wave, bp_avg, name=filt.upper())

        # populate contents lists with relevant information
        tablenames.append(filter_name)
        observatories.append("HST")
        instruments.append("WFPC2")
        names.append(newfilt.name)
        norms.append(newfilt.norm.value)
        cwaves.append(newfilt.cl.value)
        pwaves.append(newfilt.lpivot.value)
        comments.append("avg of 1, 2, 3, 4")

    # Loop through ACS WFC filters
    for filt in acs_wfc:

        # define wfc1, wfc2 modes
        mode_1 = "acs, wfc1, " + filt
        mode_2 = "acs, wfc2, " + filt

        # pull bandpasses from stsynphot for the two uvis modes
        bp_1 = stsyn.band(mode_1)
        bp_2 = stsyn.band(mode_2)

        # extract the wavelength array
        wave = bp_1.waveset

        # compute the average bandpass between uvis1 and uvis2
        bp_avg = np.average([bp_1(wave), bp_2(wave)], axis=0)

        # define the filter name
        filter_name = "HST_ACS_WFC_" + filt.upper()

        # build array of wavelength and throughput
        arr = np.array(
            list(zip(wave.value.astype(np.float64), bp_avg.astype(np.float64))),
            dtype=[("WAVELENGTH", "float64"), ("THROUGHPUT", "float64")],
        )

        # append dataset to the hdf5 filters group
        f.create_dataset(filter_name, data=arr)

        # generate filter instance to compute relevant info
        newfilt = phot.Filter(wave, bp_avg, name=filt.upper())

        # populate contents lists with relevant information
        tablenames.append(filter_name)
        observatories.append("HST")
        instruments.append("ACS_WFC")
        names.append(newfilt.name)
        norms.append(newfilt.norm.value)
        cwaves.append(newfilt.cl.value)
        pwaves.append(newfilt.lpivot.value)
        comments.append("avg of wfc1 and wfc2")

    # Loop through ACS SBC filters
    for filt in acs_sbc:

        # define ir mode
        mode = "acs, sbc, " + filt

        # pull bandpasses from stsynphot for the two uvis modes
        bp = stsyn.band(mode)

        # extract the wavelength array
        wave = bp.waveset

        # define the filter name
        filter_name = "HST_ACS_SBC_" + filt.upper()

        # build array of wavelength and throughput
        arr = np.array(
            list(zip(wave.value.astype(np.float64), bp(wave).astype(np.float64))),
            dtype=[("WAVELENGTH", "float64"), ("THROUGHPUT", "float64")],
        )

        # append dataset to the hdf5 filters group
        f.create_dataset(filter_name, data=arr)

        # generate filter instance to compute relevant info
        newfilt = phot.Filter(wave, bp(wave), name=filt.upper())

        # populate contents lists with relevant information
        tablenames.append(filter_name)
        observatories.append("HST")
        instruments.append("ACS")
        names.append(newfilt.name)
        norms.append(newfilt.norm.value)
        cwaves.append(newfilt.cl.value)
        pwaves.append(newfilt.lpivot.value)
        comments.append("")

    # Loop through GALEX filters:
    for filt in galex:
        # define ir mode
        mode = "galex," + filt

        # pull bandpasses from stsynphot for the two uvis modes
        bp = stsyn.band(mode)

        # extract the wavelength array
        wave = bp.waveset

        # define the filter name
        filter_name = "GALEX_" + filt.upper()

        # build array of wavelength and throughput
        arr = np.array(
            list(zip(wave.value.astype(np.float64), bp(wave).astype(np.float64))),
            dtype=[("WAVELENGTH", "float64"), ("THROUGHPUT", "float64")],
        )

        # append dataset to the hdf5 filters group
        f.create_dataset(filter_name, data=arr)

        # generate filter instance to compute relevant info
        newfilt = phot.Filter(wave, bp(wave), name=filt.upper())

        # populate contents lists with relevant information
        tablenames.append(filter_name)
        observatories.append("GALEX")
        instruments.append("GALEX")
        names.append(newfilt.name)
        norms.append(newfilt.norm.value)
        cwaves.append(newfilt.cl.value)
        pwaves.append(newfilt.lpivot.value)
        comments.append("")

    for filt in jwst_nircam_sw:
        # mock configuration
        conf = {
            "detector": {
                "nexp": 1,
                "ngroup": 10,
                "nint": 1,
                "readout_pattern": "deep2",
                "subarray": "full"
            },
            "instrument": {
                "aperture": "sw",
                "disperser": "null",
                "filter": filt,
                "instrument": "nircam",
                "mode": "sw_imaging"
            },
        }

        # create a configured instrument
        instrument_factory = InstrumentFactory(config=conf)

        # set up your wavelengths
        pwave = np.logspace(np.log10(0.5), np.log10(7.0), 501) * u.micron

        # get the throughput of the instrument over the desired wavelength range
        eff = instrument_factory.get_total_eff(pwave.value)

        # get wavelengths in Angstroms
        wave = pwave.to(u.AA)

        # define the filter name
        filter_name = "JWST_NIRCAM_" + filt.upper()

        # build array of wavelength and throughput
        arr = np.array(
            list(zip(wave.value.astype(np.float64), eff.astype(np.float64))),
            dtype=[("WAVELENGTH", "float64"), ("THROUGHPUT", "float64")],
        )

        # append dataset to the hdf5 filters group
        f.create_dataset(filter_name, data=arr)

        # generate filter instance to compute relevant info
        newfilt = phot.Filter(wave, eff, name=filt.upper())

        # populate contents lists with relevant information
        tablenames.append(filter_name)
        observatories.append("JWST")
        instruments.append("NIRCAM")
        names.append(newfilt.name)
        norms.append(newfilt.norm.value)
        cwaves.append(newfilt.cl.value)
        pwaves.append(newfilt.lpivot.value)
        comments.append("")

    for filt in jwst_nircam_lw:
        # mock configuration
        conf = {
            "detector": {
                "nexp": 1,
                "ngroup": 10,
                "nint": 1,
                "readout_pattern": "deep2",
                "subarray": "full",
            },
            "instrument": {
                "aperture": "lw",
                "disperser": "null",
                "filter": filt,
                "instrument": "nircam",
                "mode": "lw_imaging",
            },
        }

        # create a configured instrument
        instrument_factory = InstrumentFactory(config=conf)

        # set up your wavelengths
        pwave = np.logspace(np.log10(0.5), np.log10(7.0), 501) * u.micron

        # get the throughput of the instrument over the desired wavelength range
        eff = instrument_factory.get_total_eff(pwave.value)

        # get wavelengths in Angstroms
        wave = pwave.to(u.AA)

        # define the filter name
        filter_name = "JWST_NIRCAM_" + filt.upper()

        # build array of wavelength and throughput
        arr = np.array(
            list(zip(wave.value.astype(np.float64), eff.astype(np.float64))),
            dtype=[("WAVELENGTH", "float64"), ("THROUGHPUT", "float64")],
        )

        # append dataset to the hdf5 filters group
        f.create_dataset(filter_name, data=arr)

        # generate filter instance to compute relevant info
        newfilt = phot.Filter(wave, eff, name=filt.upper())

        # populate contents lists with relevant information
        tablenames.append(filter_name)
        observatories.append("JWST")
        instruments.append("NIRCAM")
        names.append(newfilt.name)
        norms.append(newfilt.norm.value)
        cwaves.append(newfilt.cl.value)
        pwaves.append(newfilt.lpivot.value)
        comments.append("")

    for filt in jwst_niriss:
        # mock configuration
        conf = {
            "detector": {
                "nexp": 1,
                "ngroup": 10,
                "nint": 1,
                "readout_pattern": "nis",
                "subarray": "full",
            },
            "instrument": {
                "aperture": "imager",
                "disperser": "null",
                "filter": filt,
                "instrument": "niriss",
                "mode": "imaging",
            },
        }

        # create a configured instrument
        instrument_factory = InstrumentFactory(config=conf)

        # set up your wavelengths
        pwave = np.logspace(np.log10(0.5), np.log10(7.0), 501) * u.micron

        # get the throughput of the instrument over the desired wavelength range
        eff = instrument_factory.get_total_eff(pwave.value)

        # get wavelengths in Angstroms
        wave = pwave.to(u.AA)

        # define the filter name
        filter_name = "JWST_NIRISS_" + filt.upper()

        # build array of wavelength and throughput
        arr = np.array(
            list(zip(wave.value.astype(np.float64), eff.astype(np.float64))),
            dtype=[("WAVELENGTH", "float64"), ("THROUGHPUT", "float64")],
        )

        # append dataset to the hdf5 filters group
        f.create_dataset(filter_name, data=arr)

        # generate filter instance to compute relevant info
        newfilt = phot.Filter(wave, eff, name=filt.upper())

        # populate contents lists with relevant information
        tablenames.append(filter_name)
        observatories.append("JWST")
        instruments.append("NIRISS")
        names.append(newfilt.name)
        norms.append(newfilt.norm.value)
        cwaves.append(newfilt.cl.value)
        pwaves.append(newfilt.lpivot.value)
        comments.append("")

    for filt in jwst_miri:
        # mock configuration
        conf = {
            "detector": {
                "nexp": 1,
                "ngroup": 10,
                "nint": 1,
                "readout_pattern": "fastr1",
                "subarray": "full",
            },
            "dynamic_scene": True,
            "instrument": {
                "aperture": "imager",
                "filter": filt,
                "instrument": "miri",
                "mode": "imaging",
            },
        }

        # create a configured instrument
        instrument_factory = InstrumentFactory(config=conf)

        # set up your wavelengths
        pwave = np.logspace(np.log10(3.0), np.log10(40.0), 501) * u.micron

        # get the throughput of the instrument over the desired wavelength range
        eff = instrument_factory.get_total_eff(pwave.value)

        # get wavelengths in Angstroms
        wave = pwave.to(u.AA)

        # define the filter name
        filter_name = "JWST_MIRI_" + filt.upper()

        # build array of wavelength and throughput
        arr = np.array(
            list(zip(wave.value.astype(np.float64), eff.astype(np.float64))),
            dtype=[("WAVELENGTH", "float64"), ("THROUGHPUT", "float64")],
        )

        # append dataset to the hdf5 filters group
        f.create_dataset(filter_name, data=arr)

        # generate filter instance to compute relevant info
        newfilt = phot.Filter(wave, eff, name=filt.upper())

        # populate contents lists with relevant information
        tablenames.append(filter_name)
        observatories.append("JWST")
        instruments.append("MIRI")
        names.append(newfilt.name)
        norms.append(newfilt.norm.value)
        cwaves.append(newfilt.cl.value)
        pwaves.append(newfilt.lpivot.value)
        comments.append("")

    # Loop through Roman WFI filters
    for filt in roman_wfi:

        # define ir mode
        mode = "roman, wfi, " + filt

        # pull bandpasses from stsynphot for the two uvis modes
        bp = stsyn.band(mode)

        # extract the wavelength array
        wave = bp.waveset

        # define the filter name
        filter_name = "ROMAN_WFI_" + filt.upper()

        # build array of wavelength and throughput
        arr = np.array(
            list(zip(wave.value.astype(np.float64), bp(wave).astype(np.float64))),
            dtype=[("WAVELENGTH", "float64"), ("THROUGHPUT", "float64")],
        )

        # append dataset to the hdf5 filters group
        f.create_dataset(filter_name, data=arr)

        # generate filter instance to compute relevant info
        newfilt = phot.Filter(wave, bp(wave), name=filt.upper())

        # populate contents lists with relevant information
        tablenames.append(filter_name)
        observatories.append("HST")
        instruments.append("ACS")
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
                instruments,
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


def make_vega_libfile():

    # Read an updated filters.hd5 lib
    __default_filtlist__ = __ROOT__ + "/filters.hd5"
    filtlist = QTable.read(__default_filtlist__, path="content")
    filters = [cfilt.decode("UTF-8") for cfilt in filtlist["TABLENAME"].data]

    # filenames for vega info
    __default_vega__ = __ROOT__ + "/vega.hd5"
    __default_vega_old__ = __ROOT__ + "/vega_old.hd5"

    # rename the current file so we can write a new version
    os.rename(__default_vega__, __default_vega_old__)

    # Get a spectrum from the need-to-be updated vega lib
    vega_old = QTable.read(__default_vega_old__, path="spectrum")
    vl = vega_old["WAVELENGTH"].data
    vf = vega_old["FLUX"].data

    # Write out an updated vega.hd5 file
    vega = h5py.File(__default_vega__, "w")

    vega.create_dataset("spectrum", data=vega_old)
    flist = phot.load_filters(
        filters, interp=True, lamb=vl, filterLib=__default_filtlist__
    )

    fname = []
    cwave = []
    lum = []
    mag = []

    for cfilt in flist:
        fname.append(cfilt.name)
        cwave.append(cfilt.cl)
        flux = cfilt.getFlux(vl, vf)
        lum.append(flux)
        mag.append(-2.5 * np.log10(flux))

    contents = np.array(
        list(
            zip(
                fname,
                cwave,
                lum,
                mag,
            )
        ),
        dtype=[
            ("FNAME", "S30"),
            ("CWAVE", "<f8"),
            ("LUM", "<f8"),
            ("MAG", "<f8"),
        ],
    )

    vega.create_dataset("sed", data=contents)
    vega.close()


if __name__ == "__main__":  # pragma: no cover
    make_filters_libfile()
    make_vega_libfile()
