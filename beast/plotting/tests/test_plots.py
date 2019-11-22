import matplotlib.pyplot as plt

import pytest

import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.tests.helper import remote_data

from beast.plotting import plot_indiv_fit, plot_cmd, plot_cmd_with_fits, plot_filters
from beast.tests.helpers import download_rename

plt.switch_backend("agg")


@remote_data
@pytest.mark.mpl_image_compare(tolerance=25)
def test_indiv_plot():

    # download cached version of fitting results
    stats_fname_cache = download_rename("beast_example_phat_stats.fits")
    pdf1d_fname_cache = download_rename("beast_example_phat_pdf1d.fits")

    starnum = 0

    # read in the stats
    stats = Table.read(stats_fname_cache)
    # open 1D PDF file
    pdf1d_hdu = fits.open(pdf1d_fname_cache)

    filters = [
        "HST_WFC3_F275W",
        "HST_WFC3_F336W",
        "HST_ACS_WFC_F475W",
        "HST_ACS_WFC_F814W",
        "HST_WFC3_F110W",
        "HST_WFC3_F160W",
    ]
    waves = np.asarray(
        [
            2722.05531502,
            3366.00507206,
            4763.04670013,
            8087.36760191,
            11672.35909295,
            15432.7387546,
        ]
    )

    fig, ax = plt.subplots(figsize=(8, 8))

    # make the plot!
    plot_indiv_fit.plot_beast_ifit(filters, waves, stats, pdf1d_hdu, starnum)

    return fig


@remote_data
@pytest.mark.mpl_image_compare(tolerance=10)
def test_plot_cmd():

    # Download example data from phat_small
    fitsfile = download_rename("b15_4band_det_27_A.fits")

    # Plot CMD using defaults
    fig = plot_cmd.plot_cmd(fitsfile, show_plot=False)

    return fig


@remote_data
@pytest.mark.mpl_image_compare(tolerance=55)
def test_plot_cmd_with_fits():

    # Download example data from phat_small
    fitsfile = download_rename("b15_4band_det_27_A.fits")

    # Download BEAST fits to example data
    beast_fitsfile = download_rename("beast_example_phat_stats.fits")

    # Plot CMD using defaults
    fig = plot_cmd_with_fits.plot(fitsfile, beast_fitsfile)

    return fig


@remote_data
@pytest.mark.mpl_image_compare(tolerance=18)
def test_plot_filters():

    filter_names = [
        "HST_WFC3_F225W",
        "HST_WFC3_F275W",
        "HST_WFC3_F336W",
        "HST_ACS_WFC_F475W",
        "HST_ACS_WFC_F550M",
        "HST_ACS_WFC_F814W",
        "HST_WFC3_F110W",
        "HST_WFC3_F160W",
    ]

    filters = download_rename("filters.hd5")

    # Plot filters using above arguments (the defaults)
    fig = plot_filters.plot_filters(filter_names, filterLib=filters, show_plot=False)

    return fig
