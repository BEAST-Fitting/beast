import matplotlib.pyplot as plt
import pytest

from beast.plotting import plot_indiv_fit, plot_cmd, plot_cmd_with_fits, plot_filters
from beast.tests.helpers import download_rename

plt.switch_backend("agg")


@pytest.mark.remote_data
@pytest.mark.mpl_image_compare(tolerance=25)
def test_indiv_plot():

    # download cached version of fitting results
    stats_fname_cache = download_rename("phat_small/beast_example_phat_stats.fits")
    pdf1d_fname_cache = download_rename("phat_small/beast_example_phat_pdf1d.fits")

    # make the plot!
    fig = plot_indiv_fit.plot_indiv_fit(
        [stats_fname_cache, pdf1d_fname_cache], 0, plotfig=False
    )

    return fig


@pytest.mark.remote_data
@pytest.mark.mpl_image_compare(tolerance=10)
def test_plot_cmd():

    # Download example data from phat_small
    fitsfile = download_rename("phat_small/b15_4band_det_27_A.fits")

    # Plot CMD using defaults
    fig = plot_cmd.plot_cmd(fitsfile, show_plot=False)

    return fig


@pytest.mark.remote_data
@pytest.mark.mpl_image_compare(tolerance=55)
def test_plot_cmd_with_fits():

    # Download example data from phat_small
    fitsfile = download_rename("phat_small/b15_4band_det_27_A.fits")

    # Download BEAST fits to example data
    beast_fitsfile = download_rename("phat_small/beast_example_phat_stats.fits")

    # Plot CMD using defaults
    fig = plot_cmd_with_fits.plot_cmd_with_fits(fitsfile, beast_fitsfile)

    return fig


@pytest.mark.remote_data
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
