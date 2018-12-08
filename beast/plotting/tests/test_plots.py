import matplotlib.pyplot as plt

import os.path
import pytest

import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.utils.data import download_file
from astropy.tests.helper import remote_data

# from beast.plotting.plot_indiv_fit import plot_indiv_fit
from beast.plotting import plot_indiv_fit
# from beast.plotting import plot_filters

plt.switch_backend('agg')


def _download_rename(filename):
    """
    Download a file and rename it to have the right extension

    Otherwise, downloaded file will not have an extension at all
    """
    url_loc = 'http://www.stsci.edu/~kgordon/beast/'
    fname_dld = download_file('%s%s' % (url_loc, filename))
    extension = filename.split('.')[-1]
    fname = '%s.%s' % (fname_dld, extension)
    os.rename(fname_dld, fname)
    return fname


# @pytest.mark.skip(reason="awaiting resolution of pytest-mpl")
# @remote_data
@pytest.mark.mpl_image_compare
def test_indiv_plot():

    # download cached version of fitting results
    stats_fname_cache = _download_rename('beast_example_phat_stats.fits')
    pdf1d_fname_cache = _download_rename('beast_example_phat_pdf1d.fits')

    # results_dir = '../../examples/phat_small/beast_example_phat/'
    # stats_fname_cache = results_dir + 'beast_example_phat_stats.fits'
    # pdf1d_fname_cache = results_dir + 'beast_example_phat_pdf1d.fits'

    starnum = 0

    # read in the stats
    stats = Table.read(stats_fname_cache)
    # open 1D PDF file
    pdf1d_hdu = fits.open(pdf1d_fname_cache)

    filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_ACS_WFC_F475W',
               'HST_ACS_WFC_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']
    waves = np.asarray([2722.05531502, 3366.00507206, 4763.04670013,
                        8087.36760191, 11672.35909295, 15432.7387546])

    fig, ax = plt.subplots(figsize=(8, 8))

    # make the plot!
    plot_indiv_fit.plot_beast_ifit(filters, waves, stats, pdf1d_hdu, starnum)

    return fig
