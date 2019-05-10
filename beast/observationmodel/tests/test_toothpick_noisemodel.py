from astropy.tests.helper import remote_data

from ..noisemodel import generic_noisemodel as noisemodel
from ..noisemodel.absflux_covmat import hst_frac_matrix
from ...physicsmodel.grid import FileSEDGrid
from beast.tests.helpers import (download_rename, compare_hdf5)


@remote_data
def test_toothpick_noisemodel():

    # download the needed files
    asts_fname = download_rename('fake_stars_b15_27_all.hd5')
    filter_fname = download_rename('filters.hd5')
    vega_fname = download_rename('vega.hd5')
    hst_fname = download_rename('hst_whitedwarf_frac_covar.fits')
    seds_fname = download_rename('beast_example_phat_seds.grid.hd5')

    # download cached version of noisemodel on the sed grid
    noise_fname_cache = download_rename('beast_example_phat_noisemodel.hd5')

    ################
    # get the modesedgrid on which to generate the noisemodel
    modelsedgrid = FileSEDGrid(seds_fname)

    # absflux calibration covariance matrix for HST specific filters (AC)
    filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_ACS_WFC_F475W',
               'HST_ACS_WFC_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']
    absflux_a_matrix = hst_frac_matrix(filters,
                                       hst_fname=hst_fname,
                                       filterLib=filter_fname)

    # generate the AST noise model
    noise_fname = '/tmp/beast_example_phat_noisemodel.hd5'
    noisemodel.make_toothpick_noise_model(
                            noise_fname,
                            asts_fname,
                            modelsedgrid,
                            absflux_a_matrix=absflux_a_matrix,
                            vega_fname=vega_fname,
                            use_rate=False)

    # compare the new to the cached version
    compare_hdf5(noise_fname_cache, noise_fname)
