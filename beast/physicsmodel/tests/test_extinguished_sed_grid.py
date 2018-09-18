from astropy.tests.helper import remote_data

from .. import grid
from ..model_grid import make_extinguished_sed_grid
from ..dust import extinction
from beast.tests.helpers import (download_rename, compare_hdf5)


@remote_data
def test_make_extinguished_sed_grid():

    # download the needed files
    priors_fname = download_rename('beast_example_phat_spec_w_priors.grid.hd5')
    filter_fname = download_rename('filters.hd5')

    # download cached version of sed grid
    seds_fname_cache = download_rename('beast_example_phat_seds.grid.hd5')

    ################
    # generate the same extinguished SED grid from the code

    # Add in the filters
    filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_ACS_WFC_F475W',
               'HST_ACS_WFC_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']
    add_spectral_properties_kwargs = dict(filternames=filters)

    g_pspec = grid.FileSpectralGrid(priors_fname, backend='memory')

    # generate the SED grid by integrating the filter response functions
    #   effect of dust extinction applied before filter integration
    #   also computes the dust priors as weights
    seds_fname = '/tmp/beast_example_phat_sed.grid.hd5'
    seds_fname, g_seds = make_extinguished_sed_grid(
        'test',
        g_pspec,
        filters,
        seds_fname=seds_fname,
        filterLib=filter_fname,
        extLaw=extinction.Gordon16_RvFALaw(),
        av=[0.0, 10.055, 1.0],
        rv=[2.0, 6.0, 1.0],
        fA=[0.0, 1.0, 0.25],
        av_prior_model={'name': 'flat'},
        rv_prior_model={'name': 'flat'},
        fA_prior_model={'name': 'flat'},
        add_spectral_properties_kwargs=add_spectral_properties_kwargs)

    # compare the new to the cached version
    compare_hdf5(seds_fname_cache, seds_fname)
