import os.path

import numpy as np
import h5py

from astropy.utils.data import download_file
from astropy.tests.helper import remote_data

from .. import grid
from ..model_grid import make_extinguished_sed_grid
from ..dust import extinction


@remote_data
# @pytest.mark.skip(reason="temporarily disable")
def test_make_extinguished_sed_grid():

    # download the needed files
    url_loc = 'http://www.stsci.edu/~kgordon/beast/'
    priors_fname_dld = download_file('%s%s' % (url_loc,
                                'beast_example_phat_spec_w_priors.grid.hd5'))
    filter_fname_dld = download_file('%s%s' % (url_loc, 'filters.hd5'))

    # rename files to have the correct extensions
    priors_fname = '%s.hd5' % (priors_fname_dld)
    os.rename(priors_fname_dld, priors_fname)
    filter_fname = '%s.hd5' % (filter_fname_dld)
    os.rename(filter_fname_dld, filter_fname)

    # download cached version of sed grid
    filename = download_file('%s%s' % (url_loc,
                                       'beast_example_phat_seds.grid.hd5'))

    hdf_cache = h5py.File(filename, 'r')

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

    # open the hdf file with the specral grid
    hdf_new = h5py.File(seds_fname, 'r')

    # go through the file and check if it is exactly the same
    for sname in hdf_cache.keys():
        if isinstance(hdf_cache[sname], h5py.Dataset):
            cvalue = hdf_cache[sname]
            cvalue_new = hdf_new[sname]
            if cvalue.dtype.fields is None:
                np.testing.assert_equal(cvalue.value, cvalue_new.value,
                                        'testing %s' % (sname))
            else:
                for ckey in cvalue.dtype.fields.keys():
                    np.testing.assert_equal(cvalue.value[ckey],
                                            cvalue_new.value[ckey],
                                            'testing %s/%s' % (sname, ckey))


if __name__ == '__main__':

    test_make_extinguished_sed_grid()
