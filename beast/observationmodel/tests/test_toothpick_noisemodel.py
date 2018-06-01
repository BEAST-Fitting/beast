import os.path

import numpy as np
import h5py

from astropy.utils.data import download_file
from astropy.tests.helper import remote_data

from ..noisemodel import generic_noisemodel as noisemodel
from ..noisemodel.absflux_covmat import hst_frac_matrix
from ...physicsmodel.grid import FileSEDGrid


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


@remote_data
# @pytest.mark.skip(reason="temporarily disable")
def test_toothpick_noisemodel():

    # download the needed files
    asts_fname = _download_rename('fake_stars_b15_27_all.hd5')
    filter_fname = _download_rename('filters.hd5')
    vega_fname = _download_rename('vega.hd5')
    hst_fname = _download_rename('hst_whitedwarf_frac_covar.fits')
    seds_fname = _download_rename('beast_example_phat_seds.grid.hd5')

    # download cached version of noisemodel on the sed grid
    filename = _download_rename('beast_example_phat_noisemodel.hd5')

    hdf_cache = h5py.File(filename, 'r')

    ################
    # get the modesedgrid on which to generate the noisemodel
    modelsedgrid = FileSEDGrid(seds_fname)

    # absflux calibration covariance matrix for HST specific filters (AC)
    filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_ACS_WFC_F475W',
               'HST_ACS_WFC_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']
    absflux_a_matrix = hst_frac_matrix(filters,
                                       hst_fname=hst_fname,
                                       filterLib=filter_fname)

#    ast_colnames = np.array(['F275W', 'F336W', 'F475W',
#                             'F814W', 'F110W', 'F160W'])
    # generate the AST noise model
    noise_fname = '/tmp/beast_example_phat_noisemodel.hd5'
    noisemodel.make_toothpick_noise_model(
                            noise_fname,
                            asts_fname,
                            modelsedgrid,
                            absflux_a_matrix=absflux_a_matrix,
                            vega_fname=vega_fname)

    # open the hdf file with the specral grid
    hdf_new = h5py.File(noise_fname, 'r')

    # go through the file and check if it is exactly the same
    for sname in hdf_cache.keys():
        if isinstance(hdf_cache[sname], h5py.Dataset):
            cvalue = hdf_cache[sname]
            cvalue_new = hdf_new[sname]
            if cvalue.dtype.fields is None:
                np.testing.assert_allclose(cvalue.value, cvalue_new.value,
                                           err_msg='testing %s' % (sname),
                                           rtol=1e-6)
            else:
                for ckey in cvalue.dtype.fields.keys():
                    err_msg = 'testing %s/%s' % (sname, ckey)
                    np.testing.assert_allclose(cvalue.value[ckey],
                                               cvalue_new.value[ckey],
                                               err_msg=err_msg,
                                               rtol=1e-5)


if __name__ == '__main__':

    test_toothpick_noisemodel()
