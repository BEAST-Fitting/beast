import os.path

import numpy as np
import h5py

from astropy.utils.data import download_file
from astropy.tests.helper import remote_data

from ..noisemodel import generic_noisemodel as noisemodel
from ..noisemodel.absflux_covmat import hst_frac_matrix
from ...physicsmodel.grid import FileSEDGrid
from astropy.utils.data import download_file
from ...observationmodel.ast import make_ast_input_list


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
def test_pick_models():
    # download the needed files
    vega_fname = _download_rename('vega.hd5')

    # download cached version of noisemodel on the sed grid
    filename = _download_rename('beast_example_phat_noisemodel.hd5')
    filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_ACS_WFC_F475W',
               'HST_ACS_WFC_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']
    mag_cuts = [1.]
    make_ast_input_list.pick_models(filename, filters, mag_cuts, vega_fname=vega_fname)
    if len(make_ast_input_list) < 1:
        raise ValueError('fake star catalog is empty')


if __name__ == '__main__':

    test_ast_input_list()
