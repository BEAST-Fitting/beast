import os.path

import numpy as np
import h5py
import shutil

from astropy.table import Table
from astropy.utils.data import download_file
from astropy.tests.helper import remote_data

from beast.observationmodel.noisemodel import generic_noisemodel as noisemodel
from beast.observationmodel.noisemodel.absflux_covmat import hst_frac_matrix
from beast.physicsmodel.grid import FileSEDGrid
from astropy.utils.data import download_file
from beast.observationmodel.ast import make_ast_input_list
from beast.tests.helpers import (download_rename, compare_tables)


@remote_data
def test_pick_models():
    # download the needed files
    vega_fname = download_rename('vega.hd5')
    filename = download_rename('beast_example_phat_seds.grid.hd5')

    filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_ACS_WFC_F475W',
               'HST_ACS_WFC_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']
    mag_cuts = [1.]

    table_out = make_ast_input_list.pick_models(filename, filters, mag_cuts, vega_fname=vega_fname,
                                                outfile='/tmp/test_inputAST.txt', ranseed=1234)

    table_new = Table.read('/tmp/test_inputAST.txt', format='ascii')

    # download cached version of the file and compare it to new file
    cached_table_filename = download_rename('cache_inputAST.txt')
    table_cache = Table.read(cached_table_filename, format='csv', delimiter=' ')
    compare_tables(table_new, table_cache)


if __name__ == '__main__':
    test_pick_models()
