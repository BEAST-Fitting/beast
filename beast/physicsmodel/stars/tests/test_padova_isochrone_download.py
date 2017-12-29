import os.path
import filecmp

import numpy as np
import pytest

from astropy.utils.data import download_file
from astropy.tests.helper import remote_data
from astropy.table import Table

from ...model_grid import make_iso_table

@remote_data
#@pytest.mark.skip(reason="temporarily disable")
def test_padova_isochrone_download():
    # download the cached version
    url_loc = 'http://www.stsci.edu/~kgordon/beast/'
    filename = download_file('%s%s'%(url_loc,
                                     'beast_example_phat_iso.csv'))
    table_cache = Table.read(filename, format='ascii.csv', comment='#',
                             delimiter=',')

    savename = '/tmp/padova_iso.csv'
    iso_fname, oiso = make_iso_table('test', iso_fname=savename,
                                     logtmin=6.0, logtmax=10.13, dlogt=1.0,
                                     z=[0.03, 0.019, 0.008, 0.004])

    # get the new table
    table_new = Table.read(iso_fname, format='ascii.csv', comment='#',
                             delimiter=',')
    
    # compare
    assert len(table_new) == len(table_cache)

    for tcolname in table_new.colnames:
        np.testing.assert_equal(table_new[tcolname], table_cache[tcolname],
                                '%s columns not equal'%tcolname)

                       
