import os.path
import filecmp

import numpy as np
import pytest

from astropy.utils.data import download_file
from astropy.tests.helper import remote_data
from astropy.table import Table

from beast.physicsmodel.stars import isochrone

@remote_data
#@pytest.mark.skip(reason="temporarily disable")
def test_padova_isochrone_download():
    # download the cached version
    filename = download_file('http://www.stsci.edu/~kgordon/beast/beast_example_phat_iso_new.csv', cache=True)
    table_cache = Table.read(filename, format='ascii.csv', comment='#',
                             delimiter=',')
    
    # initialize isochrone
    oiso = isochrone.PadovaWeb()

    # download from the padova website
    t = oiso._get_t_isochrones(6.0, 10.13, 1.0, [0.03, 0.019, 0.008, 0.004])
    t.header['NAME'] = 'Test of Cached Isochrones'

    # save
    savename = '/tmp/padova_iso.csv'
    t.write(savename)

    # get the new table
    table_new = Table.read(savename, format='ascii.csv', comment='#',
                             delimiter=',')
    
    # compare
    assert len(table_new) == len(table_cache)

    for tcolname in table_new.colnames:
        np.testing.assert_equal(table_new[tcolname], table_cache[tcolname],
                                '%s columns not equal'%tcolname)

                       
