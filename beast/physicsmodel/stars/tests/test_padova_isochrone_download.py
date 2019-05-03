from astropy.tests.helper import remote_data
from astropy.table import Table

from ...model_grid import make_iso_table
from beast.tests.helpers import (download_rename, compare_tables)


@remote_data
def test_padova_isochrone_download():
    # download the cached version
    iso_fname_cache = download_rename('beast_example_phat_iso.csv')

    # download the file live from the website
    savename = '/tmp/padova_iso.csv'
    iso_fname, oiso = make_iso_table('test', iso_fname=savename,
                                     logtmin=6.0, logtmax=10.13, dlogt=1.0,
                                     z=[0.03, 0.019, 0.008, 0.004])

    # read the cached and new tables using astropy tables
    table_cache = Table.read(iso_fname_cache, format='ascii.csv', comment='#',
                             delimiter=',')
    table_new = Table.read(iso_fname, format='ascii.csv', comment='#',
                           delimiter=',')

    # compare
    compare_tables(table_cache, table_new)
