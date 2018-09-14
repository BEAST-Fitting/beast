# useful functions for BEAST tests
# put here instead of having in every tests
import os.path

import numpy as np

from astropy.utils.data import download_file

__all__ = ['download_rename']


def download_rename(filename):
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


def compare_tables(table_cache, table_new):
    """
    Compare two tables using astropy tables routines

    """
    assert len(table_new) == len(table_cache)

    for tcolname in table_new.colnames:
        np.testing.assert_equal(table_new[tcolname], table_cache[tcolname],
                                '%s columns not equal' % tcolname)
