# useful functions for BEAST tests
# put here instead of having in every tests
import os.path

import numpy as np
import h5py
import tempfile

from astropy.io import fits
from astropy.utils.data import download_file

__all__ = ["download_rename", "compare_tables", "compare_fits", "compare_hdf5"]


def download_rename(filename, tmpdir=""):
    """Download a file and rename it to have the right extension.
    Otherwise, downloaded file will not have an extension at all and an
    extension is needed for the BEAST.

    Parameters
    ----------
    filename : str
        name of file to download
    """
    url_loc = "http://www.stsci.edu/~kgordon/beast/"
    fname_dld = download_file("%s%s" % (url_loc, filename))
    extension = filename.split(".")[-1]
    fname = tempfile.NamedTemporaryFile(suffix=f".{extension}").name
    os.rename(fname_dld, fname)
    return fname


def compare_tables(table_cache, table_new, rtol=1e-7, otag=""):
    """
    Compare two tables using astropy tables routines.

    Parameters
    ----------
    table_cache : astropy table
    table_new : astropy table
        data for comparision.
    rtol : float, optional (default=1e-7)
        relative tolerance for np.testing.assert_allclose (default of 1e-7
        matches the default in assert_allclose)
    otag : str
        extra string for error messages
    """
    etag = ""
    if otag != "":
        etag = f"{otag}:"
    if not len(table_new) == len(table_cache):
        raise AssertionError(
            "{0}Table lengths don't match: {1} (new) and {2} (cached)".format(
                etag, len(table_new), len(table_cache)
            )
        )

    if not len(table_new.colnames) == len(table_cache.colnames):
        raise AssertionError(
            "{0} Table number of colums don't match: {1} (new) and {2} (cached)".format(
                etag, len(table_new.colnames), len(table_cache.colnames)
            )
        )

    for tcolname in table_new.colnames:
        # test numerical types for closeness
        # and other types for equality
        if table_new[tcolname].data.dtype.kind in ["f", "i"]:
            np.testing.assert_allclose(
                table_new[tcolname],
                table_cache[tcolname],
                rtol=rtol,
                err_msg=("%s columns not equal" % tcolname),
            )
        else:
            np.testing.assert_equal(
                table_new[tcolname],
                table_cache[tcolname],
                err_msg=("%s columns not equal" % tcolname),
            )


def compare_fits(fname_cache, fname_new):
    """
    Compare two FITS files.

    Parameters
    ----------
    fname_cache : str
    fname_new : type
        names to FITS files
    """
    fits_cache = fits.open(fname_cache)
    fits_new = fits.open(fname_new)

    if not len(fits_new) == len(fits_cache):
        raise AssertionError(
            "FITS data lengths don't match: {0} (new) and {1} (cached)".format(
                len(fits_new), len(fits_cache)
            )
        )

    for k in range(1, len(fits_new)):
        qname = fits_new[k].header["EXTNAME"]
        np.testing.assert_allclose(
            fits_new[k].data,
            fits_cache[qname].data,
            err_msg=("%s FITS extension not equal" % qname),
        )

    fits_cache.close()
    fits_new.close()


def compare_hdf5(fname_cache, fname_new, ctype=None):
    """
    Compare two hdf files.

    Parameters
    ----------
    fname_cache : str
    fname_new : type
        names to hdf5 files

    ctype : str
        if set, string to identify the type of data being tested
    """
    hdf_cache = h5py.File(fname_cache, "r")
    hdf_new = h5py.File(fname_new, "r")

    # go through the file and check if it is exactly the same

    for sname in hdf_cache.keys():
        if isinstance(hdf_cache[sname], h5py.Dataset):
            cvalue = hdf_cache[sname]
            cvalue_new = hdf_new[sname]
            if ctype is not None:
                osname = "%s/%s" % (ctype, sname)
            else:
                osname = sname
            if cvalue.dtype.fields is None:
                np.testing.assert_allclose(
                    cvalue[()],
                    cvalue_new[()],
                    err_msg="testing %s" % (osname),
                    rtol=2e-6,
                )
            else:
                for ckey in cvalue.dtype.fields.keys():
                    err_msg = "testing %s/%s" % (osname, ckey)
                    np.testing.assert_allclose(
                        cvalue[()][ckey],
                        cvalue_new[()][ckey],
                        err_msg=err_msg,
                        rtol=1e-5,
                    )
