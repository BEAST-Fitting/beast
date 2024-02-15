# script to download the BEAST library files
from shutil import copyfile

from astropy.utils.data import download_file

from beast.config import libs_server, libs, __ROOT__


def _download_rename(filename, url_loc, local_loc):
    """
    Download a file and rename it to have correct name and location
    """
    fname_dld = download_file("%s%s" % (url_loc, filename))
    copyfile(fname_dld, local_loc + filename)
    return filename


def get_libfiles(vega_filters_only=False):
    """
    Download all the library files needed by the BEAST
    """
    for ckey, clib in libs.items():
        if ((not vega_filters_only) or
            (vega_filters_only & (ckey in ["vega", "filters"]))):
            _download_rename(clib, libs_server, __ROOT__)


if __name__ == "__main__":  # pragma: no cover
    get_libfiles()
